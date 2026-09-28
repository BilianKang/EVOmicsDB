#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg)) {
  script_dir <- dirname(normalizePath(gsub("~\\+~", " ", sub("^--file=", "", script_arg[1]))))
  project_lib <- normalizePath(file.path(script_dir, "..", "local_R_libs"),
                               mustWork = FALSE)
  if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
} else {
  script_dir <- getwd()
}
source(file.path(script_dir, "evomics_explore_utils.R"))

# ============================================================
# go_network.r
# GO 分层网络 (mRNA + Protein)
# 三层结构: TopCore → SmallCore → GO Term
# 支持组合: mRNA + Protein
# ============================================================

library(clusterProfiler)
library(dplyr)
library(scales)
library(visNetwork)
library(org.Hs.eg.db)
library(jsonlite)
library(argparse)

parser <- ArgumentParser(description = "Build GO hierarchical network (mRNA + Protein)")
# 基础参数
parser$add_argument("--mrna_dataset", type = "character", required = TRUE, help = "Path to mRNA dataset file")
parser$add_argument("--protein_dataset", type = "character", required = TRUE, help = "Path to Protein dataset file")
parser$add_argument("--pvalue_threshold", type = "double", default = 0.05, help = "P-value cutoff")
parser$add_argument("--feature_pvalue_threshold", type = "double", default = NULL,
                    help = "Feature-level raw/adjusted P cutoff")
parser$add_argument("--term_fdr_threshold", type = "double", default = NULL,
                    help = "GO term-level adjusted FDR cutoff")
parser$add_argument("--p_adjust_method", type = "character", default = "BH", help = "P-value adjustment method")
parser$add_argument("--use_padj_for_features", type = "logical", default = TRUE)
parser$add_argument("--mrna_universe_file", type = "character", default = NULL,
                    help = "One detected mRNA gene symbol per line")
parser$add_argument("--protein_universe_file", type = "character", default = NULL,
                    help = "One detected protein gene symbol per line")
parser$add_argument("--common_only", type = "logical", default = FALSE,
                    help = "If TRUE, retain only terms enriched in both assays")
parser$add_argument("--ont", type = "character", default = "ALL", help = "GO ontology: ALL/BP/MF/CC")
parser$add_argument("--output_html", type = "character", required = TRUE, help = "Output HTML file path")
parser$add_argument("--output_csv", type = "character", required = TRUE, help = "Output CSV file path")
parser$add_argument("--output_png", type = "character", default = NULL)
parser$add_argument("--output_pdf", type = "character", default = NULL)
parser$add_argument("--status_json", type = "character", default = NULL,
                    help = "Structured success/empty/error status output")
parser$add_argument("--static_top_n", type = "integer", default = 30,
                    help = "Maximum GO terms in static PNG/PDF; HTML and CSV retain all terms.")
# Titles
parser$add_argument("--main_title", type = "character", default = "GO Hierarchical Network", help = "Main title")
parser$add_argument("--main_title_size", type = "integer", default = 14, help = "Main title font size")
parser$add_argument("--axis_text_size", type = "integer", default = 10, help = "Node label font size")
# Dot (节点颜色)
parser$add_argument("--BP_fill_color", type = "character", default = "#00007F", help = "BP SmallCore fill color")
parser$add_argument("--CC_fill_color", type = "character", default = "#B22222", help = "CC SmallCore fill color")
parser$add_argument("--MF_fill_color", type = "character", default = "#3B4CC0", help = "MF SmallCore fill color")
# GO Term 节点颜色渐变 (按 Count 映射)
parser$add_argument("--low_color", type = "character", default = "#deebf7", help = "Low Count color")
parser$add_argument("--mid_color", type = "character", default = "#9ecae1", help = "Mid Count color")
parser$add_argument("--high_color", type = "character", default = "#3182bd", help = "High Count color")
parser$add_argument("--border_color", type = "character", default = "black", help = "Border color")
parser$add_argument("--border_size", type = "double", default = 1.0, help = "Border size")
# TopCore 颜色
parser$add_argument("--rna_core_color", type = "character", default = "#9467bd", help = "RNA_core color")
parser$add_argument("--protein_core_color", type = "character", default = "#2ca02c", help = "Protein_core color")
parser$add_argument("--common_core_color", type = "character", default = "#d62728", help = "Common_core color")
# Legend
parser$add_argument("--is_show_legend", type = "character", default = "TRUE", help = "Show legend (TRUE/FALSE)")
parser$add_argument("--legend_title_size", type = "integer", default = 10, help = "Legend title font size")
parser$add_argument("--legend_text_size", type = "integer", default = 9, help = "Legend text font size")
parser$add_argument("--legend_position", type = "character", default = "right", help = "Legend position")
# Export settings
parser$add_argument("--is_show_border", type = "character", default = "TRUE", help = "Show border (TRUE/FALSE)")
parser$add_argument("--width", type = "integer", default = 1200, help = "Output width (px)")
parser$add_argument("--height", type = "integer", default = 800, help = "Output height (px)")
args <- parser$parse_args()

if (!is.null(args$status_json) && file.exists(args$status_json)) unlink(args$status_json)
resolve_go_threshold <- function(value, legacy) {
  if (is.null(value) || !length(value) || is.na(value)) as.numeric(legacy) else as.numeric(value)
}
feature_pvalue_threshold <- resolve_go_threshold(args$feature_pvalue_threshold, args$pvalue_threshold)
term_fdr_threshold <- resolve_go_threshold(args$term_fdr_threshold, args$pvalue_threshold)
if (!is.finite(feature_pvalue_threshold) || feature_pvalue_threshold <= 0 || feature_pvalue_threshold > 1) {
  stop("feature_pvalue_threshold must be > 0 and <= 1")
}
if (!is.finite(term_fdr_threshold) || term_fdr_threshold <= 0 || term_fdr_threshold > 1) {
  stop("term_fdr_threshold must be > 0 and <= 1")
}
write_go_status <- function(status, code, message_text, counts = list(), artifacts = list()) {
  evomics_status_write(args$status_json, "go_network", status, code, message_text,
                       counts, artifacts)
}
go_empty <- function(code, message_text, counts = list(), artifacts = list()) {
  # Keep the empty-result download contract uniform across early exits.  A
  # later no_display_terms branch may already have written the complete CSV;
  # early no-significant/no-mapped branches receive a small typed export and
  # an explanatory HTML page instead of returning no artifacts at all.
  if (!file.exists(args$output_csv)) {
    empty_terms <- data.frame(
      term = character(), Description = character(), Ontology = character(),
      source = character(), pvalue = numeric(), p.adjust = numeric(),
      stringsAsFactors = FALSE
    )
    utils::write.table(empty_terms, file = args$output_csv, sep = "\t",
                       quote = FALSE, row.names = FALSE)
  }
  if (!file.exists(args$output_html)) {
    html <- paste0(
      "<!doctype html><html><head><meta charset='utf-8'><title>No GO terms</title>",
      "</head><body><h2>No GO enrichment result</h2><p>",
      gsub("[<&>]", "", as.character(message_text)), "</p></body></html>"
    )
    writeLines(html, args$output_html, useBytes = TRUE)
  }
  artifacts <- c(list(html = args$output_html, csv = args$output_csv), artifacts)
  evomics_status_empty(args$status_json, "go_network", code, message_text,
                       counts, artifacts = artifacts)
}

# Inline htmlwidget assets so a valid GO result does not depend on a system
# Pandoc installation. This is artifact packaging only; the GO enrichment,
# hierarchy, and displayed terms are unchanged.
save_portable_widget <- function(widget, output_html) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE) ||
      !requireNamespace("base64enc", quietly = TRUE)) {
    stop("The server is missing the htmlwidgets/base64enc export dependency.")
  }
  temp_root <- tempfile("go_widget_")
  dir.create(temp_root, recursive = TRUE)
  on.exit(unlink(temp_root, recursive = TRUE, force = TRUE), add = TRUE)
  temp_html <- file.path(temp_root, "widget.html")
  htmlwidgets::saveWidget(widget, file = temp_html, selfcontained = FALSE, libdir = "lib")
  html <- readLines(temp_html, warn = FALSE, encoding = "UTF-8")
  inline_asset <- function(line, attribute, mime_type) {
    pattern <- paste0(attribute, "=\"([^\"]+)\"")
    match <- regexec(pattern, line, perl = TRUE)
    values <- regmatches(line, match)[[1]]
    if (length(values) < 2 || grepl("^(data:|https?:|//)", values[[2]])) return(line)
    asset_path <- file.path(temp_root, utils::URLdecode(values[[2]]))
    if (!file.exists(asset_path)) return(line)
    data_uri <- base64enc::dataURI(file = asset_path, mime = mime_type)
    sub(pattern, paste0(attribute, "=\"", data_uri, "\""), line, perl = TRUE)
  }
  for (index in seq_along(html)) {
    if (grepl("<script[^>]+src=\"", html[[index]], perl = TRUE))
      html[[index]] <- inline_asset(html[[index]], "src", "application/javascript")
    if (grepl("<link[^>]+href=\"", html[[index]], perl = TRUE))
      html[[index]] <- inline_asset(html[[index]], "href", "text/css")
  }
  writeLines(html, output_html, useBytes = TRUE)
}

# Parse boolean arguments
args$is_show_legend <- toupper(args$is_show_legend) == "TRUE"
args$is_show_border <- toupper(args$is_show_border) == "TRUE"

# ======================== 辅助函数 ========================

get_significant <- function(df, pvalue_threshold = 0.05, use_padj = TRUE) {
  col <- if (isTRUE(use_padj)) "adj.P.Val" else "P.Value"
  if (!col %in% colnames(df)) stop("Requested ", col, " is missing in differential input")
  df %>% dplyr::filter(is.finite(.data[[col]]) & .data[[col]] <= pvalue_threshold)
}

normalize_gene_symbols <- function(ids) {
  ids <- toupper(sub("\\.[0-9]+$", "", trimws(as.character(ids))))
  unique(ids[!is.na(ids) & nzchar(ids)])
}

read_universe <- function(path) {
  if (is.null(path) || !length(path)) return(NULL)
  if (length(path) > 1L || (length(path) == 1L && !file.exists(path))) {
    return(unique(trimws(as.character(path))))
  }
  if (!file.exists(path)) stop("Universe file does not exist: ", path)
  unique(trimws(readLines(path, warn = FALSE)))
}

if (is.null(args$mrna_universe_file)) {
  args$mrna_universe_file <- evomics_read_analysis_universe_from_rda(args$mrna_dataset)
}
if (is.null(args$protein_universe_file)) {
  args$protein_universe_file <- evomics_read_analysis_universe_from_rda(args$protein_dataset)
}

map_node_color_by_count <- function(nodes_df,
                                    source_color = "#F5DEB3",
                                    na_color     = "#bdbdbd",
                                    palette_cols = c(args$low_color, args$mid_color, args$high_color)) {
  count_values <- nodes_df$Count[!is.na(nodes_df$Count)]
  if (length(count_values) == 0) {
    return(nodes_df %>% dplyr::mutate(color = na_color))
  }
  
  count_range <- range(count_values, na.rm = TRUE)
  if (diff(count_range) == 0) {
    count_range <- c(count_range[1] - 1e-6, count_range[2] + 1e-6)
  }
  
  col_fun <- scales::col_numeric(
    palette = palette_cols,
    domain  = count_range
  )
  
  nodes_df %>%
    dplyr::mutate(
      color = dplyr::case_when(
        type == "Source"  ~ source_color,
        is.na(Count)      ~ na_color,
        TRUE              ~ col_fun(Count)
      )
    )
}

# ======================== 读取RDA函数 ========================

read_deg_from_rda <- function(rda_file) {
  env <- new.env()
  load(rda_file, envir = env)
  if (!exists("deg", envir = env)) {
    stop(paste("deg object not found in", rda_file))
  }
  df <- env$deg
  df$ID <- rownames(df)
  return(df)
}

# ======================== 主逻辑 ========================

cat("Building GO hierarchical network: mRNA + Protein\n")
cat("  Ontology:", args$ont, "\n")
cat("  Feature P cutoff:", feature_pvalue_threshold,
    " Term FDR cutoff:", term_fdr_threshold,
    " Feature P type:", ifelse(isTRUE(args$use_padj_for_features), "adjusted", "raw"), "\n")

# 读取数据 (从rda文件)
mrna_df <- read_deg_from_rda(args$mrna_dataset)
protein_df <- read_deg_from_rda(args$protein_dataset)

mrna_sig <- get_significant(mrna_df, feature_pvalue_threshold, args$use_padj_for_features)
protein_sig <- get_significant(protein_df, feature_pvalue_threshold, args$use_padj_for_features)

cat("  Significant mRNA:", nrow(mrna_sig), "\n")
cat("  Significant Protein:", nrow(protein_sig), "\n")

if (nrow(mrna_sig) == 0 || nrow(protein_sig) == 0) {
  go_empty("no_significant_features", "GO has no features passing the selected feature cutoff.",
           counts = list(mrna_features = nrow(mrna_sig), protein_features = nrow(protein_sig)))
}

# 基因 ID 映射
rna_genes <- mrna_sig$ID
rna_entrez <- bitr(rna_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
if (nrow(rna_entrez) == 0) {
  go_empty("no_mapped_features", "mRNA significant features could not be mapped to ENTREZID.",
           counts = list(mrna_features = nrow(mrna_sig), protein_features = nrow(protein_sig)))
}
rna_unmapped_count <- length(setdiff(
  normalize_gene_symbols(rna_genes), normalize_gene_symbols(rna_entrez$SYMBOL)
))

protein_genes <- protein_sig$ID
protein_entrez <- bitr(protein_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
if (nrow(protein_entrez) == 0) {
  go_empty("no_mapped_features", "Protein significant features could not be mapped to ENTREZID.",
           counts = list(mrna_features = nrow(mrna_sig), protein_features = nrow(protein_sig)))
}
protein_unmapped_count <- length(setdiff(
  normalize_gene_symbols(protein_genes), normalize_gene_symbols(protein_entrez$SYMBOL)
))
mrna_universe <- read_universe(args$mrna_universe_file)
protein_universe <- read_universe(args$protein_universe_file)
mrna_universe_entrez <- if (is.null(mrna_universe)) NULL else
  unique(bitr(mrna_universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)
protein_universe_entrez <- if (is.null(protein_universe)) NULL else
  unique(bitr(protein_universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)
if (is.null(mrna_universe_entrez) || is.null(protein_universe_entrez)) {
  warning("No assay-specific universe supplied for one or both omics; results are exploratory.")
}

# GO 富集 (ont = ALL)
cat("  Running mRNA GO enrichment...\n")
rna_go <- enrichGO(
  gene          = rna_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = args$ont,
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  pAdjustMethod = args$p_adjust_method,
  readable      = TRUE
  , universe    = mrna_universe_entrez
)
rna_go_df <- as.data.frame(rna_go)
if (nrow(rna_go_df) == 0) go_empty("no_enriched_terms", "mRNA GO enrichment returned no terms.")

cat("  Running Protein GO enrichment...\n")
protein_go <- enrichGO(
  gene          = protein_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = args$ont,
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  pAdjustMethod = args$p_adjust_method,
  readable      = TRUE
  , universe    = protein_universe_entrez
)
protein_go_df <- as.data.frame(protein_go)
if (nrow(protein_go_df) == 0) go_empty("no_enriched_terms", "Protein GO enrichment returned no terms.")

cat("  mRNA GO terms:", nrow(rna_go_df), "\n")
cat("  Protein GO terms:", nrow(protein_go_df), "\n")

# Classification follows assay-specific significance, independently of whether
# a term was returned in both complete ORA result tables. Joint Fisher/BH
# remains descriptive and is corrected over all common-tested terms.
classify_go_terms <- function(go_table, fdr_threshold, adjust_method = "BH") {
  go_table <- go_table %>% dplyr::mutate(
    tested_RNA = !is.na(Count_RNA),
    tested_Protein = !is.na(Count_Protein),
    tested_in_both = tested_RNA & tested_Protein,
    significant_RNA = is.finite(p.adjust_RNA) & p.adjust_RNA <= fdr_threshold,
    significant_Protein = is.finite(p.adjust_Protein) & p.adjust_Protein <= fdr_threshold,
    significant_in_both = significant_RNA & significant_Protein,
    SourceType = dplyr::case_when(
      significant_in_both ~ "Common",
      significant_RNA ~ "RNA",
      significant_Protein ~ "Protein",
      TRUE ~ "Not significant"
    ),
    Count = dplyr::case_when(
      SourceType == "Common" ~ Count_RNA + Count_Protein,
      SourceType == "RNA" ~ Count_RNA,
      SourceType == "Protein" ~ Count_Protein,
      TRUE ~ NA_real_
    ),
    joint_p = ifelse(tested_in_both,
      stats::pchisq(-2 * (log(pmax(pvalue_RNA, .Machine$double.xmin)) +
                         log(pmax(pvalue_Protein, .Machine$double.xmin))),
                   df = 4, lower.tail = FALSE), NA_real_)
  )
  go_table$joint_FDR <- NA_real_
  both <- go_table$tested_in_both
  go_table$joint_FDR[both] <- p.adjust(go_table$joint_p[both], method = adjust_method)
  go_table %>% dplyr::mutate(
    display_FDR = dplyr::case_when(
      SourceType == "RNA" ~ p.adjust_RNA,
      SourceType == "Protein" ~ p.adjust_Protein,
      SourceType == "Common" ~ joint_FDR,
      TRUE ~ NA_real_
    ),
    p_RNA = pvalue_RNA, FDR_RNA = p.adjust_RNA,
    p_Protein = pvalue_Protein, FDR_Protein = p.adjust_Protein
  )
}

go_combined <- dplyr::full_join(
  rna_go_df     %>% dplyr::select(ID, Description, Count, pvalue, p.adjust, ONTOLOGY),
  protein_go_df %>% dplyr::select(ID, Description, Count, pvalue, p.adjust, ONTOLOGY),
  by = c("ID", "Description", "ONTOLOGY"),
  suffix = c("_RNA", "_Protein")
) %>% classify_go_terms(term_fdr_threshold, args$p_adjust_method)
display_keep <- if (isTRUE(args$common_only)) {
  go_combined$significant_in_both
} else {
  go_combined$significant_RNA | go_combined$significant_Protein
}
if (!is.finite(args$static_top_n) || args$static_top_n < 1L)
  stop("static_top_n must be a positive integer")
# A reproducible quantity cap is separate from statistical eligibility.
static_ids <- go_combined[display_keep, , drop = FALSE] %>%
  dplyr::arrange(display_FDR, ID) %>%
  dplyr::slice_head(n = args$static_top_n) %>% dplyr::pull(ID)
go_export <- go_combined
go_export$html_displayed <- display_keep
go_export$static_displayed <- go_export$ID %in% static_ids
go_export$static_selection_rule <- paste0("Among HTML-eligible terms: ascending display_FDR, then GO ID; top ", args$static_top_n)
go_export$feature_pvalue_threshold <- feature_pvalue_threshold
go_export$term_fdr_threshold <- term_fdr_threshold
go_export$use_padj_for_features <- isTRUE(args$use_padj_for_features)
go_export$feature_pvalue_type <- ifelse(isTRUE(args$use_padj_for_features), "adjusted", "raw")
go_export$threshold_contract <- "feature_pvalue_threshold=feature filter; term_fdr_threshold=GO term FDR"
# Always preserve the complete tested GO table, including when no term passes
# the display gate.  This is the artifact returned with an empty status.
write.table(go_export, file = args$output_csv, sep = "\t", row.names = FALSE, quote = FALSE)

go_combined <- go_combined[display_keep, , drop = FALSE]
go_combined <- go_combined[!is.na(go_combined$SourceType), ]
if (nrow(go_combined) == 0) {
  go_empty("no_display_terms", "GO enrichment completed, but no term passed the selected display FDR cutoff.",
           counts = list(tested_terms = nrow(go_export), displayed_terms = 0L),
           artifacts = list(csv = args$output_csv))
}
go_combined$unmapped_features_count_RNA <- rna_unmapped_count
go_combined$unmapped_features_count_Protein <- protein_unmapped_count
go_combined$unmapped_features_count <- rna_unmapped_count + protein_unmapped_count

cat("  Combined GO terms:", nrow(go_combined), "\n")

cat("  Combined GO results saved:", args$output_csv, "\n")

# ======================== 构建分层网络 ========================

# 4) 分层网络节点：TopCore / SmallCore / GO Term

# TopCore coordinates reproduce the Figure 6D reading order.
top_core <- data.frame(
  name = c("Protein_core", "Common_core", "RNA_core"),
  label = c("Protein", "Common", "RNA"),
  title_text = c("Protein core", "Common GO core", "RNA core"),
  type = "TopCore", Count = NA_real_,
  x = c(-360, 0, 360), y = 0,
  stringsAsFactors = FALSE
)

# Create only ontology cores that actually contain significant terms.  This
# prevents the empty/orphan BP/CC/MF nodes produced by the former 3 x 3 grid.
small_core <- go_combined %>%
  dplyr::count(SourceType, ONTOLOGY, name = "term_count") %>%
  dplyr::mutate(
    Top = dplyr::case_when(
      SourceType == "RNA" ~ "RNA_core",
      SourceType == "Protein" ~ "Protein_core",
      TRUE ~ "Common_core"
    ),
    name = paste(Top, ONTOLOGY, sep = "_"),
    label = paste0(SourceType, " core ", ONTOLOGY, "\n(n = ", term_count, ")"),
    title_text = paste0(SourceType, " ", ONTOLOGY, " core; ", term_count, " GO terms"),
    type = "SmallCore", Count = NA_real_,
    x = dplyr::case_when(
      name == "Protein_core_BP" ~ -360,
      name == "Protein_core_CC" ~ -650,
      name == "Protein_core_MF" ~ -360,
      name == "RNA_core_CC" ~ 360,
      name == "RNA_core_MF" ~ 650,
      name == "RNA_core_BP" ~ 360,
      TRUE ~ 0
    ),
    y = dplyr::case_when(
      name == "Protein_core_BP" ~ -300,
      name == "Protein_core_MF" ~ 300,
      name == "RNA_core_CC" ~ -300,
      name == "RNA_core_BP" ~ 300,
      name == "Common_core_MF" ~ -190,
      name == "Common_core_CC" ~ 190,
      name == "Common_core_BP" ~ -260,
      TRUE ~ 0
    )
  ) %>%
  dplyr::select(name, label, title_text, type, Count, x, y, ONTOLOGY, Top, term_count)

# GO term nodes use stable GO IDs internally and fixed coordinates between
# their assay core and ontology core.  Only shared terms are labelled in the
# canvas; every term remains available through hover and the CSV table.
go_nodes <- go_combined %>%
  dplyr::mutate(
    name = paste0("term_", ID),
    label = ifelse(SourceType == "Common", Description, ""),
    title_text = paste0(Description, "<br>GO ID: ", ID,
      "<br>Class: ", SourceType, "<br>Ontology: ", ONTOLOGY,
      "<br>Gene count: ", Count),
    type = "GO Term",
    TopCore = dplyr::case_when(
      SourceType == "RNA" ~ "RNA_core",
      SourceType == "Protein" ~ "Protein_core",
      TRUE ~ "Common_core"
    ),
    SmallCore = paste(TopCore, ONTOLOGY, sep = "_")
  ) %>%
  dplyr::group_by(SmallCore) %>%
  dplyr::arrange(display_FDR, ID, .by_group = TRUE) %>%
  dplyr::mutate(cluster_index = dplyr::row_number(), cluster_n = dplyr::n()) %>%
  dplyr::ungroup()

position_cluster <- function(df) {
  small_idx <- match(df$SmallCore[1], small_core$name)
  top_idx <- match(df$TopCore[1], top_core$name)
  sx <- small_core$x[small_idx]; sy <- small_core$y[small_idx]
  tx <- top_core$x[top_idx]; ty <- top_core$y[top_idx]
  dx <- sx - tx; dy <- sy - ty
  distance <- sqrt(dx^2 + dy^2)
  if (!is.finite(distance) || distance == 0) distance <- 1
  ux <- dx / distance; uy <- dy / distance
  px <- -uy; py <- ux
  columns <- max(1L, ceiling(sqrt(nrow(df))))
  rows <- ceiling(nrow(df) / columns)
  col_index <- (seq_len(nrow(df)) - 1L) %% columns + 1L
  row_index <- (seq_len(nrow(df)) - 1L) %/% columns + 1L
  perpendicular <- if (columns == 1L) 0 else seq(-120, 120, length.out = columns)[col_index]
  longitudinal <- if (rows == 1L) 0 else seq(-75, 75, length.out = rows)[row_index]
  df$x <- tx + 0.55 * dx + px * perpendicular + ux * longitudinal
  df$y <- ty + 0.55 * dy + py * perpendicular + uy * longitudinal
  df
}
go_nodes <- go_nodes %>%
  dplyr::group_split(SmallCore) %>%
  lapply(position_cluster) %>%
  dplyr::bind_rows()

# 5) 汇总所有节点，设置颜色和大小
nodes_vis_go <- bind_rows(
  top_core %>% dplyr::select(name, label, title_text, type, Count, x, y),
  small_core %>% dplyr::select(name, label, title_text, type, Count, x, y),
  go_nodes %>% dplyr::select(name, label, title_text, type, Count, x, y)
) %>%
  dplyr::mutate(
    id    = name,
    title = title_text,
    fixed = TRUE,
    physics = FALSE,
    size = dplyr::case_when(
      type == "TopCore"   ~ 25,
      type == "SmallCore" ~ 12,
      type == "GO Term"   ~ ifelse(is.na(Count), 10, log2(Count + 1) * 5),
      TRUE                ~ 10
    )
  )

# 1) 先给 TopCore / SmallCore 设好颜色
nodes_vis_go <- nodes_vis_go %>%
  dplyr::mutate(
    color = dplyr::case_when(
      # TopCore 颜色
      type == "TopCore" & name == "RNA_core"     ~ args$rna_core_color,
      type == "TopCore" & name == "Protein_core" ~ args$protein_core_color,
      type == "TopCore" & name == "Common_core"  ~ args$common_core_color,
      # Figure 6D uses the same orange anchor for all ontology cores.
      type == "SmallCore" ~ "#F5A623",
      TRUE                                       ~ NA_character_
    )
  )

# 2) 对 GO Term 节点，用 Count 做蓝色渐变映射
go_term_nodes <- nodes_vis_go %>%
  dplyr::filter(type == "GO Term") %>%
  dplyr::select(name, type, Count)

if (nrow(go_term_nodes) > 0) {
  go_term_colored <- map_node_color_by_count(
    nodes_df     = go_term_nodes,
    source_color = "#F5DEB3",
    na_color     = "#bdbdbd"
  ) %>%
    dplyr::select(name, color_byCount = color)
  
  # 把 GO Term 的颜色更新回去
  nodes_vis_go <- nodes_vis_go %>%
    dplyr::left_join(go_term_colored, by = "name") %>%
    dplyr::mutate(
      color = dplyr::if_else(
        type == "GO Term" & !is.na(color_byCount),
        color_byCount,
        color
      )
    ) %>%
    dplyr::select(-color_byCount)
} else {
  # 没有 GO Term 的情况，确保 color 非 NA
  nodes_vis_go <- nodes_vis_go %>%
    dplyr::mutate(
      color = dplyr::coalesce(color, "grey70")
    )
}

# 添加边框
nodes_vis_go <- nodes_vis_go %>%
  dplyr::mutate(
    borderWidth = if (args$is_show_border) args$border_size else 0,
    color.border = args$border_color,
    color.background = color
  )

# 6) 构建边
edges_go <- bind_rows(
  # TopCore → GO Term (直接连接)
  go_combined %>%
    dplyr::mutate(
      name = paste0("term_", ID),
      from = dplyr::case_when(
        SourceType == "RNA"     ~ "RNA_core",
        SourceType == "Protein" ~ "Protein_core",
        SourceType == "Common"  ~ "Common_core"
      ),
      to = name
    ) %>%
    dplyr::select(from, to) %>%
    dplyr::mutate(color = "#FF8FA3", width = 1.4),
  # SmallCore → GO Term
  go_nodes %>%
    dplyr::mutate(from = SmallCore, to = name) %>%
    dplyr::select(from, to) %>%
    dplyr::mutate(color = "#8C8C8C", width = 1.2)
) %>%
  distinct()

cat("  Total nodes:", nrow(nodes_vis_go), "\n")
cat("  Total edges:", nrow(edges_go), "\n")

# 7) 画 visNetwork
go_network <- visNetwork(nodes_vis_go, edges_go, main = args$main_title,
                         width = paste0(args$width, "px"),
                         height = paste0(args$height, "px")) %>%
  visNodes(
    shape = "dot",
    font  = list(face = "arial", size = args$axis_text_size)
  ) %>%
  visEdges(smooth = FALSE) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visPhysics(enabled = FALSE)

# 添加图例
if (args$is_show_legend) {
  legend_nodes <- data.frame(
    label = c("RNA core", "Protein core", "Common core", "Ontology core", "GO Term"),
    shape = c(rep("dot", 5)),
    color = c(args$rna_core_color, args$protein_core_color, args$common_core_color,
              "#F5A623", "#3182bd"),
    font.size = args$legend_text_size,
    stringsAsFactors = FALSE
  )
  go_network <- go_network %>%
    visLegend(
      addNodes = legend_nodes,
      useGroups = FALSE,
      position = args$legend_position,
      main = list(text = "Legend", style = paste0("font-size:", args$legend_title_size, "px;"))
    )
}

# 保存网络
save_portable_widget(go_network, args$output_html)
cat("  Network output file:", args$output_html, "\n")

if (!is.null(args$output_png) || !is.null(args$output_pdf)) {
  suppressPackageStartupMessages({
    library(igraph)
    library(ggraph)
    library(ggplot2)
  })
  static_terms <- paste0("term_", static_ids)
  term_edges <- edges_go %>% filter(to %in% static_terms)
  used_static <- unique(c(static_terms, term_edges$from))
  ancestor_edges <- edges_go %>% filter(to %in% used_static, !to %in% static_terms)
  edges_static <- bind_rows(term_edges, ancestor_edges) %>% distinct()
  used_static <- unique(c(edges_static$from, edges_static$to))
  static_core_counts <- term_edges %>% dplyr::count(from, name = "shown_terms")
  vertices_static <- nodes_vis_go %>%
    filter(id %in% used_static) %>%
    dplyr::left_join(static_core_counts, by = c("id" = "from")) %>%
    dplyr::mutate(label = ifelse(type == "SmallCore",
      paste0(sub("\\n.*$", "", label), "\n(shown ", shown_terms, ")"), label)) %>%
    transmute(name = id, label = label, type = type, color = color,
              node_size = pmax(2.5, sqrt(pmax(size, 1)))) %>%
    distinct(name, .keep_all = TRUE)
  graph_static <- igraph::graph_from_data_frame(
    edges_static[, c("from", "to")], directed = TRUE,
    vertices = vertices_static
  )
  set.seed(123)
  p_static <- ggraph(graph_static, layout = "fr") +
    geom_edge_link(arrow = grid::arrow(length = grid::unit(2, "mm")),
                   end_cap = circle(2.2, "mm"), colour = "grey70",
                   alpha = 0.55) +
    geom_node_point(aes(size = node_size, fill = I(color), shape = type),
                    colour = "black", stroke = 0.3) +
    geom_node_text(aes(label = label), repel = TRUE, size = 2.7,
                   max.overlaps = Inf) +
    scale_shape_manual(values = c(TopCore = 23, SmallCore = 22,
                                  `GO Term` = 21),
                       labels = c(`GO Term` = "Term", SmallCore = "Ontology", TopCore = "Assay"),
                       name = "Node") +
    scale_size_identity() +
    labs(title = args$main_title,
         subtitle = paste0("Showing ", length(static_ids), " of ", nrow(go_combined),
                           " significant GO terms (ranked by display FDR)"),
         caption = "Shared: significant in both assays. Fisher/BH spans all common-tested terms; single-assay terms use that assay's FDR.") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.margin = margin(8, 12, 8, 12),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA))
  if (!is.null(args$output_png)) {
    ggsave(args$output_png, p_static, width = 11, height = 8, dpi = 300,
           bg = "white")
  }
  if (!is.null(args$output_pdf)) {
    ggsave(args$output_pdf, p_static, width = 11, height = 8,
           device = cairo_pdf, bg = "white")
  }
}

# 输出统计信息
result <- list(
  go_term_count   = nrow(go_combined),
  tested_term_count = nrow(go_export),
  static_term_count = length(static_ids),
  significant_rna_count = sum(go_export$significant_RNA),
  significant_protein_count = sum(go_export$significant_Protein),
  node_count      = nrow(nodes_vis_go),
  edge_count      = nrow(edges_go),
  common_count    = sum(go_combined$SourceType == "Common"),
  rna_only        = sum(go_combined$SourceType == "RNA"),
  protein_only    = sum(go_combined$SourceType == "Protein"),
  bp_count        = sum(go_combined$ONTOLOGY == "BP"),
  cc_count        = sum(go_combined$ONTOLOGY == "CC"),
  mf_count        = sum(go_combined$ONTOLOGY == "MF")
)

write_go_status("success", "completed", "GO network generated successfully",
                counts = result,
                artifacts = list(html = args$output_html, csv = args$output_csv))
cat(toJSON(result, auto_unbox = TRUE), "\n")

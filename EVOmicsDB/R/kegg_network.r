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
source(file.path(script_dir, "mirna_resolution_utils.R"))
source(file.path(script_dir, "kegg_joint_snapshot_utils.R"))

# ============================================================
# kegg_network.r
# KEGG联合网络
# 支持组合:
#   mRNA+Protein          — GSEA + 方向性 ORA 双重验证
#   mRNA+Metabolite       — gene ORA + compound ORA (Fisher joint)
#   Protein+Metabolite    — gene ORA + compound ORA (Fisher joint)
#   miRNA+Metabolite      — miRNA-target gene ORA + compound ORA
#   lncRNA+Metabolite     — lnc-related gene ORA + compound ORA
# ============================================================

library(clusterProfiler)
library(fgsea)
library(dplyr)
library(tidyr)
library(scales)
library(visNetwork)
library(KEGGREST)
library(org.Hs.eg.db)
library(stringr)
library(jsonlite)
library(argparse)

# --- Get script path for loading database ---
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grep(file_arg, args)])
  return(normalizePath(script_path))
}

parser <- ArgumentParser(description = "Build KEGG joint network for multi-omics data")
parser$add_argument("--dataset1", type = "character", required = TRUE, help = "Path to dataset1 file")
parser$add_argument("--dataset2", type = "character", required = TRUE, help = "Path to dataset2 file")
parser$add_argument("--type1", type = "character", required = TRUE, help = "Type of dataset1: mRNA/miRNA/Protein/Metabolite/lncRNA")
parser$add_argument("--type2", type = "character", required = TRUE, help = "Type of dataset2: mRNA/miRNA/Protein/Metabolite/lncRNA")
parser$add_argument("--pvalue_threshold", type = "double", default = 0.05,
                    help = "Deprecated legacy cutoff; used when semantic thresholds are omitted")
parser$add_argument("--feature_pvalue_threshold", type = "double", default = NULL,
                    help = "Feature-level raw/adjusted P cutoff")
parser$add_argument("--pathway_fdr_threshold", type = "double", default = NULL,
                    help = "Pathway-level FDR cutoff")
parser$add_argument("--p_adjust_method", type = "character", default = "BH", help = "P-value adjustment method")
parser$add_argument("--use_padj_for_features", type = "logical", default = TRUE)
parser$add_argument("--dataset1_universe_file", type = "character", default = NULL,
                    help = "One detected feature ID per line")
parser$add_argument("--dataset2_universe_file", type = "character", default = NULL,
                    help = "One detected feature ID per line")
parser$add_argument("--metabolite_mapping_file", type = "character", default = NULL,
                    help = "Curated TSV/CSV with Metabolite and KEGG_ID columns")
parser$add_argument("--organism", type = "character", default = "hsa", help = "Organism code")
parser$add_argument("--kegg_mapping_rda", type = "character", default = NULL,
                    help = "Optional local MSigDB KEGG mapping RDA for offline enrichment")
parser$add_argument("--kegg_joint_snapshot_dir", type = "character", default = NULL,
                    help = "Frozen human gene/compound KEGG annotation snapshot; required for metabolite pairs")
parser$add_argument("--lncrna_target_db", type = "character", default = NULL,
                    help = "Frozen lncRNA target SQLite database")
parser$add_argument("--max_targets_per_mirna", type = "integer", default = 15, help = "Legacy display parameter; metabolite-pair ORA uses all evidence-qualified targets")
parser$add_argument("--mirna_evidence", type = "character", default = "predicted",
                    choices = c("validated", "predicted", "all"))
parser$add_argument("--mirna_target_db", type = "character", default = NULL)
parser$add_argument("--mirna_resolution", type = "character", default = NULL,
                    choices = c("mirbase_mature", "mirbase_precursor"))
parser$add_argument("--output_html", type = "character", required = TRUE, help = "Output HTML file path")
parser$add_argument("--output_csv", type = "character", required = TRUE, help = "Output CSV file path")
parser$add_argument("--output_png", type = "character", default = NULL)
parser$add_argument("--output_pdf", type = "character", default = NULL)
parser$add_argument("--status_json", type = "character", default = NULL,
                    help = "Structured success/empty/error status output")
# Titles
parser$add_argument("--main_title", type = "character", default = "", help = "Optional main title; resolved by pair when empty")
parser$add_argument("--main_title_size", type = "integer", default = 14, help = "Main title font size")
parser$add_argument("--axis_text_size", type = "integer", default = 10, help = "Node label font size")
# Dot
parser$add_argument("--fill_colors", type = "character", default = "#F5DEB3,#9ecae1,#3182bd", help = "Fill colors (Source,Low,High)")
parser$add_argument("--border_color", type = "character", default = "black", help = "Border color")
parser$add_argument("--border_size", type = "double", default = 1.0, help = "Border size")
# Legend
parser$add_argument("--is_show_legend", type = "character", default = "TRUE", help = "Show legend (TRUE/FALSE)")
parser$add_argument("--legend_title_size", type = "integer", default = 10, help = "Legend title font size")
parser$add_argument("--legend_text_size", type = "integer", default = 9, help = "Legend label font size")
parser$add_argument("--legend_position", type = "character", default = "right", help = "Legend position")
# Export settings
parser$add_argument("--is_show_border", type = "character", default = "TRUE", help = "Show border (TRUE/FALSE)")
parser$add_argument("--width", type = "integer", default = 1200, help = "Output width (px)")
parser$add_argument("--height", type = "integer", default = 800, help = "Output height (px)")
args <- parser$parse_args()

if (!is.null(args$status_json) && file.exists(args$status_json)) unlink(args$status_json)
write_kegg_status <- function(status, code, message_text, counts = list(), artifacts = list()) {
  evomics_status_write(args$status_json, "kegg_network", status, code, message_text,
                       counts, artifacts)
}
kegg_empty <- function(code, message_text, counts = list(), artifacts = list()) {
  evomics_status_empty(args$status_json, "kegg_network", code, message_text,
                       counts, artifacts)
}

resolve_threshold <- function(value, legacy) {
  if (is.null(value) || !length(value) || is.na(value)) as.numeric(legacy) else as.numeric(value)
}
feature_pvalue_threshold <- resolve_threshold(args$feature_pvalue_threshold, args$pvalue_threshold)
pathway_fdr_threshold <- resolve_threshold(args$pathway_fdr_threshold, args$pvalue_threshold)
if (!is.finite(feature_pvalue_threshold) || feature_pvalue_threshold <= 0 || feature_pvalue_threshold > 1) {
  stop("feature_pvalue_threshold must be > 0 and <= 1")
}
if (!is.finite(pathway_fdr_threshold) || pathway_fdr_threshold <= 0 || pathway_fdr_threshold > 1) {
  stop("pathway_fdr_threshold must be > 0 and <= 1")
}

production_combo <- kegg_joint_pair_key(args$type1,args$type2)
SUPPORTED_COMBOS <- c(
  "mRNA_Protein",
  "Metabolite_mRNA", "Metabolite_Protein",
  "Metabolite_miRNA", "Metabolite_lncRNA"
)
if (!production_combo %in% SUPPORTED_COMBOS) {
  stop(
    "Unsupported KEGG combination: ", args$type1, " + ", args$type2,
    ". Supported: ", paste(SUPPORTED_COMBOS, collapse = ", ")
  )
}

if (is.null(args$kegg_joint_snapshot_dir)) {
  args$kegg_joint_snapshot_dir <- evomics_default_db_path(script_dir, "kegg_joint_snapshot")
}
if (is.null(args$lncrna_target_db)) {
  args$lncrna_target_db <- evomics_default_db_path(script_dir, "external_targets.sqlite")
}
if (is.null(args$mirna_target_db)) {
  args$mirna_target_db <- evomics_default_db_path(script_dir, "mirna_targets.sqlite")
}
if (is.null(args$kegg_mapping_rda)) {
  args$kegg_mapping_rda <- evomics_default_db_path(script_dir, "kegg_mapping.rda")
}
if (is.null(args$dataset1_universe_file)) {
  args$dataset1_universe_file <- evomics_read_analysis_universe_from_rda(args$dataset1)
}
if (is.null(args$dataset2_universe_file)) {
  args$dataset2_universe_file <- evomics_read_analysis_universe_from_rda(args$dataset2)
}

offline_kegg_map <- NULL
if (!"Metabolite" %in% c(args$type1,args$type2) && !is.null(args$kegg_mapping_rda)) {
  map_env <- new.env()
  load(args$kegg_mapping_rda, envir = map_env)
  map_names <- ls(map_env)
  if (!length(map_names)) stop("No object found in kegg_mapping_rda")
  selected_name <- if ("kegg_mapping" %in% map_names) "kegg_mapping" else map_names[1]
  offline_kegg_map <- as.data.frame(map_env[[selected_name]])
  if (!all(c("gs_name", "gene_symbol") %in% names(offline_kegg_map))) {
    stop("Offline KEGG mapping must contain gs_name and gene_symbol")
  }
}

# Parse boolean arguments
args$is_show_legend <- toupper(args$is_show_legend) == "TRUE"
args$is_show_border <- toupper(args$is_show_border) == "TRUE"

# Parse fill colors
fill_colors_vec <- strsplit(args$fill_colors, ",")[[1]]
if (length(fill_colors_vec) < 3) {
  fill_colors_vec <- c("#F5DEB3", "#9ecae1", "#3182bd")
}

# ======================== 辅助函数 ========================

get_significant <- function(df, pvalue_threshold = 0.05, use_padj = TRUE, context = "input") {
  p_col <- if (isTRUE(use_padj)) "adj.P.Val" else "P.Value"
  if (!p_col %in% colnames(df)) stop("Requested ", p_col, " is missing in ", context)
  df %>% dplyr::filter(is.finite(.data[[p_col]]) & .data[[p_col]] <= pvalue_threshold)
}

normalize_gene_symbols <- function(ids) {
  ids <- toupper(sub("\\.[0-9]+$", "", trimws(as.character(ids))))
  ids[!is.na(ids) & nzchar(ids)]
}

# Save a self-contained htmlwidget without requiring a system Pandoc binary.
# ``visSave(..., selfcontained=TRUE)`` delegates to Pandoc and therefore turns
# an otherwise valid KEGG result into a technical 500 on minimal deployments.
# Inline the widget's generated JS/CSS assets directly instead; this changes
# only packaging of the artifact, not enrichment or network calculations.
save_portable_widget <- function(widget, output_html) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE) ||
      !requireNamespace("base64enc", quietly = TRUE)) {
    stop("The server is missing the htmlwidgets/base64enc export dependency.")
  }
  temp_root <- tempfile("kegg_widget_")
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

run_offline_gene_kegg <- function(symbols, pathway_fdr_threshold,
                                  p_adjust_method, universe_symbols = NULL) {
  term2gene <- unique(offline_kegg_map[, c("gs_name", "gene_symbol")])
  term2gene$gene_symbol <- toupper(sub(
    "\\.[0-9]+$", "", trimws(as.character(term2gene$gene_symbol))
  ))
  term2gene <- term2gene[!is.na(term2gene$gene_symbol) & nzchar(term2gene$gene_symbol), ]
  term2name <- unique(offline_kegg_map[, c("gs_name", "gs_description")])
  colnames(term2name) <- c("term", "name")
  res <- enricher(
    gene = unique(normalize_gene_symbols(symbols)),
    universe = if (is.null(universe_symbols)) NULL else unique(normalize_gene_symbols(universe_symbols)),
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = p_adjust_method
  )
  as.data.frame(res)
}

get_mirna_targets <- function(mirna_list, gene_list = NULL, limit = 10000,
                              max_targets_per_mirna = NULL) {
  # A paired measured miRNA dataset supplies an explicit resolution from the
  # router.  lncRNA+Metabolite uses miRNAs only as annotation intermediates;
  # those records are mature annotations and do not create expression nodes.
  resolution <- if (is.null(args$mirna_resolution) || !nzchar(args$mirna_resolution)) {
    "mirbase_mature"
  } else args$mirna_resolution
  mirna_ids <- mirna_list$ID
  target_ids <- if (!is.null(gene_list)) gene_list$ID else NULL

  if (!is.null(args$mirna_target_db) && nzchar(trimws(args$mirna_target_db))) {
    local_map <- evomics_read_target_map(
      args$mirna_target_db,
      evidence = args$mirna_evidence
    ) %>%
      dplyr::transmute(
        miRNA = .data$miRNA,
        target_symbol = .data$target_gene,
        score = suppressWarnings(as.numeric(.data$score))
      )
    wanted_keys <- resolve_mirna_feature_key(mirna_ids,resolution)
    local_map <- local_map[resolve_mirna_feature_key(local_map$miRNA,resolution) %in% wanted_keys,,drop=FALSE]
    local_map <- resolve_mirna_target_annotations(local_map %>% dplyr::transmute(miRNA, target_gene = target_symbol, score = score), mirna_ids, resolution)
    local_map <- local_map %>% dplyr::transmute(mature_mirna_id = miRNA, target_symbol = toupper(target_gene), score = score)
    if (!is.null(target_ids)) {
      wanted_targets <- toupper(sub("\\.[0-9]+$", "", trimws(as.character(target_ids))))
      local_map <- local_map %>%
        dplyr::filter(.data$target_symbol %in% wanted_targets)
    }
    # Local maps are complete and must never be globally truncated.  Apply
    # the cap only after resolution-aware matching, independently per measured
    # miRNA, so later/alphabetically distant miRNAs are not starved of targets.
    if (!is.null(max_targets_per_mirna)) {
      return(evomics_select_targets_per_mirna(
        local_map, max_targets_per_mirna,
        mirna_col = "mature_mirna_id", score_col = "score",
        target_col = "target_symbol"
      ))
    }
    return(local_map %>% dplyr::arrange(
      .data$mature_mirna_id,
      dplyr::desc(dplyr::coalesce(.data$score, -Inf)),
      .data$target_symbol
    ))
  }

  stop("A frozen local miRNA target database is required; live multiMiR fallback is disabled for reproducibility.")
}

## Combine evidence for pathways present in both omics.
## Fisher uses the two raw enrichment P values; the resulting joint P values
## are then adjusted together. Single-omics pathways retain their own FDR.
apply_fisher_to_common <- function(df, p1_col, p2_col, padj1_col, padj2_col,
                                   p_adjust_method = "BH", significance_threshold = 0.05) {
  required <- c("SourceType", p1_col, p2_col, padj1_col, padj2_col)
  if (!all(required %in% colnames(df))) {
    stop("Cannot calculate Fisher joint significance; missing columns: ",
         paste(setdiff(required, colnames(df)), collapse = ", "))
  }

  common <- !is.na(df$SourceType) & df$SourceType == "Common"
  df$joint_p <- NA_real_
  df$joint_FDR <- NA_real_

  if (any(common)) {
    p1 <- pmax(as.numeric(df[[p1_col]][common]), .Machine$double.xmin)
    p2 <- pmax(as.numeric(df[[p2_col]][common]), .Machine$double.xmin)
    valid <- is.finite(p1) & is.finite(p2) & p1 <= 1 & p2 <= 1
    joint <- rep(NA_real_, length(p1))
    joint[valid] <- stats::pchisq(
      -2 * (log(p1[valid]) + log(p2[valid])),
      df = 4,
      lower.tail = FALSE
    )
    df$joint_p[common] <- joint
    df$joint_FDR[common] <- stats::p.adjust(joint, method = p_adjust_method)
  }

  df$p.adjust <- ifelse(
    common,
    df$joint_FDR,
    ifelse(!is.na(df[[padj1_col]]), df[[padj1_col]], df[[padj2_col]])
  )
  # Keep the legacy generic p.adjust value for downstream compatibility, but
  # make its statistical role explicit.  Common rows contain exploratory
  # Fisher joint FDR; single-assay rows retain assay-specific pathway FDR.
  df$p.adjust_role <- ifelse(
    common,
    "exploratory_Fisher_joint_FDR",
    ifelse(!is.na(df[[padj1_col]]) | !is.na(df[[padj2_col]]),
           "assay_specific_pathway_FDR", NA_character_)
  )
  df$combined_p_exploratory <- df$joint_p
  df$combined_FDR_exploratory <- df$joint_FDR
  df$shared_significant <- common &
    is.finite(df[[padj1_col]]) & is.finite(df[[padj2_col]]) &
    df[[padj1_col]] <= significance_threshold & df[[padj2_col]] <= significance_threshold
  df
}

map_node_color_by_count <- function(nodes_df,
                                    source_color = "#F5DEB3",
                                    na_color     = "#bdbdbd",
                                    palette_cols = c("#deebf7", "#9ecae1", "#3182bd"),
                                    border_color = "black",
                                    border_size  = 1.0) {
  count_range <- range(nodes_df$Count, na.rm = TRUE)
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
      ),
      color.border = border_color,
      borderWidth  = border_size
    )
}

# ======================== 各组合KEGG网络构建 ========================

build_kegg_mrna_protein <- function(mrna_sig, protein_sig, feature_threshold, pathway_threshold, p_adjust_method,
                                    organism, mrna_universe = NULL, protein_universe = NULL) {
  cat("  [KEGG] mRNA + Protein enrichment network...\n")

  rna_genes  <- mrna_sig$ID
  if (!is.null(offline_kegg_map)) {
    if (!length(intersect(normalize_gene_symbols(rna_genes),
                          normalize_gene_symbols(offline_kegg_map$gene_symbol)))) {
      kegg_empty("no_mapped_features", "mRNA features could not be mapped to KEGG identifiers.")
    }
    rna_unmapped_count <- length(setdiff(
      normalize_gene_symbols(rna_genes),
      normalize_gene_symbols(offline_kegg_map$gene_symbol)
    ))
    rna_kegg_df <- run_offline_gene_kegg(
      rna_genes, pathway_threshold, p_adjust_method, mrna_universe
    )
  } else {
    rna_entrez <- bitr(rna_genes, fromType = "SYMBOL", toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)
    if (nrow(rna_entrez) == 0) kegg_empty("no_mapped_features", "mRNA features could not be mapped to KEGG identifiers.")
    rna_unmapped_count <- length(setdiff(
      normalize_gene_symbols(rna_genes), normalize_gene_symbols(rna_entrez$SYMBOL)
    ))
    rna_kegg <- enrichKEGG(
      gene = rna_entrez$ENTREZID, organism = organism,
      pvalueCutoff = pathway_threshold, pAdjustMethod = p_adjust_method,
      universe = mrna_universe
    )
    rna_kegg_df <- as.data.frame(rna_kegg)
  }
  if (nrow(rna_kegg_df) == 0) kegg_empty("no_enriched_pathways", "mRNA KEGG enrichment returned no pathways.")
  rna_kegg_df$Source <- "RNA"

  protein_genes  <- protein_sig$ID
  if (!is.null(offline_kegg_map)) {
    if (!length(intersect(normalize_gene_symbols(protein_genes),
                          normalize_gene_symbols(offline_kegg_map$gene_symbol)))) {
      kegg_empty("no_mapped_features", "Protein features could not be mapped to KEGG identifiers.")
    }
    protein_unmapped_count <- length(setdiff(
      normalize_gene_symbols(protein_genes),
      normalize_gene_symbols(offline_kegg_map$gene_symbol)
    ))
    protein_kegg_df <- run_offline_gene_kegg(
      protein_genes, pathway_threshold, p_adjust_method, protein_universe
    )
  } else {
    protein_entrez <- bitr(protein_genes, fromType = "SYMBOL",
                           toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    if (nrow(protein_entrez) == 0) kegg_empty("no_mapped_features", "Protein features could not be mapped to KEGG identifiers.")
    protein_unmapped_count <- length(setdiff(
      normalize_gene_symbols(protein_genes), normalize_gene_symbols(protein_entrez$SYMBOL)
    ))
    protein_kegg <- enrichKEGG(
      gene = protein_entrez$ENTREZID, organism = organism,
      pvalueCutoff = pathway_threshold, pAdjustMethod = p_adjust_method,
      universe = protein_universe
    )
    protein_kegg_df <- as.data.frame(protein_kegg)
  }
  if (nrow(protein_kegg_df) == 0) kegg_empty("no_enriched_pathways", "Protein KEGG enrichment returned no pathways.")
  protein_kegg_df$Source <- "Protein"
  
  combined <- full_join(
    rna_kegg_df     %>% dplyr::select(ID, Description, Count, pvalue, p.adjust, Source),
    protein_kegg_df %>% dplyr::select(ID, Description, Count, pvalue, p.adjust, Source),
    by = c("ID", "Description"),
    suffix = c("_RNA", "_Protein")
  ) %>%
    mutate(
      SourceType = case_when(
        !is.na(Count_RNA) & !is.na(Count_Protein) ~ "Common",
        !is.na(Count_RNA) &  is.na(Count_Protein) ~ "RNA",
        is.na(Count_RNA)  & !is.na(Count_Protein) ~ "Protein"
      ),
      Count = ifelse(SourceType == "Common", Count_RNA + Count_Protein,
                     ifelse(!is.na(Count_RNA), Count_RNA, Count_Protein))
    )
  combined <- apply_fisher_to_common(
    combined, "pvalue_RNA", "pvalue_Protein",
    "p.adjust_RNA", "p.adjust_Protein", p_adjust_method, pathway_threshold
  )
  combined$unmapped_features_count_RNA <- rna_unmapped_count
  combined$unmapped_features_count_Protein <- protein_unmapped_count
  combined$unmapped_features_count <- rna_unmapped_count + protein_unmapped_count
  
  return(combined)
}

# -------------------------------------------------------------------------
# Figure 6C-compatible mRNA + Protein integration
#
# The production Explore endpoint used to pre-filter each DEG table and run
# over-representation analysis.  Figure 6C instead ranks every modelled
# feature by its moderated t statistic and runs GSEA independently in each
# assay.  Shared pathways are defined only by the two assay-specific GSEA FDRs;
# directional ORA is retained as a secondary validation field.
# -------------------------------------------------------------------------

adjust_method_name <- function(method) {
  method <- tolower(trimws(as.character(method %||% "BH")))
  if (method == "bh") return("BH")
  if (method == "fdr") return("fdr")
  if (method == "bonferroni") return("bonferroni")
  stop("Unsupported p_adjust_method: ", method)
}

prepare_gsea_pathways <- function() {
  if (is.null(offline_kegg_map) || !nrow(offline_kegg_map)) {
    stop("A local KEGG mapping is required for reproducible GSEA.")
  }
  term2gene <- unique(data.frame(
    term = as.character(offline_kegg_map$gs_name),
    gene = normalize_gene_symbols(offline_kegg_map$gene_symbol),
    stringsAsFactors = FALSE
  ))
  term2gene <- term2gene[!is.na(term2gene$term) & nzchar(term2gene$term) &
                           !is.na(term2gene$gene) & nzchar(term2gene$gene), , drop = FALSE]
  term2name <- unique(data.frame(
    term = as.character(offline_kegg_map$gs_name),
    name = as.character(offline_kegg_map$gs_description),
    stringsAsFactors = FALSE
  ))
  term2name <- term2name[!is.na(term2name$term) & nzchar(term2name$term), , drop = FALSE]
  list(
    term2gene = term2gene,
    term2name = term2name,
    pathways = split(term2gene$gene, term2gene$term),
    description = setNames(term2name$name, term2name$term)
  )
}

run_ranked_gsea <- function(df, label, gene_sets, p_adjust_method) {
  ids <- normalize_gene_symbols(rownames(df))
  stat_col <- if ("statistic" %in% names(df)) "statistic" else
    if ("t" %in% names(df)) "t" else NULL
  if (is.null(stat_col)) stop(label, " differential table has no moderated statistic column")
  stats <- suppressWarnings(as.numeric(df[[stat_col]]))
  keep <- !is.na(ids) & nzchar(ids) & is.finite(stats) & !duplicated(ids)
  stats <- stats[keep]
  names(stats) <- ids[keep]
  stats <- sort(stats, decreasing = TRUE)
  if (length(stats) < 10L) {
    kegg_empty("insufficient_ranked_features",
               paste0(label, " has fewer than 10 finite ranked features."),
               counts = list(label = label, ranked_features = length(stats)))
  }
  set.seed(123)
  fg <- as.data.frame(fgsea::fgseaMultilevel(
    pathways = gene_sets$pathways, stats = stats,
    minSize = 10, maxSize = 500, eps = 0,
    nPermSimple = 100000
  ))
  if (!nrow(fg)) {
    kegg_empty("no_enriched_pathways",
               paste0(label, " GSEA returned no pathways after gene-set size filtering."),
               counts = list(label = label, ranked_features = length(stats)))
  }
  fg$pval <- suppressWarnings(as.numeric(fg$pval))
  fg$NES <- suppressWarnings(as.numeric(fg$NES))
  fg$padj <- stats::p.adjust(fg$pval, method = adjust_method_name(p_adjust_method))
  fg$Description <- unname(gene_sets$description[as.character(fg$pathway)])
  fg$Description[is.na(fg$Description)] <- fg$pathway[is.na(fg$Description)]
  fg$size <- as.integer(fg$size)
  fg$model <- label
  fg
}

empty_directional_ora <- function() {
  # Contract boundary: an empty ORA result is still a typed table.  Keeping
  # the same columns as the non-empty return makes bind_rows()/joins safe and
  # preserves the distinction between no pathways and a failed analysis.
  data.frame(
    pathway = character(0), direction = character(0),
    Description = character(0), Count = integer(0),
    pvalue = numeric(0), foreground_n = integer(0),
    stringsAsFactors = FALSE
  )
}

run_directional_ora <- function(df, label, gene_sets, feature_threshold, p_adjust_method) {
  ids <- normalize_gene_symbols(rownames(df))
  p_col <- if (isTRUE(args$use_padj_for_features)) "adj.P.Val" else "P.Value"
  if (!p_col %in% names(df)) stop(label, " differential table is missing requested significance column ", p_col)
  logfc <- if ("logFC" %in% names(df)) suppressWarnings(as.numeric(df$logFC)) else rep(NA_real_, nrow(df))
  pvals <- suppressWarnings(as.numeric(df[[p_col]]))
  universe <- unique(ids[!is.na(ids) & nzchar(ids)])
  observed <- lapply(c("up", "down"), function(direction) {
    direction_keep <- if (direction == "up") logfc > 0 else logfc < 0
    selected <- !is.na(pvals) & pvals <= feature_threshold & direction_keep
    genes <- unique(ids[selected & !is.na(ids) & nzchar(ids)])
    if (length(genes) < 3L) return(empty_directional_ora())
    result <- tryCatch(as.data.frame(clusterProfiler::enricher(
      gene = genes, universe = universe,
      TERM2GENE = gene_sets$term2gene,
      TERM2NAME = gene_sets$term2name,
      pvalueCutoff = 1, qvalueCutoff = 1,
      pAdjustMethod = adjust_method_name(p_adjust_method),
      minGSSize = 10, maxGSSize = 500
    )), error = function(e) data.frame())
    if (!nrow(result)) return(empty_directional_ora())
    data.frame(
      pathway = as.character(result$ID), direction = direction,
      Description = as.character(result$Description),
      Count = as.integer(result$Count), pvalue = as.numeric(result$pvalue),
      foreground_n = length(genes), stringsAsFactors = FALSE
    )
  }) %>% dplyr::bind_rows()
  required_observed <- c("pathway", "direction", "Description", "Count", "pvalue", "foreground_n")
  if (!all(required_observed %in% names(observed))) {
    observed <- empty_directional_ora()
  }
  all_rows <- expand.grid(
    pathway = names(gene_sets$pathways), direction = c("up", "down"),
    stringsAsFactors = FALSE
  )
  all_rows <- dplyr::left_join(all_rows, observed, by = c("pathway", "direction")) %>%
    dplyr::mutate(
      Description = ifelse(is.na(Description),
                           unname(gene_sets$description[pathway]), Description),
      Count = tidyr::replace_na(Count, 0L),
      pvalue = tidyr::replace_na(pvalue, 1),
      foreground_n = tidyr::replace_na(foreground_n, 0L),
      FDR_full_family = stats::p.adjust(pvalue, method = adjust_method_name(p_adjust_method)),
      omics = label
    )
  all_rows
}

build_gsea_mrna_protein <- function(mrna_df, protein_df, feature_threshold, pathway_threshold, p_adjust_method,
                                    mrna_universe = NULL, protein_universe = NULL,
                                    max_display = 12L) {
  cat("  [KEGG] Figure 6C-compatible ranked GSEA + directional integration...\n")
  gene_sets <- prepare_gsea_pathways()
  rna_gsea <- run_ranked_gsea(mrna_df, "RNA", gene_sets, p_adjust_method)
  protein_gsea <- run_ranked_gsea(protein_df, "Protein", gene_sets, p_adjust_method)
  rna_ora <- run_directional_ora(mrna_df, "RNA", gene_sets, feature_threshold, p_adjust_method)
  protein_ora <- run_directional_ora(protein_df, "Protein", gene_sets, feature_threshold, p_adjust_method)

  rna <- rna_gsea %>% dplyr::select(pathway, Description_RNA = Description,
    ES_RNA = ES, NES_RNA = NES, pval_RNA = pval, padj_RNA = padj, size_RNA = size)
  protein <- protein_gsea %>% dplyr::select(pathway, Description_Protein = Description,
    ES_Protein = ES, NES_Protein = NES, pval_Protein = pval, padj_Protein = padj,
    size_Protein = size)
  pair <- dplyr::full_join(rna, protein, by = "pathway") %>% dplyr::mutate(
    Description = dplyr::coalesce(Description_RNA, Description_Protein,
                                  unname(gene_sets$description[pathway])),
    significant_RNA = !is.na(padj_RNA) & padj_RNA <= pathway_threshold,
    significant_Protein = !is.na(padj_Protein) & padj_Protein <= pathway_threshold,
    shared_GSEA = significant_RNA & significant_Protein,
    direction_class = dplyr::case_when(
      shared_GSEA & sign(NES_RNA) == sign(NES_Protein) ~ "Concordant",
      shared_GSEA & sign(NES_RNA) != sign(NES_Protein) ~ "Discordant",
      TRUE ~ "Not shared"
    ),
    GSEA_class = dplyr::case_when(
      shared_GSEA & direction_class == "Concordant" ~ "Shared concordant",
      shared_GSEA & direction_class == "Discordant" ~ "Shared discordant",
      significant_RNA & !significant_Protein ~ "RNA-specific",
      significant_Protein & !significant_RNA ~ "Protein-specific",
      TRUE ~ "Not significant"
    )
  )

  common <- !is.na(pair$pval_RNA) & !is.na(pair$pval_Protein)
  pair$joint_p_all_tested <- NA_real_
  pair$joint_FDR_all_tested <- NA_real_
  if (any(common)) {
    jp <- stats::pchisq(-2 * (log(pmax(pair$pval_RNA[common], .Machine$double.xmin)) +
      log(pmax(pair$pval_Protein[common], .Machine$double.xmin))), df = 4, lower.tail = FALSE)
    pair$joint_p_all_tested[common] <- jp
    pair$joint_FDR_all_tested[common] <- stats::p.adjust(jp, method = adjust_method_name(p_adjust_method))
  }
  pair$expected_RNA_direction <- ifelse(pair$NES_RNA >= 0, "up", "down")
  pair$expected_Protein_direction <- ifelse(pair$NES_Protein >= 0, "up", "down")
  rna_key <- paste(rna_ora$pathway, rna_ora$direction)
  protein_key <- paste(protein_ora$pathway, protein_ora$direction)
  pair$ORA_FDR_RNA_expected <- rna_ora$FDR_full_family[
    match(paste(pair$pathway, pair$expected_RNA_direction), rna_key)]
  pair$ORA_FDR_Protein_expected <- protein_ora$FDR_full_family[
    match(paste(pair$pathway, pair$expected_Protein_direction), protein_key)]
  pair$ORA_validated <- pair$shared_GSEA &
    !is.na(pair$ORA_FDR_RNA_expected) & !is.na(pair$ORA_FDR_Protein_expected) &
    pair$ORA_FDR_RNA_expected <= pathway_threshold & pair$ORA_FDR_Protein_expected <= pathway_threshold

  select_display <- function(source, sig_col, nes_col, fdr_col, size_col) {
    z <- pair[pair[[sig_col]] %in% TRUE, , drop = FALSE]
    if (!nrow(z)) return(data.frame())
    # Preserve one quantitative node per assay.  A shared pathway therefore
    # has an RNA node and a Protein node linked in the integrated view; it is
    # not collapsed to an averaged pseudo-assay node.
    z$display_source <- source
    z$display_class <- ifelse(z$shared_GSEA, "Shared GSEA", paste0(source, "-specific"))
    z$source_NES <- z[[nes_col]]
    z$source_FDR <- z[[fdr_col]]
    z$source_size <- z[[size_col]]
    z$direction_arrow <- ifelse(z$source_NES >= 0, "↑", "↓")
    z$display_label <- paste0(z$direction_arrow, " ", z$Description,
      ifelse(z$ORA_validated, " †", ""))
    z[order(!z$shared_GSEA, z$source_FDR), , drop = FALSE][seq_len(min(max_display, nrow(z))), , drop = FALSE]
  }
  display <- dplyr::bind_rows(
    select_display("RNA", "significant_RNA", "NES_RNA", "padj_RNA", "size_RNA"),
    select_display("Protein", "significant_Protein", "NES_Protein", "padj_Protein", "size_Protein")
  )
  if (nrow(display)) {
    display <- display %>% dplyr::arrange(pathway, display_source) %>%
      dplyr::distinct(pathway, display_source, .keep_all = TRUE)
  }
  pair$displayed <- pair$pathway %in% display$pathway
  pair$display_source <- NA_character_
  displayed_rna <- unique(display$pathway[display$display_source == "RNA"])
  displayed_protein <- unique(display$pathway[display$display_source == "Protein"])
  pair$display_source[pair$pathway %in% displayed_rna & pair$pathway %in% displayed_protein] <- "RNA;Protein"
  pair$display_source[pair$pathway %in% displayed_rna & !pair$pathway %in% displayed_protein] <- "RNA"
  pair$display_source[pair$pathway %in% displayed_protein & !pair$pathway %in% displayed_rna] <- "Protein"
  attr(pair, "display_rows") <- display
  attr(pair, "gsea_counts") <- c(
    tested = nrow(pair), shared = sum(pair$shared_GSEA, na.rm = TRUE),
    displayed = length(unique(display$pathway))
  )
  pair
}

build_gsea_visnetwork <- function(combined, output_html, output_csv,
                                   main_title = "KEGG GSEA Integration",
                                   main_title_size = 14, axis_text_size = 10,
                                   fill_colors = c("#F5DEB3", "#9ecae1", "#3182bd"),
                                   border_color = "black", border_size = 1,
                                   is_show_legend = TRUE, legend_title_size = 10,
                                   legend_text_size = 9, legend_position = "right",
                                   is_show_border = TRUE, width = 1200, height = 800,
                                   output_png = NULL, output_pdf = NULL) {
  display <- attr(combined, "display_rows")
  if (is.null(display) || !nrow(display)) {
    # Preserve the complete primary GSEA table for an expected empty result;
    # the caller writes the structured empty status and does not fabricate a
    # blank network.
    utils::write.csv(combined, file = output_csv, row.names = FALSE, na = "")
    return(list(empty = TRUE, empty_code = "no_display_pathways",
                empty_message = "GSEA completed, but no pathway passed the selected FDR cutoff.",
                tested_pathways = nrow(combined), displayed_pathways = 0L,
                pathway_count = 0L, node_count = 0L, edge_count = 0L))
  }
  source_color <- fill_colors[1]
  class_colors <- c(
    "RNA-specific" = "#5B9BD5", "Protein-specific" = "#1F4E79",
    "Shared GSEA" = "#148F77"
  )
  border_width <- if (is_show_border) border_size else 0
  nes_limit <- max(4, max(abs(display$source_NES), na.rm = TRUE))
  nes_palette <- grDevices::colorRampPalette(c("#3C78B4", "#F7F7F7", "#D55E4A"))(201)
  nes_index <- pmax(1L, pmin(201L, round(scales::rescale(
    display$source_NES, to = c(1, 201), from = c(-nes_limit, nes_limit)
  ))))
  display$node_color <- nes_palette[nes_index]
  display$node_size <- pmax(13, pmin(40,
    10 + 1.2 * -log10(pmax(display$source_FDR, 1e-25))))
  display$node_id <- paste0("pathway_", display$display_source, "_", seq_len(nrow(display)))

  # Fixed bilateral coordinates mirror the formal Figure 6C grammar:
  # Protein on the left, RNA on the right, and paired nodes for shared paths.
  display$x <- NA_real_
  display$y <- NA_real_
  place_specific <- function(source, x_value) {
    idx <- which(display$display_source == source & !display$shared_GSEA)
    if (!length(idx)) return(invisible(NULL))
    display$x[idx] <<- x_value
    display$y[idx] <<- seq(-250, 250, length.out = length(idx))
  }
  place_specific("Protein", -470)
  place_specific("RNA", 470)
  shared_paths <- sort(unique(display$pathway[display$shared_GSEA]))
  if (length(shared_paths)) {
    shared_y <- seq(-105, 105, length.out = length(shared_paths))
    for (i in seq_along(shared_paths)) {
      idx <- which(display$pathway == shared_paths[i] & display$shared_GSEA)
      display$x[idx] <- ifelse(display$display_source[idx] == "Protein", -85, 85)
      display$y[idx] <- shared_y[i]
    }
  }

  source_nodes <- data.frame(
    id = c("Protein", "RNA"), label = c("Protein", "RNA"), type = "Source",
    class = "Source", shape = "dot", x = c(-250, 250), y = c(0, 0),
    fixed = TRUE, physics = FALSE, size = c(25, 25),
    color.background = c(source_color, source_color),
    color.border = c("#C6922A", "#C6922A"), borderWidth = border_width,
    title = c("Protein ranked GSEA", "RNA ranked GSEA"), stringsAsFactors = FALSE
  )
  pathway_nodes <- data.frame(
    id = display$node_id,
    label = ifelse(display$shared_GSEA, "", display$display_label),
    type = "Pathway", class = display$display_class, shape = "dot",
    x = display$x, y = display$y, fixed = TRUE, physics = FALSE,
    size = display$node_size,
    color.background = display$node_color,
    color.border = if (is_show_border) unname(class_colors[display$display_class]) else "transparent",
    borderWidth = border_width,
    title = paste0(display$Description, "<br>Assay: ", display$display_source,
      "<br>Class: ", display$GSEA_class,
      "<br>NES: ", round(display$source_NES, 3),
      "<br>GSEA FDR: ", format.pval(display$source_FDR, digits = 3),
      "<br>Directional ORA validated: ", ifelse(display$ORA_validated, "yes", "no")),
    stringsAsFactors = FALSE
  )
  label_nodes <- data.frame()
  if (length(shared_paths)) {
    label_rows <- lapply(seq_along(shared_paths), function(i) {
      idx <- which(display$pathway == shared_paths[i] & display$shared_GSEA)
      arrows <- paste0(display$direction_arrow[idx][match(c("Protein", "RNA"), display$display_source[idx])], collapse = "")
      description <- display$Description[idx[1]]
      data.frame(
        id = paste0("shared_label_", i),
        label = paste0(arrows, " ", description,
          ifelse(any(display$ORA_validated[idx]), " †", "")),
        type = "Label", class = "Shared GSEA", shape = "text",
        x = 0, y = unique(display$y[idx])[1] - 38,
        fixed = TRUE, physics = FALSE, size = 1,
        color.background = "transparent", color.border = "transparent",
        borderWidth = 0, title = description, stringsAsFactors = FALSE
      )
    })
    label_nodes <- dplyr::bind_rows(label_rows)
  }
  nodes <- dplyr::bind_rows(source_nodes, pathway_nodes, label_nodes)

  hub_edges <- data.frame(
    from = display$display_source, to = display$node_id, arrows = "",
    color = "#D9A18F", width = 1.5, stringsAsFactors = FALSE
  )
  shared_edges <- data.frame()
  if (length(shared_paths)) {
    shared_edges <- dplyr::bind_rows(lapply(shared_paths, function(path) {
      idx <- which(display$pathway == path & display$shared_GSEA)
      protein_id <- display$node_id[idx][display$display_source[idx] == "Protein"]
      rna_id <- display$node_id[idx][display$display_source[idx] == "RNA"]
      if (!length(protein_id) || !length(rna_id)) return(NULL)
      data.frame(from = protein_id[1], to = rna_id[1], arrows = "",
        color = "#2B9B78", width = 4, stringsAsFactors = FALSE)
    }))
  }
  edges <- dplyr::bind_rows(hub_edges, shared_edges)
  network <- visNetwork(nodes, edges, main = list(text = main_title,
    style = paste0("font-size:", main_title_size, "px; font-weight:bold;")),
    width = paste0(width, "px"), height = paste0(height, "px")) %>%
    visNodes(shape = "dot", font = list(face = "arial", size = axis_text_size),
      borderWidth = border_width) %>%
    visEdges(smooth = FALSE) %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visPhysics(enabled = FALSE)
  if (is_show_legend) {
    legend_nodes <- data.frame(label = c("Protein-specific", "RNA-specific", "Shared GSEA"),
      color = unname(class_colors[c("Protein-specific", "RNA-specific", "Shared GSEA")]),
      shape = "dot", font.size = legend_text_size)
    network <- network %>% visLegend(addNodes = legend_nodes, useGroups = FALSE,
      position = legend_position, main = list(text = "GSEA class",
        style = paste0("font-size:", legend_title_size, "px;")))
  }
  save_portable_widget(network, output_html)
  utils::write.csv(combined, file = output_csv, row.names = FALSE, na = "")
  cat("  Output HTML:", output_html, "\n")
  counts <- attr(combined, "gsea_counts")
  list(pathway_count = length(unique(display$pathway)),
       node_count = nrow(nodes), edge_count = nrow(edges),
       tested_pathways = unname(counts[["tested"]]),
       shared_gsea = unname(counts[["shared"]]),
       displayed_pathways = length(unique(display$pathway)),
       displayed_assay_nodes = nrow(display))
}

build_kegg_visnetwork <- function(combined, output_html, output_csv,
                                   output_png = NULL,
                                   output_pdf = NULL,
                                   main_title = "KEGG Joint Network",
                                   main_title_size = 14,
                                   axis_text_size = 10,
                                   fill_colors = c("#F5DEB3", "#9ecae1", "#3182bd"),
                                   border_color = "black",
                                   border_size = 1.0,
                                   is_show_legend = TRUE,
                                   legend_title_size = 10,
                                   legend_text_size = 9,
                                   legend_position = "right",
                                   is_show_border = TRUE,
                                   width = 1200,
                                   height = 800,
                                   pathway_fdr_threshold = 0.05) {
  # Metabolite integrations retain all common-tested pathways for the
  # exploratory Fisher/BH denominator, then apply the assay-specific display
  # rule here.  The mRNA+Protein GSEA branch has no shared_significant field
  # and is intentionally left unchanged.
  export_combined <- combined
  if ("shared_significant" %in% names(combined)) {
    combined <- combined %>% dplyr::filter(
      (SourceType == "Common" & shared_significant) |
      (SourceType != "Common" & is.finite(p.adjust) & p.adjust <= pathway_fdr_threshold)
    )
  }
  if (nrow(combined) == 0) {
    utils::write.table(export_combined, file = output_csv, sep = "\t", quote = FALSE, row.names = FALSE)
    kegg_joint_empty_artifacts(output_html,output_png,output_pdf,main_title)
    kegg_empty("no_display_pathways", "Enrichment completed, but no pathway passed the selected display FDR cutoff.",
                counts = list(tested_pathways = nrow(export_combined), displayed_pathways = 0L),
                artifacts = list(csv = output_csv,html=output_html,png=output_png,pdf=output_pdf))
  }
  
  # Extract colors
  source_color <- fill_colors[1]
  palette_cols <- c("#deebf7", fill_colors[2], fill_colors[3])
  
  nodes <- data.frame(
    name = unique(c(combined$SourceType, combined$ID)),
    display_label = c(unique(combined$SourceType),combined$Description[match(unique(combined$ID),combined$ID)]),
    type = c(
      rep("Source",  length(unique(combined$SourceType))),
      rep("Pathway", length(unique(combined$ID)))
    ),
    Count = c(
      rep(NA, length(unique(combined$SourceType))),
      combined$Count[match(unique(combined$ID), combined$ID)]
    ),
    stringsAsFactors = FALSE
  )
  
  nodes_colored <- map_node_color_by_count(
    nodes,
    source_color = source_color,
    palette_cols = palette_cols,
    border_color = border_color,
    border_size  = border_size
  )
  
  nodes_vis <- nodes_colored %>%
    mutate(
      id    = name,
      label = display_label,
      title = paste0(display_label, "<br>KEGG ID: ", name,"<br>Count: ", Count),
      size  = ifelse(is.na(Count), 10, log2(Count + 1) * 5)
    )
  
  edges <- combined %>%
    dplyr::select(from = SourceType, to = ID)
  
  # Build edge color based on border settings
  edge_color <- if (is_show_border) border_color else "transparent"
  
  network <- visNetwork(nodes_vis, edges, 
                        main = list(text = main_title, 
                                    style = paste0("font-size:", main_title_size, "px; font-weight:bold;")),
                        width = paste0(width, "px"),
                        height = paste0(height, "px")) %>%
    visNodes(shape = "dot", 
             font = list(face = "arial", size = axis_text_size),
             borderWidth = if (is_show_border) border_size else 0) %>%
    visEdges(color = list(color = edge_color, highlight = "red")) %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visLayout(randomSeed = 123) %>%
    visPhysics(stabilization = TRUE)
  
  # Add legend if enabled
  if (is_show_legend) {
    legend_nodes <- data.frame(
      label = c("Source", "Low Count", "High Count"),
      color = c(source_color, fill_colors[2], fill_colors[3]),
      shape = "dot",
      font.size = legend_text_size
    )
    network <- network %>%
      visLegend(addNodes = legend_nodes, useGroups = FALSE, position = legend_position,
                main = list(text = "Legend", style = paste0("font-size:", legend_title_size, "px;")))
  }
  
  save_portable_widget(network, output_html)
  # CSV retains the complete common-tested family (including exploratory
  # combined FDR values); the network itself uses the display-filtered rows.
  write.table(export_combined, file = output_csv, sep = "\t", quote = FALSE, row.names = FALSE)

  if (!is.null(output_png) || !is.null(output_pdf)) {
    suppressPackageStartupMessages({
      library(igraph)
      library(ggraph)
      library(ggplot2)
    })
    vertices_static <- nodes_vis %>%
      transmute(name = id, label = label, type = type, color = color,
                node_size = pmax(3, sqrt(pmax(size, 1)))) %>%
      distinct(name, .keep_all = TRUE)
    graph_static <- igraph::graph_from_data_frame(
      edges[, c("from", "to")], directed = TRUE,
      vertices = vertices_static
    )
    set.seed(123)
    p_static <- ggraph(graph_static, layout = "fr") +
      geom_edge_link(arrow = grid::arrow(length = grid::unit(2, "mm")),
                     end_cap = circle(2.3, "mm"), colour = "grey65",
                     alpha = 0.6) +
      geom_node_point(aes(size = node_size, fill = I(color), shape = type),
                      colour = "black", stroke = 0.35) +
      geom_node_text(aes(label = label), repel = TRUE, size = 2.8,
                     max.overlaps = Inf, seed = 123,
                     box.padding = 0.4, point.padding = 0.4, max.time = 5) +
      scale_x_continuous(expand = expansion(mult = 0.15)) +
      scale_y_continuous(expand = expansion(mult = 0.15)) +
      coord_cartesian(clip = "off") +
      scale_shape_manual(values = c(Source = 23, Pathway = 21)) +
      scale_size_identity() +
      labs(title = main_title,
           subtitle = "Shared: both assay FDRs pass the cutoff.\nFisher-combined statistics are exploratory; direction is annotation only.") +
      theme_void(base_size = 11) +
      theme(plot.title = element_text(face = "bold"),
            plot.margin = margin(12, 16, 12, 16),
            legend.position = "bottom",
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA))
    if (!is.null(output_png)) {
      ggsave(output_png, p_static, width = 10, height = 7, dpi = 300,
             bg = "white")
    }
    if (!is.null(output_pdf)) {
      ggsave(output_pdf, p_static, width = 10, height = 7,
             device = cairo_pdf, bg = "white")
    }
  }
  
  cat("  Output HTML:", output_html, "\n")
  
  return(list(
    pathway_count = length(unique(combined$Description)),
    node_count = nrow(nodes_vis),
    edge_count = nrow(edges)
  ))
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

cat("Building KEGG network:", args$type1, "+", args$type2, "\n")
cat("  Feature P cutoff:", feature_pvalue_threshold,
    " Pathway FDR cutoff:", pathway_fdr_threshold,
    " Feature P type:", ifelse(isTRUE(args$use_padj_for_features), "adjusted", "raw"), "\n")

# Load lncRNA target database if lncRNA is involved
master_lnc_targets <- NULL
combo_check <- sort(c(args$type1, args$type2))
if ("lncRNA" %in% combo_check) {
  lnc_db_path <- args$lncrna_target_db
  legacy_csv <- evomics_default_db_path(script_dir, "master_lncRNA_target_database.csv")
  if (!is.null(lnc_db_path)) {
    master_lnc_targets <- evomics_read_lnc_target_map(lnc_db_path)
    cat("  lncRNA SQLite database loaded:", nrow(master_lnc_targets), "records\n")
  } else if (!is.null(legacy_csv)) {
    master_lnc_targets <- read.csv(legacy_csv, stringsAsFactors = FALSE)
    cat("  Legacy lncRNA CSV loaded:", nrow(master_lnc_targets), "records\n")
  } else {
    stop("lncRNA target database not found (tried external_targets.sqlite and master_lncRNA_target_database.csv)")
  }
}

df1 <- read_deg_from_rda(args$dataset1)
df2 <- read_deg_from_rda(args$dataset2)

read_gene_universe <- function(path, omics_type, offline_mode = FALSE) {
  if (is.null(path) || !length(path) || !omics_type %in% c("mRNA", "Protein")) return(NULL)
  ids <- if (length(path) > 1L || (length(path) == 1L && !file.exists(path))) {
    unique(trimws(as.character(path)))
  } else {
    if (!file.exists(path)) stop("Universe file does not exist: ", path)
    unique(trimws(readLines(path, warn = FALSE)))
  }
  ids <- normalize_gene_symbols(ids)
  if (isTRUE(offline_mode)) return(unique(ids))
  mapped <- bitr(ids, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  unique(mapped$ENTREZID)
}
if ("Metabolite" %in% c(args$type1,args$type2)) {
  universe1 <- args$dataset1_universe_file
  universe2 <- args$dataset2_universe_file
} else {
universe1 <- read_gene_universe(
  args$dataset1_universe_file, args$type1,
  offline_mode = !is.null(offline_kegg_map)
)
universe2 <- read_gene_universe(
  args$dataset2_universe_file, args$type2,
  offline_mode = !is.null(offline_kegg_map)
)

}

sig1 <- get_significant(df1, feature_pvalue_threshold, args$use_padj_for_features)
sig2 <- get_significant(df2, feature_pvalue_threshold, args$use_padj_for_features)

cat("  Significant", args$type1, ":", nrow(sig1), "\n")
cat("  Significant", args$type2, ":", nrow(sig2), "\n")

combo <- sort(c(args$type1, args$type2))
combo_str <- kegg_joint_pair_key(args$type1,args$type2)

combined <- NULL
use_gsea_integration <- combo_str == "mRNA_Protein"

if (combo_str == "mRNA_Protein") {
  if (args$type1 == "mRNA") {
    combined <- build_gsea_mrna_protein(df1, df2, feature_pvalue_threshold, pathway_fdr_threshold, args$p_adjust_method,
                                        universe1, universe2)
  } else {
    combined <- build_gsea_mrna_protein(df2, df1, feature_pvalue_threshold, pathway_fdr_threshold, args$p_adjust_method,
                                        universe2, universe1)
  }
 } else if ("Metabolite" %in% combo) {
  if (args$type1 == "Metabolite") {
    combined <- build_kegg_metabolite_pair(df2,df1,args$type2,args$dataset2_universe_file,args$dataset1_universe_file,master_lnc_targets)
  } else {
    combined <- build_kegg_metabolite_pair(df1,df2,args$type1,args$dataset1_universe_file,args$dataset2_universe_file,master_lnc_targets)
  }
} else {
  stop(paste("Unsupported combination:", args$type1, "+", args$type2))
}

# Keep the two threshold meanings auditable in every exported KEGG table.
combined$feature_pvalue_threshold <- feature_pvalue_threshold
combined$pathway_fdr_threshold <- pathway_fdr_threshold
combined$use_padj_for_features <- isTRUE(args$use_padj_for_features)
combined$feature_pvalue_type <- ifelse(isTRUE(args$use_padj_for_features), "adjusted", "raw")
combined$threshold_contract <- "feature_pvalue_threshold=feature filter; pathway_fdr_threshold=pathway FDR"

# Pair-aware defensive fallback for direct R invocation.  Python normally
# resolves this title before calling R, but an empty title must never label a
# metabolite ORA/Fisher network as ranked GSEA.
resolved_main_title <- if (nzchar(trimws(args$main_title))) {
  args$main_title
} else if (use_gsea_integration) {
  "KEGG GSEA Integration"
} else {
  "KEGG Pathway Network"
}

if (use_gsea_integration) {
  result <- build_gsea_visnetwork(
    combined, args$output_html, args$output_csv,
    main_title = resolved_main_title,
    main_title_size = args$main_title_size, axis_text_size = args$axis_text_size,
    fill_colors = fill_colors_vec, border_color = args$border_color,
    border_size = args$border_size, is_show_legend = args$is_show_legend,
    legend_title_size = args$legend_title_size, legend_text_size = args$legend_text_size,
    legend_position = args$legend_position, is_show_border = args$is_show_border,
    width = args$width, height = args$height,
    output_png = args$output_png, output_pdf = args$output_pdf
  )
} else {
  result <- build_kegg_visnetwork(
    combined, args$output_html, args$output_csv,
    output_png = args$output_png, output_pdf = args$output_pdf,
    main_title = resolved_main_title, main_title_size = args$main_title_size,
    axis_text_size = args$axis_text_size, fill_colors = fill_colors_vec,
    border_color = args$border_color, border_size = args$border_size,
    is_show_legend = args$is_show_legend, legend_title_size = args$legend_title_size,
    legend_text_size = args$legend_text_size, legend_position = args$legend_position,
    is_show_border = args$is_show_border, width = args$width, height = args$height,
    pathway_fdr_threshold = pathway_fdr_threshold
  )
}
if (isTRUE(result$empty)) {
  write_kegg_status("empty", result$empty_code %||% "no_display_pathways",
                    result$empty_message %||% "No pathway passed the selected display threshold.",
                    counts = result[setdiff(names(result), c("empty", "empty_code", "empty_message"))],
                    artifacts = list(csv = args$output_csv))
  quit(save = "no", status = 0L)
}
write_kegg_status("success", "completed", "KEGG network generated successfully",
                  counts = result, artifacts = list(html = args$output_html, csv = args$output_csv))
cat(toJSON(result, auto_unbox = TRUE), "\n")

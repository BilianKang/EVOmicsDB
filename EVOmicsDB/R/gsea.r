#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(patchwork)
  library(enrichplot)
  library(argparse)
  library(AnnotationDbi)
})

# =========================================================
# Standard GSEA for gene/protein-level human data (final web-pipeline version)
# - Uses all mapped features ranked by a continuous statistic
# - Does not pre-filter by p value or logFC
# - Uses MSigDB C2:CP:KEGG_LEGACY gene sets by default
# =========================================================

parser <- ArgumentParser(description = "Generate a standard GSEA enrichment plot using all ranked genes.")
parser$add_argument("--input_rda", type = "character", required = TRUE,
                    help = "Input RDA containing a 'deg' differential-result table.")
parser$add_argument("--output_png", type = "character", required = TRUE,
                    help = "Output PNG file for GSEA plot.")
parser$add_argument("--output_csv", type = "character", required = TRUE,
                    help = "Output CSV file for GSEA results.")
parser$add_argument("--output_pdf", type = "character", default = NULL,
                    help = "Optional output PDF file for GSEA plot.")
parser$add_argument("--output_result_rds", type = "character", default = NULL,
                    help = "Optional complete GSEA result for reproducible figure redraws.")
parser$add_argument("--output_summary_csv", type = "character", default = NULL,
                    help = "Optional output CSV file recording mapping and GSEA audit information.")

# Compatibility arguments retained for the web pipeline; not used for pre-filtering.
parser$add_argument("--is_use_padj", type = "logical", default = TRUE,
                    help = "Compatibility only. Not used for standard GSEA pre-filtering.")
parser$add_argument("--log2fc_threshold", type = "numeric", default = 1,
                    help = "Compatibility only. Not used for standard GSEA pre-filtering.")
parser$add_argument("--pvalue_threshold", type = "numeric", default = 0.05,
                    help = "Compatibility only. Not used for standard GSEA pre-filtering.")

# Optional data-type guard. Keep 'auto' for backward compatibility.
parser$add_argument("--omics_type", type = "character", default = "auto",
                    help = "auto / transcriptomics / proteomics / phosphoproteomics / metabolomics / mirna / lncrna. Gene-set GSEA is only appropriate for gene/protein-level inputs or target-projected miRNA inputs.")
parser$add_argument("--kegg_mapping_rda", type = "character", default = NULL,
                    help = "Optional local MSigDB KEGG mapping RDA. EVOMICS_KEGG_MAPPING_RDA is used when omitted.")
parser$add_argument("--id_type", type = "character", default = "auto",
                    help = "Identifier namespace: auto, SYMBOL, ENTREZID, ENSEMBL, ENSEMBLPROT or UNIPROT. Auto prioritizes canonical preprocessing metadata.")
parser$add_argument("--seed", type = "integer", default = 123,
                    help = "Random seed for reproducible GSEA permutations.")
parser$add_argument("--min_gs_size", type = "integer", default = 10,
                    help = "Minimum gene-set size.")
parser$add_argument("--max_gs_size", type = "integer", default = 500,
                    help = "Maximum gene-set size.")

# Plot settings
parser$add_argument("--top_n", type = "integer", default = 5,
                    help = "Number of top pathways to display. Default = 5.")
parser$add_argument("--top_strategy", type = "character", default = "overall",
                    help = "balanced = show both positive and negative NES when available; overall = top pathways by adjusted P value.")
parser$add_argument("--main_title", type = "character", default = "Enriched Pathways",
                    help = "Main title for the plot.")
parser$add_argument("--x_axis_title", type = "character", default = "Gene Rank",
                    help = "Title for the bottom X-axis.")
parser$add_argument("--y_axis_title", type = "character", default = "Enrichment Score",
                    help = "Title for the ES-curve Y-axis.")
parser$add_argument("--main_title_size", type = "integer", default = 10,
                    help = "Font size for the main title.")
parser$add_argument("--axis_title_size", type = "integer", default = 7,
                    help = "Font size for axis titles.")
parser$add_argument("--axis_text_size", type = "integer", default = 7,
                    help = "Font size for axis tick labels.")
parser$add_argument("--line_color", type = "character", nargs = "+",
                    default = c("#E64B35", "#4DBBD5", "#8C8C8C", "#91BD91", "#F39C12"),
                    help = "Line colors.")
parser$add_argument("--line_type", type = "character", default = "solid",
                    help = "Line type.")
parser$add_argument("--line_opacity", type = "double", default = 0.8,
                    help = "Line opacity [0,1].")
parser$add_argument("--line_thickness", type = "double", default = 1.0,
                    help = "Line thickness.")
parser$add_argument("--is_show_legend", type = "logical", default = TRUE,
                    help = "Show legend.")
parser$add_argument("--legend_title", type = "character", default = "Expression",
                    help = "Legend title.")
parser$add_argument("--legend_title_size", type = "integer", default = 7,
                    help = "Legend title size.")
parser$add_argument("--legend_text_size", type = "integer", default = 7,
                    help = "Legend text size.")
parser$add_argument("--legend_position", type = "character", default = "right",
                    help = "Legend position.")
parser$add_argument("--width", type = "double", default = 2500,
                    help = "Plot width in px for PNG.")
parser$add_argument("--height", type = "double", default = 1350,
                    help = "Plot height in px for PNG.")

args <- parser$parse_args()

# =========================================================
# Utility functions
# =========================================================

ensure_parent_directory <- function(path) {
  parent <- dirname(normalizePath(path, mustWork = FALSE))
  if (!dir.exists(parent)) {
    ok <- dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(parent)) {
      stop("Unable to create output directory: ", parent)
    }
  }
  invisible(NULL)
}

validate_positive <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop("--", name, " must be one positive finite value.")
  }
}

validate_color <- function(x, name) {
  valid <- vapply(x, function(value) {
    tryCatch({
      grDevices::col2rgb(value)
      TRUE
    }, error = function(error) FALSE)
  }, logical(1))
  if (!all(valid)) {
    stop(
      "--", name, " contains invalid R color(s): ",
      paste(x[!valid], collapse = ", ")
    )
  }
}

get_script_path <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  idx <- grep(file_arg, a)
  if (length(idx) == 0) return(normalizePath("."))
  normalizePath(gsub("~\\+~", " ", sub(file_arg, "", a[idx[1]])))
}

make_error_plot <- function(
  msg,
  output_png,
  width,
  height,
  output_pdf = NULL
) {
  p <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5, label = msg,
      size = 4, color = "#555555"
    ) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))
  ensure_parent_directory(output_png)
  ggsave(filename = output_png, plot = p, width = width, height = height,
         units = "px", device = "png", dpi = 300, bg = "white")
  if (!is.null(output_pdf) && nzchar(trimws(output_pdf))) {
    ensure_parent_directory(output_pdf)
    ggsave(
      filename = output_pdf,
      plot = p,
      width = width / 300,
      height = height / 300,
      units = "in",
      device = if (isTRUE(capabilities("cairo"))) cairo_pdf else pdf,
      bg = "white"
    )
  }
  invisible(p)
}

normalize_id_type <- function(id_type) {
  id_type <- toupper(trimws(as.character(id_type)))
  aliases <- c(GENE_SYMBOL = "SYMBOL", HGNC = "SYMBOL", HGNC_SYMBOL = "SYMBOL",
               ENTREZ = "ENTREZID", NCBI_GENE = "ENTREZID", ENSEMBL_GENE = "ENSEMBL",
               ENSEMBL_PROTEIN = "ENSEMBLPROT", MIRBASE = "MIRNA")
  if (length(id_type) != 1L || is.na(id_type) || !nzchar(id_type)) return("AUTO")
  if (id_type %in% names(aliases)) unname(aliases[id_type]) else id_type
}

clean_feature_ids <- function(ids, id_type = "SYMBOL") {
  ids <- trimws(as.character(ids))
  id_type <- normalize_id_type(id_type)
  if (id_type %in% c("ENSEMBL", "ENSEMBLPROT")) {
    ids <- sub("\\.[0-9]+$", "", ids)
  } else if (id_type == "UNIPROT") {
    ids <- sub("-[0-9]+$", "", ids)
  }
  # SYMBOL/ENTREZID identifiers are not truncated: H1-3 and IGKV2-30
  # are valid symbols, not UniProt isoforms.
  ids
}

detect_id_type <- function(ids) {
  ids <- trimws(as.character(ids))
  ids_nonempty <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids_nonempty) == 0) return("UNKNOWN")
  if (mean(grepl("^ENSG[0-9]{11}(\\.[0-9]+)?$", ids_nonempty)) > 0.8) {
    "ENSEMBL"
  } else if (mean(grepl("^ENSP[0-9]{11}(\\.[0-9]+)?$", ids_nonempty)) > 0.8) {
    "ENSEMBLPROT"
  } else if (mean(grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$", ids_nonempty)) > 0.8) {
    "ENSEMBL"
  } else if (mean(grepl("^ENS[A-Z]*P[0-9]+(\\.[0-9]+)?$", ids_nonempty)) > 0.8) {
    "ENSEMBLPROT"
  } else if (mean(grepl("^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?$", ids_nonempty)) > 0.5) {
    "UNIPROT"
  } else if (mean(grepl("^[0-9]+$", ids_nonempty)) > 0.8) {
    "ENTREZID"
  } else if (mean(grepl("^(hsa-)?(mir|let)-", ids_nonempty, ignore.case = TRUE)) > 0.5) {
    "MIRNA"
  } else {
    "SYMBOL"
  }
}

choose_rank_metric <- function(deg) {
  # Schema v2 stores the model test statistic in ``statistic``.  Older input
  # formats are retained as fallbacks.  A present-but-all-NA column must not
  # mask a usable legacy statistic.
  candidates <- c("statistic", "t", "stat", "logFC", "log2FoldChange", "avg_log2FC", "avg_logFC")
  for (candidate in candidates) {
    if (!candidate %in% colnames(deg)) next
    values <- suppressWarnings(as.numeric(deg[[candidate]]))
    if (any(is.finite(values))) return(candidate)
  }
  stop("No finite ranking statistic found. Expected Schema v2 'statistic' or a supported legacy effect/statistic column.")
}

standardize_gsea_result <- function(result = NULL) {
  schema <- list(
    ID = character(),
    Description = character(),
    setSize = integer(),
    enrichmentScore = numeric(),
    NES = numeric(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    rank = integer(),
    leading_edge = character(),
    core_enrichment = character()
  )
  n <- if (is.null(result)) 0L else nrow(result)
  output <- lapply(names(schema), function(column_name) {
    if (!is.null(result) && column_name %in% names(result)) {
      value <- result[[column_name]]
      if (is.character(schema[[column_name]])) {
        return(as.character(value))
      }
      if (is.integer(schema[[column_name]])) {
        return(as.integer(value))
      }
      return(as.numeric(value))
    }
    if (is.character(schema[[column_name]])) {
      return(rep(NA_character_, n))
    }
    if (is.integer(schema[[column_name]])) {
      return(rep(NA_integer_, n))
    }
    rep(NA_real_, n)
  })
  names(output) <- names(schema)
  as.data.frame(output, stringsAsFactors = FALSE, check.names = FALSE)
}

write_empty_result <- function(output_csv) {
  ensure_parent_directory(output_csv)
  write.csv(
    standardize_gsea_result(),
    file = output_csv,
    row.names = FALSE,
    quote = TRUE
  )
}

write_summary <- function(path, summary_list) {
  if (!is.null(path) && nzchar(trimws(path))) {
    ensure_parent_directory(path)
    write.csv(as.data.frame(summary_list), file = path, row.names = FALSE, quote = TRUE)
  }
}

select_top_pathways <- function(result_df, top_n, strategy = "balanced") {
  result_df <- result_df %>%
    filter(!is.na(p.adjust), !is.na(NES), setSize > 0) %>%
    arrange(p.adjust, pvalue, desc(abs(NES)), ID)

  if (nrow(result_df) == 0) return(character())

  if (strategy == "balanced") {
    n_up <- ceiling(top_n / 2)
    n_down <- floor(top_n / 2)
    top_up <- result_df %>%
      filter(NES > 0) %>%
      arrange(p.adjust, pvalue, desc(NES), ID) %>%
      head(n_up)
    top_down <- result_df %>%
      filter(NES < 0) %>%
      arrange(p.adjust, pvalue, NES, ID) %>%
      head(n_down)
    top_ids <- c(top_up$ID, top_down$ID)
    if (length(top_ids) < top_n) {
      add_ids <- result_df %>% filter(!(ID %in% top_ids)) %>% head(top_n - length(top_ids)) %>% pull(ID)
      top_ids <- c(top_ids, add_ids)
    }
    unique(top_ids)
  } else {
    head(result_df$ID, top_n)
  }
}

# =========================================================
# Validate parameters
# =========================================================

if (!file.exists(args$input_rda)) {
  stop("Input RDA does not exist: ", args$input_rda)
}
if (!nzchar(trimws(args$output_png)) || !nzchar(trimws(args$output_csv))) {
  stop("--output_png and --output_csv cannot be empty.")
}
if (length(args$top_n) != 1L || is.na(args$top_n) || args$top_n < 1L) {
  stop("--top_n must be a positive integer.")
}
if (length(args$seed) != 1L || is.na(args$seed)) {
  stop("--seed must be one integer.")
}
if (length(args$min_gs_size) != 1L ||
    is.na(args$min_gs_size) ||
    args$min_gs_size < 1L) {
  stop("--min_gs_size must be a positive integer.")
}
if (length(args$max_gs_size) != 1L ||
    is.na(args$max_gs_size) ||
    args$max_gs_size < args$min_gs_size) {
  stop("--max_gs_size must be an integer no smaller than --min_gs_size.")
}
validate_positive(args$main_title_size, "main_title_size")
validate_positive(args$axis_title_size, "axis_title_size")
validate_positive(args$axis_text_size, "axis_text_size")
validate_positive(args$legend_title_size, "legend_title_size")
validate_positive(args$legend_text_size, "legend_text_size")
validate_positive(args$line_thickness, "line_thickness")
validate_positive(args$width, "width")
validate_positive(args$height, "height")
if (length(args$line_opacity) != 1L ||
    is.na(args$line_opacity) ||
    !is.finite(args$line_opacity) ||
    args$line_opacity < 0 ||
    args$line_opacity > 1) {
  stop("--line_opacity must be between 0 and 1.")
}
if (length(args$line_color) == 0L) {
  stop("--line_color must contain at least one color.")
}
validate_color(args$line_color, "line_color")
args$line_type <- tolower(trimws(args$line_type))
supported_line_types <- c(
  "blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash"
)
if (!args$line_type %in% supported_line_types) {
  stop(
    "--line_type must be one of: ",
    paste(supported_line_types, collapse = ", "),
    "."
  )
}
args$legend_position <- tolower(trimws(args$legend_position))
if (!args$legend_position %in% c(
  "right", "left", "top", "bottom", "inside", "none"
)) {
  stop("--legend_position must be right, left, top, bottom, inside, or none.")
}
args$top_strategy <- tolower(trimws(args$top_strategy))
if (!args$top_strategy %in% c("overall", "balanced")) {
  stop("--top_strategy must be overall or balanced.")
}

# =========================================================
# Load input objects
# =========================================================

input_environment <- new.env(parent = emptyenv())
load(args$input_rda, envir = input_environment)
if (!exists("deg", envir = input_environment, inherits = FALSE)) {
  stop("The 'deg' object does not exist in the provided RDA file.")
}
deg <- as.data.frame(input_environment$deg, stringsAsFactors = FALSE)
if (nrow(deg) == 0L) stop("The 'deg' object contains no features.")
if (is.null(rownames(deg)) || anyNA(rownames(deg)) || any(!nzchar(trimws(rownames(deg))))) {
  stop("The 'deg' object must have non-empty feature IDs as row names.")
}
if (anyDuplicated(rownames(deg))) {
  stop("The 'deg' object contains duplicated feature IDs in its row names.")
}

input_gene_number <- nrow(deg)
row_ids_raw <- rownames(deg)
id_type_requested <- normalize_id_type(args$id_type)
id_type_source <- "explicit_argument"
if (id_type_requested == "AUTO") {
  metadata_id_type <- if (exists("preprocessing", envir = input_environment, inherits = FALSE)) {
    normalize_id_type(input_environment$preprocessing$feature_id_type_output)
  } else "AUTO"
  if (metadata_id_type %in% c("SYMBOL", "ENTREZID", "ENSEMBL", "ENSEMBLPROT", "UNIPROT", "MIRNA")) {
    id_type_detected <- metadata_id_type
    id_type_source <- "canonical_preprocessing"
  } else {
    id_type_detected <- detect_id_type(row_ids_raw)
    id_type_source <- "non_destructive_detection"
  }
} else {
  if (!id_type_requested %in% c("SYMBOL", "ENTREZID", "ENSEMBL", "ENSEMBLPROT", "UNIPROT", "MIRNA"))
    stop("Unsupported --id_type: ", args$id_type)
  id_type_detected <- id_type_requested
}
row_ids_clean <- clean_feature_ids(row_ids_raw, id_type_detected)
message(paste("Resolved ID type:", id_type_detected, "source:", id_type_source))

omics_type <- tolower(args$omics_type)
if (omics_type %in% c("metabolomics", "metabolome")) {
  msg <- "Error: gene-set GSEA is not appropriate for metabolomics features. Use metabolite pathway enrichment instead."
  warning(msg)
  make_error_plot(msg, args$output_png, args$width, args$height, args$output_pdf)
  write_empty_result(args$output_csv)
  write_summary(args$output_summary_csv, list(
    status = "failed", reason = msg, detected_id_type = id_type_detected,
    input_feature_number = input_gene_number, mapped_unique_entrez = 0
  ))
  quit(save = "no", status = 0)
}

if (id_type_detected == "MIRNA" && !(omics_type %in% c("target_projected_mirna", "mirna_target", "target_projected"))) {
  msg <- "Error: direct gene-set GSEA is not appropriate for raw miRNA IDs. First project miRNAs to target genes, then run gene-level GSEA."
  warning(msg)
  make_error_plot(msg, args$output_png, args$width, args$height, args$output_pdf)
  write_empty_result(args$output_csv)
  write_summary(args$output_summary_csv, list(
    status = "failed", reason = msg, detected_id_type = id_type_detected,
    input_feature_number = input_gene_number, mapped_unique_entrez = 0
  ))
  quit(save = "no", status = 0)
}

# A target-projected miRNA enrichment RDA retains the original miRNA row names
# but also carries the mapped target ENTREZID column.  Use that authoritative
# target identifier for GSEA; attempting to send the raw miRNA names to
# org.Hs.eg.db would otherwise fail even though target projection succeeded.
is_target_projected <- omics_type %in% c(
  "target_projected_mirna", "mirna_target", "target_projected"
)
if (is_target_projected && "ENTREZID" %in% colnames(deg)) {
  projected_entrez <- trimws(as.character(deg$ENTREZID))
  if (any(!is.na(projected_entrez) & nzchar(projected_entrez))) {
    row_ids_clean <- projected_entrez
    id_type_detected <- "ENTREZID"
    id_type_source <- "target_projected_entrez_column"
    message("Using target-projected ENTREZID values for miRNA GSEA.")
  }
}

# Attach cleaned IDs for mapping.
deg$.raw_id <- row_ids_raw
if (!exists("row_ids_clean")) row_ids_clean <- clean_feature_ids(row_ids_raw, id_type_detected)
deg$.clean_id <- row_ids_clean
id_type <- id_type_detected
if (id_type == "UNKNOWN") id_type <- "SYMBOL"
deg[[id_type]] <- row_ids_clean

# =========================================================
# Build ranked gene list using all mapped features
# =========================================================

rank_stat <- choose_rank_metric(deg)
statistic_types <- if ("statistic_type" %in% colnames(deg)) {
  unique(tolower(trimws(as.character(deg$statistic_type))))
} else character()
rank_transform <- "none"
if (identical(rank_stat, "statistic") && any(statistic_types == "qlf_f")) {
  # edgeR QLF is an unsigned F statistic.  For a one-coefficient case-control
  # contrast, signed sqrt(F) restores the direction from logFC and has the same
  # ordering interpretation as an absolute t-like statistic.
  rank_transform <- "signed_sqrt_QLF_F_by_logFC"
}
message(paste("Ranking statistic for GSEA:", rank_stat, "transform:", rank_transform))

deg2 <- deg %>%
  mutate(
    .raw_rank_metric = suppressWarnings(as.numeric(.data[[rank_stat]])),
    .rank_logFC = if ("logFC" %in% colnames(.)) suppressWarnings(as.numeric(.data[["logFC"]])) else NA_real_,
    rank_metric = if (rank_transform == "signed_sqrt_QLF_F_by_logFC") {
      sign(.rank_logFC) * sqrt(pmax(.raw_rank_metric, 0))
    } else {
      .raw_rank_metric
    },
    P.Value = if ("P.Value" %in% colnames(.)) suppressWarnings(as.numeric(.data[["P.Value"]])) else NA_real_,
    adj.P.Val = if ("adj.P.Val" %in% colnames(.)) suppressWarnings(as.numeric(.data[["adj.P.Val"]])) else NA_real_
  ) %>%
  filter(!is.na(rank_metric), is.finite(rank_metric))
features_with_finite_rank <- nrow(deg2)

if (nrow(deg2) == 0) {
  msg <- paste0("Error: no finite numeric values in ranking statistic: ", rank_stat, ".")
  warning(msg)
  make_error_plot(msg, args$output_png, args$width, args$height, args$output_pdf)
  write_empty_result(args$output_csv)
  write_summary(args$output_summary_csv, list(
    status = "failed", reason = msg, detected_id_type = id_type_detected,
    ranking_statistic = rank_stat, ranking_transform = rank_transform,
    input_feature_number = input_gene_number,
    mapped_unique_entrez = 0
  ))
  quit(save = "no", status = 0)
}

supported_keytypes <- AnnotationDbi::keytypes(org.Hs.eg.db)
fromType <- id_type
if (!(fromType %in% supported_keytypes)) {
  warning(paste0("Detected ID type '", id_type, "' is not supported by org.Hs.eg.db. Falling back to SYMBOL mapping."))
  fromType <- "SYMBOL"
  id_type <- "SYMBOL"
  deg2[[id_type]] <- deg2$.clean_id
}

# Map IDs to ENTREZID. If already ENTREZID, keep a direct mapping.
if (fromType == "ENTREZID") {
  gene_df <- deg2 %>%
    transmute(ENTREZID = as.character(.data[[id_type]])) %>%
    filter(!is.na(ENTREZID), nzchar(ENTREZID)) %>%
    distinct()
  deg2$ENTREZID <- as.character(deg2[[id_type]])
} else {
  suppressWarnings({
    gene_df <- bitr(
      unique(deg2[[id_type]]),
      fromType = fromType,
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db,
      drop = TRUE
    )
  })
  gene_df <- gene_df %>% distinct(.data[[id_type]], ENTREZID)
  deg2 <- left_join(deg2, gene_df, by = id_type)
}

mapped_feature_number <- sum(!is.na(deg2$ENTREZID))

# De-duplicate ENTREZID by keeping the entry with the largest absolute ranking metric.
deg2 <- deg2 %>%
  filter(!is.na(ENTREZID), nzchar(as.character(ENTREZID))) %>%
  mutate(ENTREZID = as.character(ENTREZID)) %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(rank_metric), n = 1, with_ties = FALSE) %>%
  ungroup()

geneList <- deg2$rank_metric
names(geneList) <- deg2$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!duplicated(names(geneList))]

message(paste("geneList size (mapped & unique ENTREZID):", length(geneList)))

if (length(geneList) < 10) {
  msg <- paste0("Error: geneList too small (", length(geneList), "). ID mapping likely failed or the input is not gene/protein-level data.")
  warning(msg)
  make_error_plot(msg, args$output_png, args$width, args$height, args$output_pdf)
  write_empty_result(args$output_csv)
  write_summary(args$output_summary_csv, list(
    status = "failed", reason = msg, detected_id_type = id_type_detected,
    mapping_from_type = fromType, ranking_statistic = rank_stat,
    ranking_transform = rank_transform,
    input_feature_number = input_gene_number, mapped_feature_number = mapped_feature_number,
    mapped_unique_entrez = length(geneList), gene_set_collection = "MSigDB C2:CP:KEGG_LEGACY"
  ))
  quit(save = "no", status = 0)
}

# =========================================================
# Load the local MSigDB KEGG_LEGACY mapping and run GSEA
# =========================================================

script_dir <- dirname(get_script_path())
mapping_env_path <- Sys.getenv("EVOMICS_KEGG_MAPPING_RDA", unset = "")
mapping_candidates <- unique(c(
  args$kegg_mapping_rda,
  if (nzchar(mapping_env_path)) mapping_env_path else NULL,
  file.path(script_dir, "..", "local_kegg_db", "kegg_mapping.rda"),
  file.path(script_dir, "..", "..", "..", "db", "rda", "kegg_mapping.rda"),
  file.path(script_dir, "..", "msig_kegg.rda")
))
mapping_candidates <- mapping_candidates[
  !is.na(mapping_candidates) & nzchar(trimws(mapping_candidates))
]

explicit_mapping <- !is.null(args$kegg_mapping_rda) &&
  nzchar(trimws(args$kegg_mapping_rda))
if (explicit_mapping && !file.exists(args$kegg_mapping_rda)) {
  stop("The explicitly supplied --kegg_mapping_rda does not exist: ",
       args$kegg_mapping_rda)
}
if (nzchar(mapping_env_path) && !file.exists(mapping_env_path)) {
  stop("EVOMICS_KEGG_MAPPING_RDA points to a missing file: ",
       mapping_env_path)
}

existing_mapping <- mapping_candidates[file.exists(mapping_candidates)]
gene_set_source <- NA_character_
gene_set_version <- NA_character_

if (length(existing_mapping) > 0L) {
  mapping_path <- normalizePath(existing_mapping[1], mustWork = TRUE)
  mapping_environment <- new.env(parent = emptyenv())
  loaded_names <- load(mapping_path, envir = mapping_environment)
  mapping_objects <- loaded_names[
    vapply(
      loaded_names,
      function(object_name) is.data.frame(mapping_environment[[object_name]]),
      logical(1)
    )
  ]
  if (length(mapping_objects) == 0L) {
    stop("Local KEGG mapping RDA contains no data.frame: ", mapping_path)
  }

  preferred_objects <- c("kegg_mapping", "msig_kegg")
  preferred_matches <- preferred_objects[preferred_objects %in% mapping_objects]
  selected_object <- if (length(preferred_matches) > 0L) {
    preferred_matches[1]
  } else {
    mapping_objects[1]
  }
  local_mapping <- mapping_environment[[selected_object]]

  term_column <- intersect(c("gs_name", "term", "ID"), names(local_mapping))[1]
  gene_column <- intersect(
    c("entrez_id", "ncbi_gene", "gene", "ENTREZID"),
    names(local_mapping)
  )[1]
  description_column <- intersect(
    c("gs_description", "Description", "term_name"),
    names(local_mapping)
  )[1]
  if (is.na(term_column) || is.na(gene_column)) {
    stop(
      "Local KEGG mapping must contain a pathway column ",
      "(gs_name/term/ID) and an Entrez column ",
      "(entrez_id/ncbi_gene/gene/ENTREZID)."
    )
  }

  kegg_term2gene <- data.frame(
    term = as.character(local_mapping[[term_column]]),
    gene = as.character(local_mapping[[gene_column]]),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(term), nzchar(term), !is.na(gene), nzchar(gene)) %>%
    distinct()

  kegg_term2name <- NULL
  if (!is.na(description_column)) {
    kegg_term2name <- data.frame(
      term = as.character(local_mapping[[term_column]]),
      name = as.character(local_mapping[[description_column]]),
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(term), nzchar(term), !is.na(name), nzchar(name)) %>%
      distinct(term, .keep_all = TRUE)
  }

  gene_set_source <- mapping_path
  if ("kegg_mapping_metadata" %in% loaded_names) {
    metadata <- mapping_environment$kegg_mapping_metadata
    gene_set_version <- paste(
      c(
        metadata$database_version,
        paste0("msigdbr ", metadata$msigdbr_version)
      ),
      collapse = "; "
    )
  } else if ("version" %in% loaded_names) {
    gene_set_version <- paste(mapping_environment$version, collapse = "; ")
  } else if ("metadata" %in% loaded_names) {
    gene_set_version <- paste(capture.output(str(mapping_environment$metadata)),
                              collapse = " ")
  }
  message("Using local KEGG gene-set mapping: ", mapping_path)
} else {
  warning(
    "No local KEGG mapping was found. Falling back to msigdbr; ",
    "deploy a pinned mapping and pass --kegg_mapping_rda for reproducible web use."
  )
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop(
      "No local KEGG mapping was found and package 'msigdbr' is unavailable. ",
      "Provide --kegg_mapping_rda or set EVOMICS_KEGG_MAPPING_RDA."
    )
  }
  msig_kegg <- tryCatch(
    msigdbr::msigdbr(
      species = "Homo sapiens",
      collection = "C2",
      subcollection = "CP:KEGG_LEGACY"
    ),
    error = function(error) {
      stop("Unable to obtain KEGG gene sets from msigdbr: ",
           conditionMessage(error))
    }
  )
  kegg_term2gene <- msig_kegg %>%
    transmute(term = as.character(gs_name), gene = as.character(ncbi_gene)) %>%
    filter(!is.na(term), nzchar(term), !is.na(gene), nzchar(gene)) %>%
    distinct()
  kegg_term2name <- msig_kegg %>%
    transmute(term = as.character(gs_name), name = as.character(gs_description)) %>%
    filter(!is.na(term), nzchar(term), !is.na(name), nzchar(name)) %>%
    distinct(term, .keep_all = TRUE)
  gene_set_source <- "msigdbr runtime fallback"
  gene_set_version <- as.character(utils::packageVersion("msigdbr"))
}

if (nrow(kegg_term2gene) == 0L) {
  stop("The KEGG TERM2GENE mapping is empty after validation.")
}

set.seed(args$seed)
GSEA_result <- tryCatch(
  clusterProfiler::GSEA(
    geneList = geneList,
    TERM2GENE = kegg_term2gene,
    TERM2NAME = kegg_term2name,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = args$min_gs_size,
    maxGSSize = args$max_gs_size,
    eps = 0,
    seed = TRUE,
    verbose = FALSE
  ),
  error = function(error) {
    stop("GSEA failed: ", conditionMessage(error))
  }
)

if (is.null(GSEA_result) || nrow(as.data.frame(GSEA_result)) == 0L) {
  msg <- "No KEGG pathway passed the gene-set size and mapping requirements."
  warning(msg)
  make_error_plot(
    msg, args$output_png, args$width, args$height, args$output_pdf
  )
  write_empty_result(args$output_csv)
  write_summary(args$output_summary_csv, list(
    status = "empty",
    reason = msg,
    detected_id_type = id_type_detected,
    mapping_from_type = fromType,
    ranking_statistic = rank_stat,
    ranking_transform = rank_transform,
    input_feature_number = input_gene_number,
    features_with_finite_rank = features_with_finite_rank,
    mapped_feature_number = mapped_feature_number,
    mapped_unique_entrez = length(geneList),
    gene_set_collection = "MSigDB C2:CP:KEGG_LEGACY",
    gene_set_source = gene_set_source,
    gene_set_version = gene_set_version,
    seed = args$seed,
    minGSSize = args$min_gs_size,
    maxGSSize = args$max_gs_size
  ))
  quit(save = "no", status = 0)
}

# =========================================================
# Plot GSEA curves
# =========================================================

draw_gsea_multilayer_patchwork <- function(gsea_result, args) {
  result_df <- as.data.frame(gsea_result@result)
  valid_results_df <- result_df %>%
    filter(!is.na(p.adjust), !is.na(NES), setSize > 0)

  top_ids <- select_top_pathways(valid_results_df, args$top_n, args$top_strategy)

  if (length(top_ids) == 0) {
    error_message <- "Error: no valid GSEA pathways are available for plotting."
    make_error_plot(
      error_message, args$output_png, args$width, args$height, args$output_pdf
    )
    message(error_message)
    return(invisible(NULL))
  }

  gene_list <- sort(gsea_result@geneList, decreasing = TRUE)

  es_data <- tryCatch(
    purrr::map_dfr(
      top_ids,
      ~{
        gsdata <- enrichplot:::gsInfo(gsea_result, geneSetID = .x)
        gsdata$pathway <- .x
        gsdata
      }
    ),
    error = function(error) {
      stop(
        "Unable to calculate running enrichment scores with the installed ",
        "enrichplot version: ", conditionMessage(error)
      )
    }
  )

  if (nrow(es_data) == 0) {
    error_message <- "Error: no running-score data could be generated for the selected pathways."
    make_error_plot(
      error_message, args$output_png, args$width, args$height, args$output_pdf
    )
    return(invisible(NULL))
  }

  pathway_order <- unique(es_data$pathway)
  plot_colors <- rep(args$line_color, length.out = length(pathway_order))
  color_values <- setNames(plot_colors, pathway_order)
  description_lookup <- setNames(
    as.character(valid_results_df$Description),
    as.character(valid_results_df$ID)
  )
  label_values <- description_lookup[pathway_order]
  missing_labels <- is.na(label_values) | !nzchar(label_values)
  label_values[missing_labels] <- stringr::str_to_sentence(
    gsub("_", " ", sub("^KEGG_", "", pathway_order[missing_labels]))
  )
  pathway_labels <- setNames(label_values, pathway_order)
  legend_inside <- isTRUE(args$is_show_legend) &&
    identical(args$legend_position, "inside")
  legend_ncol <- if (args$legend_position %in% c("top", "bottom")) 2 else 1

  es_plot <- ggplot(es_data, aes(x = x, y = runningScore, color = pathway)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey75") +
    geom_line(linewidth = args$line_thickness, alpha = args$line_opacity, linetype = args$line_type) +
    scale_color_manual(values = color_values, labels = pathway_labels, name = args$legend_title) +
    labs(title = args$main_title, y = args$y_axis_title, x = NULL) +
    theme_classic(base_size = args$axis_text_size, base_family = "sans") +
    theme(
      plot.title = element_text(
        size = args$main_title_size, face = "bold", hjust = 0.5,
        margin = margin(b = 6)
      ),
      axis.title.y = element_text(size = args$axis_title_size),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = args$axis_text_size),
      legend.position = if (!args$is_show_legend) "none" else if (legend_inside) "inside" else args$legend_position,
      legend.position.inside = c(0.73, 0.72),
      legend.justification = c(0, 0.5),
      legend.background = element_rect(fill = scales::alpha("white", 0.82), colour = NA),
      legend.title = element_text(size = args$legend_title_size, face = "bold"),
      legend.text = element_text(size = args$legend_text_size),
      legend.key.width = grid::unit(0.8, "cm"),
      legend.spacing.y = grid::unit(0.05, "cm"),
      plot.margin = margin(8, 12, 2, 12)
    ) +
    guides(color = guide_legend(ncol = legend_ncol, byrow = TRUE))

  tick_data <- es_data %>%
    filter(position == 1) %>%
    mutate(tick_y = as.numeric(factor(pathway, levels = unique(es_data$pathway))))

  tick_plot <- ggplot(tick_data, aes(x = x, y = tick_y, color = pathway)) +
    geom_linerange(aes(ymin = tick_y - 0.4, ymax = tick_y + 0.4), linewidth = 0.55) +
    scale_color_manual(values = color_values, guide = "none") +
    theme_void() +
    theme(legend.position = "none", axis.text.y = element_blank())

  rank_df <- data.frame(x = seq_along(gene_list), value = as.numeric(gene_list))
  rank_plot <- ggplot(rank_df, aes(x = x, y = value)) +
    geom_col(fill = "grey65", width = 1.0) +
    labs(x = args$x_axis_title, y = paste0("Ranked metric (", rank_stat, ")")) +
    theme_classic(base_size = args$axis_text_size, base_family = "sans") +
    theme(
      axis.title = element_text(size = args$axis_title_size),
      axis.text = element_text(size = args$axis_text_size),
      plot.margin = margin(2, 12, 8, 12)
    )

  final_plot <- es_plot / tick_plot / rank_plot + plot_layout(heights = c(3.2, 0.5, 1.5))

  ensure_parent_directory(args$output_png)
  ggsave(
    filename = args$output_png, plot = final_plot,
    width = args$width, height = args$height,
    units = "px", device = "png", dpi = 300, bg = "white"
  )

  if (!is.null(args$output_pdf) && nzchar(trimws(args$output_pdf))) {
    ensure_parent_directory(args$output_pdf)
    ggsave(
      filename = args$output_pdf, plot = final_plot,
      width = args$width / 300, height = args$height / 300,
      units = "in",
      device = if (isTRUE(capabilities("cairo"))) cairo_pdf else pdf,
      bg = "white"
    )
  }

  invisible(list(plot = final_plot, top_ids = top_ids))
}

if (!is.null(args$output_result_rds) && nzchar(trimws(args$output_result_rds))) {
  ensure_parent_directory(args$output_result_rds)
  saveRDS(GSEA_result, args$output_result_rds)
}
draw_gsea_multilayer_patchwork(GSEA_result, args)

ensure_parent_directory(args$output_csv)
write.csv(
  standardize_gsea_result(as.data.frame(GSEA_result@result)),
  file = args$output_csv,
  row.names = FALSE,
  quote = TRUE
)
write_summary(args$output_summary_csv, list(
  status = "success",
  detected_id_type = id_type_detected,
  identifier_type_source = id_type_source,
  mapping_from_type = fromType,
  ranking_statistic = rank_stat,
  ranking_transform = rank_transform,
  input_feature_number = input_gene_number,
  features_with_finite_rank = features_with_finite_rank,
  mapped_feature_number = mapped_feature_number,
  mapped_unique_entrez = length(geneList),
  gene_set_collection = "MSigDB C2:CP:KEGG_LEGACY",
  gene_set_source = gene_set_source,
  gene_set_version = gene_set_version,
  gsea_pvalue_cutoff = 1,
  minGSSize = args$min_gs_size,
  maxGSSize = args$max_gs_size,
  seed = args$seed,
  comparison_direction = "C minus N when differential_analysis.r is used",
  top_strategy = args$top_strategy,
  top_n = args$top_n
))

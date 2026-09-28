#!/usr/bin/env Rscript
# =============================================================================
# EVOmicsDB differential analysis -- two-group production workflow
# Schema: 2.0.0
#
# Supported comparison:
#   Cancer/case versus non-cancer/control only.
#   All effects are defined as case - control (default C - N).
#
# Supported analysis profiles:
#   rnaseq_raw_counts
#   transcriptomics_normalized
#   microarray
#   proteomics_label_free
#   proteomics_tmt
#   protein_spectral_counts
#   metabolomics_abundance
#
# The script intentionally does not support multi-group, time-series, dose-
# response, or interaction analyses. Unsupported profiles stop with an
# explicit error instead of falling back to numerical guessing. Metabolomics
# is supported as a separate continuous-abundance/Wilcoxon profile; it is not
# routed through the RNA-seq count branch.
#
# Output RDA contract:
#   exp_raw, exp_model, exp_visual, exp (legacy alias of exp_visual), deg,
#   detection_results, feature_metadata, sample_metadata,
#   sample_metadata_all, group_list, group_labels, design, preprocessing,
#   qc, removed_samples, removed_features, voom_weights, model_data.
# =============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(digest)
  library(edgeR)
  library(limma)
})

SCHEMA_VERSION <- "2.0.0"
SCRIPT_VERSION <- "EVOmicsDB_differential_analysis_1.1.0"

## ------------------------------------------------------------------ utilities
stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
warnf <- function(fmt, ...) warning(sprintf(fmt, ...), call. = FALSE)

parse_bool <- function(x, name) {
  if (is.logical(x) && length(x) == 1L && !is.na(x)) return(x)
  key <- tolower(trimws(as.character(x)))
  if (key %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (key %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stopf("--%s must be TRUE or FALSE.", name)
}

normalize_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

validate_probability <- function(x, name, include_zero = TRUE) {
  lower_ok <- if (include_zero) x >= 0 else x > 0
  if (length(x) != 1L || !is.finite(x) || !lower_ok || x > 1) {
    stopf(
      "--%s must be %s.", name,
      if (include_zero) "between 0 and 1" else "> 0 and <= 1"
    )
  }
}

validate_nonnegative <- function(x, name) {
  if (length(x) != 1L || !is.finite(x) || x < 0) {
    stopf("--%s must be a finite non-negative number.", name)
  }
}

validate_positive_integer <- function(x, name) {
  if (length(x) != 1L || is.na(x) || x < 1L || x != as.integer(x)) {
    stopf("--%s must be a positive integer.", name)
  }
}

ensure_parent_directory <- function(path) {
  parent <- dirname(normalizePath(path, mustWork = FALSE))
  if (!dir.exists(parent)) {
    ok <- dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(parent)) stopf("Unable to create directory: %s", parent)
  }
  invisible(NULL)
}

derive_output_path <- function(output_rda, suffix) {
  stem <- if (grepl("\\.[Rr][Dd][Aa]$", output_rda)) {
    sub("\\.[Rr][Dd][Aa]$", "", output_rda)
  } else {
    output_rda
  }
  paste0(stem, suffix)
}

read_tabular_file <- function(path) {
  if (!file.exists(path)) stopf("File does not exist: %s", path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required to read Excel metadata files.", call. = FALSE)
    }
    return(as.data.frame(readxl::read_excel(path, sheet = 1L), check.names = FALSE))
  }
  sep <- if (ext == "csv") "," else "\t"
  read.table(
    path, sep = sep, header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, quote = "\"", comment.char = "",
    na.strings = c("NA", "NaN", "")
  )
}

looks_log2 <- function(mat) {
  x <- as.numeric(mat)
  x <- x[is.finite(x)]
  if (length(x) < 30L) return(NA)
  q <- stats::quantile(x, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  isTRUE(q[[2L]] <= 50 && q[[1L]] >= -20)
}

warn_scale_only <- function(mat, declared_scale) {
  diagnostic <- looks_log2(mat)
  if (is.na(diagnostic) || declared_scale %in% c("count", "auto")) return(invisible(NULL))
  if (declared_scale == "log2" && !diagnostic) {
    warnf("input_scale=log2, but the matrix does not look log2-like; verify Table S1 configuration.")
  }
  if (declared_scale == "linear" && diagnostic) {
    warnf("input_scale=linear, but the matrix looks log2-like; verify configuration to avoid double transformation.")
  }
  invisible(NULL)
}

assert_matrix <- function(mat, context) {
  if (!is.matrix(mat) || !is.numeric(mat)) stopf("%s is not a numeric matrix.", context)
  if (nrow(mat) == 0L) stopf("No features remain after %s.", context)
  if (ncol(mat) < 2L) stopf("Fewer than two samples remain after %s.", context)
  if (any(!is.finite(mat) & !is.na(mat))) stopf("%s contains infinite values.", context)
  invisible(NULL)
}

required_package <- function(package, purpose) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stopf("Package '%s' is required for %s. Install it in the server renv environment.", package, purpose)
  }
}

## -------------------------------------------------------------- CLI arguments
parser <- ArgumentParser(
  description = paste(
    "EVOmicsDB two-group differential analysis.",
    "Contrast is Cancer/case minus Non-cancer/control."
  )
)

# Core I/O and legacy-compatible names.
parser$add_argument("--counts_file", required = TRUE)
parser$add_argument("--output_rda", required = TRUE)
parser$add_argument("--output_csv", default = NULL)
parser$add_argument("--output_all_csv", default = NULL)
parser$add_argument("--output_detection_csv", default = NULL)
parser$add_argument("--is_output_csv", default = "TRUE")

# Profile selection. analysis_profile is authoritative; quantification_method
# is retained for old backend calls and for audit.
parser$add_argument(
  "--analysis_profile",
  default = "auto",
  choices = c(
    "auto", "rnaseq_raw_counts", "rnaseq_fractional_counts", "transcriptomics_normalized", "microarray",
    "proteomics_label_free", "proteomics_tmt", "protein_spectral_counts",
    "metabolomics_abundance"
  )
)
parser$add_argument("--quantification_method", default = "")
parser$add_argument("--omics_layer", default = "")
parser$add_argument("--molecule_type", default = "")
parser$add_argument("--detection_method", default = "")

# Dataset-level preprocessing contract.
parser$add_argument(
  "--input_scale", default = "auto",
  choices = c("auto", "count", "linear", "log2")
)
parser$add_argument("--normalization_applied", default = "unknown")
parser$add_argument(
  "--normalization_todo", default = "auto",
  choices = c("auto", "none", "median_center", "quantile", "irs", "tmmwsp")
)
parser$add_argument("--transform_applied", default = "unknown")
parser$add_argument(
  "--transform_todo", default = "auto",
  choices = c("auto", "none", "log2", "log2p1")
)
parser$add_argument(
  "--missing_value_encoding", default = "NA",
  help = "Audit label, e.g. NA|blank|zero|none. Input reader always recognizes NA/blank."
)
parser$add_argument("--zero_is_missing", default = "FALSE")

# Exactly two groups. Metadata is preferred; C/N prefixes are a legacy fallback.
parser$add_argument("--group_file", default = NULL)
parser$add_argument("--sample_column", default = "sample")
parser$add_argument("--group_column", default = "group")
parser$add_argument("--case_code", default = "C")
parser$add_argument("--control_code", default = "N")
parser$add_argument("--case_prefix", default = "^C")
parser$add_argument("--control_prefix", default = "^N")
parser$add_argument("--case_display", default = "Cancer")
parser$add_argument("--control_display", default = "Non-cancer")
parser$add_argument("--paired_design", default = "FALSE")
parser$add_argument("--pair_column", default = "pair_id")
parser$add_argument("--batch_column", default = "")
parser$add_argument("--plex_column", default = "")
parser$add_argument("--visual_remove_batch", default = "TRUE")

# Feature/sample filtering.
parser$add_argument(
  "--feature_filter_mode", default = "auto",
  choices = c("auto", "none", "mean_detection", "group_valid", "filterByExpr")
)
parser$add_argument("--sample_max_missing_prop", type = "numeric", default = 1.0)
parser$add_argument("--min_valid_prop_per_group", type = "numeric", default = 0.5)
parser$add_argument("--min_valid_n_per_group", type = "integer", default = 3L)
parser$add_argument("--min_residual_df", type = "integer", default = 2L)
parser$add_argument("--min_mean", type = "numeric", default = 1)
parser$add_argument("--min_detect_prop", type = "numeric", default = 0.5)

# Retained legacy arguments. max_na is used only when explicitly supplied by an
# old call; new calls should use group-valid filtering.
parser$add_argument(
  "--max_na", type = "numeric", default = -1,
  help = "Deprecated legacy field. Values >=0 are audit-only in schema 2.0."
)
parser$add_argument("--apply_detect_filter", default = "FALSE")
parser$add_argument("--allow_rounding", default = "FALSE")

# Model and visual missingness are deliberately separated.
parser$add_argument(
  "--protein_model", default = "limma_filtered",
  choices = c("limma_filtered")
)
parser$add_argument(
  "--model_imputation", default = "none",
  help = "Schema 2.0 production value is 'none'; model imputation is prohibited."
)
parser$add_argument(
  "--visual_imputation", default = "auto",
  choices = c("auto", "none", "row_median", "minprob", "knn", "zero")
)
parser$add_argument("--minprob_q", type = "numeric", default = 0.01)
parser$add_argument("--minprob_tune_sigma", type = "numeric", default = 0.3)
parser$add_argument("--detection_shift", default = "TRUE")
parser$add_argument("--seed", type = "integer", default = 123L)

# TMT and microarray-specific configuration.
parser$add_argument("--tmt_channel_normalized", default = "FALSE")
parser$add_argument("--probe_map_file", default = NULL)
parser$add_argument("--probe_id_column", default = "probe_id")
parser$add_argument("--gene_id_column", default = "gene_id")
parser$add_argument(
  "--probe_collapse_method", default = "median",
  choices = c("median", "mean", "max_variance")
)

# Identifier/audit contract.
parser$add_argument("--feature_id_type", default = "symbol")
parser$add_argument("--output_feature_id_type", default = "")
parser$add_argument("--organism", default = "Homo sapiens")
parser$add_argument("--organism_taxid", default = "9606")

# User-facing export thresholds. They classify the per-request CSV only; the
# cached RDA stores model statistics independently of these display choices.
parser$add_argument("--is_use_padj", default = "TRUE")
parser$add_argument("--log2fc", type = "numeric", default = 1)
parser$add_argument("--pvalue", type = "numeric", default = 0.05)

args <- parser$parse_args()

is_output_csv <- parse_bool(args$is_output_csv, "is_output_csv")
zero_is_missing <- parse_bool(args$zero_is_missing, "zero_is_missing")
paired_design <- parse_bool(args$paired_design, "paired_design")
visual_remove_batch <- parse_bool(args$visual_remove_batch, "visual_remove_batch")
allow_rounding <- parse_bool(args$allow_rounding, "allow_rounding")
apply_detect_filter <- parse_bool(args$apply_detect_filter, "apply_detect_filter")
detection_shift_enabled <- parse_bool(args$detection_shift, "detection_shift")
tmt_channel_normalized <- parse_bool(args$tmt_channel_normalized, "tmt_channel_normalized")
is_use_padj <- parse_bool(args$is_use_padj, "is_use_padj")
if (args$model_imputation != "none") {
  stop("model_imputation must be 'none'; imputation is restricted to exp_visual.", call. = FALSE)
}

validate_probability(args$sample_max_missing_prop, "sample_max_missing_prop")
validate_probability(args$min_valid_prop_per_group, "min_valid_prop_per_group")
validate_probability(args$min_detect_prop, "min_detect_prop")
validate_probability(args$minprob_q, "minprob_q", include_zero = FALSE)
validate_probability(args$pvalue, "pvalue", include_zero = FALSE)
validate_positive_integer(args$min_valid_n_per_group, "min_valid_n_per_group")
validate_positive_integer(args$min_residual_df, "min_residual_df")
validate_nonnegative(args$min_mean, "min_mean")
validate_nonnegative(args$minprob_tune_sigma, "minprob_tune_sigma")
validate_nonnegative(args$log2fc, "log2fc")
if (args$max_na >= 0) validate_probability(args$max_na, "max_na")
if (args$case_code == args$control_code) stop("case_code and control_code must differ.", call. = FALSE)
set.seed(args$seed)

## ------------------------------------------------------------- profile mapping
quant_key <- normalize_key(args$quantification_method)
quant_alias <- c(
  raw_count = "raw_counts", counts = "raw_counts", count = "raw_counts",
  log2normalized = "log2_normalized", normalized_log2 = "log2_normalized",
  lfq_ibaq = "lfq", ibaq = "lfq", label_free = "lfq",
  spectral_count = "spectral_counts", spc = "spectral_counts",
  reporter_ion_intensity = "tmt",
  # Common label-free/DIA software and metric labels used in public datasets.
  dia_nn = "dia", diann = "dia", spectronaut = "dia", swath = "dia",
  dia_label_free_intensity = "dia",
  maxlfq = "lfq", lfq_intensity = "lfq", top3 = "lfq",
  prm = "lfq", mrm = "lfq", srm = "lfq",
  intensity = "lfq", normalized_intensity = "lfq"
)
if (quant_key %in% names(quant_alias)) quant_key <- unname(quant_alias[[quant_key]])

profile_from_legacy <- function(qkey, omics_layer, detection_method) {
  # Molecule-layer information takes precedence over the ambiguous legacy
  # label "relative intensity", which is also used by label-free proteomics.
  if (grepl("metabol", omics_layer, ignore.case = TRUE)) return("metabolomics_abundance")
  if (qkey == "raw_counts") return("rnaseq_raw_counts")
  if (qkey %in% c("spectral_counts", "spectral_counting")) return("protein_spectral_counts")
  if (qkey %in% c("tmt", "itraq")) return("proteomics_tmt")
  if (qkey %in% c("lfq", "dia", "relative_intensity")) return("proteomics_label_free")
  if (qkey %in% c("log2_normalized", "fpkm", "tpm", "rpm")) {
    if (grepl("microarray", detection_method, ignore.case = TRUE)) return("microarray")
    return("transcriptomics_normalized")
  }
  stop(
    "analysis_profile=auto could not be resolved from quantification_method. ",
    "Add Analysis Profile to Table S1 and pass --analysis_profile explicitly.",
    call. = FALSE
  )
}

analysis_profile <- if (args$analysis_profile == "auto") {
  profile_from_legacy(quant_key, args$omics_layer, args$detection_method)
} else {
  args$analysis_profile
}

supported_profiles <- c(
  "rnaseq_raw_counts", "rnaseq_fractional_counts", "transcriptomics_normalized", "microarray",
  "proteomics_label_free", "proteomics_tmt", "protein_spectral_counts",
  "metabolomics_abundance"
)
if (!analysis_profile %in% supported_profiles) {
  stopf(
    "Unsupported analysis_profile '%s'. This production version supports only: %s.",
    analysis_profile, paste(supported_profiles, collapse = ", ")
  )
}

resolve_input_scale <- function(profile, requested, qkey) {
  if (requested != "auto") return(list(value = requested, source = "explicit"))
  value <- switch(
    profile,
    rnaseq_raw_counts = "count",
    rnaseq_fractional_counts = "count",
    protein_spectral_counts = "count",
    transcriptomics_normalized = if (qkey == "log2_normalized") "log2" else "linear",
    microarray = "log2",
    proteomics_label_free = if (qkey == "log2_normalized") "log2" else "linear",
    proteomics_tmt = if (qkey == "log2_normalized") "log2" else "linear",
    metabolomics_abundance = "linear"
  )
  list(value = value, source = "derived_from_profile_and_quantification_method")
}
scale_resolution <- resolve_input_scale(analysis_profile, args$input_scale, quant_key)
input_scale <- scale_resolution$value
# Rank-sum/signed-rank testing does not adjust covariates. Never claim otherwise.
if (analysis_profile == "metabolomics_abundance" &&
    (nzchar(trimws(args$batch_column)) || nzchar(trimws(args$plex_column)))) {
  stop("The metabolomics Wilcoxon profile does not support batch/plex covariate adjustment. Remove those declarations or use a separately specified adjusted model.", call. = FALSE)
}
if (analysis_profile %in% c("rnaseq_raw_counts", "rnaseq_fractional_counts", "protein_spectral_counts") && input_scale != "count") {
  stopf("%s requires input_scale=count.", analysis_profile)
}
if (!analysis_profile %in% c("rnaseq_raw_counts", "rnaseq_fractional_counts", "protein_spectral_counts") &&
    !input_scale %in% c("linear", "log2")) {
  stopf("%s requires input_scale=linear or log2.", analysis_profile)
}

## ---------------------------------------------------------------- input matrix
if (!file.exists(args$counts_file)) stopf("Input matrix does not exist: %s", args$counts_file)
input_sha256 <- digest::digest(file = args$counts_file, algo = "sha256")
input_table <- read.table(
  args$counts_file, sep = "\t", header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE, quote = "", comment.char = "",
  na.strings = c("NA", "NaN", "")
)
if (ncol(input_table) < 3L) {
  stop("Input must contain one feature-ID column and at least two sample columns.", call. = FALSE)
}

feature_ids <- trimws(as.character(input_table[[1L]]))
if (anyNA(feature_ids) || any(feature_ids == "")) stop("Feature identifiers cannot be empty.", call. = FALSE)
if (anyDuplicated(feature_ids)) {
  duplicate_ids <- unique(feature_ids[duplicated(feature_ids)])
  stopf(
    "Duplicate feature IDs are not allowed before profile-specific aggregation. Examples: %s",
    paste(utils::head(duplicate_ids, 10L), collapse = ", ")
  )
}

sample_input <- input_table[-1L]
sample_names <- colnames(sample_input)
if (anyNA(sample_names) || any(sample_names == "") || anyDuplicated(sample_names)) {
  stop("Sample column names must be non-empty and unique.", call. = FALSE)
}
original_na <- is.na(sample_input)
numeric_matrix <- suppressWarnings(matrix(
  as.numeric(as.matrix(sample_input)), nrow = nrow(sample_input),
  ncol = ncol(sample_input), dimnames = list(feature_ids, sample_names)
))
introduced_na <- is.na(numeric_matrix) & !original_na
if (any(introduced_na)) {
  bad_columns <- unique(colnames(numeric_matrix)[col(introduced_na)[introduced_na]])
  stopf("Non-numeric values found in sample columns: %s", paste(bad_columns, collapse = ", "))
}
assert_matrix(numeric_matrix, "input reading")
exp_raw <- numeric_matrix
all_samples <- colnames(exp_raw)
all_features <- rownames(exp_raw)

## ---------------------------------------------------------- grouping metadata
map_group_values <- function(values) {
  raw <- trimws(as.character(values))
  lower <- tolower(raw)
  case_aliases <- unique(tolower(c(args$case_code, args$case_display, "case", "cancer", "tumor", "patient")))
  control_aliases <- unique(tolower(c(
    args$control_code, args$control_display, "control", "non-cancer", "noncancer", "normal", "healthy"
  )))
  mapped <- rep(NA_character_, length(raw))
  mapped[lower %in% case_aliases] <- "case"
  mapped[lower %in% control_aliases] <- "control"
  mapped
}

if (!is.null(args$group_file) && nzchar(args$group_file)) {
  metadata_input <- read_tabular_file(args$group_file)
  required_columns <- c(args$sample_column, args$group_column)
  missing_columns <- setdiff(required_columns, colnames(metadata_input))
  if (length(missing_columns)) {
    stopf("group_file is missing columns: %s", paste(missing_columns, collapse = ", "))
  }
  metadata_input[[args$sample_column]] <- trimws(as.character(metadata_input[[args$sample_column]]))
  if (anyDuplicated(metadata_input[[args$sample_column]])) stop("group_file contains duplicate sample IDs.", call. = FALSE)
  index <- match(all_samples, metadata_input[[args$sample_column]])
  if (anyNA(index)) {
    stopf("group_file lacks expression samples: %s", paste(all_samples[is.na(index)], collapse = ", "))
  }
  sample_metadata_all <- metadata_input[index, , drop = FALSE]
  rownames(sample_metadata_all) <- all_samples
  canonical_group <- map_group_values(sample_metadata_all[[args$group_column]])
  if (anyNA(canonical_group)) {
    bad <- unique(sample_metadata_all[[args$group_column]][is.na(canonical_group)])
    stopf(
      "Exactly two groups are supported. Unrecognized group value(s): %s",
      paste(bad, collapse = ", ")
    )
  }
  grouping_source <- "group_file"
} else {
  case_hit <- grepl(args$case_prefix, all_samples, perl = TRUE, ignore.case = TRUE)
  control_hit <- grepl(args$control_prefix, all_samples, perl = TRUE, ignore.case = TRUE)
  if (any(case_hit & control_hit)) {
    stopf("Sample prefix rules overlap for: %s", paste(all_samples[case_hit & control_hit], collapse = ", "))
  }
  unassigned <- !case_hit & !control_hit
  if (any(unassigned)) {
    stopf(
      "Samples could not be assigned by case/control prefixes: %s. Provide --group_file or correct Table S1 group configuration.",
      paste(all_samples[unassigned], collapse = ", ")
    )
  }
  canonical_group <- ifelse(case_hit, "case", "control")
  sample_metadata_all <- data.frame(
    sample = all_samples,
    group = ifelse(canonical_group == "case", args$case_code, args$control_code),
    stringsAsFactors = FALSE,
    row.names = all_samples,
    check.names = FALSE
  )
  grouping_source <- "prefix_fallback"
}
sample_metadata_all$canonical_group <- canonical_group
# Include zero-coded missing abundance before sample QC; keep exp_raw intact.
sample_missing_mask <- is.na(exp_raw)
if (zero_is_missing && !analysis_profile %in% c("rnaseq_raw_counts", "rnaseq_fractional_counts", "protein_spectral_counts")) {
  sample_missing_mask <- sample_missing_mask | (!is.na(exp_raw) & exp_raw == 0)
}
sample_metadata_all$sample_missing_fraction <- colMeans(sample_missing_mask)

keep_samples <- sample_metadata_all$sample_missing_fraction <= args$sample_max_missing_prop
removed_samples <- all_samples[!keep_samples]
exp_work <- exp_raw[, keep_samples, drop = FALSE]
sample_metadata <- sample_metadata_all[keep_samples, , drop = FALSE]
sample_metadata$retained_for_model <- TRUE
sample_metadata_all$retained_for_model <- keep_samples
assert_matrix(exp_work, "sample missingness filtering")

group_list <- factor(sample_metadata$canonical_group, levels = c("control", "case"))
names(group_list) <- colnames(exp_work)
group_counts <- table(group_list)
if (length(group_counts) != 2L || any(group_counts == 0L)) {
  stopf(
    "Exactly two non-empty groups are required after sample QC. Counts: %s",
    paste(names(group_counts), as.integer(group_counts), collapse = ", ")
  )
}
if (any(group_counts < 2L)) {
  stopf(
    "At least two samples per group are required. Counts: %s",
    paste(names(group_counts), as.integer(group_counts), collapse = ", ")
  )
}

## ---------------------------------------------------------- design construction
design_data <- data.frame(group = group_list, row.names = colnames(exp_work))
design_terms <- character()
batch_values <- NULL
plex_values <- NULL
pair_values <- NULL

if (nzchar(args$batch_column)) {
  if (!args$batch_column %in% colnames(sample_metadata)) {
    stopf("batch_column '%s' is absent from group_file.", args$batch_column)
  }
  batch_values <- factor(sample_metadata[[args$batch_column]])
  if (anyNA(batch_values)) stop("Batch values cannot be missing.", call. = FALSE)
  if (nlevels(batch_values) > 1L) {
    design_data$batch <- batch_values
    design_terms <- c(design_terms, "batch")
  } else {
    warnf("batch_column has one level and is omitted from design.")
  }
}

if (nzchar(args$plex_column)) {
  if (!args$plex_column %in% colnames(sample_metadata)) {
    stopf("plex_column '%s' is absent from group_file.", args$plex_column)
  }
  plex_values <- factor(sample_metadata[[args$plex_column]])
  if (anyNA(plex_values)) stop("Plex values cannot be missing.", call. = FALSE)
  if (nlevels(plex_values) > 1L && args$plex_column != args$batch_column) {
    design_data$plex <- plex_values
    design_terms <- c(design_terms, "plex")
  }
}

if (paired_design) {
  if (is.null(args$group_file) || !nzchar(args$group_file)) {
    stop("paired_design=TRUE requires group_file with a pair column.", call. = FALSE)
  }
  if (!args$pair_column %in% colnames(sample_metadata)) {
    stopf("pair_column '%s' is absent from group_file.", args$pair_column)
  }
  pair_values <- factor(sample_metadata[[args$pair_column]])
  if (anyNA(pair_values)) stop("Pair IDs cannot be missing.", call. = FALSE)
  pair_table <- table(pair_values, group_list)
  if (any(pair_table != 1L)) {
    stop(
      "Paired design requires exactly one case and one control sample per pair after sample QC.",
      call. = FALSE
    )
  }
  design_data$pair <- pair_values
  design_terms <- c(design_terms, "pair")
}

design_terms <- c(design_terms, "group")
design_formula <- stats::reformulate(design_terms)
design <- stats::model.matrix(design_formula, data = design_data)
rownames(design) <- rownames(design_data)
if (qr(design)$rank < ncol(design)) {
  stop(
    "Design matrix is not full rank. Group is confounded with batch/pair, or metadata are redundant.",
    call. = FALSE
  )
}
group_coef <- grep("^groupcase$", colnames(design), value = TRUE)
if (length(group_coef) != 1L) stop("Internal error: case-control coefficient was not created.", call. = FALSE)

group_labels <- list(
  case_code = args$case_code,
  control_code = args$control_code,
  case_display = args$case_display,
  control_display = args$control_display,
  contrast = paste0(args$case_code, "_minus_", args$control_code),
  effect_definition = "case_minus_control"
)

## ---------------------------------------------------- preprocessing resolution
normalization_applied_key <- normalize_key(args$normalization_applied)
normalization_todo <- args$normalization_todo
if (normalization_todo == "auto") {
  already_normalized <- !normalization_applied_key %in% c("", "none", "unknown", "not_applied")
  normalization_todo <- switch(
    analysis_profile,
    rnaseq_raw_counts = "none",
    rnaseq_fractional_counts = "none",
    protein_spectral_counts = "none",
    transcriptomics_normalized = "none",
    microarray = "none",
    proteomics_label_free = if (already_normalized) "none" else "median_center",
    proteomics_tmt = if (already_normalized || tmt_channel_normalized) "none" else "median_center",
    metabolomics_abundance = if (already_normalized) "none" else "median_center"
  )
}
if (normalization_todo == "tmmwsp" && analysis_profile != "rnaseq_raw_counts") {
  stop("normalization_todo=tmmwsp is only valid for rnaseq_raw_counts.", call. = FALSE)
}

transform_todo <- args$transform_todo
if (transform_todo == "auto") {
  transform_todo <- if (input_scale %in% c("count", "log2")) {
    "none"
  } else if (analysis_profile %in% c("proteomics_label_free", "proteomics_tmt", "metabolomics_abundance")) {
    if (zero_is_missing || !any(exp_work == 0, na.rm = TRUE)) "log2" else "log2p1"
  } else {
    "log2p1"
  }
}
if (input_scale == "log2" && transform_todo != "none") {
  stop("A log2 input must use transform_todo=none to prevent double transformation.", call. = FALSE)
}
if (input_scale == "linear" && transform_todo == "none") {
  stop(
    "Continuous linear input requires transform_todo=log2 or log2p1 so that logFC has a log2 interpretation.",
    call. = FALSE
  )
}

feature_filter_mode <- args$feature_filter_mode
if (feature_filter_mode == "auto") {
  feature_filter_mode <- switch(
    analysis_profile,
    rnaseq_raw_counts = "filterByExpr",
    rnaseq_fractional_counts = "filterByExpr",
    protein_spectral_counts = "filterByExpr",
    transcriptomics_normalized = if (apply_detect_filter) "mean_detection" else "none",
    microarray = "none",
    proteomics_label_free = "group_valid",
    proteomics_tmt = "group_valid",
    metabolomics_abundance = "group_valid"
  )
}
if (analysis_profile %in% c("rnaseq_raw_counts", "rnaseq_fractional_counts", "protein_spectral_counts") &&
    feature_filter_mode != "filterByExpr") {
  stopf("%s requires feature_filter_mode=filterByExpr.", analysis_profile)
}

if (args$max_na >= 0) {
  warnf(
    "Legacy --max_na was supplied. Schema 2.0 uses group-aware filtering; max_na is recorded but does not override explicit profile rules."
  )
}

## --------------------------------------------------------- feature processing
removed_feature_log <- data.frame(
  feature_id = character(), reason = character(), stringsAsFactors = FALSE
)
log_removed <- function(ids, reason) {
  if (!length(ids)) return(invisible(NULL))
  removed_feature_log <<- rbind(
    removed_feature_log,
    data.frame(feature_id = as.character(ids), reason = reason, stringsAsFactors = FALSE)
  )
  invisible(NULL)
}

apply_continuous_transform <- function(mat, method) {
  warn_scale_only(mat, input_scale)
  if (method == "none") return(mat)
  if (any(mat < 0, na.rm = TRUE)) {
    stop("Linear matrix contains negative values; verify input_scale and transform settings.", call. = FALSE)
  }
  if (method == "log2") {
    if (any(mat <= 0, na.rm = TRUE)) {
      stop(
        "transform_todo=log2 requires strictly positive observed values. Set zero_is_missing=TRUE when zero denotes missing, or use log2p1.",
        call. = FALSE
      )
    }
    return(log2(mat))
  }
  if (method == "log2p1") return(log2(mat + 1))
  stopf("Unsupported transformation: %s", method)
}

median_center_log2 <- function(mat) {
  medians <- apply(mat, 2L, stats::median, na.rm = TRUE)
  if (any(!is.finite(medians))) stop("Median centering failed: a sample has no finite observations.", call. = FALSE)
  target <- stats::median(medians)
  sweep(mat, 2L, medians, FUN = "-") + target
}

apply_irs_log2 <- function(mat, plex) {
  if (is.null(plex)) stop("normalization_todo=irs requires plex_column in group_file.", call. = FALSE)
  plex <- droplevels(factor(plex))
  if (nlevels(plex) < 2L) stop("IRS requires at least two plexes.", call. = FALSE)
  # Protein-specific Internal Reference Scaling on the log2 scale:
  # compute log2 sums per protein/plex, then align each plex sum to the
  # geometric-mean target. This requires summarized reporter intensities.
  plex_log_sums <- vapply(levels(plex), function(level) {
    block <- mat[, plex == level, drop = FALSE]
    apply(block, 1L, function(x) {
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      m <- max(x)
      m + log2(sum(2^(x - m)))
    })
  }, numeric(nrow(mat)))
  if (is.null(dim(plex_log_sums))) plex_log_sums <- matrix(plex_log_sums, ncol = nlevels(plex))
  rownames(plex_log_sums) <- rownames(mat)
  colnames(plex_log_sums) <- levels(plex)
  target <- rowMeans(plex_log_sums, na.rm = TRUE)
  target[!is.finite(target)] <- NA_real_
  result <- mat
  for (level in levels(plex)) {
    shift <- target - plex_log_sums[, level]
    result[, plex == level] <- sweep(result[, plex == level, drop = FALSE], 1L, shift, FUN = "+")
  }
  result
}

collapse_microarray_probes <- function(mat) {
  if (is.null(args$probe_map_file) || !nzchar(args$probe_map_file)) {
    if (normalize_key(args$feature_id_type) == "probe") {
      stop("feature_id_type=probe requires --probe_map_file for gene-level analysis.", call. = FALSE)
    }
    metadata <- data.frame(
      feature_id = rownames(mat), source_features = rownames(mat),
      stringsAsFactors = FALSE
    )
    return(list(matrix = mat, metadata = metadata, output_id_type = args$feature_id_type))
  }
  mapping <- read_tabular_file(args$probe_map_file)
  required <- c(args$probe_id_column, args$gene_id_column)
  if (!all(required %in% colnames(mapping))) {
    stopf("probe_map_file must contain: %s", paste(required, collapse = ", "))
  }
  mapping <- mapping[, required, drop = FALSE]
  colnames(mapping) <- c("probe_id", "gene_id")
  mapping$probe_id <- trimws(as.character(mapping$probe_id))
  mapping$gene_id <- trimws(as.character(mapping$gene_id))
  mapping <- mapping[
    mapping$probe_id %in% rownames(mat) & !is.na(mapping$gene_id) & mapping$gene_id != "",
    , drop = FALSE
  ]
  if (!nrow(mapping)) stop("No expression probes mapped to genes.", call. = FALSE)
  gene_ids <- unique(mapping$gene_id)
  collapsed <- matrix(
    NA_real_, nrow = length(gene_ids), ncol = ncol(mat),
    dimnames = list(gene_ids, colnames(mat))
  )
  source_features <- character(length(gene_ids))
  for (i in seq_along(gene_ids)) {
    probes <- unique(mapping$probe_id[mapping$gene_id == gene_ids[[i]]])
    block <- mat[probes, , drop = FALSE]
    source_features[[i]] <- paste(probes, collapse = ";")
    if (args$probe_collapse_method == "median") {
      collapsed[i, ] <- apply(block, 2L, stats::median, na.rm = TRUE)
    } else if (args$probe_collapse_method == "mean") {
      collapsed[i, ] <- colMeans(block, na.rm = TRUE)
    } else {
      variances <- apply(block, 1L, stats::var, na.rm = TRUE)
      variances[!is.finite(variances)] <- -Inf
      selected <- if (all(variances == -Inf)) 1L else which.max(variances)
      collapsed[i, ] <- block[selected, ]
    }
  }
  collapsed[!is.finite(collapsed)] <- NA_real_
  unmapped <- setdiff(rownames(mat), unique(mapping$probe_id))
  log_removed(unmapped, "microarray_unmapped_probe")
  metadata <- data.frame(
    feature_id = gene_ids, source_features = source_features,
    stringsAsFactors = FALSE
  )
  output_type <- if (nzchar(args$output_feature_id_type)) args$output_feature_id_type else "symbol"
  list(matrix = collapsed, metadata = metadata, output_id_type = output_type)
}

compute_detection_results <- function(mat, group, pair = NULL) {
  detected <- !is.na(mat)
  case_index <- group == "case"
  control_index <- group == "control"
  result <- data.frame(
    feature_id = rownames(mat),
    detected_case = rowSums(detected[, case_index, drop = FALSE]),
    total_case = sum(case_index),
    detected_control = rowSums(detected[, control_index, drop = FALSE]),
    total_control = sum(control_index),
    stringsAsFactors = FALSE
  )
  result$detection_rate_case <- result$detected_case / result$total_case
  result$detection_rate_control <- result$detected_control / result$total_control
  result$detection_rate_difference <- result$detection_rate_case - result$detection_rate_control
  result$odds_ratio <- NA_real_
  result$statistic <- NA_real_
  result$statistic_type <- if (is.null(pair)) "fisher_exact" else "exact_mcnemar_binomial"
  result$P.Value <- NA_real_

  if (is.null(pair)) {
    for (i in seq_len(nrow(mat))) {
      tab <- matrix(c(
        result$detected_case[[i]], result$total_case[[i]] - result$detected_case[[i]],
        result$detected_control[[i]], result$total_control[[i]] - result$detected_control[[i]]
      ), nrow = 2L, byrow = TRUE)
      test <- suppressWarnings(stats::fisher.test(tab))
      result$odds_ratio[[i]] <- unname(test$estimate)
      result$P.Value[[i]] <- test$p.value
    }
  } else {
    pair <- factor(pair)
    pair_levels <- levels(pair)
    for (i in seq_len(nrow(mat))) {
      case_detect <- vapply(pair_levels, function(p) {
        detected[i, pair == p & case_index][[1L]]
      }, logical(1))
      control_detect <- vapply(pair_levels, function(p) {
        detected[i, pair == p & control_index][[1L]]
      }, logical(1))
      b <- sum(case_detect & !control_detect)
      c <- sum(!case_detect & control_detect)
      discordant <- b + c
      result$statistic[[i]] <- b - c
      result$P.Value[[i]] <- if (discordant == 0L) 1 else {
        stats::binom.test(b, discordant, p = 0.5, alternative = "two.sided")$p.value
      }
      result$odds_ratio[[i]] <- if (c == 0L) {
        if (b == 0L) NA_real_ else Inf
      } else {
        b / c
      }
    }
  }
  result$adj.P.Val <- stats::p.adjust(result$P.Value, method = "BH")
  result$direction <- ifelse(
    result$detection_rate_difference > 0, "Case_more_detected",
    ifelse(result$detection_rate_difference < 0, "Control_more_detected", "No_difference")
  )
  result
}

## ------------------------------------------------------------- model dispatch
exp_model <- NULL
exp_visual <- NULL
feature_metadata <- NULL
detection_results <- data.frame()
voom_weights <- NULL
model_data <- list()
normalization_performed <- "none"
transformation_performed <- transform_todo
model_method <- NULL
statistic_type <- NULL
effect_size_type <- "model_log2_coefficient"
ebayes_method <- NULL
visual_imputation_used <- "none"
model_imputation_used <- "none"

is_count_profile <- analysis_profile %in% c("rnaseq_raw_counts", "rnaseq_fractional_counts", "protein_spectral_counts")
if (is_count_profile) {
  if (anyNA(exp_work)) stopf("%s requires a complete count matrix without NA.", analysis_profile)
  if (any(exp_work < 0)) stopf("%s contains negative counts.", analysis_profile)
  noninteger <- abs(exp_work - round(exp_work)) > 1e-6
  if (analysis_profile != "rnaseq_fractional_counts" && any(noninteger)) {
    if (allow_rounding) {
      warnf("Rounding non-integer counts because allow_rounding=TRUE was explicitly supplied.")
      exp_work <- round(exp_work)
    } else {
      stopf("%s contains non-integer values; verify Table S1 Analysis Profile.", analysis_profile)
    }
  }
  dge <- edgeR::DGEList(counts = exp_work, group = group_list)
  keep <- edgeR::filterByExpr(dge, design = design)
  log_removed(rownames(exp_work)[!keep], "filterByExpr")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  if (!nrow(dge$counts)) stop("No count features remain after filterByExpr.", call. = FALSE)
  count_norm_method <- if (
    analysis_profile == "rnaseq_raw_counts" && normalization_todo == "tmmwsp"
  ) "TMMwsp" else "TMM"
  dge <- edgeR::calcNormFactors(dge, method = count_norm_method)
  normalization_performed <- count_norm_method

  if (analysis_profile %in% c("rnaseq_raw_counts", "rnaseq_fractional_counts")) {
    use_quality_weights <- normalization_todo == "tmmwsp"
    voom_object <- if (use_quality_weights) {
      limma::voomWithQualityWeights(dge, design = design, plot = FALSE)
    } else {
      limma::voom(dge, design = design, plot = FALSE)
    }
    fit <- limma::lmFit(voom_object, design)
    estimable <- is.finite(fit$coefficients[, group_coef]) & fit$df.residual >= args$min_residual_df
    if (!all(estimable)) {
      log_removed(rownames(fit$coefficients)[!estimable], "non_estimable_or_low_residual_df")
      fit <- fit[estimable, ]
      voom_object$E <- voom_object$E[estimable, , drop = FALSE]
      voom_object$weights <- voom_object$weights[estimable, , drop = FALSE]
    }
    if (!nrow(fit$coefficients)) stop("No RNA-seq features are estimable.", call. = FALSE)
    fit <- limma::eBayes(fit, robust = TRUE)
    ebayes_method <- "robust"
    table <- limma::topTable(fit, coef = group_coef, number = Inf, sort.by = "none", adjust.method = "BH")
    exp_model <- voom_object$E
    voom_weights <- voom_object$weights
    exp_visual <- exp_model
    model_method <- if (use_quality_weights) {
      "edgeR_TMMwsp_voomWithQualityWeights_limma"
    } else {
      "edgeR_TMM_voom_limma"
    }
    statistic_type <- "moderated_t"
    deg_core <- data.frame(
      feature_id = rownames(table), logFC = table$logFC,
      AveExpr = table$AveExpr, statistic = table$t,
      P.Value = table$P.Value, adj.P.Val = table$adj.P.Val,
      stringsAsFactors = FALSE
    )
    feature_metadata <- data.frame(
      feature_id = rownames(exp_model), source_features = rownames(exp_model),
      stringsAsFactors = FALSE
    )
    model_data <- list(
      counts = dge$counts, norm_factors = dge$samples$norm.factors,
      effective_library_sizes = dge$samples$lib.size * dge$samples$norm.factors
    )
  } else {
    dge <- edgeR::estimateDisp(dge, design = design, robust = TRUE)
    qlfit <- edgeR::glmQLFit(dge, design = design, robust = TRUE)
    test <- edgeR::glmQLFTest(qlfit, coef = match(group_coef, colnames(design)))
    table <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
    exp_model <- edgeR::cpm(dge, log = TRUE, prior.count = 0.5)
    exp_visual <- exp_model
    model_method <- "edgeR_TMM_QLF"
    statistic_type <- "QLF_F"
    ebayes_method <- "not_applicable"
    deg_core <- data.frame(
      feature_id = rownames(table), logFC = table$logFC,
      AveExpr = table$logCPM, statistic = table$F,
      P.Value = table$PValue, adj.P.Val = table$FDR,
      stringsAsFactors = FALSE
    )
    feature_metadata <- data.frame(
      feature_id = rownames(exp_model), source_features = rownames(exp_model),
      stringsAsFactors = FALSE
    )
    model_data <- list(
      counts = dge$counts, norm_factors = dge$samples$norm.factors,
      dispersions = dge$tagwise.dispersion
    )
  }
} else {
  if (zero_is_missing) exp_work[exp_work == 0] <- NA_real_
  detection_input <- exp_work
  case_index <- group_list == "case"
  control_index <- group_list == "control"

  if (detection_shift_enabled && analysis_profile %in% c("proteomics_label_free", "metabolomics_abundance")) {
    detection_results <- compute_detection_results(
      detection_input, group_list,
      pair = if (paired_design) pair_values else NULL
    )
  }

  # Profile-specific feature filtering is calculated before imputation and
  # before transformation, preserving the original missingness mask.
  if (feature_filter_mode == "group_valid") {
    case_n <- sum(group_list == "case")
    control_n <- sum(group_list == "control")
    if (case_n < args$min_valid_n_per_group || control_n < args$min_valid_n_per_group) {
      stopf(
        "group_valid filtering requires at least %d samples in each group; observed case=%d, control=%d.",
        args$min_valid_n_per_group, case_n, control_n
      )
    }
    valid_case <- rowSums(!is.na(exp_work[, group_list == "case", drop = FALSE]))
    valid_control <- rowSums(!is.na(exp_work[, group_list == "control", drop = FALSE]))
    required_case <- max(args$min_valid_n_per_group, ceiling(case_n * args$min_valid_prop_per_group))
    required_control <- max(args$min_valid_n_per_group, ceiling(control_n * args$min_valid_prop_per_group))
    keep <- valid_case >= required_case & valid_control >= required_control
    log_removed(rownames(exp_work)[!keep], "group_valid_filter")
    exp_work <- exp_work[keep, , drop = FALSE]
  } else if (feature_filter_mode == "mean_detection") {
    if (input_scale != "linear") {
      stop("mean_detection filtering is valid only for non-negative linear abundance input.", call. = FALSE)
    }
    feature_mean <- rowMeans(exp_work, na.rm = TRUE)
    detection_prop <- rowMeans(exp_work > 0, na.rm = TRUE)
    keep <- is.finite(feature_mean) & is.finite(detection_prop) &
      feature_mean > args$min_mean & detection_prop >= args$min_detect_prop
    log_removed(rownames(exp_work)[!keep], "mean_detection_filter")
    exp_work <- exp_work[keep, , drop = FALSE]
  } else if (feature_filter_mode == "none") {
    keep <- rowSums(is.finite(exp_work)) > 0L
    log_removed(rownames(exp_work)[!keep], "all_missing")
    exp_work <- exp_work[keep, , drop = FALSE]
  } else {
    stopf("feature_filter_mode=%s is invalid for continuous profile %s.", feature_filter_mode, analysis_profile)
  }
  assert_matrix(exp_work, "continuous feature filtering")

  exp_transformed <- apply_continuous_transform(exp_work, transform_todo)
  if (normalization_todo == "none") {
    exp_normalized <- exp_transformed
  } else if (normalization_todo == "median_center") {
    exp_normalized <- median_center_log2(exp_transformed)
    normalization_performed <- "log2_column_median_centering"
  } else if (normalization_todo == "quantile") {
    exp_normalized <- limma::normalizeBetweenArrays(exp_transformed, method = "quantile")
    normalization_performed <- "quantile"
  } else if (normalization_todo == "irs") {
    if (analysis_profile != "proteomics_tmt") stop("IRS is supported only for proteomics_tmt.", call. = FALSE)
    exp_normalized <- apply_irs_log2(exp_transformed, plex_values)
    normalization_performed <- "protein_specific_IRS"
  } else {
    stopf("Unsupported normalization_todo: %s", normalization_todo)
  }
  if (normalization_todo == "none") normalization_performed <- "none"

  if (analysis_profile == "microarray") {
    collapse <- collapse_microarray_probes(exp_normalized)
    exp_normalized <- collapse$matrix
    feature_metadata <- collapse$metadata
    final_feature_id_type <- collapse$output_id_type
  } else {
    feature_metadata <- data.frame(
      feature_id = rownames(exp_normalized), source_features = rownames(exp_normalized),
      stringsAsFactors = FALSE
    )
    final_feature_id_type <- if (nzchar(args$output_feature_id_type)) {
      args$output_feature_id_type
    } else {
      args$feature_id_type
    }
  }

  if (analysis_profile == "metabolomics_abundance") {
    # Metabolite peak intensities are continuous and generally non-Gaussian.
    # Use a two-sided Wilcoxon rank-sum test (or paired signed-rank test when
    # explicitly configured), with BH correction. The effect is the median
    # log2 abundance difference (case minus control), which retains a clear
    # direction for volcano plots and signed GSEA-style ranking.
    safe_wilcox <- function(case_values, control_values, paired = FALSE) {
      case_values <- as.numeric(case_values)
      control_values <- as.numeric(control_values)
      if (paired) {
        ok <- is.finite(case_values) & is.finite(control_values)
        if (sum(ok) < 2L) return(list(statistic = NA_real_, p = NA_real_))
        test <- suppressWarnings(stats::wilcox.test(
          case_values[ok], control_values[ok], paired = TRUE,
          exact = FALSE, alternative = "two.sided"
        ))
      } else {
        case_values <- case_values[is.finite(case_values)]
        control_values <- control_values[is.finite(control_values)]
        if (length(case_values) < 2L || length(control_values) < 2L) {
          return(list(statistic = NA_real_, p = NA_real_))
        }
        test <- suppressWarnings(stats::wilcox.test(
          case_values, control_values, paired = FALSE,
          exact = FALSE, alternative = "two.sided"
        ))
      }
      list(statistic = unname(test$statistic), p = unname(test$p.value))
    }

    if (paired_design) {
      pair_levels <- levels(pair_values)
      case_by_pair <- vapply(pair_levels, function(pair) {
        idx <- pair_values == pair & group_list == "case"
        if (sum(idx) != 1L) NA_real_ else which(idx)[1L]
      }, integer(1))
      control_by_pair <- vapply(pair_levels, function(pair) {
        idx <- pair_values == pair & group_list == "control"
        if (sum(idx) != 1L) NA_real_ else which(idx)[1L]
      }, integer(1))
    }

    metabolite_rows <- lapply(seq_len(nrow(exp_normalized)), function(i) {
      case_values <- exp_normalized[i, case_index]
      control_values <- exp_normalized[i, control_index]
      if (paired_design) {
        test <- safe_wilcox(
          exp_normalized[i, case_by_pair], exp_normalized[i, control_by_pair], paired = TRUE
        )
        effect <- stats::median(
          exp_normalized[i, case_by_pair] - exp_normalized[i, control_by_pair], na.rm = TRUE
        )
      } else {
        test <- safe_wilcox(case_values, control_values, paired = FALSE)
        effect <- stats::median(case_values, na.rm = TRUE) - stats::median(control_values, na.rm = TRUE)
      }
      if (!is.finite(effect)) effect <- NA_real_
      signed_stat <- if (is.finite(test$p) && is.finite(effect) && effect != 0) {
        sign(effect) * stats::qnorm(max(test$p / 2, .Machine$double.xmin), lower.tail = FALSE)
      } else {
        NA_real_
      }
      data.frame(
        feature_id = rownames(exp_normalized)[i],
        logFC = effect,
        AveExpr = mean(c(case_values, control_values), na.rm = TRUE),
        statistic = signed_stat,
        W = test$statistic,
        P.Value = test$p,
        stringsAsFactors = FALSE
      )
    })
    table <- do.call(rbind, metabolite_rows)
    table$adj.P.Val <- stats::p.adjust(table$P.Value, method = "BH")
    rownames(table) <- table$feature_id
    exp_model <- exp_normalized
    model_method <- if (paired_design) "wilcoxon_signed_rank" else "wilcoxon_rank_sum"
    statistic_type <- "signed_wilcoxon_z"
    effect_size_type <- "median_log2_difference"
    ebayes_method <- "not_applicable"
    deg_core <- data.frame(
      feature_id = table$feature_id, logFC = table$logFC,
      AveExpr = table$AveExpr, statistic = table$statistic,
      P.Value = table$P.Value, adj.P.Val = table$adj.P.Val,
      stringsAsFactors = FALSE
    )
  } else {
    fit <- limma::lmFit(exp_normalized, design)
    estimable <- is.finite(fit$coefficients[, group_coef]) & fit$df.residual >= args$min_residual_df
    if (!all(estimable)) {
      log_removed(rownames(fit$coefficients)[!estimable], "non_estimable_or_low_residual_df")
      fit <- fit[estimable, ]
      exp_normalized <- exp_normalized[estimable, , drop = FALSE]
      feature_metadata <- feature_metadata[
        match(rownames(exp_normalized), feature_metadata$feature_id), , drop = FALSE
      ]
    }
    if (!nrow(fit$coefficients)) stop("No continuous features are estimable.", call. = FALSE)
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
    ebayes_method <- "trend+robust"
    table <- limma::topTable(fit, coef = group_coef, number = Inf, sort.by = "none", adjust.method = "BH")
    exp_model <- exp_normalized
    model_method <- "limma_filtered"
    statistic_type <- "moderated_t"
    deg_core <- data.frame(
      feature_id = rownames(table), logFC = table$logFC,
      AveExpr = table$AveExpr, statistic = table$t,
      P.Value = table$P.Value, adj.P.Val = table$adj.P.Val,
      stringsAsFactors = FALSE
    )
  }

  # Visual imputation is isolated from exp_model and never affects statistics.
  visual_method <- args$visual_imputation
  if (visual_method == "auto") {
    visual_method <- if (!anyNA(exp_model)) {
      "none"
    } else if (analysis_profile %in% c("proteomics_label_free", "metabolomics_abundance")) {
      "minprob"
    } else {
      "row_median"
    }
  }
  exp_visual <- exp_model
  if (anyNA(exp_visual)) {
    if (visual_method == "none") {
      stop("exp_visual still contains NA; choose a visual_imputation method.", call. = FALSE)
    } else if (visual_method == "row_median") {
      row_medians <- apply(exp_visual, 1L, stats::median, na.rm = TRUE)
      if (any(!is.finite(row_medians))) stop("Row-median visual imputation encountered an all-missing feature.", call. = FALSE)
      missing <- which(is.na(exp_visual), arr.ind = TRUE)
      exp_visual[missing] <- row_medians[missing[, "row"]]
    } else if (visual_method == "minprob") {
      required_package("imputeLCMD", "MinProb visual imputation")
      set.seed(args$seed)
      invisible(capture.output({
        exp_visual <- imputeLCMD::impute.MinProb(
          exp_visual, q = args$minprob_q, tune.sigma = args$minprob_tune_sigma
        )
      }, type = "output"))
    } else if (visual_method == "knn") {
      required_package("impute", "KNN visual imputation")
      exp_visual <- impute::impute.knn(exp_visual)$data
    } else if (visual_method == "zero") {
      exp_visual[is.na(exp_visual)] <- 0
    }
  }
  visual_imputation_used <- visual_method
}

assert_matrix(exp_model, "model matrix construction")
assert_matrix(exp_visual, "visual matrix construction")
if (anyNA(exp_visual)) stop("Internal error: exp_visual contains NA after visual preprocessing.", call. = FALSE)
if (!identical(colnames(exp_model), rownames(design))) stop("exp_model/design sample order mismatch.", call. = FALSE)
if (!identical(colnames(exp_visual), rownames(design))) stop("exp_visual/design sample order mismatch.", call. = FALSE)

# Batch removal is visualization-only. The model was fit to exp_model and the
# unmodified design above, preventing double correction.
visual_batch_correction <- "none"
visual_batch_1 <- if (!is.null(batch_values) && nlevels(batch_values) > 1L) batch_values else NULL
visual_batch_2 <- if (
  !is.null(plex_values) && nlevels(plex_values) > 1L && args$plex_column != args$batch_column
) plex_values else NULL
if (visual_remove_batch && (!is.null(visual_batch_1) || !is.null(visual_batch_2))) {
  preserve_group <- stats::model.matrix(~ group_list)
  if (!is.null(visual_batch_1) && !is.null(visual_batch_2)) {
    exp_visual <- limma::removeBatchEffect(
      exp_visual, batch = visual_batch_1, batch2 = visual_batch_2,
      design = preserve_group
    )
    visual_batch_correction <- paste0(
      "removeBatchEffect(", args$batch_column, ",", args$plex_column, ")"
    )
  } else {
    selected_batch <- if (!is.null(visual_batch_1)) visual_batch_1 else visual_batch_2
    selected_name <- if (!is.null(visual_batch_1)) args$batch_column else args$plex_column
    exp_visual <- limma::removeBatchEffect(
      exp_visual, batch = selected_batch, design = preserve_group
    )
    visual_batch_correction <- paste0("removeBatchEffect(", selected_name, ")")
  }
}

## ---------------------------------------------------------- standardize result
deg_core <- deg_core[match(rownames(exp_model), deg_core$feature_id), , drop = FALSE]
if (anyNA(deg_core$feature_id)) stop("Internal error: DEG rows do not align with exp_model.", call. = FALSE)

case_index <- group_list == "case"
control_index <- group_list == "control"
case_summary <- apply(exp_model[, case_index, drop = FALSE], 1L, stats::median, na.rm = TRUE)
control_summary <- apply(exp_model[, control_index, drop = FALSE], 1L, stats::median, na.rm = TRUE)
n_valid_case <- rowSums(is.finite(exp_model[, case_index, drop = FALSE]))
n_valid_control <- rowSums(is.finite(exp_model[, control_index, drop = FALSE]))

final_feature_id_type <- if (exists("final_feature_id_type")) final_feature_id_type else {
  if (nzchar(args$output_feature_id_type)) args$output_feature_id_type else args$feature_id_type
}

deg <- data.frame(
  feature_id = deg_core$feature_id,
  feature_id_type = final_feature_id_type,
  organism = args$organism,
  organism_taxid = as.character(args$organism_taxid),
  logFC = deg_core$logFC,
  effect_size_type = effect_size_type,
  AveExpr = deg_core$AveExpr,
  statistic = deg_core$statistic,
  statistic_type = statistic_type,
  P.Value = deg_core$P.Value,
  adj.P.Val = deg_core$adj.P.Val,
  n_valid_case = n_valid_case[deg_core$feature_id],
  n_valid_control = n_valid_control[deg_core$feature_id],
  case_median = case_summary[deg_core$feature_id],
  control_median = control_summary[deg_core$feature_id],
  analysis_profile = analysis_profile,
  contrast = paste0(args$case_code, "-", args$control_code),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

p_column <- if (is_use_padj) "adj.P.Val" else "P.Value"
deg$significant <- is.finite(deg[[p_column]]) &
  deg[[p_column]] <= args$pvalue & abs(deg$logFC) >= args$log2fc
deg$direction <- ifelse(
  deg$significant & deg$logFC > 0, "Up",
  ifelse(deg$significant & deg$logFC < 0, "Down", "NS")
)
deg$flag_extreme_logFC <- abs(deg$logFC) > 10
rownames(deg) <- deg$feature_id
deg <- deg[order(deg$P.Value, na.last = TRUE), , drop = FALSE]

feature_metadata <- feature_metadata[
  match(rownames(exp_model), feature_metadata$feature_id), , drop = FALSE
]
feature_metadata$feature_id_type <- final_feature_id_type
feature_metadata$organism <- args$organism
feature_metadata$organism_taxid <- as.character(args$organism_taxid)

removed_feature_log <- unique(removed_feature_log)
removed_features <- unique(removed_feature_log$feature_id)

## ------------------------------------------------------------- audit and save
output_csv <- if (!is.null(args$output_csv) && nzchar(args$output_csv)) {
  args$output_csv
} else {
  derive_output_path(args$output_rda, "_significant.csv")
}
output_all_csv <- if (!is.null(args$output_all_csv) && nzchar(args$output_all_csv)) {
  args$output_all_csv
} else {
  derive_output_path(args$output_rda, "_all.csv")
}
output_detection_csv <- if (!is.null(args$output_detection_csv) && nzchar(args$output_detection_csv)) {
  args$output_detection_csv
} else {
  derive_output_path(args$output_rda, "_detection.csv")
}

preprocessing <- list(
  schema_version = SCHEMA_VERSION,
  script_version = SCRIPT_VERSION,
  analysis_profile = analysis_profile,
  supported_comparison = "two_group_cancer_vs_non_cancer",
  contrast = paste0(args$case_code, "_minus_", args$control_code),
  effect_definition = "case_minus_control",
  grouping_source = grouping_source,
  group_labels = group_labels,
  design_formula = paste(deparse(design_formula), collapse = ""),
  group_coefficient = group_coef,
  paired_design = paired_design,
  pair_column = if (paired_design) args$pair_column else NULL,
  batch_column = if (nzchar(args$batch_column)) args$batch_column else NULL,
  plex_column = if (nzchar(args$plex_column)) args$plex_column else NULL,
  quantification_method = args$quantification_method,
  normalized_quantification_key = quant_key,
  omics_layer = args$omics_layer,
  molecule_type = args$molecule_type,
  detection_method = args$detection_method,
  feature_id_type_input = args$feature_id_type,
  feature_id_type_output = final_feature_id_type,
  organism = args$organism,
  organism_taxid = as.character(args$organism_taxid),
  input_scale_requested = args$input_scale,
  input_scale_resolved = input_scale,
  input_scale_source = scale_resolution$source,
  normalization_applied = args$normalization_applied,
  normalization_todo = normalization_todo,
  normalization_performed = normalization_performed,
  transform_applied = args$transform_applied,
  transform_todo = transform_todo,
  transformation_performed = transformation_performed,
  missing_value_encoding = args$missing_value_encoding,
  zero_is_missing = zero_is_missing,
  feature_filter_mode = feature_filter_mode,
  model_method = model_method,
  model_imputation = model_imputation_used,
  visual_imputation = visual_imputation_used,
  visual_batch_correction = visual_batch_correction,
  ebayes_method = ebayes_method,
  statistic_type = statistic_type,
  effect_size_type = effect_size_type,
  legacy_exp_alias = "exp_visual",
  # Significance thresholds classify a particular export request and are not
  # part of the immutable model cache. The requested values are retained
  # separately for auditability, without making the RDA threshold-dependent.
  thresholds = list(
    p_column = "not_applied_to_model", pvalue = NA_real_, log2fc = NA_real_
  ),
  export_thresholds = "per-request CSV only; not stored in the model RDA",
  parameters = list(
    sample_max_missing_prop = args$sample_max_missing_prop,
    min_valid_prop_per_group = args$min_valid_prop_per_group,
    min_valid_n_per_group = args$min_valid_n_per_group,
    min_residual_df = args$min_residual_df,
    min_mean = args$min_mean,
    min_detect_prop = args$min_detect_prop,
    legacy_max_na = if (args$max_na >= 0) args$max_na else NULL,
    minprob_q = args$minprob_q,
    minprob_tune_sigma = args$minprob_tune_sigma,
    detection_shift = detection_shift_enabled,
    seed = args$seed
  ),
  input_file = normalizePath(args$counts_file, mustWork = TRUE),
  input_file_sha256 = input_sha256,
  group_file = if (!is.null(args$group_file) && nzchar(args$group_file)) {
    normalizePath(args$group_file, mustWork = TRUE)
  } else {
    NULL
  },
  group_file_sha256 = if (!is.null(args$group_file) && nzchar(args$group_file)) {
    digest::digest(file = args$group_file, algo = "sha256")
  } else {
    NULL
  },
  sample_counts = list(
    input = length(all_samples), retained = ncol(exp_model), removed = length(removed_samples),
    case = unname(sum(group_list == "case")), control = unname(sum(group_list == "control"))
  ),
  feature_counts = list(
    input = length(all_features), modelled = nrow(exp_model),
    removed = length(removed_features), detection_tested = nrow(detection_results)
  ),
  output_files = list(
    rda = args$output_rda,
    # CSV paths belong to an individual request and are intentionally not
    # persisted in the shared model cache.
    significant_csv = NULL,
    all_csv = NULL,
    detection_csv = NULL
  ),
  session = utils::sessionInfo()
)

qc <- list(
  input_dimensions = dim(exp_raw),
  model_dimensions = dim(exp_model),
  visual_dimensions = dim(exp_visual),
  input_missing_count = sum(is.na(exp_raw)),
  model_missing_count = sum(is.na(exp_model)),
  visual_missing_count = sum(is.na(exp_visual)),
  sample_missing_fraction = sample_metadata_all$sample_missing_fraction,
  removed_feature_log = removed_feature_log,
  design_rank = qr(design)$rank,
  design_columns = colnames(design),
  model_and_design_aligned = identical(colnames(exp_model), rownames(design))
)

# Keep the RDA independent of one user's significance thresholds. Downstream
# consumers must recalculate significance from P.Value/adj.P.Val and logFC.
# The request-specific classified table remains in `deg_export` for CSV output.
deg_export <- deg
deg$significant <- NA
deg$direction <- ifelse(deg$logFC > 0, "Up", ifelse(deg$logFC < 0, "Down", "NS"))

# Backward-compatible alias. New PCA/heatmap code should explicitly use
# exp_visual; ROC should use exp_model and feature-wise complete cases.
exp <- exp_visual

ensure_parent_directory(args$output_rda)
temporary_rda <- tempfile(
  pattern = paste0(basename(args$output_rda), ".tmp-"),
  tmpdir = dirname(normalizePath(args$output_rda, mustWork = FALSE))
)
on.exit(if (file.exists(temporary_rda)) unlink(temporary_rda), add = TRUE)
save(
  exp_raw, exp_model, exp_visual, exp, deg, detection_results,
  feature_metadata, sample_metadata, sample_metadata_all,
  group_list, group_labels, design, preprocessing, qc,
  removed_samples, removed_features, voom_weights, model_data,
  file = temporary_rda
)
if (!file.rename(temporary_rda, args$output_rda)) {
  stop("Could not atomically publish RDA output: ", args$output_rda, call. = FALSE)
}

if (is_output_csv) {
  ensure_parent_directory(output_all_csv)
  utils::write.csv(deg_export, output_all_csv, row.names = FALSE, quote = TRUE)
  significant_deg <- deg_export[deg_export$significant %in% TRUE, , drop = FALSE]
  ensure_parent_directory(output_csv)
  utils::write.csv(significant_deg, output_csv, row.names = FALSE, quote = TRUE)
  if (nrow(detection_results)) {
    ensure_parent_directory(output_detection_csv)
    utils::write.csv(detection_results, output_detection_csv, row.names = FALSE, quote = TRUE)
  }
}

message(
  "Differential analysis completed. Profile=", analysis_profile,
  "; contrast=", args$case_code, "-", args$control_code,
  "; samples=", ncol(exp_model),
  "; modelled features=", nrow(exp_model),
  "; significant export features=", sum(deg_export$significant %in% TRUE),
  "; schema=", SCHEMA_VERSION, "."
)

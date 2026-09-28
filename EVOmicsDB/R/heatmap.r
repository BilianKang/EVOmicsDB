#!/usr/bin/env Rscript

# EVOmicsDB heatmap
# Uses the dedicated visualization matrix saved as `exp_visual`.
# Any configured visualization-only imputation occurs upstream in differential_analysis.r.

suppressPackageStartupMessages({
  library(argparse)
  library(ComplexHeatmap)
  library(circlize)
})

## ------------------------------------------------------------------ helpers
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

validate_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop("--", name, " must be one finite value.")
  }
}

validate_positive <- function(x, name) {
  validate_scalar(x, name)
  if (x <= 0) {
    stop("--", name, " must be greater than 0.")
  }
}

validate_probability <- function(x, name) {
  validate_scalar(x, name)
  if (x <= 0 || x > 1) {
    stop("--", name, " must be greater than 0 and no greater than 1.")
  }
}

validate_color <- function(x, name) {
  ok <- tryCatch({
    grDevices::col2rgb(x)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    stop("--", name, " is not a valid R color: ", x)
  }
}

make_group <- function(sample_names) {
  prefix <- toupper(substr(sample_names, 1L, 1L))
  group <- ifelse(
    prefix == "N",
    "Normal",
    ifelse(prefix == "C", "Tumor", NA_character_)
  )
  if (anyNA(group)) {
    stop(
      "Sample columns must start with N (control) or C (cancer). Offending columns: ",
      paste(sample_names[is.na(group)], collapse = ", ")
    )
  }
  factor(group, levels = c("Normal", "Tumor"))
}

write_empty_outputs <- function(message_text, args) {
  ensure_parent_directory(args$output_png)
  grDevices::png(
    args$output_png,
    width = args$width,
    height = args$height,
    res = 300,
    bg = "white"
  )
  graphics::plot.new()
  graphics::text(0.5, 0.5, message_text, cex = 1.2)
  grDevices::dev.off()

  if (!is.null(args$output_pdf) && nzchar(trimws(args$output_pdf))) {
    ensure_parent_directory(args$output_pdf)
    grDevices::pdf(
      args$output_pdf,
      width = args$width / 300,
      height = args$height / 300
    )
    graphics::plot.new()
    graphics::text(0.5, 0.5, message_text, cex = 1.2)
    grDevices::dev.off()
  }

  if (!is.null(args$output_csv) && nzchar(trimws(args$output_csv))) {
    ensure_parent_directory(args$output_csv)
    utils::write.csv(
      data.frame(feature_id = character()),
      args$output_csv,
      row.names = FALSE,
      quote = TRUE
    )
  }
  invisible(NULL)
}

## ------------------------------------------------------------------ arguments
parser <- ArgumentParser(
  description = "Generate an EVOmicsDB heatmap from differential-analysis results."
)
parser$add_argument("--input_rda", type = "character", required = TRUE)
parser$add_argument("--output_png", type = "character", required = TRUE)
parser$add_argument("--output_csv", type = "character", default = NULL)
parser$add_argument("--output_pdf", type = "character", default = NULL)

# Differential-analysis thresholds
parser$add_argument("--is_use_padj", type = "logical", default = TRUE)
parser$add_argument("--log2fc_threshold", type = "numeric", default = 1)
parser$add_argument("--pvalue_threshold", type = "numeric", default = 0.05)

# Key parameters shown by the web UI
parser$add_argument("--top_n", type = "integer", default = 50)
parser$add_argument("--title_text", type = "character", default = "Heatmap")
parser$add_argument("--col_label_size", type = "numeric", default = 7)
parser$add_argument("--row_label_size", type = "numeric", default = 7)
parser$add_argument("--is_col_clustering", type = "logical", default = TRUE)
parser$add_argument("--is_row_clustering", type = "logical", default = TRUE)
parser$add_argument("--is_show_legend", type = "logical", default = TRUE)
parser$add_argument("--legend_title", type = "character", default = "Expression")
parser$add_argument("--normal_fill_color", type = "character", default = "#3B4CC0")
parser$add_argument("--tumor_fill_color", type = "character", default = "#B40426")
parser$add_argument("--color_scheme_low", type = "character", default = "navy")
parser$add_argument("--color_scheme_mid", type = "character", default = "white")
parser$add_argument("--color_scheme_high", type = "character", default = "red")
# Optional symmetric color cap for row Z-scores. Omission keeps the data-driven range.
parser$add_argument("--color_limit", type = "numeric", default = NULL)
parser$add_argument("--width", type = "numeric", default = 2500)
parser$add_argument("--height", type = "numeric", default = 1350)
parser$add_argument("--font_size", type = "numeric", default = 7)
parser$add_argument("--angle_col", type = "numeric", default = 45)
parser$add_argument("--show_column_names", type = "logical", default = TRUE,
                    help = "Hide dense sample labels for publication display; sample IDs remain in CSV.")

args <- parser$parse_args()

## ------------------------------------------------------------------ validation
if (!file.exists(args$input_rda)) {
  stop("Input RDA does not exist: ", args$input_rda)
}
if (!nzchar(trimws(args$output_png))) {
  stop("--output_png cannot be empty.")
}
if (length(args$top_n) != 1L || is.na(args$top_n) || args$top_n < 1L) {
  stop("--top_n must be a positive integer.")
}
validate_scalar(args$log2fc_threshold, "log2fc_threshold")
if (args$log2fc_threshold < 0) {
  stop("--log2fc_threshold must be non-negative.")
}
validate_probability(args$pvalue_threshold, "pvalue_threshold")
validate_positive(args$col_label_size, "col_label_size")
validate_positive(args$row_label_size, "row_label_size")
validate_positive(args$width, "width")
validate_positive(args$height, "height")
validate_positive(args$font_size, "font_size")
color_limit <- args$color_limit
if (is.null(color_limit)) {
  color_limit <- NA_real_
}
if (length(color_limit) != 1L ||
    (!is.na(color_limit) &&
      (!is.finite(color_limit) || color_limit <= 0))) {
  stop("--color_limit must be omitted or set to one finite value greater than 0.")
}

validate_scalar(args$angle_col, "angle_col")
angle_col <- args$angle_col %% 360

validate_color(args$normal_fill_color, "normal_fill_color")
validate_color(args$tumor_fill_color, "tumor_fill_color")
validate_color(args$color_scheme_low, "color_scheme_low")
validate_color(args$color_scheme_mid, "color_scheme_mid")
validate_color(args$color_scheme_high, "color_scheme_high")

legend_title <- trimws(args$legend_title)
if (!nzchar(legend_title)) {
  legend_title <- "Expression"
}

## ------------------------------------------------------------------ load input
command <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command, value = TRUE)
script_directory <- if (length(file_argument)) {
  # macOS Rscript may encode spaces in an absolute --file path as ~+~.
  dirname(normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_argument[[1L]]))))
} else {
  getwd()
}
source(file.path(script_directory, "evomics_rda_utils.R"))
analysis_data <- load_evomics_visual_data(args$input_rda)
exp <- analysis_data$exp_visual
deg <- analysis_data$deg
sample_group <- factor(
  ifelse(analysis_data$group == "Control", "Normal", "Tumor"),
  levels = c("Normal", "Tumor")
)
names(sample_group) <- colnames(exp)

if (nrow(exp) == 0L || ncol(exp) == 0L) {
  stop("The expression matrix is empty.")
}
if (is.null(rownames(exp)) || anyNA(rownames(exp)) || any(rownames(exp) == "")) {
  stop("The expression matrix must have non-empty feature row names.")
}
if (is.null(colnames(exp)) || anyNA(colnames(exp)) || any(colnames(exp) == "")) {
  stop("The expression matrix must have non-empty sample column names.")
}
if (anyDuplicated(rownames(exp))) {
  stop("The expression matrix contains duplicated feature identifiers.")
}
if (anyDuplicated(colnames(exp))) {
  stop("The expression matrix contains duplicated sample names.")
}
if (nrow(deg) == 0L || is.null(rownames(deg))) {
  stop("The DEG table is empty or lacks feature row names.")
}
if (anyNA(rownames(deg)) || any(rownames(deg) == "") || anyDuplicated(rownames(deg))) {
  stop("The DEG table must have non-empty, unique feature row names.")
}

p_column <- if (isTRUE(args$is_use_padj)) "adj.P.Val" else "P.Value"
required_columns <- c("logFC", p_column)
missing_columns <- setdiff(required_columns, colnames(deg))
if (length(missing_columns) > 0L) {
  stop("The DEG table is missing required columns: ", paste(missing_columns, collapse = ", "))
}

deg$logFC <- suppressWarnings(as.numeric(deg$logFC))
deg[[p_column]] <- suppressWarnings(as.numeric(deg[[p_column]]))
invalid_p <- !is.na(deg[[p_column]]) &
  (!is.finite(deg[[p_column]]) |
     deg[[p_column]] < 0 |
     deg[[p_column]] > 1)
if (any(invalid_p)) {
  stop("The selected P-value column contains values outside [0, 1] or non-finite values.")
}

## ------------------------------------------------------------------ select plottable significant features
is_significant <- !is.na(deg[[p_column]]) &
  deg[[p_column]] <= args$pvalue_threshold &
  is.finite(deg$logFC) &
  abs(deg$logFC) >= args$log2fc_threshold

significant_deg <- deg[is_significant, , drop = FALSE]
if (nrow(significant_deg) > 0L) {
  significant_deg <- significant_deg[
    order(significant_deg[[p_column]], -abs(significant_deg$logFC)),
    ,
    drop = FALSE
  ]
}

candidate_ids <- rownames(significant_deg)
present_ids <- candidate_ids[candidate_ids %in% rownames(exp)]
missing_from_exp <- setdiff(candidate_ids, present_ids)

if (length(present_ids) > 0L) {
  candidate_matrix <- exp[present_ids, , drop = FALSE]
  finite_rows <- apply(candidate_matrix, 1L, function(x) all(is.finite(x)))
  finite_ids <- rownames(candidate_matrix)[finite_rows]
  excluded_nonfinite <- rownames(candidate_matrix)[!finite_rows]

  if (length(finite_ids) > 0L) {
    finite_matrix <- candidate_matrix[finite_ids, , drop = FALSE]
    row_sd <- apply(finite_matrix, 1L, stats::sd)
    variable_rows <- is.finite(row_sd) & row_sd > 0
    plottable_ids <- rownames(finite_matrix)[variable_rows]
    excluded_zero_variance <- rownames(finite_matrix)[!variable_rows]
  } else {
    plottable_ids <- character()
    excluded_zero_variance <- character()
  }
} else {
  plottable_ids <- character()
  excluded_nonfinite <- character()
  excluded_zero_variance <- character()
}

selected_ids <- utils::head(plottable_ids, args$top_n)

if (length(selected_ids) == 0L) {
  message_text <- if (nrow(significant_deg) == 0L) {
    "No significant features to display"
  } else {
    "No complete, variable significant features to display"
  }
  warning(message_text)
  write_empty_outputs(message_text, args)
  message(
    "Heatmap completed with no plotted features. Missing from exp: ",
    length(missing_from_exp), "; non-finite: ", length(excluded_nonfinite),
    "; zero variance: ", length(excluded_zero_variance), "."
  )
  quit(save = "no", status = 0)
}

expression_matrix <- exp[selected_ids, , drop = FALSE]

# Compute the exact row z-scores plotted. NA/Inf and zero-SD rows were removed above.
scaled_matrix <- t(scale(t(expression_matrix)))
if (any(!is.finite(scaled_matrix))) {
  stop("Internal error: non-finite values remain after row z-score scaling.")
}

## ------------------------------------------------------------------ annotation and plotting
group_colors <- c(
  Normal = args$normal_fill_color,
  Tumor = args$tumor_fill_color
)

cluster_rows <- isTRUE(args$is_row_clustering) && nrow(scaled_matrix) >= 2L
cluster_cols <- isTRUE(args$is_col_clustering) && ncol(scaled_matrix) >= 2L

data_absolute_limit <- max(abs(scaled_matrix), na.rm = TRUE)
if (!is.finite(data_absolute_limit) || data_absolute_limit <= 0) {
  stop("Unable to construct a symmetric color scale from the selected features.")
}
absolute_limit <- if (is.na(color_limit)) {
  data_absolute_limit
} else {
  min(data_absolute_limit, color_limit)
}
heatmap_color_function <- circlize::colorRamp2(
  c(-absolute_limit, 0, absolute_limit),
  c(
    args$color_scheme_low,
    args$color_scheme_mid,
    args$color_scheme_high
  ),
  space = "LAB"
)
legend_ticks <- seq(-absolute_limit, absolute_limit, length.out = 5L)
legend_tick_labels <- format(
  signif(legend_ticks, digits = 3L),
  trim = TRUE,
  scientific = FALSE
)

build_heatmap <- function() {
  column_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Group = sample_group,
    col = list(Group = group_colors),
    show_legend = isTRUE(args$is_show_legend),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = args$font_size),
    annotation_legend_param = list(
      Group = list(
        title = "Group",
        title_gp = grid::gpar(
          fontsize = args$font_size,
          fontface = "bold"
        ),
        labels_gp = grid::gpar(fontsize = args$font_size)
      )
    )
  )

  ComplexHeatmap::Heatmap(
    matrix = scaled_matrix,
    name = legend_title,
    col = heatmap_color_function,
    na_col = "#F0F0F0",
    top_annotation = column_annotation,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    row_dend_reorder = FALSE,
    column_dend_reorder = FALSE,
    show_row_names = TRUE,
    show_column_names = isTRUE(args$show_column_names),
    row_names_gp = grid::gpar(fontsize = args$row_label_size),
    column_names_gp = grid::gpar(fontsize = args$col_label_size),
    column_names_rot = angle_col,
    column_title = args$title_text,
    column_title_gp = grid::gpar(
      fontsize = args$font_size + 3,
      fontface = "bold"
    ),
    show_heatmap_legend = isTRUE(args$is_show_legend),
    heatmap_legend_param = list(
      title = legend_title,
      at = legend_ticks,
      labels = legend_tick_labels,
      title_gp = grid::gpar(
        fontsize = args$font_size,
        fontface = "bold"
      ),
      labels_gp = grid::gpar(fontsize = args$font_size)
    ),
    border = FALSE,
    use_raster = FALSE
  )
}

draw_heatmap <- function() {
  ComplexHeatmap::draw(
    build_heatmap(),
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE
  )
}

ensure_parent_directory(args$output_png)
grDevices::png(
  args$output_png,
  width = args$width,
  height = args$height,
  res = 300,
  bg = "white"
)
draw_heatmap()
grDevices::dev.off()

if (!is.null(args$output_pdf) && nzchar(trimws(args$output_pdf))) {
  ensure_parent_directory(args$output_pdf)
  grDevices::pdf(
    args$output_pdf,
    width = args$width / 300,
    height = args$height / 300,
    bg = "white"
  )
  draw_heatmap()
  grDevices::dev.off()
}

if (!is.null(args$output_csv) && nzchar(trimws(args$output_csv))) {
  ensure_parent_directory(args$output_csv)
  csv_data <- data.frame(
    feature_id = rownames(scaled_matrix),
    scaled_matrix,
    check.names = FALSE
  )
  utils::write.csv(
    csv_data,
    args$output_csv,
    row.names = FALSE,
    quote = TRUE
  )
}

message(
  "Heatmap completed: ", nrow(scaled_matrix), " plotted features; ",
  length(missing_from_exp), " significant features missing from exp; ",
  length(excluded_nonfinite), " excluded for NA/Inf; ",
  length(excluded_zero_variance), " excluded for zero variance; ",
  if (cluster_rows) "row clustering enabled" else "row clustering disabled",
  "; ",
  if (cluster_cols) "column clustering enabled" else "column clustering disabled",
  "; color range = ±", format(absolute_limit, digits = 4),
  if (is.na(color_limit)) " (automatic)." else " (manual cap)."
)

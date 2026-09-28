#!/usr/bin/env Rscript

# EVOmicsDB volcano plot
# The differential-analysis convention is logFC = C - N:
#   Up   = higher in cancer
#   Down = lower in cancer

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
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

validate_nonnegative <- function(x, name, strictly_positive = FALSE) {
  validate_scalar(x, name)
  valid <- if (strictly_positive) x > 0 else x >= 0
  if (!valid) {
    stop(
      "--", name, " must be ",
      if (strictly_positive) "greater than 0." else "non-negative."
    )
  }
}

validate_probability <- function(x, name, strictly_positive = FALSE) {
  validate_scalar(x, name)
  lower_ok <- if (strictly_positive) x > 0 else x >= 0
  if (!lower_ok || x > 1) {
    stop(
      "--", name, " must be ",
      if (strictly_positive) "greater than 0 and no greater than 1."
      else "between 0 and 1."
    )
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

## ------------------------------------------------------------------ arguments
parser <- ArgumentParser(
  description = "Generate an EVOmicsDB volcano plot from a limma differential-analysis RDA."
)
parser$add_argument(
  "--input_rda",
  type = "character",
  required = TRUE,
  help = "RDA containing a 'deg' object."
)
parser$add_argument("--output_png", type = "character", required = TRUE)
parser$add_argument("--output_csv", type = "character", default = NULL)
parser$add_argument("--output_pdf", type = "character", default = NULL)

# Differential-analysis thresholds
parser$add_argument("--is_use_padj", type = "logical", default = TRUE)
parser$add_argument("--log2fc_threshold", type = "numeric", default = 1)
parser$add_argument("--pvalue_threshold", type = "numeric", default = 0.05)

# Axes
parser$add_argument("--is_x_axis_auto", type = "logical", default = TRUE)
parser$add_argument("--x_axis_low", type = "numeric", default = -1)
parser$add_argument("--x_axis_high", type = "numeric", default = 1)
parser$add_argument("--is_y_axis_auto", type = "logical", default = TRUE)
parser$add_argument("--y_axis_low", type = "numeric", default = 0)
parser$add_argument("--y_axis_high", type = "numeric", default = 10)

# Titles
parser$add_argument("--main_title", type = "character", default = "Volcano Plot")
parser$add_argument("--x_axis_title", type = "character", default = "Log2 Fold Change")
parser$add_argument("--y_axis_title", type = "character", default = "-Log10(p-value)")

# Dots
parser$add_argument("--up_fill_color", type = "character", default = "#E64B35")
parser$add_argument("--down_fill_color", type = "character", default = "#4DBBD5")
parser$add_argument("--ns_fill_color", type = "character", default = "#BEBEBE")
parser$add_argument("--shape", type = "numeric", default = 21)
parser$add_argument("--size", type = "numeric", default = 1)
parser$add_argument("--opacity", type = "numeric", default = 0.8)

# Legend and labels
parser$add_argument("--legend_title_size", type = "numeric", default = 8)
parser$add_argument("--legend_text_size", type = "numeric", default = 7)
parser$add_argument("--label_top_n_significant", type = "integer", default = 5)
parser$add_argument("--legend_position", type = "character", default = "right")

# Style and export
parser$add_argument("--is_show_border", type = "logical", default = TRUE)
parser$add_argument("--font_size", type = "numeric", default = 7)
parser$add_argument("--width", type = "numeric", default = 800)
parser$add_argument("--height", type = "numeric", default = 600)

args <- parser$parse_args()

## ------------------------------------------------------------------ validation
if (!file.exists(args$input_rda)) {
  stop("Input RDA does not exist: ", args$input_rda)
}
if (!nzchar(trimws(args$output_png))) {
  stop("--output_png cannot be empty.")
}

validate_nonnegative(args$log2fc_threshold, "log2fc_threshold")
validate_probability(args$pvalue_threshold, "pvalue_threshold", strictly_positive = TRUE)
validate_probability(args$opacity, "opacity")
validate_nonnegative(args$size, "size", strictly_positive = TRUE)
validate_nonnegative(args$font_size, "font_size", strictly_positive = TRUE)
validate_nonnegative(args$legend_title_size, "legend_title_size", strictly_positive = TRUE)
validate_nonnegative(args$legend_text_size, "legend_text_size", strictly_positive = TRUE)
validate_nonnegative(args$width, "width", strictly_positive = TRUE)
validate_nonnegative(args$height, "height", strictly_positive = TRUE)

if (
  length(args$label_top_n_significant) != 1L ||
    is.na(args$label_top_n_significant) ||
    args$label_top_n_significant < 0L
) {
  stop("--label_top_n_significant must be a non-negative integer.")
}
if (
  length(args$shape) != 1L ||
    is.na(args$shape) ||
    !is.finite(args$shape) ||
    args$shape != round(args$shape) ||
    args$shape < 0 ||
    args$shape > 25
) {
  stop("--shape must be an integer from 0 to 25.")
}

args$legend_position <- tolower(trimws(args$legend_position))
if (!args$legend_position %in% c("right", "bottom", "left", "top", "inside", "none")) {
  stop("--legend_position must be right, bottom, left, top, inside, or none.")
}
legend_inside <- identical(args$legend_position, "inside")

validate_color(args$up_fill_color, "up_fill_color")
validate_color(args$down_fill_color, "down_fill_color")
validate_color(args$ns_fill_color, "ns_fill_color")

if (!isTRUE(args$is_x_axis_auto)) {
  validate_scalar(args$x_axis_low, "x_axis_low")
  validate_scalar(args$x_axis_high, "x_axis_high")
  if (args$x_axis_low >= args$x_axis_high) {
    stop("--x_axis_low must be smaller than --x_axis_high.")
  }
}
if (!isTRUE(args$is_y_axis_auto)) {
  validate_scalar(args$y_axis_low, "y_axis_low")
  validate_scalar(args$y_axis_high, "y_axis_high")
  if (args$y_axis_low >= args$y_axis_high) {
    stop("--y_axis_low must be smaller than --y_axis_high.")
  }
}

## ------------------------------------------------------------------ load and validate DEG data
rda_environment <- new.env(parent = emptyenv())
load(args$input_rda, envir = rda_environment)
if (!exists("deg", envir = rda_environment, inherits = FALSE)) {
  stop("The input RDA does not contain a 'deg' object.")
}

res_df <- as.data.frame(rda_environment$deg)
required_columns <- c("logFC", if (isTRUE(args$is_use_padj)) "adj.P.Val" else "P.Value")
missing_columns <- setdiff(required_columns, colnames(res_df))
if (length(missing_columns) > 0L) {
  stop("The DEG table is missing required columns: ", paste(missing_columns, collapse = ", "))
}
if (nrow(res_df) == 0L) {
  stop("The DEG table is empty.")
}

gene_names <- rownames(res_df)
if (is.null(gene_names) || anyNA(gene_names) || any(gene_names == "")) {
  gene_names <- paste0("Feature_", seq_len(nrow(res_df)))
}
res_df$gene <- gene_names

res_df$logFC <- suppressWarnings(as.numeric(res_df$logFC))
p_column <- if (isTRUE(args$is_use_padj)) "adj.P.Val" else "P.Value"
res_df$pvalue_for_plot <- suppressWarnings(as.numeric(res_df[[p_column]]))

invalid_p <- !is.na(res_df$pvalue_for_plot) &
  (!is.finite(res_df$pvalue_for_plot) |
     res_df$pvalue_for_plot < 0 |
     res_df$pvalue_for_plot > 1)
if (any(invalid_p)) {
  stop("The selected P-value column contains values outside [0, 1] or non-finite values.")
}
if (!any(is.finite(res_df$logFC))) {
  stop("The DEG table contains no finite logFC values.")
}
if (!any(!is.na(res_df$pvalue_for_plot))) {
  stop("The selected P-value column contains no usable values.")
}

# Preserve the selected P value and use a clamped copy only for -log10 plotting.
res_df$plot_pvalue <- pmax(
  res_df$pvalue_for_plot,
  .Machine$double.xmin,
  na.rm = FALSE
)
res_df$neg_log10_pvalue <- -log10(res_df$plot_pvalue)

## ------------------------------------------------------------------ classify and label
is_significant <- !is.na(res_df$pvalue_for_plot) &
  res_df$pvalue_for_plot <= args$pvalue_threshold
res_df$group <- "NS"
res_df$group[
  is_significant &
    is.finite(res_df$logFC) &
    res_df$logFC > 0 &
    res_df$logFC >= args$log2fc_threshold
] <- "Up"
res_df$group[
  is_significant &
    is.finite(res_df$logFC) &
    res_df$logFC < 0 &
    -res_df$logFC >= args$log2fc_threshold
] <- "Down"
res_df$group <- factor(res_df$group, levels = c("Up", "Down", "NS"))

top_labels <- res_df[0, , drop = FALSE]
if (args$label_top_n_significant > 0L) {
  top_labels <- res_df |>
    dplyr::filter(.data$group != "NS", is.finite(.data$neg_log10_pvalue)) |>
    dplyr::arrange(.data$pvalue_for_plot, dplyr::desc(abs(.data$logFC))) |>
    dplyr::slice_head(n = args$label_top_n_significant)
}

## ------------------------------------------------------------------ plot
group_colors <- c(
  Up = args$up_fill_color,
  Down = args$down_fill_color,
  NS = args$ns_fill_color
)

resolved_y_axis_title <- args$y_axis_title
default_p_titles <- c("-Log10(p-value)", "-log10(p-value)", "-Log10(P-value)", "-log10(P-value)")
if (isTRUE(args$is_use_padj) && resolved_y_axis_title %in% default_p_titles) {
  resolved_y_axis_title <- "-Log10(adjusted P-value)"
}

p <- ggplot(
  res_df,
  aes(
    x = .data$logFC,
    y = .data$neg_log10_pvalue,
    color = .data$group,
    fill = .data$group
  )
) +
  geom_point(
    alpha = args$opacity,
    size = args$size,
    shape = as.integer(args$shape),
    na.rm = TRUE
  ) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE, guide = "none") +
  geom_vline(
    xintercept = c(-args$log2fc_threshold, args$log2fc_threshold),
    linetype = "dashed",
    color = "black",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = -log10(args$pvalue_threshold),
    linetype = "dashed",
    color = "black",
    linewidth = 0.4
  ) +
  labs(
    title = args$main_title,
    x = args$x_axis_title,
    y = resolved_y_axis_title,
    color = "Regulation",
    fill = "Regulation"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = as.integer(args$shape),
        size = max(3, args$size),
        alpha = 1,
        fill = unname(group_colors)
      )
    )
  ) +
  theme_minimal(base_size = args$font_size) +
  theme(
    plot.title = element_text(
      size = args$font_size + 3,
      hjust = 0.5,
      face = "bold"
    ),
    axis.title = element_text(size = args$font_size + 1, face = "bold"),
    axis.text = element_text(size = args$font_size),
    legend.title = element_text(size = args$legend_title_size, face = "bold"),
    legend.text = element_text(size = args$legend_text_size),
    legend.position = if (legend_inside) "inside" else args$legend_position,
    legend.position.inside = c(0.98, 0.98),
    legend.justification = if (legend_inside) c(1, 1) else "center",
    legend.background = if (legend_inside) {
      element_rect(fill = scales::alpha("white", 0.82), colour = "grey75")
    } else {
      element_blank()
    },
    panel.border = if (isTRUE(args$is_show_border)) {
      element_rect(colour = "black", fill = NA, linewidth = 0.5)
    } else {
      element_blank()
    },
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_line(colour = "grey95")
  )

if (nrow(top_labels) > 0L) {
  p <- p + ggrepel::geom_text_repel(
    data = top_labels,
    aes(label = .data$gene),
    size = max(2.5, args$font_size / 3),
    show.legend = FALSE,
    seed = 123,
    box.padding = 0.35,
    point.padding = 0.2,
    max.overlaps = Inf,
    min.segment.length = 0
  )
}

x_limits <- if (isTRUE(args$is_x_axis_auto)) {
  NULL
} else {
  c(args$x_axis_low, args$x_axis_high)
}
y_limits <- if (isTRUE(args$is_y_axis_auto)) {
  NULL
} else {
  c(args$y_axis_low, args$y_axis_high)
}
if (!is.null(x_limits) || !is.null(y_limits)) {
  p <- p + coord_cartesian(
    xlim = x_limits,
    ylim = y_limits,
    expand = TRUE
  )
}

## ------------------------------------------------------------------ outputs
ensure_parent_directory(args$output_png)
ggsave(
  filename = args$output_png,
  plot = p,
  width = args$width,
  height = args$height,
  units = "px",
  dpi = 300,
  bg = "white"
)

if (!is.null(args$output_pdf) && nzchar(trimws(args$output_pdf))) {
  ensure_parent_directory(args$output_pdf)
  ggsave(
    filename = args$output_pdf,
    plot = p,
    width = args$width / 300,
    height = args$height / 300,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
}

if (!is.null(args$output_csv) && nzchar(trimws(args$output_csv))) {
  ensure_parent_directory(args$output_csv)
  utils::write.csv(
    res_df,
    file = args$output_csv,
    row.names = FALSE,
    quote = TRUE
  )
}

message(
  "Volcano plot completed: ",
  sum(res_df$group == "Up", na.rm = TRUE), " Up, ",
  sum(res_df$group == "Down", na.rm = TRUE), " Down, ",
  sum(res_df$group == "NS", na.rm = TRUE), " NS."
)

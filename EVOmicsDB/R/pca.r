#!/usr/bin/env Rscript

# EVOmicsDB principal component analysis
# Uses the dedicated visualization matrix saved as `exp_visual` and the
# group_list stored by differential_analysis.r. Legacy ^N/^C grouping is only
# used when an old RDA lacks group_list.

required_packages <- c("argparse", "ggplot2", "FactoMineR")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them in the R environment used by the EVOmicsDB backend."
  )
}

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
  group <- factor(group, levels = c("Normal", "Tumor"))
  group_count <- table(group)
  if (any(group_count == 0L)) {
    stop("Both Normal (^N) and Tumor (^C) samples are required for PCA.")
  }
  group
}

make_confidence_ellipse <- function(data, level = 0.95, points = 200L) {
  if (nrow(data) < 3L) {
    return(NULL)
  }

  coordinates <- as.matrix(data[, c("PC1", "PC2"), drop = FALSE])
  covariance <- stats::cov(coordinates)
  if (any(!is.finite(covariance))) {
    return(NULL)
  }

  decomposition <- eigen(covariance, symmetric = TRUE)
  tolerance <- .Machine$double.eps * max(1, max(decomposition$values))
  if (any(decomposition$values <= tolerance)) {
    return(NULL)
  }

  theta <- seq(0, 2 * pi, length.out = points)
  unit_circle <- cbind(cos(theta), sin(theta))
  transform <- decomposition$vectors %*%
    diag(sqrt(decomposition$values), nrow = 2L)
  radius <- sqrt(stats::qchisq(level, df = 2L))
  ellipse <- unit_circle %*% t(transform) * radius
  ellipse <- sweep(ellipse, 2L, colMeans(coordinates), "+")

  data.frame(
    PC1 = ellipse[, 1L],
    PC2 = ellipse[, 2L],
    Group = data$Group[[1L]],
    stringsAsFactors = FALSE
  )
}

## ------------------------------------------------------------------ arguments
parser <- argparse::ArgumentParser(
  description = "Perform PCA on an EVOmicsDB differential-analysis expression matrix."
)
parser$add_argument("--input_rda", type = "character", required = TRUE)
parser$add_argument("--output_png", type = "character", required = TRUE)
parser$add_argument("--output_csv", type = "character", required = TRUE)
parser$add_argument("--output_pdf", type = "character", default = NULL)

# Key parameters shown by the web UI
parser$add_argument("--main_title", type = "character", default = "PCA Plot")
parser$add_argument("--x_axis_title", type = "character", default = "PC1")
parser$add_argument("--y_axis_title", type = "character", default = "PC2")
parser$add_argument("--axis_title_size", type = "numeric", default = 9)
parser$add_argument("--axis_text_size", type = "numeric", default = 7)
parser$add_argument("--is_show_legend", type = "logical", default = TRUE)
parser$add_argument("--legend_title", type = "character", default = "Expression")
parser$add_argument("--legend_title_size", type = "numeric", default = 9)
parser$add_argument("--legend_text_size", type = "numeric", default = 7)
parser$add_argument("--legend_position", type = "character", default = "right")
parser$add_argument("--normal_fill_color", type = "character", default = "#3B4CC0")
parser$add_argument("--tumor_fill_color", type = "character", default = "#B40426")
parser$add_argument("--shape", type = "integer", default = 21)
parser$add_argument("--size", type = "numeric", default = 3)
parser$add_argument("--opacity", type = "numeric", default = 0.8)
parser$add_argument("--is_show_border", type = "logical", default = TRUE)
parser$add_argument("--font_size", type = "numeric", default = 12)
parser$add_argument("--width", type = "numeric", default = 1800)
parser$add_argument("--height", type = "numeric", default = 1350)

args <- parser$parse_args()

## ------------------------------------------------------------------ validation
if (!file.exists(args$input_rda)) {
  stop("Input RDA does not exist: ", args$input_rda)
}
if (!nzchar(trimws(args$output_png))) {
  stop("--output_png cannot be empty.")
}
if (!nzchar(trimws(args$output_csv))) {
  stop("--output_csv cannot be empty.")
}

validate_positive(args$axis_title_size, "axis_title_size")
validate_positive(args$axis_text_size, "axis_text_size")
validate_positive(args$legend_title_size, "legend_title_size")
validate_positive(args$legend_text_size, "legend_text_size")
validate_positive(args$size, "size")
validate_positive(args$font_size, "font_size")
validate_positive(args$width, "width")
validate_positive(args$height, "height")

validate_scalar(args$opacity, "opacity")
if (args$opacity < 0 || args$opacity > 1) {
  stop("--opacity must be between 0 and 1.")
}
if (length(args$shape) != 1L ||
    is.na(args$shape) ||
    args$shape < 0L ||
    args$shape > 25L) {
  stop("--shape must be an integer from 0 to 25.")
}

legend_position <- tolower(trimws(args$legend_position))
supported_legend_positions <- c("right", "left", "top", "bottom", "inside", "none")
if (!legend_position %in% supported_legend_positions) {
  stop(
    "--legend_position must be one of: ",
    paste(supported_legend_positions, collapse = ", "),
    "."
  )
}
legend_inside <- identical(legend_position, "inside")

validate_color(args$normal_fill_color, "normal_fill_color")
validate_color(args$tumor_fill_color, "tumor_fill_color")

## ------------------------------------------------------------------ load and validate input
command <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command, value = TRUE)
script_directory <- if (length(file_argument)) {
  dirname(normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_argument[[1L]]))))
} else {
  getwd()
}
source(file.path(script_directory, "evomics_rda_utils.R"))
analysis_data <- load_evomics_visual_data(args$input_rda)
expression_matrix <- analysis_data$exp_visual

if (nrow(expression_matrix) == 0L || ncol(expression_matrix) == 0L) {
  stop("The expression matrix is empty.")
}
if (is.null(rownames(expression_matrix)) ||
    anyNA(rownames(expression_matrix)) ||
    any(rownames(expression_matrix) == "")) {
  stop("The expression matrix must have non-empty feature row names.")
}
if (is.null(colnames(expression_matrix)) ||
    anyNA(colnames(expression_matrix)) ||
    any(colnames(expression_matrix) == "")) {
  stop("The expression matrix must have non-empty sample column names.")
}
if (anyDuplicated(rownames(expression_matrix))) {
  stop("The expression matrix contains duplicated feature identifiers.")
}
if (anyDuplicated(colnames(expression_matrix))) {
  stop("The expression matrix contains duplicated sample names.")
}
if (ncol(expression_matrix) < 3L) {
  stop("PCA requires at least three samples to display PC1 and PC2.")
}

sample_group <- factor(
  ifelse(analysis_data$group == "Control", "Normal", "Tumor"),
  levels = c("Normal", "Tumor")
)
names(sample_group) <- colnames(expression_matrix)

# PCA cannot accept NA/Inf. Exclude incomplete features; do not silently impute.
finite_rows <- apply(
  expression_matrix,
  1L,
  function(x) all(is.finite(x))
)
excluded_nonfinite <- rownames(expression_matrix)[!finite_rows]
expression_matrix <- expression_matrix[finite_rows, , drop = FALSE]

if (nrow(expression_matrix) == 0L) {
  stop("No complete features remain after excluding rows containing NA/Inf.")
}

# Feature scaling requires positive variance.
feature_sd <- apply(expression_matrix, 1L, stats::sd)
variable_rows <- is.finite(feature_sd) & feature_sd > 0
excluded_zero_variance <- rownames(expression_matrix)[!variable_rows]
expression_matrix <- expression_matrix[variable_rows, , drop = FALSE]

if (nrow(expression_matrix) < 2L) {
  stop("At least two complete, variable features are required for PCA.")
}

## ------------------------------------------------------------------ PCA
# Retain the historical EVOmicsDB PCA engine and its axis orientation.
pca_fit <- FactoMineR::PCA(
  as.data.frame(t(expression_matrix), check.names = FALSE),
  scale.unit = TRUE,
  ncp = 2L,
  graph = FALSE
)

if (ncol(pca_fit$ind$coord) < 2L) {
  stop("The expression matrix does not contain two estimable principal components.")
}

score_matrix <- pca_fit$ind$coord
variance_percent <- pca_fit$eig[, 2L]
x_axis_base <- trimws(args$x_axis_title)
y_axis_base <- trimws(args$y_axis_title)
if (!nzchar(x_axis_base)) {
  x_axis_base <- "PC1"
}
if (!nzchar(y_axis_base)) {
  y_axis_base <- "PC2"
}
x_axis_label <- sprintf(
  "%s (%.2f%%)",
  x_axis_base,
  variance_percent[[1L]]
)
y_axis_label <- sprintf(
  "%s (%.2f%%)",
  y_axis_base,
  variance_percent[[2L]]
)

pca_data <- data.frame(
  sample_id = rownames(score_matrix),
  Group = sample_group,
  PC1 = score_matrix[, 1L],
  PC2 = score_matrix[, 2L],
  PC1_variance_percent = variance_percent[[1L]],
  PC2_variance_percent = variance_percent[[2L]],
  check.names = FALSE
)

group_colors <- c(
  Normal = args$normal_fill_color,
  Tumor = args$tumor_fill_color
)

ellipse_parts <- lapply(
  split(pca_data, pca_data$Group, drop = TRUE),
  make_confidence_ellipse
)
ellipse_parts <- Filter(Negate(is.null), ellipse_parts)
ellipse_data <- if (length(ellipse_parts) > 0L) {
  do.call(rbind, ellipse_parts)
} else {
  data.frame(
    PC1 = numeric(),
    PC2 = numeric(),
    Group = character()
  )
}
ellipse_data$Group <- factor(
  ellipse_data$Group,
  levels = c("Normal", "Tumor")
)

## ------------------------------------------------------------------ plotting
plot_object <- ggplot2::ggplot(
  pca_data,
  ggplot2::aes(x = PC1, y = PC2)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "black"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "black"
  )

if (nrow(ellipse_data) > 0L) {
  plot_object <- plot_object +
    ggplot2::geom_polygon(
      data = ellipse_data,
      ggplot2::aes(x = PC1, y = PC2, fill = Group, group = Group),
      inherit.aes = FALSE,
      alpha = 0.10,
      color = NA,
      show.legend = FALSE
    ) +
    ggplot2::geom_path(
      data = ellipse_data,
      ggplot2::aes(x = PC1, y = PC2, color = Group, group = Group),
      inherit.aes = FALSE,
      linewidth = 0.7,
      show.legend = FALSE
    )
}

if (args$shape >= 21L && args$shape <= 25L) {
  plot_object <- plot_object +
    ggplot2::geom_point(
      ggplot2::aes(fill = Group),
      shape = args$shape,
      size = args$size,
      alpha = args$opacity,
      color = if (isTRUE(args$is_show_border)) "black" else NA,
      stroke = if (isTRUE(args$is_show_border)) 0.45 else 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = group_colors,
      breaks = c("Normal", "Tumor"),
      drop = FALSE,
      name = args$legend_title
    ) +
    ggplot2::scale_color_manual(values = group_colors, guide = "none") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        override.aes = list(
          shape = args$shape,
          size = args$size,
          alpha = 1,
          color = if (isTRUE(args$is_show_border)) "black" else NA,
          stroke = if (isTRUE(args$is_show_border)) 0.45 else 0
        )
      )
    )
} else {
  # Shapes 0--20 do not have separate ggplot2 fill and outline aesthetics.
  # In particular, `stroke` is ignored by the solid shape 16 used by the UI.
  # Draw a slightly larger black copy underneath so Show Border has a visible,
  # consistent effect for every frontend shape, then draw the coloured point.
  if (isTRUE(args$is_show_border)) {
    border_expansion <- max(0.75, min(1.5, args$size * 0.20))
    plot_object <- plot_object +
      ggplot2::geom_point(
        shape = args$shape,
        size = args$size + border_expansion,
        alpha = args$opacity,
        color = "black",
        show.legend = FALSE
      )
  }

  plot_object <- plot_object +
    ggplot2::geom_point(
      ggplot2::aes(color = Group),
      shape = args$shape,
      size = args$size,
      alpha = args$opacity,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = group_colors,
      breaks = c("Normal", "Tumor"),
      drop = FALSE,
      name = args$legend_title
    ) +
    ggplot2::scale_fill_manual(values = group_colors, guide = "none")
}

plot_object <- plot_object +
  ggplot2::labs(
    title = args$main_title,
    x = x_axis_label,
    y = y_axis_label
  ) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw(base_size = args$font_size, base_family = "sans") +
  ggplot2::theme(
    # Keep the PCA title centered over the plotting panel, independently of
    # any legend placed to the right or elsewhere in the exported figure.
    plot.title.position = "panel",
    plot.title = ggplot2::element_text(
      size = args$font_size + 2,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = ggplot2::element_text(size = args$axis_title_size),
    axis.text = ggplot2::element_text(size = args$axis_text_size),
    legend.title = ggplot2::element_text(size = args$legend_title_size),
    legend.text = ggplot2::element_text(size = args$legend_text_size),
    legend.position = if (!isTRUE(args$is_show_legend) || identical(legend_position, "none")) {
      "none"
    } else if (legend_inside) "inside" else legend_position,
    legend.position.inside = c(0.98, 0.98),
    legend.justification = if (legend_inside) c(1, 1) else "center",
    legend.background = if (legend_inside) {
      ggplot2::element_rect(fill = scales::alpha("white", 0.82), colour = "grey75")
    } else ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(
      linewidth = 0.3,
      color = "#E5E5E5"
    ),
    panel.grid.minor = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA)
  )

## ------------------------------------------------------------------ export
ensure_parent_directory(args$output_png)
ggplot2::ggsave(
  filename = args$output_png,
  plot = plot_object,
  width = args$width,
  height = args$height,
  units = "px",
  dpi = 300,
  bg = "white"
)

if (!is.null(args$output_pdf) && nzchar(trimws(args$output_pdf))) {
  ensure_parent_directory(args$output_pdf)
  pdf_device <- if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  }
  ggplot2::ggsave(
    filename = args$output_pdf,
    plot = plot_object,
    width = args$width / 300,
    height = args$height / 300,
    units = "in",
    device = pdf_device,
    bg = "white"
  )
}

ensure_parent_directory(args$output_csv)
utils::write.csv(
  pca_data,
  file = args$output_csv,
  row.names = FALSE,
  quote = TRUE
)

message(
  "PCA completed: ", nrow(pca_data), " samples; ",
  nrow(expression_matrix), " features used; ",
  length(excluded_nonfinite), " features excluded for NA/Inf; ",
  length(excluded_zero_variance), " excluded for zero variance; ",
  "PC1 = ", format(variance_percent[[1L]], digits = 4), "%; ",
  "PC2 = ", format(variance_percent[[2L]], digits = 4), "%; ",
  if (nrow(ellipse_data) > 0L) {
    "95% confidence ellipse(s) drawn."
  } else {
    "confidence ellipses unavailable for the supplied group sizes/covariance."
  }
)

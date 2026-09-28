#!/usr/bin/env Rscript

# EVOmicsDB GO over-representation analysis: bar plot

required_packages <- c(
  "argparse", "ggplot2", "clusterProfiler", "org.Hs.eg.db"
)
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

wrap_text <- function(x, width = 55L) {
  vapply(
    as.character(x),
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1)
  )
}

standardize_go_export <- function(result = NULL) {
  schema <- list(
    ONTOLOGY = character(),
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    RichFactor = numeric(),
    FoldEnrichment = numeric(),
    zScore = numeric(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer(),
    minus_log10_adjusted_p = numeric()
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

empty_go_export <- function() {
  standardize_go_export()
}

write_empty_outputs <- function(message_text, args) {
  empty_plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = message_text,
      size = 4,
      color = "#555555"
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title = args$main_title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = args$main_title_size,
        face = "bold",
        hjust = 0.5
      ),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  ensure_parent_directory(args$output_png)
  ggplot2::ggsave(
    args$output_png,
    empty_plot,
    width = args$width,
    height = args$height,
    units = "px",
    dpi = 300,
    bg = "white"
  )
  if (!is.null(args$output_pdf) && nzchar(trimws(args$output_pdf))) {
    ensure_parent_directory(args$output_pdf)
    ggplot2::ggsave(
      args$output_pdf,
      empty_plot,
      width = args$width / 300,
      height = args$height / 300,
      units = "in",
      device = if (isTRUE(capabilities("cairo"))) {
        grDevices::cairo_pdf
      } else {
        grDevices::pdf
      },
      bg = "white"
    )
  }
  ensure_parent_directory(args$output_csv)
  utils::write.csv(
    empty_go_export(),
    args$output_csv,
    row.names = FALSE,
    quote = TRUE
  )
  invisible(NULL)
}

## ------------------------------------------------------------------ arguments
parser <- argparse::ArgumentParser(
  description = "Generate an EVOmicsDB GO enrichment bar plot."
)
parser$add_argument("--input_rda", type = "character", required = TRUE)
parser$add_argument("--output_png", type = "character", required = TRUE)
parser$add_argument("--output_csv", type = "character", required = TRUE)
parser$add_argument("--output_pdf", type = "character", default = NULL)

# Differential feature selection
parser$add_argument("--is_use_padj", type = "logical", default = TRUE)
parser$add_argument("--log2fc_threshold", type = "numeric", default = 1)
parser$add_argument("--pvalue_threshold", type = "numeric", default = 0.05)

# Web UI parameters
parser$add_argument("--top_n", type = "integer", default = 5)
parser$add_argument("--is_show_legend", type = "logical", default = TRUE)
parser$add_argument("--legend_title_size", type = "numeric", default = 7)
parser$add_argument("--legend_text_size", type = "numeric", default = 7)
parser$add_argument("--legend_position", type = "character", default = "right")
parser$add_argument("--main_title", type = "character", default = "GO Enrichment Analysis")
parser$add_argument(
  "--x_axis_title",
  type = "character",
  default = "-log10(adj. P-value)"
)
parser$add_argument("--y_axis_title", type = "character", default = "GO Term")
parser$add_argument("--main_title_size", type = "numeric", default = 10)
parser$add_argument("--axis_title_size", type = "numeric", default = 7)
parser$add_argument("--axis_text_size", type = "numeric", default = 7)
parser$add_argument("--BP_fill_color", type = "character", default = "#E64B35")
parser$add_argument("--CC_fill_color", type = "character", default = "#4DBBD5")
parser$add_argument("--MF_fill_color", type = "character", default = "#8C8C8C")
parser$add_argument("--is_show_border", type = "logical", default = TRUE)
parser$add_argument("--border_color", type = "character", default = "black")
parser$add_argument("--border_size", type = "numeric", default = 0.5)
parser$add_argument("--width", type = "numeric", default = 2500)
parser$add_argument("--height", type = "numeric", default = 1350)

# Stable production defaults not currently exposed in the UI
parser$add_argument("--enrichment_pvalue_cutoff", type = "numeric", default = 0.05)
parser$add_argument("--p_adjust_method", type = "character", default = "BH")

args <- parser$parse_args()

## ------------------------------------------------------------------ validation
if (tolower(gsub("[[:space:]_-]+", "", args$x_axis_title)) == "genecount") {
  warning(
    "Legacy X-axis title 'Gene Count' does not match the plotted quantity. ",
    "Using '-log10(adj. P-value)' instead."
  )
  args$x_axis_title <- "-log10(adj. P-value)"
}

if (!file.exists(args$input_rda)) {
  stop("Input RDA does not exist: ", args$input_rda)
}
if (!nzchar(trimws(args$output_png)) || !nzchar(trimws(args$output_csv))) {
  stop("--output_png and --output_csv cannot be empty.")
}
if (length(args$top_n) != 1L || is.na(args$top_n) || args$top_n < 1L) {
  stop("--top_n must be a positive integer.")
}
validate_scalar(args$log2fc_threshold, "log2fc_threshold")
if (args$log2fc_threshold < 0) {
  stop("--log2fc_threshold must be non-negative.")
}
validate_probability(args$pvalue_threshold, "pvalue_threshold")
validate_probability(
  args$enrichment_pvalue_cutoff,
  "enrichment_pvalue_cutoff"
)
validate_positive(args$legend_title_size, "legend_title_size")
validate_positive(args$legend_text_size, "legend_text_size")
validate_positive(args$main_title_size, "main_title_size")
validate_positive(args$axis_title_size, "axis_title_size")
validate_positive(args$axis_text_size, "axis_text_size")
validate_positive(args$width, "width")
validate_positive(args$height, "height")
validate_scalar(args$border_size, "border_size")
if (args$border_size < 0) {
  stop("--border_size must be non-negative.")
}

legend_position <- tolower(trimws(args$legend_position))
if (!legend_position %in% c("right", "left", "top", "bottom", "inside", "none")) {
  stop("--legend_position must be right, left, top, bottom, inside, or none.")
}
legend_inside <- identical(legend_position, "inside")
if (!args$p_adjust_method %in% stats::p.adjust.methods) {
  stop(
    "--p_adjust_method must be one of: ",
    paste(stats::p.adjust.methods, collapse = ", "),
    "."
  )
}
validate_color(args$BP_fill_color, "BP_fill_color")
validate_color(args$CC_fill_color, "CC_fill_color")
validate_color(args$MF_fill_color, "MF_fill_color")
validate_color(args$border_color, "border_color")

## ------------------------------------------------------------------ load and select genes
input_environment <- new.env(parent = emptyenv())
load(args$input_rda, envir = input_environment)
if (!exists("deg", envir = input_environment, inherits = FALSE)) {
  stop("The input RDA is missing the required `deg` table.")
}
deg <- as.data.frame(input_environment$deg)

p_column <- if (isTRUE(args$is_use_padj)) "adj.P.Val" else "P.Value"
required_columns <- c("logFC", p_column, "ENTREZID")
missing_columns <- setdiff(required_columns, colnames(deg))
if (length(missing_columns) > 0L) {
  stop(
    "The mapped DEG table is missing: ",
    paste(missing_columns, collapse = ", "),
    ". Run enrichment_analysis.r first."
  )
}
deg$logFC <- suppressWarnings(as.numeric(deg$logFC))
deg[[p_column]] <- suppressWarnings(as.numeric(deg[[p_column]]))
deg$ENTREZID <- trimws(as.character(deg$ENTREZID))

valid_p <- !is.na(deg[[p_column]]) &
  is.finite(deg[[p_column]]) &
  deg[[p_column]] >= 0 &
  deg[[p_column]] <= 1
is_significant <- valid_p &
  is.finite(deg$logFC) &
  abs(deg$logFC) >= args$log2fc_threshold &
  deg[[p_column]] <= args$pvalue_threshold &
  !is.na(deg$ENTREZID) &
  nzchar(deg$ENTREZID)

foreground_entrez <- unique(deg$ENTREZID[is_significant])
universe_entrez <- if (exists(
  "universe_entrez",
  envir = input_environment,
  inherits = FALSE
)) {
  unique(trimws(as.character(input_environment$universe_entrez)))
} else {
  unique(deg$ENTREZID[!is.na(deg$ENTREZID) & nzchar(deg$ENTREZID)])
}
universe_entrez <- universe_entrez[
  !is.na(universe_entrez) & nzchar(universe_entrez)
]
foreground_entrez <- intersect(foreground_entrez, universe_entrez)

if (length(foreground_entrez) == 0L) {
  write_empty_outputs("No significant mapped genes for GO enrichment", args)
  warning("No significant mapped genes for GO enrichment.")
  quit(save = "no", status = 0)
}

## ------------------------------------------------------------------ enrichment
ego <- suppressMessages(
  clusterProfiler::enrichGO(
    gene = foreground_entrez,
    universe = universe_entrez,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = args$p_adjust_method,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
)
ego_result <- as.data.frame(ego)
if (nrow(ego_result) > 0L) {
  ego_result <- ego_result[
    ego_result$ONTOLOGY %in% c("BP", "CC", "MF") &
      is.finite(ego_result$p.adjust) &
      ego_result$p.adjust < args$enrichment_pvalue_cutoff,
    ,
    drop = FALSE
  ]
}

if (nrow(ego_result) == 0L) {
  write_empty_outputs("No significant GO terms under the current criteria", args)
  warning("No significant GO terms under the current criteria.")
  quit(save = "no", status = 0)
}

ontology_order <- c("BP", "CC", "MF")
top_by_ontology <- lapply(ontology_order, function(ontology) {
  ontology_result <- ego_result[
    ego_result$ONTOLOGY == ontology,
    ,
    drop = FALSE
  ]
  ontology_result <- ontology_result[
    order(
      ontology_result$p.adjust,
      ontology_result$pvalue,
      -ontology_result$Count
    ),
    ,
    drop = FALSE
  ]
  utils::head(ontology_result, args$top_n)
})
top_by_ontology <- Filter(function(x) nrow(x) > 0L, top_by_ontology)
ego_df <- do.call(rbind, top_by_ontology)
ego_df$ONTOLOGY <- factor(ego_df$ONTOLOGY, levels = ontology_order)
ego_df$minus_log10_adjusted_p <- -log10(
  pmax(ego_df$p.adjust, .Machine$double.xmin)
)
ego_df$plot_key <- paste0(
  as.character(ego_df$ONTOLOGY),
  "__",
  seq_len(nrow(ego_df))
)
plot_label_lookup <- stats::setNames(
  wrap_text(ego_df$Description),
  ego_df$plot_key
)
ego_df$plot_key <- factor(
  ego_df$plot_key,
  levels = rev(unique(ego_df$plot_key))
)

ontology_colors <- c(
  BP = args$BP_fill_color,
  CC = args$CC_fill_color,
  MF = args$MF_fill_color
)

plot_object <- ggplot2::ggplot(
  ego_df,
  ggplot2::aes(
    x = minus_log10_adjusted_p,
    y = plot_key,
    fill = ONTOLOGY
  )
) +
  ggplot2::geom_col(
    width = 0.88,
    color = if (isTRUE(args$is_show_border)) args$border_color else NA,
    linewidth = if (isTRUE(args$is_show_border)) args$border_size else 0
  ) +
  ggplot2::facet_grid(
    ONTOLOGY ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.06))
  ) +
  ggplot2::scale_y_discrete(
    labels = plot_label_lookup,
    expand = ggplot2::expansion(add = c(0.12, 0.12))
  ) +
  ggplot2::scale_fill_manual(
    values = ontology_colors,
    breaks = c("BP", "CC", "MF"),
    drop = FALSE,
    name = "Ontology"
  ) +
  ggplot2::labs(
    title = args$main_title,
    x = args$x_axis_title,
    y = args$y_axis_title
  ) +
  ggplot2::theme_minimal(base_size = 7, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = args$main_title_size,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.title = ggplot2::element_text(
      size = args$axis_title_size,
      face = "bold",
      color = "black"
    ),
    axis.text = ggplot2::element_text(
      size = args$axis_text_size,
      color = "black"
    ),
    legend.title = ggplot2::element_text(
      size = args$legend_title_size,
      face = "bold"
    ),
    legend.text = ggplot2::element_text(size = args$legend_text_size),
    legend.position = if (!isTRUE(args$is_show_legend) || identical(legend_position, "none")) {
      "none"
    } else if (legend_inside) "inside" else legend_position,
    legend.position.inside = c(0.98, 0.98),
    legend.justification = if (legend_inside) c(1, 1) else "center",
    legend.background = if (legend_inside) {
      ggplot2::element_rect(fill = scales::alpha("white", 0.82), colour = "grey75")
    } else ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(
      size = args$axis_text_size,
      face = "bold",
      color = "black"
    ),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = if (isTRUE(args$is_show_border)) {
      ggplot2::element_rect(
        fill = NA,
        color = args$border_color,
        linewidth = args$border_size
      )
    } else {
      ggplot2::element_blank()
    },
    strip.background.y = ggplot2::element_rect(
      fill = "#F2F2F2",
      color = if (isTRUE(args$is_show_border)) args$border_color else NA,
      linewidth = if (isTRUE(args$is_show_border)) args$border_size else 0
    ),
    strip.text.y = ggplot2::element_text(
      angle = 0,
      face = "bold",
      size = args$axis_text_size,
      color = "black"
    ),
    panel.spacing.y = grid::unit(1.5, "mm"),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.margin = ggplot2::margin(8, 12, 8, 8)
  )

## ------------------------------------------------------------------ export
ensure_parent_directory(args$output_png)
ggplot2::ggsave(
  args$output_png,
  plot_object,
  width = args$width,
  height = args$height,
  units = "px",
  dpi = 300,
  bg = "white"
)
if (!is.null(args$output_pdf) && nzchar(trimws(args$output_pdf))) {
  ensure_parent_directory(args$output_pdf)
  ggplot2::ggsave(
    args$output_pdf,
    plot_object,
    width = args$width / 300,
    height = args$height / 300,
    units = "in",
    device = if (isTRUE(capabilities("cairo"))) {
      grDevices::cairo_pdf
    } else {
      grDevices::pdf
    },
    bg = "white"
  )
}
ensure_parent_directory(args$output_csv)
export_df <- ego_df
export_df$plot_key <- NULL
export_df <- standardize_go_export(export_df)
utils::write.csv(
  export_df,
  args$output_csv,
  row.names = FALSE,
  quote = TRUE
)

message(
  "GO bar plot completed: ", length(foreground_entrez),
  " foreground ENTREZ IDs; ", length(universe_entrez),
  " background IDs; ", nrow(ego_result),
  " significant GO terms; ", nrow(ego_df), " plotted."
)

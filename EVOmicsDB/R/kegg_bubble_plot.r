#!/usr/bin/env Rscript

# EVOmicsDB KEGG enrichment bubble plot.
# X = Rich Factor; size = gene count; fill = -log10(adjusted P-value).

required_packages <- c("argparse", "ggplot2", "clusterProfiler")
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

get_script_directory <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command, value = TRUE)
  if (length(file_argument) == 0L) {
    return(getwd())
  }
  dirname(normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_argument[[1L]]))))
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

wrap_text <- function(x, width = 48L) {
  vapply(
    as.character(x),
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1)
  )
}

parse_ratio <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(parts, function(value) {
    if (length(value) != 2L) {
      return(NA_real_)
    }
    numerator <- suppressWarnings(as.numeric(value[[1L]]))
    denominator <- suppressWarnings(as.numeric(value[[2L]]))
    if (!is.finite(numerator) || !is.finite(denominator) || denominator <= 0) {
      return(NA_real_)
    }
    numerator / denominator
  }, numeric(1))
}

parse_ratio_numerator <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(parts, function(value) {
    if (length(value) < 1L) {
      return(NA_real_)
    }
    suppressWarnings(as.numeric(value[[1L]]))
  }, numeric(1))
}

parse_ratio_denominator <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(parts, function(value) {
    if (length(value) != 2L) {
      return(NA_real_)
    }
    suppressWarnings(as.numeric(value[[2L]]))
  }, numeric(1))
}

standardize_kegg_result <- function(result = NULL) {
  schema <- list(
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
  if (!is.null(result) && n > 0L) {
    if (!"RichFactor" %in% names(result)) {
      result$RichFactor <- NA_real_
    }
    result$RichFactor <- suppressWarnings(as.numeric(result$RichFactor))
    missing_rich_factor <- !is.finite(result$RichFactor)
    if (any(missing_rich_factor) &&
        all(c("Count", "BgRatio") %in% names(result))) {
      background_term_size <- parse_ratio_numerator(result$BgRatio)
      calculated <- suppressWarnings(as.numeric(result$Count)) /
        background_term_size
      result$RichFactor[missing_rich_factor] <- calculated[missing_rich_factor]
    }

    if (!"FoldEnrichment" %in% names(result)) {
      result$FoldEnrichment <- NA_real_
    }
    result$FoldEnrichment <- suppressWarnings(
      as.numeric(result$FoldEnrichment)
    )
    missing_fold <- !is.finite(result$FoldEnrichment)
    if (any(missing_fold) &&
        all(c("GeneRatio", "BgRatio") %in% names(result))) {
      calculated <- parse_ratio(result$GeneRatio) / parse_ratio(result$BgRatio)
      result$FoldEnrichment[missing_fold] <- calculated[missing_fold]
    }
  }

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

resolve_local_mapping <- function(explicit_path, script_directory) {
  explicit_path <- if (is.null(explicit_path)) "" else trimws(explicit_path)
  environment_path <- trimws(
    Sys.getenv("EVOMICS_KEGG_MAPPING_RDA", unset = "")
  )
  if (nzchar(explicit_path) && !file.exists(explicit_path)) {
    stop("--kegg_mapping_rda does not exist: ", explicit_path)
  }
  if (!nzchar(explicit_path) &&
      nzchar(environment_path) &&
      !file.exists(environment_path)) {
    stop(
      "EVOMICS_KEGG_MAPPING_RDA points to a missing file: ",
      environment_path
    )
  }
  candidates <- c(
    explicit_path,
    environment_path,
    file.path(script_directory, "../local_kegg_db/kegg_mapping.rda"),
    file.path(script_directory, "../../../db/rda/kegg_mapping.rda"),
    file.path(script_directory, "../../local_kegg_db/kegg_mapping.rda")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    return(NULL)
  }
  normalizePath(existing[[1L]])
}

load_local_mapping <- function(path) {
  mapping_environment <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = mapping_environment)
  data_names <- loaded_names[vapply(
    loaded_names,
    function(name) {
      object <- mapping_environment[[name]]
      is.data.frame(object) || is.matrix(object)
    },
    logical(1)
  )]
  if (length(data_names) == 0L) {
    stop("The local KEGG mapping RDA does not contain a table.")
  }
  mapping <- as.data.frame(
    mapping_environment[[data_names[[1L]]]],
    stringsAsFactors = FALSE
  )
  lower_names <- tolower(names(mapping))
  find_column <- function(candidates, required = TRUE) {
    position <- match(tolower(candidates), lower_names)
    position <- position[!is.na(position)]
    if (length(position) == 0L) {
      if (required) {
        stop(
          "Local KEGG mapping is missing one of: ",
          paste(candidates, collapse = ", ")
        )
      }
      return(NULL)
    }
    names(mapping)[position[[1L]]]
  }
  term_column <- find_column(c("gs_name", "term", "pathway_id", "id"))
  gene_column <- find_column(
    c("entrez_id", "entrezid", "target_entrez", "gene_symbol", "symbol")
  )
  name_column <- find_column(
    c("gs_description", "description", "pathway_name", "name"),
    required = FALSE
  )
  gene_type <- if (tolower(gene_column) %in%
      c("gene_symbol", "symbol")) "SYMBOL" else "ENTREZID"

  term2gene <- data.frame(
    term = trimws(as.character(mapping[[term_column]])),
    gene = trimws(as.character(mapping[[gene_column]])),
    stringsAsFactors = FALSE
  )
  if (gene_type == "SYMBOL") {
    term2gene$gene <- toupper(term2gene$gene)
  }
  keep <- !is.na(term2gene$term) &
    nzchar(term2gene$term) &
    !is.na(term2gene$gene) &
    nzchar(term2gene$gene)
  term2gene <- unique(term2gene[keep, , drop = FALSE])
  if (nrow(term2gene) == 0L) {
    stop("The local KEGG mapping contains no usable term-gene pairs.")
  }

  if (is.null(name_column)) {
    term2name <- unique(data.frame(
      term = term2gene$term,
      name = term2gene$term,
      stringsAsFactors = FALSE
    ))
  } else {
    term2name <- unique(data.frame(
      term = trimws(as.character(mapping[[term_column]])),
      name = trimws(as.character(mapping[[name_column]])),
      stringsAsFactors = FALSE
    ))
    term2name <- term2name[
      !is.na(term2name$term) & nzchar(term2name$term) &
        !is.na(term2name$name) & nzchar(term2name$name),
      ,
      drop = FALSE
    ]
  }
  list(term2gene = term2gene, term2name = term2name, gene_type = gene_type)
}

write_empty_outputs <- function(message_text, args) {
  empty_plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0.5, y = 0.5, label = message_text,
      size = 4, color = "#555555"
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
    args$output_png, empty_plot,
    width = args$width, height = args$height,
    units = "px", dpi = 300, bg = "white"
  )
  if (!is.null(args$output_pdf) && nzchar(trimws(args$output_pdf))) {
    ensure_parent_directory(args$output_pdf)
    ggplot2::ggsave(
      args$output_pdf, empty_plot,
      width = args$width / 300, height = args$height / 300,
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
    standardize_kegg_result(),
    args$output_csv,
    row.names = FALSE,
    quote = TRUE
  )
  invisible(NULL)
}

## ------------------------------------------------------------------ arguments
parser <- argparse::ArgumentParser(
  description = "Generate an EVOmicsDB KEGG enrichment bubble plot."
)
parser$add_argument("--input_rda", type = "character", required = TRUE)
parser$add_argument("--output_png", type = "character", required = TRUE)
parser$add_argument("--output_csv", type = "character", required = TRUE)
parser$add_argument("--output_pdf", type = "character", default = NULL)
parser$add_argument("--is_use_padj", type = "logical", default = TRUE)
parser$add_argument("--log2fc_threshold", type = "numeric", default = 1)
parser$add_argument("--pvalue_threshold", type = "numeric", default = 0.05)
parser$add_argument("--top_n", type = "integer", default = 5)
parser$add_argument("--wrap_width", type = "integer", default = 48)
parser$add_argument("--kegg_mapping_rda", type = "character", default = NULL)
parser$add_argument("--is_show_legend", type = "logical", default = TRUE)
parser$add_argument("--legend_title_size", type = "numeric", default = 7)
parser$add_argument("--legend_text_size", type = "numeric", default = 7)
parser$add_argument("--legend_position", type = "character", default = "right")
parser$add_argument("--main_title", type = "character", default = "KEGG Analysis")
parser$add_argument("--x_axis_title", type = "character", default = "Rich Factor")
parser$add_argument("--y_axis_title", type = "character", default = "KEGG Pathway")
parser$add_argument("--main_title_size", type = "numeric", default = 10)
parser$add_argument("--axis_title_size", type = "numeric", default = 7)
parser$add_argument("--axis_text_size", type = "numeric", default = 7)
parser$add_argument("--low_fill_color", type = "character", default = "#4DBBD5")
parser$add_argument("--high_fill_color", type = "character", default = "#E64B35")
parser$add_argument("--border_color", type = "character", default = "black")
parser$add_argument("--border_size", type = "numeric", default = 0.5)
parser$add_argument("--is_show_border", type = "logical", default = TRUE)
parser$add_argument("--width", type = "numeric", default = 2500)
parser$add_argument("--height", type = "numeric", default = 1350)
parser$add_argument("--enrichment_pvalue_cutoff", type = "numeric", default = 0.05)
parser$add_argument("--p_adjust_method", type = "character", default = "BH")
args <- parser$parse_args()

## ------------------------------------------------------------------ validation
if (!file.exists(args$input_rda)) {
  stop("Input RDA does not exist: ", args$input_rda)
}
if (!nzchar(trimws(args$output_png)) || !nzchar(trimws(args$output_csv))) {
  stop("--output_png and --output_csv cannot be empty.")
}
if (length(args$top_n) != 1L || is.na(args$top_n) || args$top_n < 1L) {
  stop("--top_n must be a positive integer.")
}
if (length(args$wrap_width) != 1L ||
    is.na(args$wrap_width) ||
    args$wrap_width < 1L) {
  stop("--wrap_width must be a positive integer.")
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
validate_color(args$low_fill_color, "low_fill_color")
validate_color(args$high_fill_color, "high_fill_color")
validate_color(args$border_color, "border_color")

## ------------------------------------------------------------------ load/select
input_environment <- new.env(parent = emptyenv())
load(args$input_rda, envir = input_environment)
if (!exists("deg", envir = input_environment, inherits = FALSE)) {
  stop("The input RDA is missing the required `deg` table.")
}
deg <- as.data.frame(input_environment$deg)
p_column <- if (isTRUE(args$is_use_padj)) "adj.P.Val" else "P.Value"
required_columns <- c("logFC", p_column, "ENTREZID")
missing_columns <- setdiff(required_columns, names(deg))
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
  write_empty_outputs("No significant mapped genes for KEGG enrichment", args)
  warning("No significant mapped genes for KEGG enrichment.")
  quit(save = "no", status = 0)
}

## ------------------------------------------------------------------ enrichment
mapping_path <- resolve_local_mapping(
  args$kegg_mapping_rda,
  get_script_directory()
)
enrichment_source <- "KEGG online"
if (!is.null(mapping_path)) {
  local_mapping <- load_local_mapping(mapping_path)
  if (local_mapping$gene_type == "ENTREZID") {
    foreground_ids <- foreground_entrez
    universe_ids <- universe_entrez
  } else {
    if (!"SYMBOL" %in% names(deg)) {
      stop(
        "The local KEGG mapping uses gene symbols, but the mapped DEG table ",
        "does not contain a SYMBOL column."
      )
    }
    deg$SYMBOL <- toupper(trimws(as.character(deg$SYMBOL)))
    foreground_ids <- unique(deg$SYMBOL[
      is_significant & !is.na(deg$SYMBOL) & nzchar(deg$SYMBOL)
    ])
    universe_ids <- unique(deg$SYMBOL[
      deg$ENTREZID %in% universe_entrez &
        !is.na(deg$SYMBOL) &
        nzchar(deg$SYMBOL)
    ])
  }
  foreground_ids <- intersect(
    foreground_ids,
    unique(local_mapping$term2gene$gene)
  )
  universe_ids <- intersect(
    universe_ids,
    unique(local_mapping$term2gene$gene)
  )
  foreground_ids <- intersect(foreground_ids, universe_ids)
  if (length(foreground_ids) == 0L) {
    write_empty_outputs("No foreground genes matched the local KEGG mapping", args)
    warning("No foreground genes matched the local KEGG mapping.")
    quit(save = "no", status = 0)
  }
  kk <- suppressMessages(clusterProfiler::enricher(
    gene = foreground_ids,
    universe = universe_ids,
    TERM2GENE = local_mapping$term2gene,
    TERM2NAME = local_mapping$term2name,
    pAdjustMethod = args$p_adjust_method,
    pvalueCutoff = 1,
    qvalueCutoff = 1
  ))
  enrichment_source <- paste0(
    "local mapping (",
    local_mapping$gene_type,
    ")"
  )
} else {
  kk <- tryCatch(
    suppressMessages(clusterProfiler::enrichKEGG(
      gene = foreground_entrez,
      universe = universe_entrez,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pAdjustMethod = args$p_adjust_method,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )),
    error = function(error) {
      stop(
        "Online KEGG enrichment failed: ", conditionMessage(error),
        ". Deploy a local mapping and pass --kegg_mapping_rda, or set ",
        "EVOMICS_KEGG_MAPPING_RDA."
      )
    }
  )
}

kegg_result <- standardize_kegg_result(as.data.frame(kk))
eligible_foreground_counts <- unique(
  parse_ratio_denominator(kegg_result$GeneRatio)
)
eligible_foreground_counts <- eligible_foreground_counts[
  is.finite(eligible_foreground_counts)
]
eligible_background_counts <- unique(
  parse_ratio_denominator(kegg_result$BgRatio)
)
eligible_background_counts <- eligible_background_counts[
  is.finite(eligible_background_counts)
]
kegg_result <- kegg_result[
  is.finite(kegg_result$p.adjust) &
    kegg_result$p.adjust < args$enrichment_pvalue_cutoff &
    is.finite(kegg_result$RichFactor) &
    kegg_result$RichFactor >= 0,
  ,
  drop = FALSE
]
if (nrow(kegg_result) == 0L) {
  write_empty_outputs("No significant KEGG pathways under current criteria", args)
  warning("No significant KEGG pathways under current criteria.")
  quit(save = "no", status = 0)
}
kegg_result <- kegg_result[
  order(kegg_result$p.adjust, kegg_result$pvalue, -kegg_result$Count),
  ,
  drop = FALSE
]
kegg_df <- utils::head(kegg_result, args$top_n)
kegg_df$minus_log10_adjusted_p <- -log10(
  pmax(kegg_df$p.adjust, .Machine$double.xmin)
)
kegg_df$plot_label <- make.unique(
  wrap_text(kegg_df$Description, width = args$wrap_width)
)
kegg_df$plot_label <- factor(
  kegg_df$plot_label,
  levels = rev(kegg_df$plot_label)
)

## ------------------------------------------------------------------ plot/export
plot_object <- ggplot2::ggplot(
  kegg_df,
  ggplot2::aes(
    x = RichFactor,
    y = plot_label,
    fill = minus_log10_adjusted_p,
    size = Count
  )
) +
  ggplot2::geom_point(
    shape = 21,
    alpha = 0.9,
    color = if (isTRUE(args$is_show_border)) args$border_color else NA,
    stroke = if (isTRUE(args$is_show_border)) args$border_size else 0
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0.04, 0.08))
  ) +
  ggplot2::scale_size_continuous(
    range = c(3, 10),
    name = "Gene Count"
  ) +
  ggplot2::scale_fill_gradient(
    low = args$low_fill_color,
    high = args$high_fill_color,
    name = "-log10(adj. P-value)"
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
    plot.background = ggplot2::element_rect(fill = "white", color = NA)
  )

ensure_parent_directory(args$output_png)
ggplot2::ggsave(
  args$output_png, plot_object,
  width = args$width, height = args$height,
  units = "px", dpi = 300, bg = "white"
)
if (!is.null(args$output_pdf) && nzchar(trimws(args$output_pdf))) {
  ensure_parent_directory(args$output_pdf)
  ggplot2::ggsave(
    args$output_pdf, plot_object,
    width = args$width / 300, height = args$height / 300,
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
  standardize_kegg_result(kegg_df),
  args$output_csv,
  row.names = FALSE,
  quote = TRUE
)
message(
  "KEGG bubble plot completed: ", length(foreground_entrez),
  " foreground ENTREZ IDs; ", length(universe_entrez),
  " background IDs; ",
  if (length(eligible_foreground_counts) == 1L) {
    paste0(eligible_foreground_counts, " KEGG-eligible foreground; ")
  } else {
    ""
  },
  if (length(eligible_background_counts) == 1L) {
    paste0(eligible_background_counts, " KEGG-eligible background; ")
  } else {
    ""
  },
  nrow(kegg_result),
  " significant pathways; ", nrow(kegg_df),
  " plotted; source = ", enrichment_source, "."
)

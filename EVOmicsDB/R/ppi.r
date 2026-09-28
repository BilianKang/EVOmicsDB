#!/usr/bin/env Rscript

# EVOmicsDB STRING protein-association network
#
# Selection contract:
#   1. apply the requested differential P/FDR and |log2FC| thresholds;
#   2. rank the passing proteins by |log2FC|, then P/FDR and identifier;
#   3. submit up to submit_n identifiers and construct the complete STRING graph;
#   4. calculate degree on that complete graph, then draw the top display_n hubs;
#   5. retain every submitted identifier and every full-network edge in audit CSVs.
#
# The production API is pinned to STRING v12.0 for reproducibility. Both
# functional associations and physical interactions are supported explicitly.
# STRING responses are cached by a hash of the complete version-pinned,
# canonical request (the URL for GET or the equivalent endpoint/parameters for POST).

required_packages <- c("argparse", "digest", "ggraph", "ggplot2", "httr", "igraph", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them in the backend R environment before deployment."
  )
}

options(timeout = max(300, getOption("timeout", 60)))

ensure_parent_directory <- function(path) {
  parent <- dirname(normalizePath(path, mustWork = FALSE))
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE) &&
      !dir.exists(parent)) {
    stop("Unable to create output directory: ", parent)
  }
}

validate_probability <- function(value, name, lower_open = FALSE) {
  valid_lower <- if (lower_open) value > 0 else value >= 0
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      !valid_lower || value > 1) {
    boundary <- if (lower_open) "(0, 1]" else "[0, 1]"
    stop("--", name, " must be in ", boundary, ".")
  }
}

string_api_request <- function(endpoint, parameters, api_base, cache_dir = NULL,
                               retries = 3L) {
  request_base <- paste0(sub("/+$", "", api_base), "/tsv/", endpoint)
  encoded_parameters <- vapply(
    parameters,
    function(value) utils::URLencode(as.character(value), reserved = TRUE),
    character(1)
  )
  request_url <- paste0(
    request_base, "?",
    paste(names(encoded_parameters), encoded_parameters, sep = "=", collapse = "&")
  )
  cache_file <- NULL
  if (!is.null(cache_dir) && nzchar(trimws(cache_dir))) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    cache_file <- file.path(
      cache_dir,
      paste0("string_v12_", digest::digest(request_url, algo = "sha256"), ".tsv")
    )
    if (file.exists(cache_file) && file.info(cache_file)$size > 0L) {
      return(utils::read.delim(cache_file, header = TRUE, sep = "\t", quote = "",
        check.names = FALSE, stringsAsFactors = FALSE))
    }
  }
  response_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(response_file), add = TRUE)
  # Keep short requests compatible with the historical GET cache, but use
  # POST for long identifier lists so submit_n=500 cannot exceed the web
  # server request-line limit. The cache key remains the canonical endpoint
  # plus encoded parameters, so GET and POST responses are interchangeable.
  use_post <- nchar(request_url, type = "bytes") > 7500L
  status <- 1L
  last_error <- NULL
  for (attempt in seq_len(max(1L, retries))) {
    status <- if (use_post) {
      tryCatch({
        response <- httr::POST(
          request_base,
          body = lapply(parameters, as.character),
          encode = "form",
          httr::timeout(max(30, getOption("timeout", 300)))
        )
        response_status <- httr::status_code(response)
        if (response_status >= 200L && response_status < 300L) {
          writeBin(httr::content(response, as = "raw"), response_file)
          0L
        } else {
          last_error <<- paste0("HTTP ", response_status)
          response_status
        }
      }, error = function(error) {
        last_error <<- conditionMessage(error)
        1L
      })
    } else {
      tryCatch(
        utils::download.file(request_url, response_file, mode = "wb", quiet = TRUE),
        error = function(error) {
          last_error <<- conditionMessage(error)
          1L
        }
      )
    }
    if (identical(status, 0L) && file.exists(response_file) && file.info(response_file)$size > 0L) break
    if (attempt < retries) Sys.sleep(attempt)
  }
  if (length(status) != 1L || is.na(status) || status != 0L) {
    stop("STRING API request failed after ", retries, " attempts: ", last_error %||% paste("status", status), ".")
  }
  if (!file.exists(response_file) || file.info(response_file)$size == 0L) {
    return(data.frame())
  }
  result <- utils::read.delim(
    response_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!is.null(cache_file)) file.copy(response_file, cache_file, overwrite = TRUE)
  result
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x)) y else x

parser <- argparse::ArgumentParser(
  description = "Generate a version-pinned STRING protein-association network."
)
parser$add_argument("--input_rda", type = "character", required = TRUE)
parser$add_argument("--output_png", type = "character", required = TRUE)
parser$add_argument("--output_csv", type = "character", required = TRUE)
parser$add_argument("--output_nodes_csv", type = "character", default = NULL)
parser$add_argument("--output_pdf", type = "character", default = NULL)
parser$add_argument("--output_status_json", type = "character", default = NULL)
parser$add_argument("--is_use_padj", type = "logical", default = TRUE)
parser$add_argument("--log2fc_threshold", type = "numeric", default = 0.585)
parser$add_argument("--pvalue_threshold", type = "numeric", default = 0.05)
parser$add_argument("--submit_n", type = "integer", default = 300)
parser$add_argument("--display_n", type = "integer", default = 10)
parser$add_argument("--top_n", type = "integer", default = NULL,
                    help = "Deprecated alias for display_n")
parser$add_argument("--interaction_score_threshold", type = "numeric", default = 0.7)
parser$add_argument(
  "--network_type",
  type = "character",
  choices = c("functional", "physical"),
  default = "functional"
)
parser$add_argument("--width", type = "numeric", default = 1800)
parser$add_argument("--height", type = "numeric", default = 1350)
parser$add_argument("--string_cache_dir", type = "character", default = NULL)
parser$add_argument("--api_retries", type = "integer", default = 3)
args <- parser$parse_args()

if (!is.null(args$top_n)) args$display_n <- args$top_n
if (length(args$submit_n) != 1L || is.na(args$submit_n) ||
    args$submit_n < 2L || args$submit_n > 500L) {
  stop("--submit_n must be an integer between 2 and 500.")
}
if (length(args$display_n) != 1L || is.na(args$display_n) ||
    args$display_n < 2L || args$display_n > 100L) {
  stop("--display_n must be an integer between 2 and 100.")
}
if (length(args$log2fc_threshold) != 1L || is.na(args$log2fc_threshold) ||
    !is.finite(args$log2fc_threshold) || args$log2fc_threshold < 0) {
  stop("--log2fc_threshold must be finite and non-negative.")
}
validate_probability(args$pvalue_threshold, "pvalue_threshold", lower_open = TRUE)
validate_probability(args$interaction_score_threshold, "interaction_score_threshold")
if (args$interaction_score_threshold < 0.15) {
  stop("--interaction_score_threshold must be between 0.15 and 1.")
}
if (length(args$width) != 1L || length(args$height) != 1L ||
    is.na(args$width) || is.na(args$height) ||
    !is.finite(args$width) || !is.finite(args$height) ||
    args$width < 600 || args$height < 600 ||
    args$width > 10000 || args$height > 10000) {
  stop("--width and --height must each be between 600 and 10000 pixels.")
}

node_output <- args$output_nodes_csv
if (is.null(node_output) || !nzchar(trimws(node_output))) {
  node_output <- sub("\\.csv$", "_nodes.csv", args$output_csv, ignore.case = TRUE)
  if (identical(node_output, args$output_csv)) node_output <- paste0(args$output_csv, "_nodes.csv")
}
node_export_columns <- c(
  "feature_id", "logFC", "selected_p_value", "regulation", "p_value_tie_size",
  "submission_rank", "STRING_id", "preferred_name", "string_mapped",
  "string_id_collapsed", "graph_representative", "represented_by", "degree_full",
  "displayed", "organism_taxid", "string_version", "network_type", "score_threshold"
)

write_node_audit <- function(data, output_path, degree_lookup = NULL,
                             displayed_ids = character(), string_version_value = "12.0") {
  defaults <- list(
    STRING_id = NA_character_, preferred_name = NA_character_,
    string_mapped = FALSE, string_id_collapsed = FALSE,
    graph_representative = FALSE, represented_by = NA_character_,
    degree_full = NA_integer_, displayed = FALSE,
    organism_taxid = organism_taxid, string_version = string_version_value,
    network_type = args$network_type,
    score_threshold = args$interaction_score_threshold
  )
  for (column in names(defaults)) {
    if (!(column %in% colnames(data))) {
      data[[column]] <- rep(defaults[[column]], nrow(data))
    }
  }
  if (!is.null(degree_lookup)) {
    data$degree_full <- as.integer(unname(degree_lookup[data$STRING_id]))
    data$displayed <- data$graph_representative & data$STRING_id %in% displayed_ids
  }
  ensure_parent_directory(output_path)
  utils::write.csv(data[, node_export_columns, drop = FALSE], output_path,
                   row.names = FALSE, quote = TRUE)
}

write_status <- function(code, message_text, details = list()) {
  if (is.null(args$output_status_json) || !nzchar(trimws(args$output_status_json))) {
    return(invisible(NULL))
  }
  ensure_parent_directory(args$output_status_json)
  jsonlite::write_json(
    list(code = code, message = message_text, details = details),
    args$output_status_json,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
}

finish_semantic_empty <- function(code, message_text, details = list()) {
  write_status(code, message_text, details)
  message(message_text)
  quit(save = "no", status = 0L)
}

input_environment <- new.env(parent = emptyenv())
load(args$input_rda, envir = input_environment)
if (!exists("deg", envir = input_environment, inherits = FALSE)) {
  stop("The input RDA is missing the required `deg` table.")
}
deg <- as.data.frame(input_environment$deg)
preprocessing <- if (exists(
  "preprocessing",
  envir = input_environment,
  inherits = FALSE
)) {
  input_environment$preprocessing
} else {
  list()
}

organism_taxid <- suppressWarnings(as.integer(preprocessing$organism_taxid))
if (length(organism_taxid) != 1L || is.na(organism_taxid) || organism_taxid < 1L) {
  organism_taxid <- 9606L
  warning("RDA organism_taxid is unavailable; falling back to Homo sapiens (9606).")
}
group_labels <- preprocessing$group_labels
case_code <- if (is.list(group_labels) && !is.null(group_labels$case_code)) {
  as.character(group_labels$case_code)
} else {
  "case"
}
control_code <- if (is.list(group_labels) && !is.null(group_labels$control_code)) {
  as.character(group_labels$control_code)
} else {
  "control"
}
contrast_label <- paste0(case_code, " - ", control_code)

p_column <- if (isTRUE(args$is_use_padj)) "adj.P.Val" else "P.Value"
required_columns <- c("logFC", p_column)
missing_columns <- setdiff(required_columns, colnames(deg))
if (length(missing_columns)) {
  stop("The DEG table is missing: ", paste(missing_columns, collapse = ", "))
}

feature_id <- if ("feature_id" %in% colnames(deg)) {
  trimws(as.character(deg$feature_id))
} else {
  trimws(rownames(deg))
}
logfc <- suppressWarnings(as.numeric(deg$logFC))
p_value <- suppressWarnings(as.numeric(deg[[p_column]]))
significant <- nzchar(feature_id) & is.finite(logfc) & is.finite(p_value) &
  p_value <= args$pvalue_threshold & abs(logfc) >= args$log2fc_threshold
candidate <- data.frame(
  feature_id = feature_id[significant],
  logFC = logfc[significant],
  selected_p_value = p_value[significant],
  regulation = ifelse(logfc[significant] > 0, "Up", "Down"),
  stringsAsFactors = FALSE
)
candidate <- candidate[!duplicated(candidate$feature_id), , drop = FALSE]
candidate <- candidate[
  order(-abs(candidate$logFC), candidate$selected_p_value, candidate$feature_id),
  ,
  drop = FALSE
]
candidate$p_value_tie_size <- ave(
  candidate$selected_p_value,
  candidate$selected_p_value,
  FUN = length
)
submitted <- utils::head(candidate, args$submit_n)
submitted$submission_rank <- seq_len(nrow(submitted))
if (nrow(submitted) < 2L) {
  write_node_audit(submitted, node_output)
  finish_semantic_empty(
    "insufficient_significant_proteins",
    paste0(
      "Only ", nrow(submitted), " protein(s) passed ", p_column, " ≤ ",
      args$pvalue_threshold, " and |log2FC| ≥ ", args$log2fc_threshold,
      "; at least 2 are required for STRING analysis."
    ),
    list(
      significant_count = nrow(submitted),
      p_column = p_column,
      pvalue_threshold = args$pvalue_threshold,
      log2fc_threshold = args$log2fc_threshold
    )
  )
}

string_version <- "12.0"
api_base <- Sys.getenv(
  "EVOMICS_STRING_API_BASE",
  unset = "https://version-12-0.string-db.org/api"
)
caller_identity <- Sys.getenv(
  "EVOMICS_STRING_CALLER_IDENTITY",
  unset = "EVOmicsDB"
)
message(
  "Submitting ", nrow(submitted), " of ", nrow(candidate),
  " significant protein identifiers to STRING v", string_version,
  "; taxid = ", organism_taxid,
  "; network type = ", args$network_type,
  "; score threshold = ", args$interaction_score_threshold, "."
)

mapping_response <- string_api_request(
  "get_string_ids",
  list(
    identifiers = paste(submitted$feature_id, collapse = "\r"),
    species = organism_taxid,
    limit = 1,
    echo_query = 1,
    caller_identity = caller_identity
  ),
  api_base,
  cache_dir = args$string_cache_dir,
  retries = args$api_retries
)
required_mapping_columns <- c("queryItem", "stringId", "preferredName")
if (nrow(mapping_response)) {
  missing_mapping_columns <- setdiff(required_mapping_columns, colnames(mapping_response))
  if (length(missing_mapping_columns)) {
    stop(
      "STRING identifier-mapping response is missing: ",
      paste(missing_mapping_columns, collapse = ", ")
    )
  }
  mapping_response <- mapping_response[
    !duplicated(tolower(trimws(mapping_response$queryItem))),
    ,
    drop = FALSE
  ]
} else {
  mapping_response <- data.frame(
    queryItem = character(), stringId = character(), preferredName = character(),
    stringsAsFactors = FALSE
  )
}

mapping_index <- match(
  tolower(trimws(submitted$feature_id)),
  tolower(trimws(mapping_response$queryItem))
)
submitted$STRING_id <- as.character(mapping_response$stringId[mapping_index])
submitted$preferred_name <- as.character(mapping_response$preferredName[mapping_index])
submitted$string_mapped <- !is.na(submitted$STRING_id) & nzchar(submitted$STRING_id)
submitted$string_id_collapsed <- FALSE
submitted$graph_representative <- FALSE

mapped_all <- submitted[submitted$string_mapped, , drop = FALSE]
if (nrow(mapped_all)) {
  mapped_all$string_id_collapsed <- duplicated(mapped_all$STRING_id)
  submitted$string_id_collapsed[match(mapped_all$feature_id, submitted$feature_id)] <-
    mapped_all$string_id_collapsed
  mapped <- mapped_all[!mapped_all$string_id_collapsed, , drop = FALSE]
  submitted$graph_representative[
    match(mapped$feature_id, submitted$feature_id)
  ] <- TRUE
  representative_by_string <- stats::setNames(mapped$feature_id, mapped$STRING_id)
  submitted$represented_by <- unname(representative_by_string[submitted$STRING_id])
} else {
  mapped <- mapped_all
  submitted$represented_by <- NA_character_
}
collapsed_synonym_count <- sum(submitted$string_id_collapsed)
unmapped_count <- sum(!submitted$string_mapped)

if (nrow(mapped) < 2L) {
  write_node_audit(submitted, node_output, string_version_value = string_version)
  finish_semantic_empty(
    "insufficient_string_mappings",
    paste0(
      nrow(submitted), " significant protein(s) were submitted, but only ",
      nrow(mapped), " unique STRING identifier(s) were mapped (", unmapped_count,
      " unmapped; ", collapsed_synonym_count,
      " synonym mapping(s) collapsed); at least 2 unique mappings are required."
    ),
    list(
      submitted_count = nrow(submitted), unique_mapped_count = nrow(mapped),
      unmapped_count = unmapped_count,
      collapsed_synonym_count = collapsed_synonym_count
    )
  )
}

message(
  "STRING mapping audit: ", nrow(submitted), " submitted; ", nrow(mapped),
  " unique mapped; ", unmapped_count, " unmapped; ",
  collapsed_synonym_count, " synonym mapping(s) collapsed."
)

network_response <- string_api_request(
  "network",
  list(
    identifiers = paste(mapped$STRING_id, collapse = "\r"),
    species = organism_taxid,
    required_score = round(args$interaction_score_threshold * 1000),
    network_type = args$network_type,
    caller_identity = caller_identity
  ),
  api_base,
  cache_dir = args$string_cache_dir,
  retries = args$api_retries
)
interactions <- if (nrow(network_response)) {
  required_network_columns <- c("stringId_A", "stringId_B", "score")
  missing_network_columns <- setdiff(
    required_network_columns,
    colnames(network_response)
  )
  if (length(missing_network_columns)) {
    stop(
      "STRING network response is missing: ",
      paste(missing_network_columns, collapse = ", ")
    )
  }
  data.frame(
    from = as.character(network_response$stringId_A),
    to = as.character(network_response$stringId_B),
    combined_score = round(as.numeric(network_response$score) * 1000),
    stringsAsFactors = FALSE
  )
} else {
  data.frame(from = character(), to = character(), combined_score = numeric())
}
if (nrow(interactions)) {
  interactions <- interactions[
    interactions$from %in% mapped$STRING_id &
      interactions$to %in% mapped$STRING_id &
      interactions$combined_score >= round(args$interaction_score_threshold * 1000),
    ,
    drop = FALSE
  ]
  canonical_edge <- apply(
    interactions[c("from", "to")],
    1L,
    function(pair) paste(sort(pair), collapse = "|")
  )
  interactions <- interactions[!duplicated(canonical_edge), , drop = FALSE]
}

edge_frame <- if (nrow(interactions)) {
  data.frame(
    from = interactions$from,
    to = interactions$to,
    combined_score = interactions$combined_score,
    stringsAsFactors = FALSE
  )
} else {
  data.frame(from = character(), to = character(), combined_score = numeric())
}

# Calculate hub degree on the complete submitted-and-mapped network.  The plot
# is an induced subnetwork selected only after this full-network calculation.
vertices_full <- data.frame(
  name = mapped$STRING_id,
  label = mapped$feature_id,
  preferred_name = mapped$preferred_name,
  logFC = mapped$logFC,
  regulation = mapped$regulation,
  stringsAsFactors = FALSE
)
network_full <- igraph::graph_from_data_frame(
  edge_frame, directed = FALSE, vertices = vertices_full
)
mapped$degree_full <- as.integer(
  igraph::degree(network_full)[match(mapped$STRING_id, igraph::V(network_full)$name)]
)
displayed_nodes <- mapped[
  order(
    -mapped$degree_full,
    -abs(mapped$logFC),
    mapped$selected_p_value,
    mapped$feature_id
  ),
  ,
  drop = FALSE
]
displayed_nodes <- utils::head(displayed_nodes, min(args$display_n, nrow(displayed_nodes)))
displayed_ids <- displayed_nodes$STRING_id
displayed_edges <- edge_frame[
  edge_frame$from %in% displayed_ids & edge_frame$to %in% displayed_ids,
  ,
  drop = FALSE
]
vertices_display <- data.frame(
  name = displayed_nodes$STRING_id,
  label = displayed_nodes$feature_id,
  preferred_name = displayed_nodes$preferred_name,
  logFC = displayed_nodes$logFC,
  regulation = displayed_nodes$regulation,
  degree_full = displayed_nodes$degree_full,
  stringsAsFactors = FALSE
)
network_display <- igraph::graph_from_data_frame(
  displayed_edges, directed = FALSE, vertices = vertices_display
)
displayed_interaction_count <- nrow(displayed_edges)
maximum_degree <- max(displayed_nodes$degree_full, na.rm = TRUE)
degree_breaks <- unique(round(pretty(c(0, maximum_degree), n = 4)))
degree_breaks <- degree_breaks[
  degree_breaks >= 0 & degree_breaks <= maximum_degree
]
if (!length(degree_breaks)) {
  degree_breaks <- 0
}

set.seed(123)
layout_name <- if (igraph::ecount(network_display) > 0L) "fr" else "circle"
label_size <- max(2.2, min(3.2, 3.2 * sqrt(30 / nrow(vertices_display))))
network_plot <- ggraph::ggraph(network_display, layout = layout_name)
if (igraph::ecount(network_display) > 0L) {
  network_plot <- network_plot +
    ggraph::geom_edge_link(
      ggplot2::aes(width = combined_score / 1000),
      colour = "#9E9E9E",
      alpha = 0.7,
      show.legend = TRUE
    ) +
    ggraph::scale_edge_width(range = c(0.4, 1.8), name = "STRING score")
}
network_title <- if (identical(args$network_type, "physical")) {
  "STRING Physical Interaction Network"
} else {
  "STRING Functional Association Network"
}
network_plot <- network_plot +
  ggraph::geom_node_point(
    ggplot2::aes(fill = logFC, size = degree_full),
    shape = 21,
    colour = "black",
    stroke = 0.45
  ) +
  ggraph::geom_node_text(
    ggplot2::aes(label = label),
    repel = TRUE,
    size = label_size,
    family = "sans",
    box.padding = 0.35,
    point.padding = 0.15,
    max.overlaps = Inf,
    seed = 123
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#4DBBD5",
    mid = "#F2F2F2",
    high = "#E64B35",
    midpoint = 0,
    name = paste0("log2FC (", contrast_label, ")")
  ) +
  ggplot2::scale_size_continuous(
    range = c(3.5, 8.5),
    breaks = degree_breaks,
    name = "Full-network degree"
  ) +
  ggplot2::labs(
    title = network_title,
    subtitle = paste0(
      "Top ", nrow(displayed_nodes), " hub proteins among ", nrow(submitted),
      " submitted DE proteins (", nrow(candidate), " passed thresholds); ",
      displayed_interaction_count, " of ", nrow(edge_frame),
      " full-network interactions shown; STRING v",
      string_version, " score >= ",
      format(args$interaction_score_threshold, trim = TRUE)
    )
  ) +
  ggplot2::theme_void(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
    legend.position = "right",
    plot.margin = ggplot2::margin(12, 18, 12, 18)
  )

ensure_parent_directory(args$output_png)
ggplot2::ggsave(
  args$output_png,
  network_plot,
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
    network_plot,
    width = args$width / 300,
    height = args$height / 300,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}

from_index <- match(edge_frame$from, mapped$STRING_id)
to_index <- match(edge_frame$to, mapped$STRING_id)
n_edges <- nrow(edge_frame)
edge_export <- data.frame(
  from_symbol = mapped$feature_id[from_index],
  to_symbol = mapped$feature_id[to_index],
  from_STRING_id = edge_frame$from,
  to_STRING_id = edge_frame$to,
  combined_score = edge_frame$combined_score,
  score_fraction = edge_frame$combined_score / 1000,
  from_logFC = mapped$logFC[from_index],
  to_logFC = mapped$logFC[to_index],
  in_displayed_subnetwork = edge_frame$from %in% displayed_ids &
    edge_frame$to %in% displayed_ids,
  organism_taxid = rep(organism_taxid, n_edges),
  string_version = rep(string_version, n_edges),
  network_type = rep(args$network_type, n_edges),
  score_threshold = rep(args$interaction_score_threshold, n_edges),
  stringsAsFactors = FALSE
)
ensure_parent_directory(args$output_csv)
utils::write.csv(edge_export, args$output_csv, row.names = FALSE, quote = TRUE)

ensure_parent_directory(node_output)
degree_by_string <- stats::setNames(mapped$degree_full, mapped$STRING_id)
submitted$degree_full <- as.integer(unname(degree_by_string[submitted$STRING_id]))
submitted$displayed <- submitted$graph_representative & submitted$STRING_id %in% displayed_ids
submitted$organism_taxid <- organism_taxid
submitted$string_version <- string_version
submitted$network_type <- args$network_type
submitted$score_threshold <- args$interaction_score_threshold
write_node_audit(
  submitted,
  node_output,
  degree_lookup = degree_by_string,
  displayed_ids = displayed_ids,
  string_version_value = string_version
)

status_details <- list(
  significant_count = nrow(candidate),
  submitted_count = nrow(submitted),
  unique_mapped_count = nrow(mapped),
  unmapped_count = unmapped_count,
  collapsed_synonym_count = collapsed_synonym_count,
  displayed_count = nrow(displayed_nodes),
  displayed_interaction_count = displayed_interaction_count,
  retained_interaction_count = nrow(edge_export),
  p_value_tie_protein_count = sum(submitted$p_value_tie_size > 1L),
  score_threshold = args$interaction_score_threshold,
  network_type = args$network_type,
  string_version = string_version,
  cache_policy = "Canonical-request-hashed, version-pinned STRING v12.0 response cache"
)

if (igraph::ecount(network_full) == 0L) {
  no_edge_message <- paste0(
    nrow(mapped), " unique STRING-mapped protein(s) were analysed, but no ",
    args$network_type, " interaction passed score ≥ ",
    format(args$interaction_score_threshold, trim = TRUE), " (", unmapped_count,
    " submitted identifier(s) unmapped; ", collapsed_synonym_count,
    " synonym mapping(s) collapsed)."
  )
  write_status("no_interactions", no_edge_message, status_details)
  message(no_edge_message)
} else {
  success_message <- paste0(
    "PPI completed: ", nrow(submitted), " submitted; ", nrow(mapped),
    " unique mapped; ", unmapped_count, " unmapped; ",
    collapsed_synonym_count, " synonym mapping(s) collapsed; ",
    nrow(edge_export), " retained full-network interactions; ",
    displayed_interaction_count, " shown in the displayed subnetwork; ",
    nrow(displayed_nodes), " displayed hub proteins."
  )
  write_status("success", success_message, status_details)
  message(success_message)
}

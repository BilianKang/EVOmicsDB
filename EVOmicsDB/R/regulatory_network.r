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

# ============================================================
# regulatory_network.r
# miRNA调控网络 (miRNA → mRNA / Protein)
# 支持组合: miRNA+mRNA, miRNA+Protein
# ============================================================

library(dplyr)
library(scales)
library(jsonlite)
library(argparse)

parser <- ArgumentParser(description = "Build miRNA regulatory network (miRNA -> mRNA/Protein)")
parser$add_argument("--mirna_dataset", type = "character", required = TRUE, help = "Path to miRNA dataset file")
parser$add_argument("--target_dataset", type = "character", required = TRUE, help = "Path to target dataset file (mRNA/Protein)")
parser$add_argument("--target_group", type = "character", default = "mRNA", help = "Target type: mRNA or Protein")
# Use a numeric sentinel rather than R's NA here.  The argparse R package
# translates defaults through its Python parser, where an unquoted `NA` would
# be evaluated as an undefined Python name and abort every network request.
parser$add_argument("--dataset_id", type = "integer", default = -1,
                    help = "Numeric miRNA dataset ID for audit provenance")
parser$add_argument("--mirna_resolution", type = "character", default = NULL,
                    choices = c("mirbase_mature", "mirbase_precursor"),
                    help = "Configured measurement resolution; never inferred from ID strings")
parser$add_argument("--strict_arm_resolution", type = "logical", default = FALSE,
                    help = "Exclude precursor edges from hairpins with both mature arms")
parser$add_argument("--max_targets_per_mirna", type = "integer", default = 15, help = "Max targets per miRNA")
parser$add_argument("--output_html", type = "character", required = TRUE, help = "Output HTML file path")
parser$add_argument("--output_csv", type = "character", required = TRUE, help = "Output CSV file path")
# miRNA p-value parameters
parser$add_argument("--mirna_is_use_padj", type = "logical", default = TRUE, help = "Whether to use adjusted p-value for miRNA")
parser$add_argument("--mirna_pvalue", type = "double", default = 0.05, help = "P-value threshold for miRNA")
# Target p-value parameters
parser$add_argument("--target_is_use_padj", type = "logical", default = TRUE, help = "Whether to use adjusted p-value for target")
parser$add_argument("--target_pvalue", type = "double", default = 0.05, help = "P-value threshold for target")
parser$add_argument("--mirna_log2fc", type = "double", default = 0)
parser$add_argument("--target_log2fc", type = "double", default = 0)
parser$add_argument("--evidence", type = "character", default = "predicted",
                    choices = c("validated", "predicted", "all"),
                    help = "Local target-map evidence mode")
parser$add_argument("--minimum_mirdb_score", type = "double", default = 50,
                    help = "Minimum local miRDB prediction score (0-100)")
parser$add_argument("--require_inverse_direction", type = "logical", default = FALSE,
                    help = "Keep miRNA-target pairs with opposite logFC signs")
parser$add_argument("--multimir_limit", type = "integer", default = 100000)
parser$add_argument("--target_map", type = "character", default = NULL,
                    help = "Optional local CSV/TSV target map; avoids live multiMiR access")
parser$add_argument("--mirna_column", type = "character", default = "miRNA")
parser$add_argument("--target_column", type = "character", default = "target_gene")
parser$add_argument("--score_column", type = "character", default = "score")
parser$add_argument("--output_png", type = "character", default = NULL)
parser$add_argument("--output_pdf", type = "character", default = NULL)
parser$add_argument("--output_audit_csv", type = "character", default = NULL,
                    help = "Optional matching audit TSV output")
parser$add_argument("--output_feature_audit_csv", type = "character", default = NULL,
                    help = "Optional feature-level miRNA resolution audit TSV output")
parser$add_argument("--status_json", type = "character", default = NULL,
                    help = "Structured Explore execution status JSON")
# Node shape parameters (visNetwork shapes: dot, square, triangle, triangleDown, diamond, star, ellipse, box, circle, database)
parser$add_argument("--mirna_shape", type = "character", default = "diamond", help = "miRNA node shape")
parser$add_argument("--target_shape", type = "character", default = "dot", help = "Target node shape")
# logFC color mapping parameters
parser$add_argument("--low_color", type = "character", default = "#1f78b4", help = "Color for low logFC (downregulated)")
parser$add_argument("--mid_color", type = "character", default = "white", help = "Color for zero logFC")
parser$add_argument("--high_color", type = "character", default = "#e31a1c", help = "Color for high logFC (upregulated)")
args <- parser$parse_args()
if (!is.null(args$status_json) && file.exists(args$status_json)) unlink(args$status_json)
if (is.null(args$mirna_resolution) || !nzchar(args$mirna_resolution)) {
  stop("miRNA regulatory analysis requires an explicit --mirna_resolution from dataset configuration")
}
write_regulatory_status <- function(status, code, message_text, counts = list(), artifacts = list()) {
  evomics_status_write(args$status_json, "regulatory_network", status, code,
                       message_text, counts, artifacts)
}

write_empty_network_result <- function(reason, code = "no_interactions") {
  escaped <- gsub("&", "&amp;", as.character(reason), fixed = TRUE)
  escaped <- gsub("<", "&lt;", escaped, fixed = TRUE)
  escaped <- gsub(">", "&gt;", escaped, fixed = TRUE)
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<style>body{margin:0;font-family:Arial,sans-serif;background:#fff;color:#344054}",
    ".empty{max-width:760px;margin:80px auto;padding:32px;border:1px solid #d0d5dd;",
    "border-radius:10px;text-align:center;background:#f9fafb}.empty h2{color:#0b447c}",
    "</style></head><body><div class='empty'><h2>No regulatory interactions</h2><p>",
    escaped, "</p><p>Try a lower <b>miRDB score</b>, disable <b>Inverse direction only</b>, ",
    "or use a less restrictive significance threshold.</p></div></body></html>"
  )
  writeLines(html, args$output_html, useBytes = TRUE)
  empty_edges <- data.frame(
    miRNA_id = character(), miRNA_resolution = character(),
    measured_mirna_id = character(), canonical_hairpin_id = character(),
    supporting_arm = character(), arm_ambiguous = logical(),
    mature_mirna_id = character(), target_symbol = character(),
    source_type = character(), target_type = character(),
    regulatory_direction = character(),
    score = numeric(), support_count = numeric(), evidence_query = character()
  )
  write.table(
    empty_edges, file = args$output_csv, sep = "\t", quote = FALSE,
    row.names = FALSE
  )
  cat("  Empty network:", reason, "\n")
  write_regulatory_status(
    "empty", code, reason,
    counts = list(mirna_count = 0, target_count = 0, edge_count = 0),
    artifacts = list(html = args$output_html, csv = args$output_csv,
                     matching_audit = args$output_audit_csv,
                     feature_matching_audit = args$output_feature_audit_csv)
  )
  cat(jsonlite::toJSON(
    list(mirna_count = 0, target_count = 0, edge_count = 0, reason = reason,
         mirna_resolution = args$mirna_resolution),
    auto_unbox = TRUE
  ), "\n")
  quit(save = "no", status = 0)
}

# ======================== 辅助函数 ========================

get_significant <- function(df, pvalue_threshold = 0.05, is_use_padj = TRUE, context = "input") {
  requested <- if (isTRUE(is_use_padj)) "adj.P.Val" else "P.Value"
  if (!requested %in% colnames(df)) {
    stop("Requested ", requested, " is missing in ", context, "; no significance fallback is permitted.")
  }
  df %>% dplyr::filter(is.finite(.data[[requested]]) & .data[[requested]] <= pvalue_threshold)
}

get_mirna_targets <- function(mirna_list, gene_list = NULL, limit = 100000,
                              evidence = "validated") {
  if (!requireNamespace("multiMiR", quietly = TRUE)) {
    stop("multiMiR is required when --target_map is not supplied")
  }
  mirna_ids <- mirna_list$ID
  target_ids <- if (!is.null(gene_list)) gene_list$ID else NULL
  
  results <- multiMiR::get_multimir(
    mirna   = mirna_ids,
    target  = target_ids,
    table   = evidence,
    summary = TRUE,
    limit   = limit
  )
  return(results)
}

extract_mirna_targets_summary <- function(mirna_targets, evidence) {
  if (nrow(mirna_targets@summary) > 0) {
    sm <- as.data.frame(mirna_targets@summary)
    keep <- intersect(c("mature_mirna_id", "target_entrez", "target_symbol",
                        "validated", "predicted", "disease.drug"), colnames(sm))
    associations <- sm[, keep, drop = FALSE]
    if ("mature_mirna_id" %in% names(associations)) {
      names(associations)[names(associations) == "mature_mirna_id"] <- "miRNA"
    }
    if ("target_symbol" %in% names(associations)) {
      names(associations)[names(associations) == "target_symbol"] <- "target_gene"
    }
    associations$evidence_query <- evidence
    support_cols <- intersect(c("validated", "predicted"), colnames(associations))
    associations$support_count <- if (length(support_cols)) {
      rowSums(sapply(associations[support_cols], function(x) suppressWarnings(as.numeric(x))),
              na.rm = TRUE)
    } else 1
    return(associations)
  } else {
    # A successful target query with no rows is a valid biological empty
    # result, not a technical failure.  Keep a typed zero-row table so the
    # downstream resolution/audit code can emit the normal no_interactions
    # status.  Schema/package/database failures still stop before reaching
    # this branch.
    data.frame(
      miRNA = character(), target_gene = character(),
      evidence_query = character(), support_count = numeric(),
      stringsAsFactors = FALSE
    )
  }
}

filter_mirna_associations_for_network <- function(
    associations, mirna_sig, target_sig, resolution
) {
  mirna_ids  <- unique(resolve_mirna_feature_key(mirna_sig$ID, resolution))
  target_ids <- unique(toupper(sub("\\.[0-9]+$", "", target_sig$ID)))
  if (resolution == "mirbase_precursor") {
    associations <- resolve_mirna_target_annotations(
      associations, mirna_sig$ID, resolution
    )
    associations <- associations[associations$hairpin_id %in% mirna_ids, , drop = FALSE]
  } else {
    associations <- resolve_mirna_target_annotations(
      associations, mirna_sig$ID, resolution
    )
    associations <- associations[associations$miRNA_match_key %in% mirna_ids, , drop = FALSE]
  }
  df <- associations %>%
    dplyr::mutate(target_symbol = toupper(target_gene)) %>%
    dplyr::filter(target_symbol %in% target_ids)
  return(df)
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

cat("Building regulatory network: miRNA ->", args$target_group, "\n")
cat("  Configured miRNA resolution:", args$mirna_resolution, "\n")

# 读取数据 (从rda文件)
mirna_df <- read_deg_from_rda(args$mirna_dataset)
target_df <- read_deg_from_rda(args$target_dataset)
# Keep the network's significance-filtered DEG table separate from the
# annotation-coverage universe.  The latter must represent every measured
# miRNA, otherwise the audit changes when a p-value threshold changes.
mirna_feature_universe <- evomics_read_measured_universe_from_rda(args$mirna_dataset)
if (!length(mirna_feature_universe)) {
  stop("Could not recover the complete measured miRNA feature universe from the canonical RDA")
}

# 筛选显著
mirna_sig <- get_significant(mirna_df, args$mirna_pvalue, args$mirna_is_use_padj, "miRNA differential input")
target_sig <- get_significant(target_df, args$target_pvalue, args$target_is_use_padj, "target differential input")
mirna_sig <- mirna_sig %>% dplyr::filter(abs(logFC) >= args$mirna_log2fc)
target_sig <- target_sig %>% dplyr::filter(abs(logFC) >= args$target_log2fc)

cat("  Significant miRNAs:", nrow(mirna_sig), "\n")
cat("  Significant targets:", nrow(target_sig), "\n")

if (!is.null(args$target_map) && nzchar(trimws(args$target_map))) {
  cat("  Reading local miRNA target map:", args$target_map, "\n")
  map_query_ids <- if (args$mirna_resolution == "mirbase_precursor") {
    hairpins <- unique(canonical_hairpin_id(mirna_feature_universe))
    unique(c(hairpins, paste0(hairpins, "-3p"), paste0(hairpins, "-5p")))
  } else {
    paste0("hsa-", unique(normalize_mirna_exact(mirna_feature_universe)))
  }
  local_map <- evomics_read_target_map(
    args$target_map, args$mirna_column, args$target_column,
    evidence = args$evidence, minimum_score = args$minimum_mirdb_score,
    mirna_ids = map_query_ids
  )
  support <- rep(1, nrow(local_map))
  mirna_associations_all <- data.frame(
    miRNA = as.character(local_map$miRNA),
    target_gene = as.character(local_map$target_gene),
    score = suppressWarnings(as.numeric(local_map$score)),
    support_count = support,
    evidence_query = rep("local target map", nrow(local_map)),
    stringsAsFactors = FALSE
  )
} else {
  if (!nrow(mirna_sig)) {
    mirna_associations_all <- data.frame(miRNA = character(), target_gene = character(),
      score = numeric(), support_count = numeric(), evidence_query = character(),
      stringsAsFactors = FALSE)
  } else {
    cat("  Using multiMiR to retrieve targets...\n")
    mirna_targets <- get_mirna_targets(mirna_sig, target_sig,
                                       limit = args$multimir_limit,
                                       evidence = args$evidence)
    mirna_associations_all <- extract_mirna_targets_summary(mirna_targets, args$evidence)
  }
}

# 过滤
mirna_associations <- filter_mirna_associations_for_network(
  associations           = mirna_associations_all,
  mirna_sig              = mirna_sig,
  target_sig             = target_sig,
  resolution             = args$mirna_resolution
)

# Audit annotation coverage independently of significance thresholds.  This is
# deliberately calculated from the measured miRNA feature list and the full
# local library, so precursor statistics are never copied to mature arms.
matching_audit <- mirna_matching_audit(
  mirna_feature_universe, mirna_associations_all, args$mirna_resolution
)
matching_audit$dataset_id <- ifelse(is.na(args$dataset_id) || args$dataset_id < 0,
                                    basename(args$mirna_dataset),
                                    as.character(args$dataset_id))
matching_audit$mirna_resolution <- args$mirna_resolution
matching_audit <- matching_audit[, c("dataset_id", "mirna_resolution",
  "total_mirna_features", "exact_unarmed_matches",
  "unique_single_arm_hairpins", "dual_arm_hairpins",
  "no_target_library_match", "matched_total", "match_rate",
  "strict_eligible_features", "strict_eligible_rate")]
if (!is.null(args$output_audit_csv) && nzchar(args$output_audit_csv)) {
  write.table(matching_audit, file = args$output_audit_csv, sep = "\t",
              quote = FALSE, row.names = FALSE)
}
feature_audit <- mirna_feature_matching_audit(
  mirna_feature_universe, mirna_associations_all, args$mirna_resolution,
  dataset_id = if (args$dataset_id < 0) basename(args$mirna_dataset) else as.character(args$dataset_id)
)
if (!is.null(args$output_feature_audit_csv) && nzchar(args$output_feature_audit_csv)) {
  write.table(feature_audit, file = args$output_feature_audit_csv, sep = "\t",
              quote = FALSE, row.names = FALSE)
}
cat("  Matching audit:", paste(names(matching_audit), matching_audit[1, ],
  sep = "=", collapse = "; "), "\n")

if (nrow(mirna_sig) == 0 || nrow(target_sig) == 0) {
  write_empty_network_result("No significant miRNA or target features passed the selected filters.",
                             code = "no_significant_features")
}

if (args$mirna_resolution == "mirbase_precursor" && isTRUE(args$strict_arm_resolution)) {
  before_strict <- nrow(mirna_associations)
  mirna_associations <- mirna_associations %>% dplyr::filter(!arm_ambiguous)
  cat("  Strict arm resolution removed", before_strict - nrow(mirna_associations),
      "ambiguous precursor-target edges.\n")
}

if (isTRUE(args$require_inverse_direction) && nrow(mirna_associations) > 0) {
  mirna_fc <- mirna_sig %>% transmute(
    mirna_key = resolve_mirna_feature_key(ID, args$mirna_resolution), mirna_logFC = logFC
  )
  target_fc <- target_sig %>% transmute(
    target_key = toupper(sub("\\.[0-9]+$", "", ID)), target_logFC = logFC
  )
  mirna_associations <- mirna_associations %>%
    mutate(mirna_key = if (args$mirna_resolution == "mirbase_precursor") {
             hairpin_id
           } else {
             miRNA_match_key
           },
           target_key = toupper(target_symbol)) %>%
    inner_join(mirna_fc, by = "mirna_key") %>%
    inner_join(target_fc, by = "target_key") %>%
    filter(mirna_logFC * target_logFC < 0)
}

# Apply the display cap last, so significance, mapping-score and optional
# inverse-direction filters are all respected before selecting the top N.
if (nrow(mirna_associations) > 0) {
  if (!"score" %in% names(mirna_associations)) {
    mirna_associations$score <- NA_real_
  }
  if (args$evidence == "validated" && !any(is.finite(mirna_associations$score))) {
    message("Validated target map has no comparable score; max_targets_per_mirna is a display cap and ties use support_count then target ID.")
  }
  mirna_associations <- mirna_associations %>%
    dplyr::arrange(dplyr::desc(score), dplyr::desc(support_count), target_gene) %>%
    dplyr::group_by(miRNA) %>%
    dplyr::slice_head(n = args$max_targets_per_mirna) %>%
    dplyr::ungroup()
}

if (nrow(mirna_associations) == 0) {
  write_empty_network_result("No miRNA-target associations passed the selected evidence and significance filters.")
}

# Use one canonical identifier contract for both graph edges and nodes.  The
# display label remains the original matrix ID, but visNetwork joins on these
# normalized IDs.  This is essential because local target maps intentionally
# normalize e.g. ``hsa-miR-21`` to ``mir-21`` and gene versions to symbols.
normalize_mirna_id <- function(x) {
  tolower(sub("^hsa-", "", trimws(as.character(x))))
}
normalize_target_id <- function(x) {
  toupper(sub("\\.[0-9]+$", "", trimws(as.character(x))))
}
mirna_associations <- mirna_associations %>%
  dplyr::mutate(
    annotation_mirna_id = as.character(miRNA),
    miRNA_id = vapply(miRNA, function(x) {
      hit <- which(resolve_mirna_feature_key(mirna_sig$ID, args$mirna_resolution) ==
                   resolve_mirna_feature_key(x, args$mirna_resolution))
      if (length(hit)) as.character(mirna_sig$ID[hit[[1]]]) else as.character(x)
    }, character(1)),
    measured_mirna_id = normalize_mirna_id(miRNA_id),
    canonical_hairpin_id = ifelse(args$mirna_resolution == "mirbase_precursor",
                                  canonical_hairpin_id(miRNA_id), NA_character_),
    mature_mirna_id = ifelse(args$mirna_resolution == "mirbase_mature",
                             normalize_mirna_id(miRNA_id), NA_character_),
    target_symbol = normalize_target_id(target_gene),
    miRNA_resolution = args$mirna_resolution,
    resolution_mode = ifelse(args$mirna_resolution == "mirbase_precursor",
                              "hairpin_collapsed", "mature_exact"),
    source_type = "miRNA",
    target_type = args$target_group,
    regulatory_direction = paste0("miRNA_to_", tolower(args$target_group))
  ) %>%
  dplyr::filter(nzchar(measured_mirna_id), nzchar(target_symbol)) %>%
  dplyr::distinct(miRNA_id, target_symbol, .keep_all = TRUE)

cat("  Filtered associations:", nrow(mirna_associations), "\n")

# 颜色映射 (基于logFC)
range_all <- range(c(mirna_sig$logFC, target_sig$logFC), na.rm = TRUE)

col_fun <- scales::col_numeric(
  palette = c(args$low_color, args$mid_color, args$high_color),
  domain  = range_all
)

# miRNA节点
mirna_nodes_all <- mirna_sig %>%
  dplyr::select(ID, logFC) %>%
  dplyr::mutate(
    id    = normalize_mirna_id(ID),
    label = ID,
    group = "miRNA",
    color = col_fun(logFC),
    miRNA_resolution = args$mirna_resolution,
    title = paste0("miRNA: ", ID, "<br>logFC: ", round(logFC, 3),
      ifelse(args$mirna_resolution == "mirbase_precursor",
        "<br>Resolution: precursor/hairpin", "")),
    value = abs(logFC) * 5
  )

# Target节点
target_nodes_all <- target_sig %>%
  dplyr::select(ID, logFC) %>%
  dplyr::mutate(
    id    = normalize_target_id(ID),
    label = ID,
    group = args$target_group,
    color = col_fun(logFC),
    title = paste0(args$target_group, ": ", ID, "<br>logFC: ", round(logFC, 3)),
    value = abs(logFC) * 5
  )

# 构建边
edges <- mirna_associations %>%
  dplyr::mutate(
    from   = measured_mirna_id,
    to     = target_symbol,
    color  = "#425b76",
    width  = 2.4,
    title  = paste0("miRNA ", miRNA_id, " regulates ",
                    args$target_group, " ", target_symbol),
    arrows = "to"
  ) %>%
dplyr::select(from, to, color, width, title, arrows)

# 去掉孤立节点
used_ids <- unique(c(edges$from, edges$to))
mirna_nodes  <- mirna_nodes_all  %>% dplyr::filter(id %in% used_ids)
target_nodes <- target_nodes_all %>% dplyr::filter(id %in% used_ids)

all_nodes <- dplyr::bind_rows(mirna_nodes, target_nodes) %>%
  dplyr::distinct(id, .keep_all = TRUE)

if (nrow(all_nodes) == 0) {
  write_empty_network_result("All filtered nodes were isolated after identifier matching.")
}

cat("  Nodes:", nrow(all_nodes), ", Edges:", nrow(edges), "\n")

build_static_network_plot <- function() {
  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("ggraph", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("igraph, ggraph, and ggplot2 are required for the network fallback.")
  }
  graph_static <- igraph::graph_from_data_frame(
    edges[, c("from", "to")], directed = TRUE,
    vertices = all_nodes %>% transmute(
      name = id, label = label, group = group, color = color,
      node_size = pmax(2.5, sqrt(pmax(value, 1)))
    )
  )
  set.seed(123)
  ggraph::ggraph(graph_static, layout = "fr") +
    ggraph::geom_edge_link(
      # Preserve the original light, open-V arrow geometry while retaining
      # the higher-contrast colour that makes direction easier to inspect.
      arrow = grid::arrow(length = grid::unit(2.5, "mm")),
      end_cap = ggraph::circle(2.5, "mm"),
      colour = "#425b76", linewidth = 0.45, alpha = 0.78
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(size = node_size, fill = I(color), shape = group),
      colour = "black", stroke = 0.35
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = label), repel = TRUE, size = 3.2
    ) +
    ggplot2::scale_shape_manual(
      values = setNames(c(23, 21), c("miRNA", args$target_group))
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::labs(
      title = paste("miRNA–", args$target_group, "regulatory network"),
      subtitle = paste0(
        "Direction: miRNA (diamond) → ", args$target_group,
        " (circle). Node colour denotes C−N log2FC; ", nrow(all_nodes),
        " nodes and ", nrow(edges), " edges",
        ifelse(args$mirna_resolution == "mirbase_precursor",
          ". Precursor-level exploratory network: mature-miRNA annotations are collapsed to the measured hairpin; no mature-arm expression is inferred.", "")
      ),
      shape = "Node type",
      caption = paste0("Arrowheads point to regulated ", args$target_group, " target nodes.",
        ifelse(args$mirna_resolution == "mirbase_precursor",
          " Precursor abundance should not be interpreted as mature-miRNA activity.", ""))
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      plot.caption = ggplot2::element_text(face = "bold", hjust = 0.5, colour = "#425b76"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position = "bottom"
    )
}

# Prefer the interactive widget when both optional runtime dependencies are
# available.  Otherwise emit a fully self-contained static HTML result.  This
# keeps local and migrated deployments functional without visNetwork/pandoc.
has_interactive_runtime <- requireNamespace("visNetwork", quietly = TRUE) &&
  nzchar(Sys.which("pandoc"))
p_static <- NULL
if (has_interactive_runtime) {
  network <- visNetwork::visNetwork(all_nodes, edges) %>%
    visNetwork::visNodes(shape = "dot", font = list(size = 14, face = "arial")) %>%
    visNetwork::visEdges(
      arrows = "to", arrowStrikethrough = FALSE,
      color = list(color = "#425b76", highlight = "#0b5ea8"),
      width = 2.4, smooth = list(enabled = TRUE, type = "dynamic")
    ) %>%
    visNetwork::visGroups(groupname = "miRNA", shape = args$mirna_shape, size = 25) %>%
    visNetwork::visGroups(
      groupname = args$target_group, shape = args$target_shape, size = 20
    ) %>%
    visNetwork::visLegend() %>%
    visNetwork::visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      nodesIdSelection = TRUE, selectedBy = "group"
    ) %>%
    visNetwork::visLayout(randomSeed = 123) %>%
    visNetwork::visPhysics(stabilization = TRUE)
  visNetwork::visSave(network, file = args$output_html, selfcontained = TRUE)
} else {
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("base64enc is required for self-contained static network HTML.")
  }
  message("visNetwork and/or pandoc unavailable; using static HTML fallback.")
  p_static <- build_static_network_plot()
  fallback_png <- tempfile(fileext = ".png")
  on.exit(unlink(fallback_png), add = TRUE)
  ggplot2::ggsave(
    fallback_png, p_static, width = 10, height = 7.5, dpi = 180, bg = "white"
  )
  image_uri <- base64enc::dataURI(file = fallback_png, mime = "image/png")
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<style>html,body{margin:0;background:#fff;font-family:Arial,sans-serif}",
    ".wrap{padding:12px;text-align:center}.wrap img{display:block;width:100%;",
    "max-width:1800px;height:auto;margin:auto}</style></head><body>",
    "<div class='wrap'><img alt='miRNA regulatory network' src='", image_uri,
    "'></div></body></html>"
  )
  writeLines(html, args$output_html, useBytes = TRUE)
}

write.table(
  mirna_associations,
  file      = args$output_csv,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

if (!is.null(args$output_png) || !is.null(args$output_pdf)) {
  if (is.null(p_static)) p_static <- build_static_network_plot()
  if (!is.null(args$output_png)) {
    ggplot2::ggsave(args$output_png, p_static, width = 9, height = 7,
                    dpi = 300, bg = "white")
  }
  if (!is.null(args$output_pdf)) {
    ggplot2::ggsave(args$output_pdf, p_static, width = 9, height = 7,
                    device = grDevices::cairo_pdf, bg = "white")
  }
}

cat("  Output HTML:", args$output_html, "\n")
cat("  Output CSV:", args$output_csv, "\n")

# 输出统计信息
result <- list(
  mirna_count = nrow(mirna_nodes),
  target_count = nrow(target_nodes),
  edge_count = nrow(edges),
  mirna_resolution = args$mirna_resolution,
  strict_arm_resolution = isTRUE(args$strict_arm_resolution),
  matching_audit = matching_audit[1, , drop = FALSE]
)

cat(toJSON(result, auto_unbox = TRUE), "\n")
write_regulatory_status(
  "success", "completed", "miRNA regulatory network generated successfully",
  counts = result,
  artifacts = list(html = args$output_html, csv = args$output_csv,
                   matching_audit = args$output_audit_csv,
                   feature_matching_audit = args$output_feature_audit_csv,
                   png = args$output_png, pdf = args$output_pdf)
)

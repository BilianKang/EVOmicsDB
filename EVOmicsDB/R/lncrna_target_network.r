#!/usr/bin/env Rscript

# ============================================================
# lncrna_target_network.r
# lncRNA靶点交集网络
# 支持组合: lncRNA+mRNA, lncRNA+miRNA, lncRNA+Protein
# 基于本地 external_targets.sqlite（旧 master CSV 兼容）
# ============================================================

library(dplyr)
library(scales)
library(visNetwork)
library(jsonlite)
library(stringr)
library(argparse)
library(igraph)

# --- 获取脚本路径 ---
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grep(file_arg, args)])
  return(normalizePath(script_path))
}

script_dir <- dirname(get_script_path())
source(file.path(script_dir, "evomics_explore_utils.R"))
source(file.path(script_dir, "mirna_resolution_utils.R"))

# --- 命令行参数 ---
parser <- ArgumentParser(description = "Build lncRNA target intersection network")
parser$add_argument("--lncrna_dataset", type = "character", required = TRUE, help = "Path to lncRNA dataset file")
parser$add_argument("--target_dataset", type = "character", required = TRUE, help = "Path to target dataset file (mRNA/miRNA/Protein)")
parser$add_argument("--target_type", type = "character", required = TRUE, help = "Target type: mRNA/miRNA/Protein")
parser$add_argument("--mirna_resolution", type = "character", default = NULL,
                    choices = c("mirbase_mature", "mirbase_precursor"))
parser$add_argument("--min_degree", type = "integer", default = 0, help = "Display-only minimum node degree; 0 retains all")
parser$add_argument("--keep_target_classes", type = "character", default = "miRNA,gene_protein", help = "Target classes to keep (comma-separated)")
parser$add_argument("--output_html", type = "character", required = TRUE, help = "Output HTML file path")
parser$add_argument("--output_csv", type = "character", required = TRUE, help = "Output CSV file path")
# lncRNA p-value parameters
parser$add_argument("--lncrna_is_use_padj", type = "logical", default = TRUE, help = "Whether to use adjusted p-value for lncRNA")
parser$add_argument("--lncrna_pvalue", type = "double", default = 0.05, help = "P-value threshold for lncRNA")
# Target p-value parameters
parser$add_argument("--target_is_use_padj", type = "logical", default = TRUE, help = "Whether to use adjusted p-value for target")
parser$add_argument("--target_pvalue", type = "double", default = 0.05, help = "P-value threshold for target")
parser$add_argument("--lncrna_log2fc", type = "double", default = 0)
parser$add_argument("--target_log2fc", type = "double", default = 0)
parser$add_argument("--direction_rule", type = "character", default = "none",
                    choices = c("none", "same", "opposite"))
parser$add_argument("--minimum_evidence_score", type = "double", default = 0)
# Titles
parser$add_argument("--main_title", type = "character", default = "lncRNA Target Network", help = "Main title")
parser$add_argument("--main_title_size", type = "integer", default = 14, help = "Main title font size")
parser$add_argument("--axis_text_size", type = "integer", default = 10, help = "Node label font size")
# logFC color mapping parameters
parser$add_argument("--low_color", type = "character", default = "#1f78b4", help = "Color for low logFC (downregulated)")
parser$add_argument("--mid_color", type = "character", default = "white", help = "Color for zero logFC")
parser$add_argument("--high_color", type = "character", default = "#e31a1c", help = "Color for high logFC (upregulated)")
# Node shapes
parser$add_argument("--lncrna_shape", type = "character", default = "dot", help = "Shape for lncRNA nodes")
parser$add_argument("--target_shape", type = "character", default = "diamond", help = "Shape for target nodes (mRNA/miRNA/Protein)")
# Legend
parser$add_argument("--is_show_legend", type = "character", default = "TRUE", help = "Show legend (TRUE/FALSE)")
parser$add_argument("--legend_title_size", type = "integer", default = 10, help = "Legend title font size")
parser$add_argument("--legend_position", type = "character", default = "right", help = "Legend position")
# Export settings
parser$add_argument("--width", type = "integer", default = 1200, help = "Output width (px)")
parser$add_argument("--height", type = "integer", default = 800, help = "Output height (px)")
parser$add_argument("--master_lnc_target_db", type = "character", default = NULL,
                    help = "Normalized external_targets.sqlite or legacy CSV/TSV")
parser$add_argument("--evidence", type = "character", default = "all",
                    choices = c("validated", "experimental", "clip", "predicted", "all"),
                    help = "Evidence tier to include")
parser$add_argument("--source_databases", type = "character", default = NULL,
                    help = "Comma-separated source databases, or all")
parser$add_argument("--species", type = "character", default = "hsa",
                    help = "Species filter; production default is human (hsa/TaxID 9606)")
parser$add_argument("--include_derived_cerna", type = "logical", default = FALSE,
                    help = "Include lncRNA-mRNA edges derived only from ceRNA triplets")
parser$add_argument("--status_json", type = "character", default = NULL,
                    help = "Structured Explore execution status JSON")
args <- parser$parse_args()
if (!is.null(args$status_json) && file.exists(args$status_json)) unlink(args$status_json)

write_lnc_status <- function(status, code, message_text, counts = list(), artifacts = list()) {
  evomics_status_write(args$status_json, "lncrna_target_network", status, code,
                       message_text, counts, artifacts)
}

master_lnc_targets_path <- args$master_lnc_target_db
if (is.null(master_lnc_targets_path) || !nzchar(trimws(master_lnc_targets_path))) {
  master_lnc_targets_path <- evomics_default_db_path(script_dir, "external_targets.sqlite")
  if (is.null(master_lnc_targets_path)) {
    master_lnc_targets_path <- evomics_default_db_path(script_dir, "master_lncRNA_target_database.csv")
  }
}

# Parse keep_target_classes
keep_target_classes <- unlist(strsplit(args$keep_target_classes, ","))

# Parse boolean arguments
show_legend <- toupper(args$is_show_legend) == "TRUE"

# ======================== 辅助函数 ========================

# 生成错误页面 HTML
# 生成错误页面 HTML
generate_error_html <- function(output_html, error_message, title = "Analysis Failed") {
  html_content <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>%s</title>
  <style>
    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
      background-color: #f9fafb;
      color: #111827;
      display: flex;
      justify-content: center;
      align-items: center;
      min-height: 100vh;
      margin: 0;
    }
    .error-container {
      background-color: #ffffff;
      padding: 2rem;
      border: 1px solid #e5e7eb;
      border-radius: 0.5rem;
      max-width: 32rem;
      width: 100%%;
      box-shadow: 0 1px 3px 0 rgba(0, 0, 0, 0.1);
    }
    .error-title {
      font-size: 1.125rem;
      font-weight: 600;
      color: #ef4444;
      margin-bottom: 1rem;
      display: flex;
      align-items: center;
      gap: 0.5rem;
    }
    .error-message {
      background-color: #fef2f2;
      border: 1px solid #fee2e2;
      color: #991b1b;
      padding: 1rem;
      border-radius: 0.375rem;
      font-size: 0.875rem;
      line-height: 1.5;
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
      word-break: break-word;
      white-space: pre-wrap;
    }
    .suggestion {
      margin-top: 1.5rem;
      font-size: 0.875rem;
      color: #6b7280;
    }
  </style>
</head>
<body>
  <div class="error-container">
    <div class="error-title">
      <svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="10"></circle><line x1="12" y1="8" x2="12" y2="12"></line><line x1="12" y1="16" x2="12.01" y2="16"></line></svg>
      %s
    </div>
    <div class="error-message">%s</div>
    <div class="suggestion">Please adjust your parameters and try again.</div>
  </div>
</body>
</html>', title, title, error_message)
  
  writeLines(html_content, output_html)
  cat("  [Error HTML] Generated error page:", output_html, "\n")
}

# 带错误 HTML 输出的 stop 函数
stop_with_html <- function(error_message, output_html = NULL) {
  if (!is.null(output_html)) {
    generate_error_html(output_html, error_message)
  }
  stop(error_message)
}

lnc_empty <- function(code, message_text, counts = list()) {
  escaped <- gsub("&", "&amp;", as.character(message_text), fixed = TRUE)
  escaped <- gsub("<", "&lt;", escaped, fixed = TRUE)
  escaped <- gsub(">", "&gt;", escaped, fixed = TRUE)
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'><title>No lncRNA interactions</title>",
    "<style>body{font-family:Arial,sans-serif;background:#fff;color:#344054}",
    ".empty{max-width:760px;margin:80px auto;padding:32px;border:1px solid #f2c94c;",
    "border-radius:10px;background:#fffaf0}.empty h2{color:#8a5a00}</style></head>",
    "<body><div class='empty'><h2>No lncRNA target interactions</h2><p>",
    escaped, "</p></div></body></html>"
  )
  writeLines(html, args$output_html, useBytes = TRUE)
  empty_edges <- data.frame(
    from = character(), to = character(), title = character(),
    source_database = character(), evidence_type = character(), pmid = character(),
    stringsAsFactors = FALSE
  )
  write.table(empty_edges, file = args$output_csv, sep = "\t", quote = FALSE, row.names = FALSE)
  write_lnc_status(
    "empty", code, message_text, counts = counts,
    artifacts = list(html = args$output_html, csv = args$output_csv)
  )
  quit(save = "no", status = 0L)
}

get_significant <- function(df, pvalue_threshold = 0.05, is_use_padj = TRUE, context = "input") {
  requested <- if (isTRUE(is_use_padj)) "adj.P.Val" else "P.Value"
  if (!requested %in% colnames(df)) {
    stop_with_html(paste0("Requested ", requested, " is missing in ", context, "; no significance fallback is permitted."), args$output_html)
  }
  df %>% dplyr::filter(is.finite(.data[[requested]]) & .data[[requested]] <= pvalue_threshold)
}

# Save a portable htmlwidget without requiring a system Pandoc installation.
# htmlwidgets can write the dependency files itself; converting their local
# src/href values to data URIs keeps the returned HTML usable inside srcdoc and
# in downloaded standalone files.
save_portable_widget <- function(widget, output_html) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE) ||
      !requireNamespace("base64enc", quietly = TRUE)) {
    stop_with_html(
      "The server is missing the htmlwidgets/base64enc export dependency.",
      output_html
    )
  }

  temp_root <- tempfile("lncrna_widget_")
  dir.create(temp_root, recursive = TRUE)
  on.exit(unlink(temp_root, recursive = TRUE, force = TRUE), add = TRUE)
  temp_html <- file.path(temp_root, "widget.html")

  htmlwidgets::saveWidget(
    widget,
    file = temp_html,
    selfcontained = FALSE,
    libdir = "lib"
  )
  html <- readLines(temp_html, warn = FALSE, encoding = "UTF-8")

  inline_asset <- function(line, attribute, mime_type) {
    pattern <- paste0(attribute, "=\"([^\"]+)\"")
    match <- regexec(pattern, line, perl = TRUE)
    values <- regmatches(line, match)[[1]]
    if (length(values) < 2 || grepl("^(data:|https?:|//)", values[[2]])) {
      return(line)
    }
    asset_path <- file.path(temp_root, utils::URLdecode(values[[2]]))
    if (!file.exists(asset_path)) {
      return(line)
    }
    data_uri <- base64enc::dataURI(file = asset_path, mime = mime_type)
    sub(pattern, paste0(attribute, "=\"", data_uri, "\""), line, perl = TRUE)
  }

  for (index in seq_along(html)) {
    if (grepl("<script[^>]+src=\"", html[[index]], perl = TRUE)) {
      html[[index]] <- inline_asset(html[[index]], "src", "application/javascript")
    }
    if (grepl("<link[^>]+href=\"", html[[index]], perl = TRUE)) {
      html[[index]] <- inline_asset(html[[index]], "href", "text/css")
    }
  }
  writeLines(html, output_html, useBytes = TRUE)
}

# ======================== 主逻辑 ========================

cat("[lncRNA Network] lncRNA +", args$target_type, " regulatory network (showing intersection targets only)...\n")

# ======================== 读取RDA函数 ========================

read_deg_from_rda <- function(rda_file, output_html = NULL) {
  env <- new.env()
  load(rda_file, envir = env)
  if (!exists("deg", envir = env)) {
    stop_with_html(paste("deg object not found in", rda_file), output_html)
  }
  df <- env$deg
  df$ID <- rownames(df)
  return(df)
}

# 加载统一外部靶标数据库；旧版 master CSV 仍可通过同一参数传入
if (is.null(master_lnc_targets_path) || !nzchar(trimws(master_lnc_targets_path)) ||
    !file.exists(master_lnc_targets_path)) {
  stop_with_html(paste("lncRNA target database not found at:", master_lnc_targets_path), args$output_html)
}
source_databases <- NULL
if (!is.null(args$source_databases) && nzchar(trimws(args$source_databases)) &&
    tolower(trimws(args$source_databases)) != "all") {
  source_databases <- trimws(unlist(strsplit(args$source_databases, ",", fixed = TRUE)))
}
master_lnc_targets <- evomics_read_lnc_target_map(
  master_lnc_targets_path, target_type = args$target_type,
  evidence = args$evidence, species = args$species,
  source_databases = source_databases,
  include_derived_cerna = args$include_derived_cerna
)
cat("  Local target index loaded successfully, total", nrow(master_lnc_targets), "records\n")

# 读取差异表达数据 (从rda文件)
lnc_all <- read_deg_from_rda(args$lncrna_dataset, args$output_html)
other_all <- read_deg_from_rda(args$target_dataset, args$output_html)

# 按 P 值过滤
lnc_sig <- get_significant(lnc_all, args$lncrna_pvalue, args$lncrna_is_use_padj, "lncRNA differential input")
other_sig <- get_significant(other_all, args$target_pvalue, args$target_is_use_padj, "target differential input")
lnc_sig <- lnc_sig %>% dplyr::filter(abs(logFC) >= args$lncrna_log2fc)
other_sig <- other_sig %>% dplyr::filter(abs(logFC) >= args$target_log2fc)

## 1. 取差异 lnc 列名
if ("lncRNA" %in% colnames(lnc_sig)) {
  lnc_id_col <- "lncRNA"
} else if ("ID" %in% colnames(lnc_sig)) {
  lnc_id_col <- "ID"
    } else {
  stop_with_html("Cannot find lncRNA or ID column in lnc_sig, please check column names.", args$output_html)
}
lnc_list <- unique(lnc_sig[[lnc_id_col]])
lnc_list <- lnc_list[!is.na(lnc_list) & lnc_list != ""]

if (length(lnc_list) == 0) {
  metric <- if (args$lncrna_is_use_padj) "adjusted P" else "raw P"
  minimum_value <- if (args$lncrna_is_use_padj && "adj.P.Val" %in% names(lnc_all)) {
    suppressWarnings(min(lnc_all$adj.P.Val, na.rm = TRUE))
  } else {
    suppressWarnings(min(lnc_all$P.Value, na.rm = TRUE))
  }
  lnc_empty(
    "no_significant_features",
    paste0(
      "No lncRNA passed ", metric, " ≤ ", args$lncrna_pvalue,
      ". The minimum available ", metric, " is ", signif(minimum_value, 4), ". ",
      if (args$lncrna_is_use_padj) {
        "For an exploratory network, disable 'Use adjusted P for lncRNA'; keep it enabled for FDR-controlled inference."
      } else {
        "Consider a less stringent threshold only if it is scientifically justified."
      }
    ), counts = list(lncrna_features = 0)
  )
}
cat("  Number of differentially expressed lncRNA:", length(lnc_list), "\n")

## 2. 另一组学 ID
if ("ID" %in% colnames(other_sig)) {
  other_id_col <- "ID"
} else if (args$target_type == "miRNA" && "miRNA" %in% colnames(other_sig)) {
  other_id_col <- "miRNA"
} else {
  stop_with_html("Cannot find ID (or miRNA) column in other_sig, please check.", args$output_html)
}
other_ids <- unique(other_sig[[other_id_col]])
other_ids <- other_ids[!is.na(other_ids) & other_ids != ""]
cat("  Number of differentially expressed", args$target_type, ":", length(other_ids), "\n")

## 3. 从 master_lnc_targets 中筛选差异 lnc 的靶，并保留证据字段
required_master_cols <- c("lncRNA", "target", "target_type")
if (!all(required_master_cols %in% colnames(master_lnc_targets))) {
  stop_with_html(paste("Master target map is missing:",
                       paste(setdiff(required_master_cols, colnames(master_lnc_targets)), collapse = ", ")),
                 args$output_html)
}
master_lnc_targets <- master_lnc_targets %>%
  dplyr::mutate(
    lnc_key = toupper(sub("\\.[0-9]+$", "", trimws(lncRNA))),
    target_key = toupper(sub("\\.[0-9]+$", "", trimws(target)))
  )
if (args$target_type == "miRNA") {
  if (is.null(args$mirna_resolution) || !nzchar(args$mirna_resolution)) {
    stop_with_html("lncRNA + miRNA matching requires an explicit --mirna_resolution from dataset configuration", args$output_html)
  }
  master_lnc_targets$target_key <- resolve_mirna_feature_key(
    master_lnc_targets$target, args$mirna_resolution
  )
}
if ("evidence_score" %in% colnames(master_lnc_targets)) {
  master_lnc_targets <- master_lnc_targets %>%
    dplyr::filter(is.na(evidence_score) | evidence_score >= args$minimum_evidence_score)
}
lnc_keys <- toupper(sub("\\.[0-9]+$", "", lnc_list))
network_df0 <- master_lnc_targets %>%
  dplyr::filter(lnc_key %in% lnc_keys) %>%
  dplyr::select(dplyr::any_of(c("lncRNA", "target", "target_type", "source",
                                "source_database", "source_version", "species",
                                "evidence_type", "evidence_score", "pmid", "direction",
                                "relation_type", "metadata_json",
                                "lnc_key", "target_key"))) %>%
  distinct()

# Align lncRNA→miRNA annotations to the measured miRNA resolution.  For
# precursor datasets, 3p/5p records collapse to one measured hairpin key and
# never become separate expression nodes.
if (args$target_type == "miRNA" && nrow(network_df0) > 0) {
  network_df0 <- network_df0 %>%
    group_by(lnc_key, target_key) %>%
    summarise(
      across(everything(), ~ dplyr::first(.x)),
      .groups = "drop"
    )
}

if (nrow(network_df0) == 0) {
  lnc_empty("no_mapped_features", "No targets found for these lncRNAs in master_lnc_targets, skipping.",
            counts = list(lncrna_features = length(lnc_list), mapped_edges = 0))
}
cat("  Number of lnc-target relationships in master table:", nrow(network_df0), "\n")

## 4. 靶标类型分类
network_df <- network_df0 %>%
  mutate(
    class = dplyr::case_when(
      stringr::str_detect(target_type, stringr::regex("miRNA", ignore_case = TRUE))             ~ "miRNA",
      stringr::str_detect(target_type, stringr::regex("protein|pcg|gene|tf|mrna", ignore_case = TRUE)) ~ "gene_protein",
      stringr::str_detect(target_type, stringr::regex("RNA",   ignore_case = TRUE))             ~ "RNA",
      stringr::str_detect(target_type, stringr::regex("DNA",   ignore_case = TRUE))             ~ "DNA_TF",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::filter(class %in% keep_target_classes)

if (nrow(network_df) == 0) {
  lnc_empty("no_mapped_features", "No targets remain after filtering by keep_target_classes, skipping.",
            counts = list(lncrna_features = length(lnc_list), mapped_edges = 0))
}
cat("  Number of lnc-target relationships after type filtering:", nrow(network_df), "\n")

## 5. 标准化 ID 后与差异 omics 取交集
other_keys <- if (args$target_type == "miRNA") {
  resolve_mirna_feature_key(other_ids, args$mirna_resolution)
} else {
  toupper(sub("\\.[0-9]+$", "", trimws(other_ids)))
}
network_df_int <- network_df %>%
  dplyr::filter(target_key %in% other_keys)

if (args$direction_rule != "none" && nrow(network_df_int) > 0) {
  lfc <- lnc_sig %>% transmute(
    lnc_key = toupper(sub("\\.[0-9]+$", "", .data[[lnc_id_col]])),
    lnc_logFC = logFC
  )
  ofc <- other_sig %>% transmute(
    target_key = if (args$target_type == "miRNA")
      resolve_mirna_feature_key(.data[[other_id_col]], args$mirna_resolution)
    else toupper(sub("\\.[0-9]+$", "", .data[[other_id_col]])),
    target_logFC = logFC
  )
  network_df_int <- network_df_int %>%
    inner_join(lfc, by = "lnc_key") %>%
    inner_join(ofc, by = "target_key") %>%
    filter(if (args$direction_rule == "same")
      lnc_logFC * target_logFC > 0 else lnc_logFC * target_logFC < 0)
}

if (nrow(network_df_int) == 0) {
  lnc_empty(
    "no_interactions",
    paste0(
      "The selected lncRNAs have mapped targets, but none overlap the filtered ",
      args$target_type, " features. ",
      if (args$target_is_use_padj) {
        "For an exploratory network, disable 'Use adjusted P for Target' or review the threshold."
      } else {
        "Review the target P-value threshold and target-database coverage."
      }
    ), counts = list(lncrna_features = length(lnc_list), eligible_edges = 0)
  )
}
cat("  lnc targets ∩ differentially expressed", args$target_type, " size:", nrow(network_df_int), "\n")

## ========= 统一 logFC 范围用于颜色映射 =========
if (!"logFC" %in% colnames(lnc_sig)) {
  stop_with_html("lnc_sig is missing logFC column, cannot perform logFC mapping.", args$output_html)
}
if (!"logFC" %in% colnames(other_sig)) {
  stop_with_html("other_sig is missing logFC column, cannot perform logFC mapping.", args$output_html)
}

range_all <- range(c(lnc_sig$logFC, other_sig$logFC), na.rm = TRUE)

## 构建颜色渐变函数（使用用户指定的颜色）
color_palette <- c(args$low_color, args$mid_color, args$high_color)

## 为了后面 join，先整理出 id-logFC 表
lnc_logfc_tbl <- lnc_sig %>%
  dplyr::transmute(
    id = toupper(sub("\\.[0-9]+$", "", trimws(.data[[lnc_id_col]]))),
    logFC = logFC
  ) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(logFC = stats::median(logFC, na.rm = TRUE), .groups = "drop")

other_logfc_tbl <- other_sig %>%
  dplyr::transmute(
    id = if (args$target_type == "miRNA") {
      resolve_mirna_feature_key(.data[[other_id_col]], args$mirna_resolution)
    } else {
      toupper(sub("\\.[0-9]+$", "", trimws(.data[[other_id_col]])))
    },
    logFC = logFC
  ) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(logFC = stats::median(logFC, na.rm = TRUE), .groups = "drop")

## ========= 构建 nodes & edges =========

# lncRNA 节点
nodes_lnc <- data.frame(
  id    = lnc_keys,
  label = lnc_keys,
  group = "lncRNA",
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(lnc_logfc_tbl, by = "id") %>%
  dplyr::mutate(
    color = ifelse(
      is.na(logFC),
      "#bdbdbd",  ## Gray color for nodes without logFC
      scales::col_numeric(
        palette = color_palette,
        domain  = range_all
      )(logFC)
    ),
    value = ifelse(
      is.na(logFC),
      5,                  ## Medium size for nodes without logFC
      abs(logFC) * 5
    ),
    title = ifelse(
      is.na(logFC),
      paste0("lncRNA: ", id, "<br>logFC: NA"),
      paste0("lncRNA: ", id, "<br>logFC: ", round(logFC, 3))
    )
  )

# Target 节点 - 统一使用 target_type 作为 group（简化分组）
nodes_target <- network_df_int %>%
  dplyr::select(id = target_key) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    label = id,
    group = args$target_type  # 使用 target_type（mRNA/miRNA/Protein）作为组名
  ) %>%
  dplyr::left_join(other_logfc_tbl, by = "id") %>%
      dplyr::mutate(
    color = ifelse(
      is.na(logFC),
      "#bdbdbd",
      scales::col_numeric(
        palette = color_palette,
        domain  = range_all
      )(logFC)
    ),
    value = ifelse(
      is.na(logFC),
      5,
      abs(logFC) * 5
    ),
    title = ifelse(
      is.na(logFC),
      paste0(args$target_type, ": ", id, "<br>logFC: NA"),
      paste0(args$target_type, ": ", id, "<br>logFC: ", round(logFC, 3))
    )
  )

nodes <- dplyr::bind_rows(nodes_lnc, nodes_target) %>%
  dplyr::distinct(id, .keep_all = TRUE)

edges <- network_df_int %>%
  dplyr::transmute(
    from = lnc_key, to = target_key,
    title = paste0("Source: ", ifelse(is.na(source), "unknown", source),
                   "<br>Evidence: ", ifelse(is.na(evidence_type), "unknown", evidence_type),
                   ifelse(is.na(pmid) | !nzchar(pmid), "", paste0("<br>PMID: ", pmid))),
    source_database = source,
    evidence_type = evidence_type,
    pmid = pmid,
    mirna_resolution = if (args$target_type == "miRNA") args$mirna_resolution else NA_character_,
    resolution_mode = if (args$target_type == "miRNA" && args$mirna_resolution == "mirbase_precursor") "hairpin_collapsed" else if (args$target_type == "miRNA") "mature_exact" else NA_character_
  ) %>% dplyr::distinct()

## ========= 度过滤 =========
g <- igraph::graph_from_data_frame(d = edges,
                                   directed = TRUE,
                                   vertices = nodes)

deg <- igraph::degree(g, mode = "all")
keep_ids <- names(deg[deg >= args$min_degree])

if (length(keep_ids) == 0) {
  cat("  [Warning] All node degrees < min_degree, consider lowering min_degree. Not drawing network.\n")
  lnc_empty("no_interactions",
            paste0("All node degrees are less than min_degree (", args$min_degree,
                   "), cannot construct network. Please try lowering the min_degree parameter."),
            counts = list(eligible_edges = nrow(network_df_int), displayed_edges = 0))
}

nodes_filt <- nodes %>%
  dplyr::filter(id %in% keep_ids)

edges_filt <- edges %>%
  dplyr::filter(from %in% keep_ids, to %in% keep_ids) %>%
  dplyr::mutate(displayed = TRUE)

if (nrow(edges_filt) == 0) {
  lnc_empty(
    "no_interactions",
    paste0(
      "No complete lncRNA-target edge remains at minimum node degree ",
      args$min_degree, ". Lower Minimum node degree to 1."
    ), counts = list(eligible_edges = nrow(network_df_int), displayed_edges = 0)
  )
}

cat("  After degree filtering: nodes", nrow(nodes_filt), ", edges", nrow(edges_filt), "\n")

## ========= 构建网络图 =========
net <- visNetwork::visNetwork(nodes_filt, edges_filt,
                      main = list(text = args$main_title, 
                                  style = paste0("font-size:", args$main_title_size, "px;font-weight:bold;")),
                      width = paste0(args$width, "px"),
                      height = paste0(args$height, "px")) %>%
  visNetwork::visNodes(
    font    = list(face = "arial", size = args$axis_text_size),
    scaling = list(min = 5, max = 40)   ## Let value have more noticeable size differences
  ) %>%
  visNetwork::visEdges(smooth = FALSE) %>%
  ## Only define two groups: lncRNA and target_type (using user-specified shapes)
  visNetwork::visGroups(groupname = "lncRNA", shape = args$lncrna_shape) %>%
  visNetwork::visGroups(groupname = args$target_type, shape = args$target_shape) %>%
  visNetwork::visOptions(
    highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
    nodesIdSelection = TRUE,
    selectedBy       = "group"
  ) %>%
  visNetwork::visLayout(randomSeed = 123)

# 根据设置添加或隐藏图例
if (show_legend) {
  net <- net %>%
    visNetwork::visLegend(
      useGroups = TRUE,
      position = if (args$legend_position %in% c("left", "right")) args$legend_position else "right",
      main = list(text = "Legend", style = paste0("font-size:", args$legend_title_size, "px;")),
      ncol = 1,
      width = 0.15
    )
}

## ========= 保存输出 =========
save_portable_widget(net, args$output_html)
write.table(
  edges_filt,
  file      = args$output_csv,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("  Output file:", args$output_html, "\n")
cat("  Output file:", args$output_csv, "\n")

# 输出统计信息
result <- list(
  lncrna_count = sum(nodes_filt$group == "lncRNA"),
  target_count = sum(nodes_filt$group != "lncRNA"),
  edge_count   = nrow(edges_filt),
  full_eligible_edge_count = nrow(network_df_int),
  displayed_edge_count = nrow(edges_filt),
  min_degree = args$min_degree,
  mirna_resolution = if (args$target_type == "miRNA") args$mirna_resolution else NA_character_,
  resolution_mode = if (args$target_type == "miRNA" && args$mirna_resolution == "mirbase_precursor") "hairpin_collapsed" else if (args$target_type == "miRNA") "mature_exact" else NA_character_,
  target_type  = args$target_type,
  evidence = args$evidence,
  source_databases = if (is.null(source_databases)) "all" else paste(source_databases, collapse = ",")
)

cat(toJSON(result, auto_unbox = TRUE), "\n")
write_lnc_status(
  "success", "completed", "lncRNA target network generated successfully",
  counts = result,
  artifacts = list(html = args$output_html, csv = args$output_csv)
)

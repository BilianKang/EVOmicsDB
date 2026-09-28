#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(
  description = "Export Schema v2 model/group/statistic objects for the search index."
)
parser$add_argument("--input_rda", required = TRUE)
parser$add_argument("--output_model_tsv", required = TRUE)
parser$add_argument("--output_zscore_tsv", required = TRUE)
parser$add_argument("--output_stats_tsv", required = TRUE)
parser$add_argument("--output_groups_tsv", required = TRUE)
args <- parser$parse_args()

env <- new.env(parent = emptyenv())
load(args$input_rda, envir = env)
required <- c("exp_model", "group_list", "deg")
missing <- required[!vapply(required, exists, logical(1), envir = env, inherits = FALSE)]
if (length(missing)) stop("Schema v2 RDA is missing: ", paste(missing, collapse = ", "))

exp_model <- as.matrix(env$exp_model)
storage.mode(exp_model) <- "numeric"
if (is.null(rownames(exp_model)) || is.null(colnames(exp_model))) {
  stop("exp_model must have feature and sample names.")
}
if (anyDuplicated(rownames(exp_model)) || anyDuplicated(colnames(exp_model))) {
  stop("exp_model feature/sample names must be unique.")
}

groups <- as.character(env$group_list)
group_names <- names(env$group_list)
if (is.null(group_names) || !length(group_names)) {
  if (length(groups) != ncol(exp_model)) stop("Unnamed group_list does not align with exp_model.")
  group_names <- colnames(exp_model)
}
group_index <- match(colnames(exp_model), group_names)
if (anyNA(group_index)) stop("group_list does not cover every exp_model sample.")
groups <- groups[group_index]
if (!setequal(unique(groups), c("control", "case"))) {
  stop("group_list must contain canonical control/case labels.")
}

control_index <- groups == "control"
control_mean <- rowMeans(exp_model[, control_index, drop = FALSE], na.rm = TRUE)
control_sd <- apply(exp_model[, control_index, drop = FALSE], 1L, stats::sd, na.rm = TRUE)
control_mean[!is.finite(control_mean)] <- NA_real_
control_sd[!is.finite(control_sd) | control_sd <= sqrt(.Machine$double.eps)] <- NA_real_
zscore <- sweep(exp_model, 1L, control_mean, "-")
zscore <- sweep(zscore, 1L, control_sd, "/")

deg <- as.data.frame(env$deg, stringsAsFactors = FALSE)
if (!all(c("logFC", "P.Value", "adj.P.Val") %in% colnames(deg))) {
  stop("deg must contain logFC, P.Value and adj.P.Val.")
}
deg_index <- match(rownames(exp_model), rownames(deg))
if (anyNA(deg_index)) stop("deg does not cover every exp_model feature.")
stats <- data.frame(
  ID = rownames(exp_model),
  logFC = as.numeric(deg$logFC[deg_index]),
  P.Value = as.numeric(deg$P.Value[deg_index]),
  adj.P.Val = as.numeric(deg$adj.P.Val[deg_index]),
  statistic = if ("statistic" %in% colnames(deg)) as.numeric(deg$statistic[deg_index]) else NA_real_,
  statistic_type = if ("statistic_type" %in% colnames(deg)) as.character(deg$statistic_type[deg_index]) else NA_character_,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write.table(data.frame(ID = rownames(exp_model), exp_model, check.names = FALSE),
            args$output_model_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(data.frame(ID = rownames(zscore), zscore, check.names = FALSE),
            args$output_zscore_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(stats, args$output_stats_tsv, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
write.table(data.frame(sample = colnames(exp_model), group = groups),
            args$output_groups_tsv, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")

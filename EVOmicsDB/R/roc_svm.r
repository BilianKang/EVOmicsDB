#!/usr/bin/env Rscript
# Exploratory linear-SVM feature ranking.  [FULL REVISED SCRIPT]
# This script does NOT build or validate a diagnostic model.
#
# CHANGELOG vs original
#  - is_use_padj / legacy-prefix switch normalised via as_flag().
#  - allow_legacy_prefix_groups threaded into load_roc_input().
#  - weight_scale column added: weights are on the STANDARDIZED feature scale
#    (svm(scale = TRUE)), so they are not directly comparable to logFC.
#  - [STANDALONE-1] constant/non-finite candidate columns are rejected before fit.
#  - [STANDALONE-2] the fitted e1071 scale object and weight vector are verified
#    before the standardized-z-score audit label is emitted.
#  - cost fixed at 1 and recorded (kept out of the user surface on purpose to
#    avoid re-introducing tuning degrees of freedom into an exploratory ranking).
#  - audit columns added (threshold_comparator, cohort small-sample fields).

suppressPackageStartupMessages({ library(e1071); library(argparse) })

parser <- ArgumentParser(description = "Exploratory linear-SVM feature ranking; this script does not build or validate a diagnostic model.")
parser$add_argument("--input_rda", type = "character", required = TRUE)
parser$add_argument("--utils_r", type = "character", required = TRUE)
parser$add_argument("--output_csv", type = "character", required = TRUE)
parser$add_argument("--is_use_padj", type = "character", default = "TRUE")
parser$add_argument("--log2fc_threshold", type = "double", default = 1)
parser$add_argument("--pvalue_threshold", type = "double", default = 0.05)
parser$add_argument("--allow_legacy_prefix_groups", type = "character", default = "TRUE")
parser$add_argument("--group_mapping_case_label", type = "character", default = NULL)
parser$add_argument("--group_mapping_control_label", type = "character", default = NULL)
parser$add_argument("--group_mapping_source", type = "character", default = NULL)
parser$add_argument("--seed", type = "integer", default = 123)
args <- parser$parse_args()

source(args$utils_r)
args$is_use_padj <- as_flag(args$is_use_padj, "is_use_padj")
allow_legacy <- as_flag(args$allow_legacy_prefix_groups, "allow_legacy_prefix_groups")
group_mapping <- make_authoritative_group_mapping(
  args$group_mapping_case_label, args$group_mapping_control_label, args$group_mapping_source
)
set.seed(args$seed)
input <- load_roc_input(args$input_rda, allow_legacy_prefix = allow_legacy,
                         group_mapping = group_mapping,
                         group_mapping_source = args$group_mapping_source)
candidates <- candidate_features(input$deg, input$exp, args$is_use_padj,
                                 args$pvalue_threshold, args$log2fc_threshold,
                                 require_complete = TRUE, min_features = 2L)

x <- t(input$exp[candidates$genes, , drop = FALSE])
y <- input$group
if (ncol(x) < 2L) stop("At least two complete candidate features are required for SVM ranking.")

# Do not let e1071 silently disable scaling for a constant column.  This is a
# data-not-evaluable condition, distinct from a package/runtime failure.
if (any(!is.finite(x))) stop("SVM input is not evaluable: non-finite candidate expression values remain after filtering.")
feature_variance <- vapply(seq_len(ncol(x)), function(j) {
  value <- stats::var(x[, j])
  is.finite(value) && value > 0
}, logical(1))
if (any(!feature_variance)) {
  bad <- colnames(x)[!feature_variance]
  stop("SVM input is not evaluable: candidate feature(s) have zero or non-finite variance: ",
       paste(bad, collapse = ", "))
}

svm_cost <- 1  # fixed; recorded below for reproducibility
model <- tryCatch(
  e1071::svm(x = x, y = y, kernel = "linear", scale = TRUE, cost = svm_cost),
  error = function(e) stop("SVM technical failure during standardized fit: ", conditionMessage(e))
)
if (is.null(model$coefs) || is.null(model$SV)) stop("Unable to extract linear SVM feature weights.")
# With scale = TRUE, model$SV holds standardized support vectors.  Verify that
# the fitted object actually contains finite, non-zero scale parameters before
# claiming a standardized feature scale.
scale_center <- model$x.scale[["scaled:center"]]
scale_scale <- model$x.scale[["scaled:scale"]]
if (is.null(scale_center) || is.null(scale_scale) ||
    length(scale_center) != ncol(x) || length(scale_scale) != ncol(x) ||
    any(!is.finite(scale_center)) || any(!is.finite(scale_scale)) ||
    any(scale_scale <= 0)) {
  stop("SVM technical failure: fitted e1071 model did not retain valid standardization parameters.")
}
weights <- as.vector(t(model$coefs) %*% model$SV)
names(weights) <- colnames(x)
if (length(weights) != ncol(x) || any(!is.finite(weights))) {
  stop("SVM technical failure: fitted standardized feature weights are missing or non-finite.")
}

ranked <- data.frame(feature = names(weights), importance = abs(weights), stringsAsFactors = FALSE)
ranked <- ranked[order(-ranked$importance, ranked$feature), , drop = FALSE]
ranked$rank <- seq_len(nrow(ranked))
ranked$logFC <- input$deg[ranked$feature, "logFC"]
ranked$p_value <- input$deg[ranked$feature, if (args$is_use_padj) "adj.P.Val" else "P.Value"]
ranked$ranking_method <- "linear_svm_weight"
ranked$weight_scale <- "standardized_z_score"
ranked$analysis_level <- "exploratory_feature_ranking"
ranked$validation_status <- "not_a_predictive_model"
ranked$significance_metric <- candidates$metric
ranked$significance_threshold <- candidates$p_threshold
ranked$log2fc_threshold <- candidates$log2fc_threshold
ranked$threshold_comparator <- candidates$threshold_comparator
ranked$candidate_count_before_missing_filter <- candidates$candidate_before_missing_filter
ranked$excluded_missing_features <- candidates$excluded_missing
ranked$n_control <- input$n_control
ranked$n_tumor <- input$n_tumor
ranked$cohort_smallest_group_n <- input$smallest_group_n
ranked$small_sample_warning <- input$small_sample_warning
ranked$expression_matrix <- input$expression_source
ranked$group_source <- input$group_source
ranked$seed <- args$seed
ranked$svm_cost <- svm_cost
ranked$ranking_conditional_on <- "differential-expression pre-filter using the same cohort"
ranked$interpretation_note <- paste0(
  "Exploratory same-cohort |w| ranking on standardized features; no diagnostic-model ",
  "performance was estimated. Linear-SVM weights are unstable under strong feature ",
  "collinearity and p >> n.")

# Keep legacy names consumed by the current front end, while exposing explicit audit columns.
out <- ranked[, c("feature", "importance", "rank", setdiff(colnames(ranked), c("feature", "importance", "rank"))), drop = FALSE]
# Legacy `Gene` kept for the current front end; `Feature` is the correct name for a
# multi-omics resource (protein/miRNA/lncRNA/metabolite are not genes).
out$Feature <- out$feature
colnames(out)[colnames(out) == "feature"] <- "Gene"
colnames(out)[colnames(out) == "importance"] <- "Importance"
out <- out[, c("Gene", "Feature", setdiff(colnames(out), c("Gene", "Feature"))), drop = FALSE]
write.csv(out, args$output_csv, row.names = FALSE, quote = TRUE)

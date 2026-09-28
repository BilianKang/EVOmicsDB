#!/usr/bin/env Rscript
# Exploratory random-forest feature ranking.  [FULL REVISED SCRIPT]
# This script does NOT build or validate a diagnostic model.
#
# CHANGELOG vs original
#  - is_use_padj / legacy-prefix switch normalised via as_flag().
#  - allow_legacy_prefix_groups threaded into load_roc_input().
#  - audit columns added (threshold_comparator, cohort small-sample fields,
#    ranking-conditional-on note).
#  - ntree fixed-default 500 kept and recorded; class-balancing intentionally
#    NOT exposed (keeps the exploratory-ranking contract free of extra method
#    degrees of freedom).

suppressPackageStartupMessages({ library(randomForest); library(argparse) })

parser <- ArgumentParser(description = "Exploratory random-forest feature ranking; this script does not build or validate a diagnostic model.")
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
parser$add_argument("--ntree", type = "integer", default = 500,
                    help = "Compatibility argument; the freeze contract fixes ntree at 500.")
args <- parser$parse_args()

source(args$utils_r)
args$is_use_padj <- as_flag(args$is_use_padj, "is_use_padj")
allow_legacy <- as_flag(args$allow_legacy_prefix_groups, "allow_legacy_prefix_groups")
group_mapping <- make_authoritative_group_mapping(
  args$group_mapping_case_label, args$group_mapping_control_label, args$group_mapping_source
)
if (!is.finite(args$ntree) || args$ntree != 500L) {
  stop("The ROC freeze contract fixes random-forest ntree at 500; tuning ntree is not permitted.")
}
set.seed(args$seed)
input <- load_roc_input(args$input_rda, allow_legacy_prefix = allow_legacy,
                         group_mapping = group_mapping,
                         group_mapping_source = args$group_mapping_source)
candidates <- candidate_features(input$deg, input$exp, args$is_use_padj,
                                 args$pvalue_threshold, args$log2fc_threshold,
                                 require_complete = TRUE, min_features = 2L)

x <- t(input$exp[candidates$genes, , drop = FALSE])
y <- input$group
if (ncol(x) < 2L) stop("At least two complete candidate features are required for random-forest ranking.")

model <- randomForest::randomForest(x = x, y = y, importance = TRUE, ntree = args$ntree)
importance_matrix <- randomForest::importance(model, type = 1)
if (is.null(dim(importance_matrix))) importance_matrix <- matrix(importance_matrix, ncol = 1L, dimnames = list(names(importance_matrix), "MeanDecreaseAccuracy"))
# [FAIL-CLOSED] do not silently substitute a different importance column across
# package versions; a changed ranking metric must be an explicit error.
if (!"MeanDecreaseAccuracy" %in% colnames(importance_matrix)) {
  stop("MeanDecreaseAccuracy (permutation importance, type = 1) is unavailable; ",
       "ranking method was not substituted. Check the randomForest version/output.")
}
importance_column <- "MeanDecreaseAccuracy"
scores <- importance_matrix[, importance_column]

ranked <- data.frame(feature = names(scores), importance = as.numeric(scores), stringsAsFactors = FALSE)
ranked <- ranked[order(-ranked$importance, ranked$feature), , drop = FALSE]
ranked$rank <- seq_len(nrow(ranked))
ranked$logFC <- input$deg[ranked$feature, "logFC"]
ranked$p_value <- input$deg[ranked$feature, if (args$is_use_padj) "adj.P.Val" else "P.Value"]
ranked$ranking_method <- paste0("random_forest_", importance_column)
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
ranked$ntree <- args$ntree
ranked$ranking_conditional_on <- "differential-expression pre-filter using the same cohort"
ranked$interpretation_note <- "Exploratory same-cohort permutation-importance ranking; no diagnostic-model performance was estimated."

out <- ranked[, c("feature", "importance", "rank", setdiff(colnames(ranked), c("feature", "importance", "rank"))), drop = FALSE]
# Legacy `Gene` kept for the current front end; `Feature` is the correct name for a
# multi-omics resource (protein/miRNA/lncRNA/metabolite are not genes).
out$Feature <- out$feature
colnames(out)[colnames(out) == "feature"] <- "Gene"
colnames(out)[colnames(out) == "importance"] <- "Importance"
out <- out[, c("Gene", "Feature", setdiff(colnames(out), c("Gene", "Feature"))), drop = FALSE]
write.csv(out, args$output_csv, row.names = FALSE, quote = TRUE)

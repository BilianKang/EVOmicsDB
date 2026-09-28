#!/usr/bin/env Rscript
# Exploratory single-marker ROC.  [REVISED]
#
# CHANGELOG
#  [FIX-3]  --legend_position default is now valid; legacy "bottomright" accepted.
#  [FIX-6]  explicit --custom_gene now fails closed on any missing feature.
#  [NEW-1]  pROC direction is pre-specified from deg logFC (no max(AUC,1-AUC)).
#  [NEW-2]  Top-N pool can optionally be restricted to DE candidates.
#  [NEW-3]  selection_status / selection_bias_note audit columns.
#  [NEW-4]  small_sample_warning propagated to the summary table.
#  [NEW-5]  geom_line instead of geom_step -> curve matches the trapezoidal AUC.
#  [NEW-6]  ribbon drawn UNDER the curve; auto-disabled for multi-marker plots.
#  [NEW-7]  annotate() text size converted from pt to mm (was ~2.8x too large).
#  [NEW-8]  colour palette must cover all curves.
#  [NEW-9]  character(0) --custom_gene no longer bypasses top_n validation.
#  [NEW-10] each marker's ROC is computed once, not twice.
#  [NEW-11] cutoff labels no longer clipped at the right edge.
#  [STANDALONE-1] known data-not-evaluable checks are separated from pROC faults.
#  [STANDALONE-2] unknown pROC failures remain technical R errors.

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
  library(argparse)
})

parser <- ArgumentParser(description = "Generate exploratory single-marker ROC curves from a differential-analysis RDA.")
parser$add_argument("--input_rda", type = "character", required = TRUE)
parser$add_argument("--utils_r", type = "character", required = TRUE)
parser$add_argument("--output_png", type = "character", required = TRUE)
parser$add_argument("--output_csv", type = "character", required = TRUE)
parser$add_argument("--output_summary_csv", type = "character", required = TRUE)
parser$add_argument("--output_pdf", type = "character", default = NULL)
parser$add_argument("--top_n", type = "integer", default = 10)
parser$add_argument("--custom_gene", type = "character", nargs = "*", default = NULL)
parser$add_argument("--top_n_pool", type = "character", default = "all_features",
                    help = "all_features | deg_candidates")
parser$add_argument("--is_use_padj", type = "character", default = "TRUE")
parser$add_argument("--pvalue_threshold", type = "double", default = 0.05)
parser$add_argument("--log2fc_threshold", type = "double", default = 1)
parser$add_argument("--direction_mode", type = "character", default = "prespecified",
                    help = "prespecified (fixed from deg logFC) | auto (pROC data-driven)")
parser$add_argument("--allow_auto_direction_fallback", type = "character", default = "FALSE",
                    help = "In prespecified mode, whether a marker with no reliable logFC direction may fall back to pROC auto. Default FALSE (fail closed).")
parser$add_argument("--allow_legacy_prefix_groups", type = "character", default = "TRUE")
parser$add_argument("--group_mapping_case_label", type = "character", default = NULL,
                    help = "Explicit canonical case label; requires mapping source.")
parser$add_argument("--group_mapping_control_label", type = "character", default = NULL,
                    help = "Explicit canonical control label; requires mapping source.")
parser$add_argument("--group_mapping_source", type = "character", default = NULL,
                    help = "Authoritative provenance source for the group mapping.")
parser$add_argument("--line_color", type = "character", nargs = "+",
                    default = c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
                                "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"))
parser$add_argument("--line_type", type = "character", default = "solid")
parser$add_argument("--line_opacity", type = "double", default = 1.0)
parser$add_argument("--line_thickness", type = "double", default = 1.0)
parser$add_argument("--main_title", type = "character", default = "ROC Curve")
parser$add_argument("--x_axis_title", type = "character", default = "1 - Specificity (False Positive Rate)")
parser$add_argument("--y_axis_title", type = "character", default = "Sensitivity (True Positive Rate)")
parser$add_argument("--main_title_size", type = "double", default = 16)
parser$add_argument("--axis_title_size", type = "double", default = 14)
parser$add_argument("--axis_text_size", type = "double", default = 12)
parser$add_argument("--is_show_legend", type = "character", default = "TRUE")
parser$add_argument("--legend_text_size", type = "double", default = 10)
parser$add_argument("--legend_position", type = "character", default = "inside")
parser$add_argument("--is_show_diag_line", type = "character", default = "TRUE")
parser$add_argument("--diag_line_color", type = "character", default = "grey50")
parser$add_argument("--diag_line_type", type = "character", default = "dashed")
parser$add_argument("--diag_line_thickness", type = "double", default = 0.6)
parser$add_argument("--is_show_auc", type = "character", default = "TRUE")
parser$add_argument("--auc_text_opacity", type = "double", default = 0.1,
                    help = "Area fill opacity; this does not control AUC text opacity.")
parser$add_argument("--is_show_cutoff", type = "character", default = "TRUE")
parser$add_argument("--cutoff_text_size", type = "double", default = 10)
parser$add_argument("--cutoff_point_size", type = "double", default = 2.0)
parser$add_argument("--width", type = "double", default = 2400)
parser$add_argument("--height", type = "double", default = 1800)
parser$add_argument("--dpi", type = "double", default = 300)
args <- parser$parse_args()

source(args$utils_r)

# ---- normalise flags once (fixes NA-valued logicals) -------------------------
args$is_show_legend   <- as_flag(args$is_show_legend, "is_show_legend")
args$is_show_diag_line<- as_flag(args$is_show_diag_line, "is_show_diag_line")
args$is_show_auc      <- as_flag(args$is_show_auc, "is_show_auc")
args$is_show_cutoff   <- as_flag(args$is_show_cutoff, "is_show_cutoff")
args$is_use_padj      <- as_flag(args$is_use_padj, "is_use_padj")
allow_legacy          <- as_flag(args$allow_legacy_prefix_groups, "allow_legacy_prefix_groups")

supported_line_types <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
args$line_type <- tolower(trimws(args$line_type))
args$diag_line_type <- tolower(trimws(args$diag_line_type))
if (!args$line_type %in% supported_line_types) stop("Unsupported --line_type.")
if (!args$diag_line_type %in% supported_line_types) stop("Unsupported --diag_line_type.")
legend_cfg <- normalize_legend_position(args$legend_position)

args$top_n_pool <- tolower(trimws(args$top_n_pool))
if (!args$top_n_pool %in% c("all_features", "deg_candidates")) stop("Unsupported --top_n_pool.")
args$direction_mode <- tolower(trimws(args$direction_mode))
if (!args$direction_mode %in% c("prespecified", "auto")) stop("Unsupported --direction_mode.")
allow_auto_fallback <- as_flag(args$allow_auto_direction_fallback, "allow_auto_direction_fallback")
group_mapping <- make_authoritative_group_mapping(
  args$group_mapping_case_label, args$group_mapping_control_label, args$group_mapping_source
)
if (!is.finite(args$dpi) || args$dpi <= 0) stop("dpi must be positive.")

validate_style(args)

# [NEW-9] single source of truth for the custom-marker branch
has_custom <- !is.null(args$custom_gene) && length(args$custom_gene) > 0L
if (!has_custom && (!is.finite(args$top_n) || args$top_n < 1L)) {
  stop("top_n must be at least 1 when custom_gene is not provided.")
}

input <- load_roc_input(args$input_rda, allow_legacy_prefix = allow_legacy,
                         group_mapping = group_mapping,
                         group_mapping_source = args$group_mapping_source)
exp <- input$exp
group <- input$group
deg <- input$deg

# ---- ROC construction -------------------------------------------------------

# Returns either a result list, or a list(reason = <chr>) so the caller can decide
# whether an omission should be silent (Top-N) or fail-closed (user markers).
make_marker_roc <- function(gene) {
  values <- as.numeric(exp[gene, ])
  keep <- is.finite(values)
  group_keep <- droplevels(group[keep])
  if (nlevels(group_keep) < 2L) return(list(reason = "only one group after missing-value removal"))
  if (length(group_keep) < 4L || any(table(group_keep) < 2L)) {
    return(list(reason = "fewer than two usable samples in at least one group"))
  }
  if (length(unique(values[keep])) < 2L) {
    return(list(reason = "constant expression values after missing-value removal"))
  }

  dir_info <- if (args$direction_mode == "prespecified") {
    roc_direction_for(gene, deg, allow_auto_fallback = allow_auto_fallback)
  } else {
    list(direction = "auto", source = "data_driven_auto", ok = TRUE)
  }
  if (!isTRUE(dir_info$ok)) {
    return(list(reason = "no reliable logFC-based direction and auto fallback is disabled"))
  }

  # Test-only fault injection is explicit and never active in normal runs.  It
  # lets the standalone harness prove that unknown pROC failures are technical
  # errors rather than silently skipped Top-N markers.
  if (identical(Sys.getenv("ROC_STANDALONE_INJECT_PROC_ERROR", ""), "1")) {
    stop("pROC technical failure (injected standalone test fault)")
  }
  obj <- try(pROC::roc(group_keep, values[keep],
                       levels = c("Control", "Tumor"),
                       direction = dir_info$direction, quiet = TRUE), silent = TRUE)
  if (inherits(obj, "try-error")) {
    cond <- attr(obj, "condition")
    detail <- if (!is.null(cond)) conditionMessage(cond) else as.character(obj)
    stop("pROC technical failure for feature ", gene, ": ", detail)
  }
  ss <- small_sample_audit(group_keep)
  list(roc = obj, keep = keep, values = values[keep], group = group_keep,
       direction_source = dir_info$source,
       smallest_group_n = ss$smallest_group_n,
       small_sample_warning = ss$small_sample_warning,
       small_sample_note = ss$small_sample_note,
       auc = as.numeric(pROC::auc(obj)))
}

is_valid_marker <- function(z) is.null(z$reason)

# ---- marker selection -------------------------------------------------------

# DE-selection audit fields (populated only when Top-N uses deg_candidates)
sel_metric <- NA_character_; sel_p <- NA_real_; sel_lfc <- NA_real_; sel_cmp <- NA_character_

if (has_custom) {
  genes <- unique(args$custom_gene)
  # [FIX-6] fail closed on features absent from the matrix
  missing_genes <- setdiff(genes, rownames(exp))
  if (length(missing_genes) > 0L) {
    stop("Selected feature(s) are absent from the expression matrix: ",
         paste(missing_genes, collapse = ", "))
  }
  cached <- setNames(lapply(genes, make_marker_roc), genes)
  # [P2-2] user-specified markers must NOT be silently dropped when data are
  # insufficient; report each one with the reason and fail closed.
  invalid <- cached[!vapply(cached, is_valid_marker, logical(1))]
  if (length(invalid) > 0L) {
    detail <- paste(sprintf("%s (%s)", names(invalid),
                            vapply(invalid, function(z) z$reason, character(1))),
                    collapse = "; ")
    stop("Selected feature(s) cannot be evaluated: ", detail,
         ". Remove them or choose different markers.")
  }
  marker_results <- cached
  genes_to_plot <- names(marker_results)
  selection_method <- "User-selected markers"
  selection_status <- "user_specified"
  selection_bias_note <- "Markers were specified by the user; ROC was evaluated on the same cohort."
  pool_size <- NA_integer_

} else {
  pool <- rownames(exp)
  if (args$top_n_pool == "deg_candidates") {
    # [P2-3] single-marker contract: no complete-matrix requirement, >=1 candidate ok
    cand <- candidate_features(deg, exp, args$is_use_padj,
                               args$pvalue_threshold, args$log2fc_threshold,
                               require_complete = FALSE, min_features = 1L)
    pool <- cand$genes
    sel_metric <- cand$metric; sel_p <- cand$p_threshold
    sel_lfc <- cand$log2fc_threshold; sel_cmp <- cand$threshold_comparator
  }
  # [NEW-10] compute each ROC once; Top-N MAY silently skip unusable features
  cached <- setNames(lapply(pool, make_marker_roc), pool)
  ok <- vapply(cached, is_valid_marker, logical(1))
  if (!any(ok)) stop("No features have sufficient complete measurements for ROC analysis.")
  cached <- cached[ok]
  ranking <- data.frame(gene = names(cached),
                        apparent_auc = vapply(cached, function(z) z$auc, numeric(1)),
                        stringsAsFactors = FALSE)
  ranking <- ranking[order(-ranking$apparent_auc, ranking$gene), , drop = FALSE]
  genes_to_plot <- head(ranking$gene, min(args$top_n, nrow(ranking)))
  marker_results <- cached[genes_to_plot]
  selection_method <- "Top markers ranked by same-cohort apparent AUC"
  selection_status <- "same_cohort_auc_selected"
  selection_bias_note <- paste0(
    "Marker selection and ROC evaluation used the same cohort. The reported ",
    "confidence interval, cutoff, sensitivity and specificity are apparent, ",
    "same-cohort estimates and are NOT adjusted for feature selection.")
  pool_size <- nrow(ranking)
}

genes_to_plot <- names(marker_results)
if (length(genes_to_plot) == 0L) {
  stop("Selected features do not have sufficient complete measurements for ROC analysis.")
}

# [NEW-8] a colour must exist for every curve
if (length(args$line_color) < length(genes_to_plot)) {
  stop("Only ", length(args$line_color), " colour(s) supplied for ",
       length(genes_to_plot), " curves; supply at least one colour per marker.")
}
colors <- args$line_color[seq_along(genes_to_plot)]

# ---- assemble tables --------------------------------------------------------

curve_rows <- list(); summary_rows <- list(); labels <- character(length(genes_to_plot))
for (i in seq_along(genes_to_plot)) {
  gene <- genes_to_plot[[i]]; result <- marker_results[[i]]; roc_obj <- result$roc
  auc_value <- result$auc
  ci <- safe_ci_auc(roc_obj)
  cut_info <- first_youden(roc_obj); cut <- cut_info$row
  tumor_median <- stats::median(result$values[result$group == "Tumor"], na.rm = TRUE)
  control_median <- stats::median(result$values[result$group == "Control"], na.rm = TRUE)
  tumor_direction <- ifelse(tumor_median >= control_median, "Higher in Tumor", "Lower in Tumor")
  labels[[i]] <- if (args$is_show_auc) sprintf("%s (AUC=%.3f)", gene, auc_value) else gene

  curve_rows[[i]] <- data.frame(
    FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities,
    marker = labels[[i]], original_marker = gene, apparent_auc = auc_value,
    stringsAsFactors = FALSE)

  summary_rows[[i]] <- data.frame(
    marker = gene,
    selection_method = selection_method,
    selection_status = selection_status,
    selection_bias_note = selection_bias_note,
    selection_pool_size = pool_size,
    top_n_pool = if (has_custom) NA_character_ else args$top_n_pool,
    selection_significance_metric = sel_metric,
    selection_pvalue_threshold = sel_p,
    selection_log2fc_threshold = sel_lfc,
    selection_threshold_comparator = sel_cmp,
    validation_status = "exploratory_same_cohort",
    apparent_auc = auc_value,
    auc_ci_low = ci[["low"]], auc_ci_high = ci[["high"]],
    auc_ci_type = "naive_DeLong_same_cohort_not_selection_adjusted",
    auc_ci_method = ci[["method"]],
    auc_ci_status = ci[["status"]],
    auc_ci_note = ci[["note"]],
    generalization_performance_estimated = FALSE,
    roc_direction = roc_obj$direction,
    roc_direction_source = result$direction_source,
    direction_in_tumor = tumor_direction,
    cutoff = as.numeric(cut$threshold),
    cutoff_type = "same_cohort_youden",
    cutoff_tied_solutions = cut_info$n_ties,
    sensitivity = as.numeric(cut$sensitivity),
    specificity = as.numeric(cut$specificity),
    n_used = sum(result$keep),
    n_control = sum(result$group == "Control"),
    n_tumor = sum(result$group == "Tumor"),
    n_excluded_missing = sum(!result$keep),
    smallest_group_n = result$smallest_group_n,
    small_sample_warning = result$small_sample_warning,
    small_sample_note = result$small_sample_note,
    cohort_smallest_group_n = input$smallest_group_n,
    expression_matrix = input$expression_source,
    group_source = input$group_source,
    interpretation_note = "Exploratory same-cohort ROC; not independently validated.",
    stringsAsFactors = FALSE)
}

roc_df <- do.call(rbind, curve_rows)
roc_df <- roc_df[order(roc_df$marker, roc_df$FPR, roc_df$TPR), , drop = FALSE]
roc_df$marker <- factor(roc_df$marker, levels = labels)
summary_df <- do.call(rbind, summary_rows)

# [P1-2] "common sample set" must compare the actual samples used, not just their
# counts: two markers can both use n=80 while using different 80 samples.
reference_keep <- marker_results[[1]]$keep
common_sample_set <- all(vapply(marker_results,
                                function(z) identical(z$keep, reference_keep),
                                logical(1)))
summary_df$common_sample_set <- common_sample_set

# ---- plot -------------------------------------------------------------------

# [NEW-6] ribbon fill is only meaningful for a single curve
show_fill <- args$is_show_auc && length(genes_to_plot) == 1L

layers <- list()
if (args$is_show_diag_line) {
  layers <- c(layers, list(geom_abline(slope = 1, intercept = 0,
                                       linetype = args$diag_line_type,
                                       color = args$diag_line_color,
                                       linewidth = args$diag_line_thickness)))
}
if (show_fill) {
  layers <- c(layers, list(
    geom_ribbon(aes(ymin = 0, ymax = TPR, fill = marker),
                alpha = args$auc_text_opacity, color = NA),
    scale_fill_manual(values = colors, guide = "none")))
}
# [NEW-5] geom_line matches pROC's trapezoidal AUC; geom_step did not.
layers <- c(layers, list(
  geom_line(linewidth = args$line_thickness, linetype = args$line_type,
            alpha = args$line_opacity)))

gg <- ggplot(roc_df, aes(x = FPR, y = TPR, color = marker, group = marker)) +
  layers +
  scale_color_manual(values = colors, name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0), limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(expand = expansion(mult = 0), limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # [PLOT-FIX] small symmetric padding so a perfect ROC (edges at x=0 / y=1)
  # sits INSIDE the panel border instead of being painted over by it (Fig 5G).
  coord_equal(xlim = c(-0.015, 1.015), ylim = c(-0.015, 1.015), clip = "off") +
  labs(x = args$x_axis_title, y = args$y_axis_title, title = args$main_title,
       subtitle = "Exploratory same-cohort ROC; performance is not independently validated.") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(size = args$main_title_size, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = max(8, args$axis_text_size - 1), hjust = 0.5),
    axis.title = element_text(size = args$axis_title_size, face = "bold"),
    axis.text = element_text(size = args$axis_text_size),
    legend.text = element_text(size = args$legend_text_size),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(8, 18, 8, 8),
    legend.position = if (args$is_show_legend) legend_cfg$position else "none")

if (args$is_show_legend && legend_cfg$position == "inside") {
  gg <- gg + theme(
    legend.position.inside = legend_cfg$coord,
    legend.justification = legend_cfg$justification,
    legend.background = element_rect(fill = scales::alpha("white", 0.82), colour = "grey75"))
}

if (args$is_show_cutoff) {
  for (i in seq_len(nrow(summary_df))) {
    px <- 1 - summary_df$specificity[[i]]
    py <- summary_df$sensitivity[[i]]
    # [NEW-11] flip the label inward near the right edge instead of clipping it
    right_side <- px > 0.75
    gg <- gg +
      annotate("point", x = px, y = py, color = colors[[i]], size = args$cutoff_point_size) +
      annotate("text",
               x = if (right_side) px - 0.02 else px + 0.02,
               y = max(0, py - 0.03),
               label = format(round(summary_df$cutoff[[i]], 3), trim = TRUE),
               color = colors[[i]],
               # [NEW-7] annotate() size is mm, element_text() size is pt
               size = args$cutoff_text_size / ggplot2::.pt,
               hjust = if (right_side) 1 else 0)
  }
}

ggsave(args$output_png, gg, width = args$width, height = args$height,
       units = "px", device = "png", dpi = args$dpi)
if (!is.null(args$output_pdf) && nzchar(trimws(args$output_pdf))) {
  ggsave(args$output_pdf, gg, width = args$width / args$dpi,
         height = args$height / args$dpi, units = "in", device = cairo_pdf)
}
write.csv(roc_df, args$output_csv, row.names = FALSE)
write.csv(summary_df, args$output_summary_csv, row.names = FALSE)

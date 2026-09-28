#!/usr/bin/env Rscript
# Exploratory user-defined logistic-panel ROC.  [FULL REVISED SCRIPT]
#
# CHANGELOG vs original
#  [P1-1] first_youden() now returns list(row, n_ties); adapted here.
#  [P1b]  logical CLI args normalised with as_flag() (no NA-valued flags).
#  [P2-1] small-sample audit computed on the complete-case panel (y), not the
#         whole-cohort counts.
#  [P2-a] geom_line instead of geom_step (matches trapezoidal AUC).
#  [P2-b] ribbon drawn UNDER the curve.
#  [P2-c] annotate() text size converted pt -> mm.
#  [P2-d] DeLong CI semantics made explicit (naive, in-sample fitted scores).
#  [P3]   `warnings` local renamed to `model_warnings` (was shadowing base::warnings).
#  [order] source(utils_r) moved BEFORE any validation that needs its helpers.

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
  library(argparse)
})

parser <- ArgumentParser(description = "Generate an exploratory, user-defined logistic-panel ROC curve.")
parser$add_argument("--input_rda", type = "character", required = TRUE)
parser$add_argument("--utils_r", type = "character", required = TRUE)
parser$add_argument("--output_png", type = "character", required = TRUE)
parser$add_argument("--output_csv", type = "character", required = TRUE)
parser$add_argument("--output_summary_csv", type = "character", required = TRUE)
parser$add_argument("--output_coefficients_csv", type = "character", required = TRUE)
parser$add_argument("--output_pdf", type = "character", default = NULL)
parser$add_argument("--custom_gene", type = "character", nargs = "+", required = TRUE)
parser$add_argument("--allow_legacy_prefix_groups", type = "character", default = "TRUE")
parser$add_argument("--group_mapping_case_label", type = "character", default = NULL)
parser$add_argument("--group_mapping_control_label", type = "character", default = NULL)
parser$add_argument("--group_mapping_source", type = "character", default = NULL)
parser$add_argument("--line_color", type = "character", default = "#1f77b4")
parser$add_argument("--line_type", type = "character", default = "solid")
parser$add_argument("--line_opacity", type = "double", default = 1.0)
parser$add_argument("--line_thickness", type = "double", default = 1.0)
parser$add_argument("--main_title", type = "character", default = "Logistic ROC Plot")
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
parser$add_argument("--auc_text_opacity", type = "double", default = 0.1)
parser$add_argument("--is_show_cutoff", type = "character", default = "TRUE")
parser$add_argument("--cutoff_text_size", type = "double", default = 10)
parser$add_argument("--cutoff_point_size", type = "double", default = 2.0)
parser$add_argument("--width", type = "double", default = 2400)
parser$add_argument("--height", type = "double", default = 1800)
parser$add_argument("--dpi", type = "double", default = 300)
args <- parser$parse_args()

# [order] helpers must exist before we validate anything that uses them.
source(args$utils_r)

args$is_show_legend    <- as_flag(args$is_show_legend,    "is_show_legend")
args$is_show_diag_line <- as_flag(args$is_show_diag_line, "is_show_diag_line")
args$is_show_auc       <- as_flag(args$is_show_auc,       "is_show_auc")
args$is_show_cutoff    <- as_flag(args$is_show_cutoff,    "is_show_cutoff")
allow_legacy           <- as_flag(args$allow_legacy_prefix_groups, "allow_legacy_prefix_groups")
group_mapping <- make_authoritative_group_mapping(
  args$group_mapping_case_label, args$group_mapping_control_label, args$group_mapping_source
)

supported_line_types <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
args$line_type <- tolower(trimws(args$line_type))
args$diag_line_type <- tolower(trimws(args$diag_line_type))
if (!args$line_type %in% supported_line_types) stop("Unsupported --line_type.")
if (!args$diag_line_type %in% supported_line_types) stop("Unsupported --diag_line_type.")
legend_cfg <- normalize_legend_position(args$legend_position)
if (!is.finite(args$dpi) || args$dpi <= 0) stop("dpi must be positive.")

validate_style(args)
input <- load_roc_input(args$input_rda, allow_legacy_prefix = allow_legacy,
                         group_mapping = group_mapping,
                         group_mapping_source = args$group_mapping_source)
exp <- input$exp
group <- input$group
genes <- unique(args$custom_gene)
missing_genes <- setdiff(genes, rownames(exp))
if (length(missing_genes) > 0L) stop("Selected feature(s) are absent from the expression matrix: ", paste(missing_genes, collapse = ", "))
if (length(genes) < 1L) stop("At least one feature must be selected.")

x <- t(exp[genes, , drop = FALSE])
keep <- apply(x, 1L, function(z) all(is.finite(z))) & !is.na(group)
x <- x[keep, , drop = FALSE]
y <- factor(group[keep], levels = c("Control", "Tumor"))
if (nrow(x) < 4L || any(table(y) < 2L)) stop("The selected panel has insufficient complete samples in one or both groups.")
if (ncol(x) >= nrow(x)) stop("Number of selected features must be smaller than the number of complete samples.")
if (any(apply(x, 2L, stats::sd) == 0)) {
  flat <- genes[apply(x, 2L, stats::sd) == 0]
  stop("Selected feature(s) have zero variance after complete-case filtering: ", paste(flat, collapse = ", "))
}

# [P2-1] small-sample audit on the complete-case panel, not the whole cohort.
panel_ss <- small_sample_audit(y)

# Make syntactically safe model columns, while preserving original IDs in all outputs.
safe_names <- make.unique(make.names(colnames(x)))
colnames(x) <- safe_names
model_data <- data.frame(Tumor = as.integer(y == "Tumor"), x, check.names = FALSE)
model_warnings <- character()
model <- withCallingHandlers(
  stats::glm(Tumor ~ ., data = model_data, family = stats::binomial()),
  warning = function(w) {
    model_warnings <<- c(model_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
if (!isTRUE(model$converged)) model_warnings <- c(model_warnings, "Logistic regression did not converge.")
if (any(!is.finite(stats::coef(model)))) model_warnings <- c(model_warnings, "Non-finite coefficient detected; complete or quasi-complete separation is possible.")

prob <- as.numeric(stats::predict(model, type = "response"))
roc_obj <- pROC::roc(response = model_data$Tumor, predictor = prob, levels = c(0, 1), direction = "<", quiet = TRUE)
auc_value <- as.numeric(pROC::auc(roc_obj)); ci <- safe_ci_auc(roc_obj)
# [P1-1] adapt to first_youden() -> list(row, n_ties)
cut_info <- first_youden(roc_obj); cut <- cut_info$row
# [EPV-FIX] strict EPV uses the event (Tumor) count; the conservative
# minority-class ratio drives the stability warning.
n_tumor_evt <- sum(model_data$Tumor == 1)
n_control_evt <- sum(model_data$Tumor == 0)
events_per_variable <- n_tumor_evt / length(genes)
minority_class_per_variable <- min(n_tumor_evt, n_control_evt) / length(genes)
if (minority_class_per_variable < 10) model_warnings <- c(model_warnings, sprintf("Minority-class observations per variable = %.2f (<10); estimates may be unstable.", minority_class_per_variable))
if (min(n_tumor_evt, n_control_evt) < 10L) model_warnings <- c(model_warnings, "Smallest group has fewer than 10 samples; this exploratory model should not be interpreted as validated.")
warning_text <- if (length(model_warnings)) paste(unique(model_warnings), collapse = " | ") else ""

curve_df <- data.frame(
  FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities,
  model = sprintf("Logistic Regression (Apparent AUC=%.3f)", auc_value),
  apparent_auc = auc_value, validation_status = "exploratory_same_cohort", stringsAsFactors = FALSE
)
curve_df <- curve_df[order(curve_df$FPR, curve_df$TPR), , drop = FALSE]
summary_df <- data.frame(
  analysis_type = "user_defined_logistic_panel",
  validation_status = "exploratory_same_cohort",
  selection_status = "user_specified_panel",
  selection_bias_note = paste0(
    "Apparent AUC from in-sample fitted probabilities. The DeLong interval treats ",
    "the fitted scores as a fixed predictor and does NOT account for model-fitting ",
    "uncertainty, optimism, or variable selection."),
  auc_label = "Apparent AUC",
  apparent_auc = auc_value, auc_ci_low = ci[["low"]], auc_ci_high = ci[["high"]],
  auc_ci_type = "naive_DeLong_on_same_cohort_fitted_scores",
  auc_ci_method = ci[["method"]],
  auc_ci_status = ci[["status"]],
  auc_ci_note = ci[["note"]],
  generalization_performance_estimated = FALSE,
  cutoff = as.numeric(cut$threshold),
  cutoff_type = "same_cohort_youden",
  cutoff_tied_solutions = cut_info$n_ties,
  sensitivity = as.numeric(cut$sensitivity), specificity = as.numeric(cut$specificity),
  panel_size = length(genes), panel_features = paste(genes, collapse = ";"),
  n_total_in_rda = ncol(exp), n_used_complete_case = nrow(x), n_excluded_missing = sum(!keep),
  n_control = sum(y == "Control"), n_tumor = sum(y == "Tumor"),
  events_per_variable = events_per_variable,
  minority_class_per_variable = minority_class_per_variable,
  smallest_group_n = panel_ss$smallest_group_n,
  small_sample_warning = panel_ss$small_sample_warning,
  small_sample_note = panel_ss$small_sample_note,
  cohort_smallest_group_n = input$smallest_group_n,
  expression_matrix = input$expression_source, group_source = input$group_source,
  interpretation_note = "User-defined exploratory panel. Candidate/panel selection and evaluation may use the same cohort; apparent AUC is not independently validated.",
  model_warning = warning_text, stringsAsFactors = FALSE
)

sm <- summary(model)$coefficients
coef_df <- data.frame(
  term = rownames(sm), estimate = sm[, "Estimate"], odds_ratio = exp(sm[, "Estimate"]),
  std_error = sm[, "Std. Error"], z_value = sm[, "z value"], p_value = sm[, "Pr(>|z|)"], stringsAsFactors = FALSE
)
coef_df$feature <- ifelse(coef_df$term == "(Intercept)", "(Intercept)", genes[match(coef_df$term, safe_names)])
coef_df <- coef_df[, c("feature", "term", "estimate", "odds_ratio", "std_error", "z_value", "p_value")]

# ---- plot (ribbon under curve; geom_line; mm text size) ---------------------
layers <- list()
if (args$is_show_diag_line) {
  layers <- c(layers, list(geom_abline(slope = 1, intercept = 0,
                                       linetype = args$diag_line_type,
                                       color = args$diag_line_color,
                                       linewidth = args$diag_line_thickness)))
}
if (args$is_show_auc) {
  layers <- c(layers, list(
    geom_ribbon(aes(ymin = 0, ymax = TPR, fill = model),
                alpha = args$auc_text_opacity, color = NA),
    scale_fill_manual(values = args$line_color, guide = "none")))
}
layers <- c(layers, list(
  geom_line(linewidth = args$line_thickness, linetype = args$line_type,
            alpha = args$line_opacity)))

gg <- ggplot(curve_df, aes(x = FPR, y = TPR, color = model, group = model)) +
  layers +
  scale_color_manual(values = args$line_color, name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0), limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(expand = expansion(mult = 0), limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  # [PLOT-FIX] padding so a perfect ROC is not hidden under the panel border.
  coord_equal(xlim = c(-0.015, 1.015), ylim = c(-0.015, 1.015), clip = "off") +
  labs(x = args$x_axis_title, y = args$y_axis_title, title = args$main_title,
       subtitle = "Exploratory user-defined panel; apparent same-cohort AUC, not independently validated.") +
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
  px <- 1 - summary_df$specificity[[1]]
  py <- summary_df$sensitivity[[1]]
  right_side <- px > 0.75
  gg <- gg +
    annotate("point", x = px, y = py, color = args$line_color, size = args$cutoff_point_size) +
    annotate("text",
             x = if (right_side) px - 0.02 else px + 0.02,
             y = max(0, py - 0.03),
             label = format(round(summary_df$cutoff[[1]], 3), trim = TRUE),
             color = args$line_color,
             size = args$cutoff_text_size / ggplot2::.pt,
             hjust = if (right_side) 1 else 0)
}

ggsave(args$output_png, gg, width = args$width, height = args$height,
       units = "px", device = "png", dpi = args$dpi)
if (!is.null(args$output_pdf) && nzchar(trimws(args$output_pdf))) {
  ggsave(args$output_pdf, gg, width = args$width / args$dpi,
         height = args$height / args$dpi, units = "in", device = cairo_pdf)
}
write.csv(curve_df, args$output_csv, row.names = FALSE)
write.csv(summary_df, args$output_summary_csv, row.names = FALSE)
write.csv(coef_df, args$output_coefficients_csv, row.names = FALSE)

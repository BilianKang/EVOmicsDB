# Shared utilities for the exploratory ROC module.  [REVISED]
#
# Contract:
# - `exp_model` is the measured/model-scale matrix saved by differential analysis.
# - `group_list` is the authoritative group vector when present.
#   * If group_list is NAMED, names MUST cover all expression sample IDs.
#     Positional fallback is NEVER used for a named group_list (fail closed).
# - Legacy RDA files may only contain `exp`; in that case a STRICT C/N sample-ID
#   pattern is required (e.g. C1, C_001, N-12). Loose "^C"/"^N" matching is gone.
#
# CHANGELOG vs previous version
#  [FIX-1] resolve_roc_group(): named group_list can no longer fall back to
#          positional matching; duplicate names rejected.
#  [FIX-2] legacy prefix inference tightened to LEGACY_GROUP_PATTERN.
#  [FIX-3] canonical_group(): ambiguous labels (patient/case/...) now stop with an
#          explicit message instead of being silently mapped to Tumor.
#          >>> ACTION REQUIRED: audit unique(group_list) across all production RDAs
#          >>> before finalising CONTROL_LABELS / TUMOR_LABELS below.
#  [FIX-4] normalize_legend_position(): accepts legacy pROC-style values.
#  [FIX-5] as_flag(): logical CLI args can no longer become NA.
#  [FIX-6] small-sample audit fields returned instead of silent pass.
#  [FIX-7] candidate_features() reports its comparator explicitly.
#  [FIX-8] dead/broken csv_value() removed.

# ---------------------------------------------------------------- label tables

CONTROL_LABELS <- c("n", "normal", "control", "controls", "healthy",
                    "noncancer", "non-cancer")
TUMOR_LABELS   <- c("c", "cancer", "tumor", "tumour")
# Semantically ambiguous in a cancer-vs-healthy resource: a benign patient is
# still a "patient", and a "case" is study-design dependent. Fail closed.
AMBIGUOUS_LABELS <- c("patient", "patients", "case", "cases", "benign",
                      "adjacent", "para-cancerous", "paracancerous")

# Strict legacy sample-ID grammar. C1 / C_001 / C-12 / N.3 ...
# NOTE: confirm against real legacy sample IDs before deploying.
LEGACY_GROUP_PATTERN <- "^([CN])[._-]?[0-9]+$"

SMALL_GROUP_WARN_N <- 10L   # below this -> warning flag, analysis still allowed
MIN_GROUP_N        <- 2L    # below this -> hard stop

# ---------------------------------------------------------------- CLI helpers

as_flag <- function(x, label) {
  if (length(x) != 1L) stop(label, " must be a single TRUE/FALSE value.")
  if (is.logical(x) && !is.na(x)) return(x)
  key <- tolower(trimws(as.character(x)))
  if (key %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (key %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(label, " must be TRUE or FALSE; received: ", as.character(x))
}

LEGEND_POSITIONS <- c("top", "right", "bottom", "left", "inside", "none")

# Legacy pROC-style positions are mapped to inside + justification so that
# existing frontend/API payloads keep working.
LEGEND_LEGACY_MAP <- list(
  bottomright = c("right", "bottom"),
  bottomleft  = c("left",  "bottom"),
  topright    = c("right", "top"),
  topleft     = c("left",  "top"),
  center      = c("center", "center")
)

normalize_legend_position <- function(value) {
  key <- tolower(trimws(value))
  if (key %in% names(LEGEND_LEGACY_MAP)) {
    just <- LEGEND_LEGACY_MAP[[key]]
    coord <- c(
      switch(just[1], left = 0.02, center = 0.5, right = 0.98),
      switch(just[2], bottom = 0.02, center = 0.5, top = 0.98)
    )
    return(list(position = "inside", justification = just,
                coord = coord, requested = key))
  }
  if (!key %in% LEGEND_POSITIONS) {
    stop("Unsupported --legend_position: ", value, ". Allowed: ",
         paste(c(LEGEND_POSITIONS, names(LEGEND_LEGACY_MAP)), collapse = ", "))
  }
  list(position = key,
       justification = c("right", "bottom"),
       coord = c(0.98, 0.02),
       requested = key)
}

# ---------------------------------------------------------------- matrix / group

as_numeric_matrix <- function(x, object_name) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(object_name, " must have feature IDs as row names and sample IDs as column names.")
  }
  if (anyDuplicated(rownames(x))) stop(object_name, " contains duplicate feature IDs.")
  if (anyDuplicated(colnames(x))) stop(object_name, " contains duplicate sample IDs.")
  x
}

validate_group_mapping <- function(group_mapping) {
  if (is.null(group_mapping)) return(NULL)
  if (is.data.frame(group_mapping) || is.matrix(group_mapping) || is.list(group_mapping)) {
    stop("authoritative group_mapping must be a named atomic vector.")
  }
  values <- trimws(as.character(group_mapping))
  labels <- names(group_mapping)
  if (is.null(labels) || any(!nzchar(trimws(labels)))) {
    stop("authoritative group_mapping must have non-empty names.")
  }
  keys <- tolower(trimws(labels))
  if (anyDuplicated(keys)) stop("authoritative group_mapping contains duplicate labels.")
  if (any(!values %in% c("Control", "Tumor"))) {
    stop("authoritative group_mapping values must be exactly Control or Tumor.")
  }
  names(values) <- keys
  values
}

make_authoritative_group_mapping <- function(case_label = NULL,
                                             control_label = NULL,
                                             mapping_source = NULL) {
  case_label <- if (is.null(case_label)) "" else trimws(as.character(case_label))
  control_label <- if (is.null(control_label)) "" else trimws(as.character(control_label))
  if (!nzchar(case_label) && !nzchar(control_label)) return(NULL)
  if (!nzchar(case_label) || !nzchar(control_label)) {
    stop("Both case_label and control_label are required for authoritative mapping.")
  }
  if (is.null(mapping_source) || !nzchar(trimws(as.character(mapping_source)))) {
    stop("mapping_source is required for authoritative group mapping.")
  }
  mapping <- c("Tumor", "Control")
  names(mapping) <- c(case_label, control_label)
  validate_group_mapping(mapping)
}

canonical_group <- function(x, group_mapping = NULL) {
  raw <- trimws(as.character(x))
  key <- tolower(raw)

  mapping <- validate_group_mapping(group_mapping)
  if (!is.null(mapping)) {
    out <- unname(mapping[key])
    if (anyNA(out)) {
      bad <- unique(raw[is.na(out)])
      stop("group labels are not covered by the authoritative group_mapping: ",
           paste(bad, collapse = ", "))
    }
    return(out)
  }

  hit_ambiguous <- key %in% AMBIGUOUS_LABELS
  if (any(hit_ambiguous)) {
    stop("group_list contains semantically ambiguous labels that cannot be mapped ",
         "to Control/Tumor without curation: ",
         paste(sort(unique(raw[hit_ambiguous])), collapse = ", "),
         ". Curate the dataset's group encoding before enabling ROC analysis.")
  }

  out <- rep(NA_character_, length(key))
  out[key %in% CONTROL_LABELS] <- "Control"
  out[key %in% TUMOR_LABELS]   <- "Tumor"
  out
}

resolve_roc_group <- function(exp, env, allow_legacy_prefix = TRUE,
                              group_mapping = NULL,
                              group_mapping_source = NULL) {
  samples <- colnames(exp)
  mapping <- validate_group_mapping(group_mapping)

  if (exists("group_list", envir = env, inherits = FALSE)) {
    candidate <- get("group_list", envir = env, inherits = FALSE)

    if (is.data.frame(candidate) || is.matrix(candidate) || is.list(candidate)) {
      stop("group_list must be an atomic vector or factor, not a ",
           class(candidate)[1], ".")
    }

    nm <- names(candidate)
    if (!is.null(nm)) {
      # [FIX-1] a named group_list is matched BY NAME or not at all.
      if (anyDuplicated(nm)) {
        stop("group_list contains duplicate sample names: ",
             paste(unique(nm[duplicated(nm)]), collapse = ", "))
      }
      missing_names <- setdiff(samples, nm)
      if (length(missing_names) > 0L) {
        stop("Named group_list does not cover all expression sample IDs; ",
             "positional fallback is not permitted. Missing: ",
             paste(head(missing_names, 20L), collapse = ", "),
             if (length(missing_names) > 20L)
               sprintf(" (and %d more)", length(missing_names) - 20L) else "")
      }
      candidate <- candidate[samples]
    } else {
      if (length(candidate) != length(samples)) {
        stop("Unnamed group_list length (", length(candidate),
             ") does not match the number of expression matrix columns (",
             length(samples), ").")
      }
    }

    group <- canonical_group(candidate, group_mapping = mapping)
    source <- if (is.null(nm)) "RDA_group_list_positional" else "RDA_group_list_named"

  } else {
    if (!isTRUE(allow_legacy_prefix)) {
      stop("RDA has no group_list and legacy sample-ID group inference is disabled.")
    }
    # [FIX-2] strict pattern; "Control01" / "CRC01" no longer silently become Tumor.
    m <- regmatches(samples, regexec(LEGACY_GROUP_PATTERN, samples))
    prefix <- vapply(m, function(z) if (length(z) == 2L) toupper(z[2]) else NA_character_,
                     character(1))
    if (anyNA(prefix)) {
      bad <- samples[is.na(prefix)]
      stop("Legacy RDA has no group_list and sample IDs do not follow the strict ",
           "C<n>/N<n> convention required for group inference. Unmatched sample IDs: ",
           paste(head(bad, 20L), collapse = ", "),
           if (length(bad) > 20L) sprintf(" (and %d more)", length(bad) - 20L) else "")
    }
    group <- ifelse(prefix == "N", "Control", "Tumor")
    source <- "legacy_strict_C_N_prefix"
  }

  if (anyNA(group)) {
    stop("Group labels contain unsupported values; expected one of {",
         paste(c(CONTROL_LABELS, TUMOR_LABELS), collapse = ", "), "}.")
  }

  group <- factor(group, levels = c("Control", "Tumor"))
  counts <- table(group)
  if (any(counts < MIN_GROUP_N)) {
    stop("ROC analysis requires at least ", MIN_GROUP_N,
         " samples in each group; observed: ",
         paste(names(counts), counts, collapse = ", "))
  }

  smallest <- as.integer(min(counts))
  list(
    group = group,
    source = source,
    group_mapping_source = if (is.null(mapping)) "none_generic" else
      if (is.null(group_mapping_source) || !nzchar(trimws(group_mapping_source)))
        "explicit_authoritative_mapping" else as.character(group_mapping_source),
    counts = counts,
    smallest_group_n = smallest,
    small_sample_warning = smallest < SMALL_GROUP_WARN_N,
    small_sample_note = if (smallest < SMALL_GROUP_WARN_N)
      sprintf(paste0("Smallest group n = %d (< %d). ROC estimates are unstable and ",
                     "should be read descriptively, not as diagnostic performance."),
              smallest, SMALL_GROUP_WARN_N) else ""
  )
}

load_roc_input <- function(input_rda, allow_legacy_prefix = TRUE,
                           group_mapping = NULL,
                           group_mapping_source = NULL) {
  env <- new.env(parent = emptyenv())
  load(input_rda, envir = env)
  if (!exists("deg", envir = env, inherits = FALSE)) {
    stop("The 'deg' object does not exist in the provided RDA file.")
  }
  matrix_name <- if (exists("exp_model", envir = env, inherits = FALSE)) "exp_model" else "exp"
  if (!exists(matrix_name, envir = env, inherits = FALSE)) {
    stop("The RDA must contain exp_model (preferred) or legacy exp.")
  }
  exp <- as_numeric_matrix(get(matrix_name, envir = env, inherits = FALSE), matrix_name)
  deg <- as.data.frame(get("deg", envir = env, inherits = FALSE))
  if (is.null(rownames(deg))) stop("deg must have feature IDs as row names.")
  groups <- resolve_roc_group(
    exp, env,
    allow_legacy_prefix = allow_legacy_prefix,
    group_mapping = group_mapping,
    group_mapping_source = group_mapping_source
  )
  list(
    exp = exp,
    deg = deg,
    group = groups$group,
    group_source = groups$source,
    group_mapping_source = groups$group_mapping_source,
    n_control = unname(groups$counts[["Control"]]),
    n_tumor = unname(groups$counts[["Tumor"]]),
    smallest_group_n = groups$smallest_group_n,
    small_sample_warning = groups$small_sample_warning,
    small_sample_note = groups$small_sample_note,
    expression_source = matrix_name
  )
}

# ---------------------------------------------------------------- validation

check_probability <- function(x, label) {
  if (!is.finite(x) || x <= 0 || x >= 1) stop(label, " must be strictly between 0 and 1.")
}

validate_style <- function(args) {
  if (!is.finite(args$line_opacity) || args$line_opacity < 0 || args$line_opacity > 1) {
    stop("line_opacity must be between 0 and 1.")
  }
  if (!is.finite(args$auc_text_opacity) || args$auc_text_opacity < 0 || args$auc_text_opacity > 1) {
    stop("auc_text_opacity must be between 0 and 1.")
  }
  if (args$line_thickness <= 0 || args$diag_line_thickness <= 0 ||
      args$cutoff_text_size <= 0 || args$cutoff_point_size <= 0 ||
      args$width <= 0 || args$height <= 0) {
    stop("Plot sizes and dimensions must be positive.")
  }
}

# ---------------------------------------------------------------- DE candidates

# [P2-3] require_complete / min_features make this reusable for both contracts:
#   RF / SVM         -> require_complete = TRUE,  min_features = 2  (need a full matrix)
#   single-marker ROC-> require_complete = FALSE, min_features = 1  (per-marker missing
#                       filtering is done later by make_marker_roc()).
candidate_features <- function(deg, exp, is_use_padj, p_cut, logfc_cut,
                               require_complete = TRUE, min_features = 2L) {
  check_probability(p_cut, "pvalue_threshold")
  if (!is.finite(logfc_cut) || logfc_cut < 0) stop("log2fc_threshold must be non-negative.")
  p_column <- if (isTRUE(is_use_padj)) "adj.P.Val" else "P.Value"
  required <- c("logFC", p_column)
  missing <- setdiff(required, colnames(deg))
  if (length(missing) > 0) stop("deg is missing required column(s): ", paste(missing, collapse = ", "))

  common <- intersect(rownames(deg), rownames(exp))
  if (length(common) == 0) stop("No feature IDs overlap between deg and expression matrix.")
  tbl <- deg[common, , drop = FALSE]
  p_values <- as.numeric(tbl[[p_column]])
  logfc <- as.numeric(tbl[["logFC"]])
  # Comparator is strict on both sides; frontend/Methods text must read
  # "P < threshold" and "|log2FC| > threshold" to stay consistent.
  selected <- common[is.finite(p_values) & is.finite(logfc) &
                       p_values < p_cut & abs(logfc) > logfc_cut]
  if (length(selected) == 0) {
    stop("No candidate features meet the requested differential-analysis thresholds; thresholds were not relaxed.")
  }

  if (isTRUE(require_complete)) {
    genes <- selected[rowSums(is.finite(exp[selected, , drop = FALSE])) == ncol(exp)]
    excluded_missing <- length(selected) - length(genes)
  } else {
    genes <- selected
    excluded_missing <- 0L
  }
  if (length(genes) < min_features) {
    stop("Fewer than ", min_features,
         if (isTRUE(require_complete)) " complete" else "",
         " candidate feature(s) remain after differential filtering",
         if (isTRUE(require_complete)) " and missing-value exclusion." else ".")
  }
  list(
    genes = genes,
    metric = p_column,
    p_threshold = p_cut,
    log2fc_threshold = logfc_cut,
    threshold_comparator = sprintf("%s < %g & abs(logFC) > %g", p_column, p_cut, logfc_cut),
    candidate_before_missing_filter = length(selected),
    excluded_missing = excluded_missing
  )
}

# [P2-1 helper] marker/panel-specific small-sample audit computed on the samples
# actually used for THIS ROC, not on the whole-cohort counts.
small_sample_audit <- function(group_used) {
  counts <- table(droplevels(as.factor(group_used)))
  smallest <- if (length(counts)) as.integer(min(counts)) else 0L
  warn <- smallest < SMALL_GROUP_WARN_N
  list(
    smallest_group_n = smallest,
    small_sample_warning = warn,
    small_sample_note = if (warn)
      sprintf(paste0("Smallest group n = %d (< %d) for the samples used in this ",
                     "analysis. Estimates are unstable and should be read descriptively."),
              smallest, SMALL_GROUP_WARN_N) else ""
  )
}

# ---------------------------------------------------------------- ROC helpers

safe_ci_auc <- function(roc_obj) {
  ci <- try(pROC::ci.auc(roc_obj, method = "delong"), silent = TRUE)
  if (inherits(ci, "try-error")) {
    cond <- attr(ci, "condition")
    note <- if (!is.null(cond)) conditionMessage(cond) else as.character(ci)
    return(list(low = NA_real_, high = NA_real_, method = "DeLong",
                status = "unavailable", note = trimws(gsub("\\s+", " ", note))))
  }
  low <- suppressWarnings(as.numeric(ci[1]))
  high <- suppressWarnings(as.numeric(ci[3]))
  # ci.auc can return non-finite bounds WITHOUT throwing (e.g. degenerate small
  # samples); that is still an unavailable CI.
  if (!all(is.finite(c(low, high)))) {
    return(list(low = low, high = high, method = "DeLong",
                status = "unavailable",
                note = "DeLong interval returned non-finite bounds (likely degenerate/very small sample)."))
  }
  list(low = low, high = high, method = "DeLong", status = "available", note = "")
}

# Direction is FIXED from the differential-analysis logFC sign whenever available,
# so that AUC is NOT max(AUC, 1 - AUC).
#
# PRECONDITION (must be audited before deployment): canonical deg$logFC is defined
# as Cancer - Control, i.e. positive logFC == higher in Tumor. If that contract is
# not guaranteed across all ROC-eligible datasets, do NOT enable prespecified mode.
#
# The source label is deliberately "fixed_from_same_cohort_DE_logFC": the direction
# is fixed before ROC, but logFC comes from the SAME cohort, so this is not external
# prespecification.
#
# allow_auto_fallback = FALSE means: if a marker has no reliable logFC-based
# direction, it is refused rather than silently reverting to pROC "auto".
roc_direction_for <- function(gene, deg, allow_auto_fallback = FALSE) {
  no_dir <- is.null(deg) || !"logFC" %in% colnames(deg) || !gene %in% rownames(deg)
  lfc <- if (no_dir) NA_real_ else suppressWarnings(as.numeric(deg[gene, "logFC"]))
  if (no_dir || !is.finite(lfc) || lfc == 0) {
    if (isTRUE(allow_auto_fallback)) {
      return(list(direction = "auto", source = "data_driven_auto", ok = TRUE))
    }
    return(list(direction = NA_character_, source = "no_reliable_logFC_direction", ok = FALSE))
  }
  list(direction = if (lfc > 0) "<" else ">",
       source = "fixed_from_same_cohort_DE_logFC", ok = TRUE)
}

first_youden <- function(roc_obj) {
  best <- pROC::coords(
    roc_obj, x = "best", best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE
  )
  best <- as.data.frame(best)
  list(row = best[1, , drop = FALSE], n_ties = nrow(best))
}

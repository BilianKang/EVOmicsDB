# Auditable miRNA resolution helpers shared by Explore regulatory analyses.
# They align annotation resolution to the measured feature resolution.  They
# never infer mature-arm expression from a precursor measurement.

normalize_mirna_exact <- function(x) {
  x <- tolower(trimws(as.character(x)))
  sub("^hsa-", "", x, perl = TRUE)
}

mirna_arm <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(grepl("-3p$", x, perl = TRUE), "3p",
    ifelse(grepl("-5p$", x, perl = TRUE), "5p", "unarmed"))
}

canonical_hairpin_id <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- sub("^hsa-", "", x, perl = TRUE)
  # Only an explicit terminal arm suffix is removed.  No family-level or
  # fuzzy matching is performed, so mir-10a and mir-10b remain distinct.
  x <- sub("-(3p|5p)$", "", x, perl = TRUE)
  paste0("hsa-", x)
}

resolve_mirna_feature_key <- function(x, resolution) {
  resolution <- match.arg(resolution, c("mirbase_mature", "mirbase_precursor"))
  if (resolution == "mirbase_precursor") canonical_hairpin_id(x) else normalize_mirna_exact(x)
}

collapse_targets_to_hairpin <- function(associations) {
  required <- c("miRNA", "target_gene")
  if (!all(required %in% names(associations))) {
    stop("miRNA target annotations must contain: ", paste(required, collapse = ", "))
  }
  if (!nrow(associations)) {
    associations$hairpin_id <- character()
    associations$supporting_arm <- character()
    associations$arm_ambiguous <- logical()
    return(associations)
  }
  annotations <- associations
  annotations$hairpin_id <- canonical_hairpin_id(annotations$miRNA)
  annotations$annotation_arm <- mirna_arm(annotations$miRNA)
  annotations$target_gene <- toupper(sub("\\.[0-9]+$", "", trimws(as.character(annotations$target_gene))))

  # Vectorized grouping preserves first-occurrence order and the original
  # collapse contract while supporting whole-assay target backgrounds.
  row_keys <- paste(annotations$hairpin_id, annotations$target_gene, sep = "\t")
  keys <- unique(row_keys)
  group_index <- match(row_keys, keys)
  first <- match(keys, row_keys)
  out <- annotations[first, , drop = FALSE]
  has_3p <- keys %in% row_keys[annotations$annotation_arm == "3p"]
  has_5p <- keys %in% row_keys[annotations$annotation_arm == "5p"]
  support <- ifelse(has_3p & has_5p, "both", ifelse(has_3p, "3p", ifelse(has_5p, "5p", "unarmed")))
  hairpin_3p <- unique(annotations$hairpin_id[annotations$annotation_arm == "3p"])
  hairpin_5p <- unique(annotations$hairpin_id[annotations$annotation_arm == "5p"])
  out$miRNA <- out$hairpin_id
  out$supporting_arm <- support
  out$annotation_arm <- support
  out$arm_ambiguous <- out$hairpin_id %in% intersect(hairpin_3p, hairpin_5p)
  out$score <- NA_real_
  if ("score" %in% names(annotations)) {
    values <- suppressWarnings(as.numeric(annotations$score))
    grouped <- split(values, factor(group_index, levels = seq_along(keys)))
    out$score <- vapply(grouped, function(x) {
      if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
    }, numeric(1))
  }
  out$support_count <- tabulate(group_index, nbins = length(keys))
  if ("support_count" %in% names(annotations)) {
    values <- suppressWarnings(as.numeric(annotations$support_count))
    values[is.na(values)] <- 0
    out$support_count <- as.numeric(rowsum(values, group_index, reorder = FALSE))
  }
  out
}

resolve_mirna_target_annotations <- function(associations, feature_ids, resolution) {
  resolution <- match.arg(resolution, c("mirbase_mature", "mirbase_precursor"))
  if (!nrow(associations)) return(associations)
  if (resolution == "mirbase_precursor") {
    associations <- collapse_targets_to_hairpin(associations)
    feature_keys <- unique(canonical_hairpin_id(feature_ids))
    return(associations[associations$hairpin_id %in% feature_keys, , drop = FALSE])
  }
  associations$miRNA_match_key <- normalize_mirna_exact(associations$miRNA)
  feature_keys <- unique(normalize_mirna_exact(feature_ids))
  associations <- associations[associations$miRNA_match_key %in% feature_keys, , drop = FALSE]
  if (nrow(associations)) {
    associations$hairpin_id <- NA_character_
    associations$supporting_arm <- mirna_arm(associations$miRNA)
    associations$arm_ambiguous <- FALSE
  }
  associations
}

# Feature-level, mutually exclusive audit.  This is deliberately separate from
# pair-level collapse: a hairpin may have an exact/unarmed annotation for one
# target and arm annotations for another target, but the feature audit must
# retain all provenance flags and assign exactly one class.
mirna_feature_matching_audit <- function(feature_ids, associations, resolution,
                                         dataset_id = NA_character_) {
  resolution <- match.arg(resolution, c("mirbase_mature", "mirbase_precursor"))
  features <- unique(as.character(feature_ids))
  if (!length(features)) return(data.frame())
  if (is.null(associations) || !nrow(associations)) {
    return(data.frame(dataset_id = dataset_id, original_mirna_id = features,
      canonical_hairpin_id = if (resolution == "mirbase_precursor") canonical_hairpin_id(features) else NA_character_,
      has_exact_unarmed = FALSE, has_3p_annotation = FALSE, has_5p_annotation = FALSE,
      arm_ambiguous = FALSE, resolution_class = "unmatched", matched = FALSE,
      strict_eligible = FALSE, stringsAsFactors = FALSE))
  }
  ann <- associations
  ann$.key <- if (resolution == "mirbase_precursor") canonical_hairpin_id(ann$miRNA) else normalize_mirna_exact(ann$miRNA)
  ann$.arm <- mirna_arm(ann$miRNA)
  keys <- if (resolution == "mirbase_precursor") canonical_hairpin_id(features) else normalize_mirna_exact(features)
  out <- lapply(seq_along(features), function(i) {
    sub <- ann[ann$.key == keys[[i]], , drop = FALSE]
    arms <- unique(sub$.arm[sub$.arm %in% c("3p", "5p")])
    exact_unarmed <- any(sub$.arm == "unarmed")
    matched <- nrow(sub) > 0
    ambiguous <- all(c("3p", "5p") %in% arms)
    cls <- if (!matched) "unmatched" else if (exact_unarmed) "exact_unarmed" else if (identical(sort(arms), "3p")) "unique_3p" else if (identical(sort(arms), "5p")) "unique_5p" else if (ambiguous) "dual_arm" else "unresolved"
    data.frame(dataset_id = dataset_id, original_mirna_id = features[[i]],
      canonical_hairpin_id = if (resolution == "mirbase_precursor") keys[[i]] else NA_character_,
      has_exact_unarmed = exact_unarmed, has_3p_annotation = "3p" %in% arms,
      has_5p_annotation = "5p" %in% arms, arm_ambiguous = ambiguous,
      resolution_class = cls, matched = matched,
      strict_eligible = matched && (!identical(resolution, "mirbase_precursor") || !ambiguous),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

mirna_matching_audit <- function(feature_ids, associations, resolution) {
  resolution <- match.arg(resolution, c("mirbase_mature", "mirbase_precursor"))
  features <- unique(as.character(feature_ids))
  total <- length(features)
  feature_audit <- mirna_feature_matching_audit(features, associations, resolution)
  if (!total || !nrow(feature_audit)) {
    return(data.frame(total_mirna_features = total, exact_unarmed_matches = 0L,
      unique_single_arm_hairpins = 0L, dual_arm_hairpins = 0L,
      no_target_library_match = total, matched_total = 0L,
      match_rate = 0, strict_eligible_features = if (resolution == "mirbase_precursor") 0L else NA_integer_,
      strict_eligible_rate = if (resolution == "mirbase_precursor") 0 else NA_real_,
      stringsAsFactors = FALSE))
  }
  cls <- feature_audit$resolution_class
  strict <- feature_audit$strict_eligible
  data.frame(total_mirna_features = total,
    exact_unarmed_matches = sum(cls == "exact_unarmed"),
    unique_single_arm_hairpins = sum(cls %in% c("unique_3p", "unique_5p")),
    dual_arm_hairpins = sum(cls == "dual_arm"),
    no_target_library_match = sum(cls == "unmatched"),
    matched_total = sum(feature_audit$matched), match_rate = sum(feature_audit$matched) / total,
    strict_eligible_features = if (resolution == "mirbase_precursor") sum(strict) else NA_integer_,
    strict_eligible_rate = if (resolution == "mirbase_precursor") sum(strict) / total else NA_real_,
    stringsAsFactors = FALSE)
}

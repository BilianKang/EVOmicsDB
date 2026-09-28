# Shared helpers for the Explore network scripts.
# The file deliberately uses only base R + DBI/RSQLite when an SQLite target
# database is requested, so web jobs do not need live multiMiR access.

evomics_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) return(getwd())
  dirname(normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg[1])), mustWork = FALSE))
}

evomics_read_target_map <- function(path, mirna_column = "miRNA",
                                    target_column = "target_gene",
                                    species = "hsa", evidence = "all",
                                    minimum_score = 0, mirna_ids = NULL) {
  if (is.null(path) || !nzchar(trimws(path))) {
    stop("A local miRNA target database/map is required for this analysis.")
  }
  path <- normalizePath(path, mustWork = FALSE)
  if (!file.exists(path)) stop("Target database/map does not exist: ", path)
  evidence <- tolower(trimws(evidence %||% "all"))
  if (!evidence %in% c("all", "validated", "predicted")) {
    stop("Unsupported miRNA evidence mode: ", evidence)
  }
  minimum_score <- suppressWarnings(as.numeric(minimum_score %||% 0))
  if (!is.finite(minimum_score) || minimum_score < 0 || minimum_score > 100) {
    stop("minimum_score must be between 0 and 100.")
  }

  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("sqlite", "db", "sqlite3")) {
    if (!requireNamespace("DBI", quietly = TRUE) ||
        !requireNamespace("RSQLite", quietly = TRUE)) {
      stop("DBI and RSQLite are required to read the local target database.")
    }
    con <- DBI::dbConnect(RSQLite::SQLite(), path)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    tables <- DBI::dbListTables(con)
    table_name <- if ("web_target_map" %in% tables) "web_target_map" else
      if ("validated_targets" %in% tables) "validated_targets" else NULL
    if (is.null(table_name)) stop("SQLite target database has no supported target table.")
    fields <- DBI::dbListFields(con, table_name)
    mirna_field <- if ("mirna_id" %in% fields) "mirna_id" else
      if ("mature_mirna_id" %in% fields) "mature_mirna_id" else NULL
    target_field <- if ("target_symbol" %in% fields) "target_symbol" else
      if ("gene_symbol" %in% fields) "gene_symbol" else NULL
    if (is.null(mirna_field) || is.null(target_field)) {
      stop("SQLite target table must contain miRNA and target-symbol columns.")
    }
    score_select <- if ("mirdb_score" %in% fields) {
      paste0(", ", DBI::dbQuoteIdentifier(con, "mirdb_score"), " AS score")
    } else {
      ", NULL AS score"
    }
    sql <- sprintf(
      "SELECT %s AS miRNA, %s AS target_gene%s FROM %s",
      DBI::dbQuoteIdentifier(con, mirna_field),
      DBI::dbQuoteIdentifier(con, target_field),
      score_select,
      DBI::dbQuoteIdentifier(con, table_name)
    )
    where <- character()
    if ("species" %in% fields) {
      where <- c(where, paste0("species = ", DBI::dbQuoteString(con, species)))
    }
    if (evidence == "validated" && "is_validated" %in% fields) {
      where <- c(where, "coalesce(is_validated, 0) = 1")
    }
    if (evidence == "predicted" && "predicted_database_count" %in% fields) {
      where <- c(where, "coalesce(predicted_database_count, 0) > 0")
    }
    if (minimum_score > 0 && "mirdb_score" %in% fields) {
      where <- c(where, paste0("coalesce(mirdb_score, 0) >= ", minimum_score))
    }
    if (!is.null(mirna_ids) && length(mirna_ids)) {
      ids <- unique(tolower(trimws(as.character(mirna_ids))))
      ids <- ids[nzchar(ids)]
      if (length(ids)) {
        quoted <- paste(DBI::dbQuoteString(con, ids), collapse = ", ")
        # Use the resolved schema field here.  Some target indexes expose the
        # identifier as `mature_mirna_id` rather than `mirna_id`; hard-coding
        # the latter would make an otherwise valid map fail at runtime.
        where <- c(where, paste0(
          "lower(", DBI::dbQuoteIdentifier(con, mirna_field), ") IN (",
          quoted, ")"
        ))
      }
    }
    if (length(where)) sql <- paste0(sql, " WHERE ", paste(where, collapse = " AND "))
    raw <- DBI::dbGetQuery(con, sql)
  } else {
    sep <- if (ext == "csv") "," else "\t"
    raw <- utils::read.table(path, header = TRUE, sep = sep, quote = "\"",
                             comment.char = "", check.names = FALSE,
                             stringsAsFactors = FALSE)
    required <- c(mirna_column, target_column)
    if (!all(required %in% names(raw))) {
      stop("Target map is missing columns: ", paste(setdiff(required, names(raw)), collapse = ", "))
    }
    score <- if ("score" %in% names(raw)) {
      suppressWarnings(as.numeric(raw[["score"]]))
    } else {
      rep(NA_real_, nrow(raw))
    }
    if (minimum_score > 0 && any(is.finite(score))) {
      keep_score <- is.finite(score) & score >= minimum_score
      raw <- raw[keep_score, , drop = FALSE]
      score <- score[keep_score]
    }
    raw <- data.frame(miRNA = raw[[mirna_column]], target_gene = raw[[target_column]],
                      score = score, stringsAsFactors = FALSE)
  }
  if (!"score" %in% names(raw)) raw$score <- NA_real_
  raw <- raw[!is.na(raw$miRNA) & !is.na(raw$target_gene),
             c("miRNA", "target_gene", "score"), drop = FALSE]
  raw$miRNA <- tolower(sub("^hsa-", "", trimws(as.character(raw$miRNA))))
  raw$target_gene <- toupper(sub("\\.[0-9]+$", "", trimws(as.character(raw$target_gene))))
  raw <- raw[nzchar(raw$miRNA) & nzchar(raw$target_gene), , drop = FALSE]
  if (!ext %in% c("sqlite", "db", "sqlite3")) {
    score_order <- ifelse(is.finite(raw$score), raw$score, -Inf)
    raw <- raw[order(raw$miRNA, raw$target_gene, -score_order), , drop = FALSE]
    raw <- raw[!duplicated(paste(raw$miRNA, raw$target_gene, sep = "\t")), , drop = FALSE]
  }
  raw
}

# Read the normalized lncRNA interaction index.  The SQLite projection keeps
# source/evidence/PMID fields so the web result remains auditable.  A legacy
# master_lncRNA_target_database.csv is still accepted for backward
# compatibility with older deployments.
evomics_read_lnc_target_map <- function(path, target_type = NULL,
                                        evidence = "all", species = "hsa",
                                        source_databases = NULL,
                                        include_derived_cerna = FALSE) {
  if (is.null(path) || !nzchar(trimws(path))) {
    stop("A local lncRNA target database/map is required for this analysis.")
  }
  path <- normalizePath(path, mustWork = FALSE)
  if (!file.exists(path)) stop("lncRNA target database/map does not exist: ", path)
  evidence <- tolower(trimws(evidence %||% "all"))
  if (!evidence %in% c("all", "validated", "experimental", "clip", "predicted")) {
    stop("Unsupported evidence mode: ", evidence)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("sqlite", "db", "sqlite3")) {
    if (!requireNamespace("DBI", quietly = TRUE) ||
        !requireNamespace("RSQLite", quietly = TRUE)) {
      stop("DBI and RSQLite are required to read the local lncRNA target database.")
    }
    con <- DBI::dbConnect(RSQLite::SQLite(), path)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    if (!"interaction_map" %in% DBI::dbListTables(con)) {
      stop("SQLite lncRNA database has no interaction_map table.")
    }
    where <- c("(source_type = 'lncRNA' OR target_type = 'lncRNA')")
    params <- list()
    if (!is.null(species) && nzchar(species) && tolower(species) != "all") {
      if (tolower(species) %in% c("hsa", "human", "homo sapiens", "9606")) {
        where <- c(where, paste0(
          "lower(species) IN ('hsa','human','homo sapiens',",
          "'homo sapiens (human)','9606')"
        ))
      } else {
        where <- c(where, "lower(species) = lower(?)")
        params <- c(params, list(species))
      }
    }
    if (!is.null(target_type) && nzchar(target_type)) {
      # Target type is applied after orientation, but this predicate avoids
      # loading unrelated protein/RNA classes from the large NPInter index.
      target_types <- switch(tolower(target_type),
        mirna = "miRNA", mrna = "mRNA", protein = "Protein",
        target_type)
      where <- c(where, "(target_type = ? OR source_type = ?)")
      params <- c(params, list(target_types, target_types))
    }
    if (evidence != "all") {
      pattern <- switch(evidence,
        validated = "%valid%", experimental = "%experimental%",
        clip = "%clip%", predicted = "%predict%")
      where <- c(where, "lower(evidence_type) LIKE lower(?)")
      params <- c(params, list(pattern))
    }
    if (!is.null(source_databases) && length(source_databases)) {
      marks <- paste(rep("?", length(source_databases)), collapse = ",")
      where <- c(where, paste0("source_database IN (", marks, ")"))
      params <- c(params, as.list(source_databases))
    }
    if (!isTRUE(include_derived_cerna)) {
      where <- c(where, "lower(coalesce(metadata_json, '')) NOT LIKE '%derived_from%cerna_triplet%'")
    }
    sql <- paste0("SELECT source_database, source_version, species, source_id,
      source_type, target_id, target_type, relation_type, evidence_type,
      evidence_score, pmid, direction, metadata_json FROM interaction_map WHERE ",
      paste(where, collapse = " AND "))
    raw <- if (length(params)) {
      DBI::dbGetQuery(con, sql, params = params)
    } else {
      DBI::dbGetQuery(con, sql)
    }
    if (!nrow(raw)) return(data.frame(lncRNA=character(), target=character(),
      target_type=character(), source=character(), evidence_type=character(),
      evidence_score=numeric(), pmid=character(), stringsAsFactors = FALSE))
    # Orient reverse miRNA--lncRNA records (e.g. ENCORI) so the network has a
    # stable lncRNA -> target direction.
    source_is_lnc <- tolower(raw$source_type) == "lncrna"
    out <- data.frame(
      lncRNA = ifelse(source_is_lnc, raw$source_id, raw$target_id),
      target = ifelse(source_is_lnc, raw$target_id, raw$source_id),
      target_type = ifelse(source_is_lnc, raw$target_type, raw$source_type),
      source = raw$source_database,
      source_database = raw$source_database,
      source_version = raw$source_version,
      species = raw$species,
      evidence_type = raw$evidence_type,
      evidence_score = raw$evidence_score,
      pmid = raw$pmid,
      direction = raw$direction,
      relation_type = raw$relation_type,
      metadata_json = raw$metadata_json,
      stringsAsFactors = FALSE
    )
  } else {
    sep <- if (ext == "csv") "," else "\t"
    raw <- utils::read.table(path, header = TRUE, sep = sep, quote = "\"",
      comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
    required <- c("lncRNA", "target", "target_type")
    if (!all(required %in% names(raw))) {
      stop("Legacy lncRNA map is missing columns: ",
           paste(setdiff(required, names(raw)), collapse = ", "))
    }
    out <- raw
    out$source <- if ("source" %in% names(out)) out$source else "legacy_master_csv"
    out$evidence_type <- if ("evidence_type" %in% names(out)) out$evidence_type else "unspecified"
    out$evidence_score <- if ("evidence_score" %in% names(out)) out$evidence_score else NA_real_
    out$pmid <- if ("pmid" %in% names(out)) out$pmid else ""
  }
  out$target_type <- vapply(out$target_type, function(x) {
    x <- tolower(trimws(as.character(x)))
    if (grepl("mir", x)) "miRNA" else if (grepl("protein|peptide", x)) "Protein" else if (grepl("pcg|gene|tf|mrna", x)) "mRNA" else if (grepl("rna", x)) "RNA" else "Other"
  }, character(1))
  out$lncRNA <- toupper(sub("\\.[0-9]+$", "", trimws(as.character(out$lncRNA))))
  out$target <- ifelse(out$target_type == "miRNA",
    tolower(sub("^hsa-", "", trimws(as.character(out$target)))),
    toupper(sub("\\.[0-9]+$", "", trimws(as.character(out$target)))))
  out <- out[nzchar(out$lncRNA) & nzchar(out$target), , drop = FALSE]
  if (!is.null(target_type) && nzchar(target_type)) {
    want <- ifelse(tolower(target_type) == "mirna", "miRNA",
      ifelse(tolower(target_type) == "mrna", "mRNA",
        ifelse(tolower(target_type) == "protein", "Protein", target_type)))
    out <- out[out$target_type == want, , drop = FALSE]
  }
  unique(out)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

# Every standalone Explore entry point writes the same machine-readable
# execution contract.  A missing status file is never interpreted as a
# successful run by the Python dispatcher.
evomics_status_write <- function(path, analysis, status, code, message_text,
                                 counts = list(), artifacts = list()) {
  if (is.null(path) || !length(path) || !nzchar(trimws(path))) return(invisible(NULL))
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to write Explore status files")
  }
  parent <- dirname(normalizePath(path, mustWork = FALSE))
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    list(analysis = analysis, status = status, code = code,
         message = message_text, counts = counts, artifacts = artifacts),
    path, auto_unbox = TRUE, pretty = TRUE, null = "null"
  )
  invisible(path)
}

evomics_status_empty <- function(path, analysis, code, message_text,
                                 counts = list(), artifacts = list()) {
  evomics_status_write(path, analysis, "empty", code, message_text, counts, artifacts)
  message(message_text)
  quit(save = "no", status = 0L)
}

evomics_read_measured_universe_from_rda <- function(path) {
  env <- new.env(parent = emptyenv())
  tryCatch(load(path, envir = env), error = function(e) {
    stop("Could not load canonical RDA for measured universe: ", conditionMessage(e))
  })
  if (!exists("exp_raw", envir = env, inherits = FALSE)) {
    stop("Canonical RDA is missing exp_raw; measured-universe audit cannot fall back to filtered objects")
  }
  measured <- env$exp_raw
  ids <- if (is.null(dim(measured)) || is.null(rownames(measured))) character() else
    unique(trimws(as.character(rownames(measured))))
  ids <- ids[nzchar(ids) & !is.na(ids)]
  if (!length(ids)) stop("Canonical RDA exp_raw has no measured feature rownames")
  ids
}

# Schema v2's statistical model is fit on exp_model after feature filtering;
# deg is aligned to exp_model before it is saved.  GO/KEGG therefore continue
# to use this tested/enrichment-eligible universe, independently of the strict
# exp_raw measured-universe audit above.
evomics_read_analysis_universe_from_rda <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  candidate <- NULL
  if (exists("exp_model", envir = env, inherits = FALSE)) candidate <- env$exp_model
  if (is.null(candidate) && exists("exp", envir = env, inherits = FALSE)) candidate <- env$exp
  if (!is.null(candidate) && !is.null(rownames(candidate))) {
    ids <- unique(trimws(as.character(rownames(candidate))))
    ids <- ids[nzchar(ids) & !is.na(ids)]
    if (length(ids)) return(ids)
  }
  stop("Canonical RDA has no tested/enrichment-eligible feature universe")
}

# Backward-compatible name for existing GO/KEGG callers.  It deliberately
# retains the tested-universe contract; measured audits must call the strict
# helper above explicitly.
evomics_read_universe_from_rda <- function(path) {
  evomics_read_analysis_universe_from_rda(path)
}

# Resolve the canonical Schema v2 differential cache for a dataset.  Cache
# filenames are fingerprinted and multiple historical fingerprints may exist,
# so filename order/mtime is not an authoritative selector.  The cache's
# recorded input SHA-256 must match the configured dataset input.
evomics_sha256_file <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NA_character_)
  commands <- c("shasum", "sha256sum")
  for (cmd in commands) {
    if (nzchar(Sys.which(cmd))) {
      argv <- if (identical(cmd, "shasum")) c("-a", "256", path) else c(path)
      out <- tryCatch(system2(cmd, argv, stdout = TRUE, stderr = FALSE),
                      error = function(e) character())
      if (length(out)) {
        hash <- sub("\\s+.*$", "", out[[1]])
        if (grepl("^[0-9a-fA-F]{64}$", hash)) return(tolower(hash))
      }
    }
  }
  NA_character_
}

evomics_resolve_canonical_dataset_rda <- function(dataset_id, backend_root,
                                                   cache_dir = file.path(backend_root, "db", "cache")) {
  dataset_id <- as.character(as.integer(dataset_id))
  registry_path <- file.path(backend_root, "config", "differential_analysis.json")
  if (!file.exists(registry_path)) stop("Differential registry not found: ", registry_path)
  registry <- jsonlite::fromJSON(registry_path, simplifyDataFrame = TRUE)
  entries <- if (is.data.frame(registry)) registry else registry
  if (is.list(registry) && !is.null(registry$datasets)) entries <- registry$datasets
  if (is.data.frame(entries)) {
    row <- entries[as.character(entries$dataset_id) == dataset_id, , drop = FALSE]
  } else {
    row <- NULL
  }
  if (is.null(row) || !nrow(row) || !"data_file" %in% names(row)) {
    stop("No canonical data_file configured for dataset ", dataset_id)
  }
  configured <- as.character(row$data_file[[1]])
  data_path <- if (grepl("^/", configured)) configured else file.path(backend_root, configured)
  data_path <- normalizePath(data_path, mustWork = FALSE)
  if (!file.exists(data_path)) stop("Canonical dataset input does not exist: ", data_path)
  expected_hash <- evomics_sha256_file(data_path)
  candidates <- Sys.glob(file.path(cache_dir, paste0("differential_analysis_v2_", dataset_id, "_*.rda")))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) stop("No differential cache found for dataset ", dataset_id)
  exact_path_matches <- character()
  hash_matches <- character()
  for (candidate in candidates) {
    env <- new.env(parent = emptyenv())
    ok <- tryCatch({ load(candidate, envir = env); TRUE }, error = function(e) FALSE)
    if (!ok || !exists("preprocessing", envir = env, inherits = FALSE)) next
    prep <- env$preprocessing
    hash <- if (is.list(prep) && !is.null(prep$input_file_sha256))
      tolower(as.character(prep$input_file_sha256[[1]])) else NA_character_
    input_file <- if (is.list(prep) && !is.null(prep$input_file))
      normalizePath(as.character(prep$input_file[[1]]), mustWork = FALSE) else ""
    if (identical(input_file, data_path)) exact_path_matches <- c(exact_path_matches, candidate)
    if (!is.na(expected_hash) && identical(hash, expected_hash)) hash_matches <- c(hash_matches, candidate)
  }
  if (!is.na(expected_hash)) {
    matches <- unique(hash_matches)
    if (length(matches) == 0L) {
      stop("No differential cache matches configured input SHA-256 for dataset ", dataset_id)
    }
    if (length(matches) > 1L) {
      stop("Canonical differential cache is ambiguous for dataset ", dataset_id,
           "; more than one cache matches configured input SHA-256")
    }
  } else {
    # Only an actual platform inability to calculate SHA-256 permits this
    # legacy fallback.  Exact normalized paths are required; basename-only
    # matching is intentionally forbidden.
    matches <- unique(exact_path_matches)
    if (length(matches) == 0L) {
      stop("Unable to calculate SHA-256 and no normalized exact-path legacy cache matches dataset ", dataset_id)
    }
    if (length(matches) > 1L) {
      stop("Exact-path legacy differential cache is ambiguous for dataset ", dataset_id)
    }
  }
  normalizePath(matches[[1]], mustWork = TRUE)
}

# Apply an evidence-ranked cap independently within every measured miRNA.
# This is intentionally base-R so it can be regression-tested without loading
# the full clusterProfiler/KEGG stack.
evomics_select_targets_per_mirna <- function(associations, max_targets_per_mirna,
                                             mirna_col = "mature_mirna_id",
                                             score_col = "score",
                                             target_col = "target_symbol") {
  if (!nrow(associations)) return(associations)
  limit <- suppressWarnings(as.integer(max_targets_per_mirna))
  if (!is.finite(limit) || limit < 1L) stop("max_targets_per_mirna must be a positive integer")
  if (!all(c(mirna_col, score_col, target_col) %in% names(associations))) {
    stop("Target associations are missing required ranking columns")
  }
  score <- suppressWarnings(as.numeric(associations[[score_col]]))
  score[!is.finite(score)] <- -Inf
  ids <- as.character(associations[[mirna_col]])
  targets <- as.character(associations[[target_col]])
  ord <- order(ids, -score, targets, na.last = TRUE)
  ranked <- associations[ord, , drop = FALSE]
  ranked <- ranked[!duplicated(ranked[, c(mirna_col, target_col), drop = FALSE]), , drop = FALSE]
  keep <- ave(seq_len(nrow(ranked)), ranked[[mirna_col]], FUN = function(i) seq_along(i) <= limit)
  ranked[as.logical(keep), , drop = FALSE]
}

validate_mirna_audit <- function(summary, feature_audit) {
  required <- c("total_mirna_features", "exact_unarmed_matches",
    "unique_single_arm_hairpins", "dual_arm_hairpins", "no_target_library_match",
    "matched_total")
  if (!all(required %in% names(summary))) return(list(ok = FALSE, reason = "summary columns missing"))
  total <- as.integer(summary$total_mirna_features[[1]])
  exact <- as.integer(summary$exact_unarmed_matches[[1]])
  single <- as.integer(summary$unique_single_arm_hairpins[[1]])
  dual <- as.integer(summary$dual_arm_hairpins[[1]])
  unmatched <- as.integer(summary$no_target_library_match[[1]])
  matched <- as.integer(summary$matched_total[[1]])
  checks <- c(total_partition = total == exact + single + dual + unmatched,
    matched_partition = matched == total - unmatched,
    matched_classes = matched == exact + single + dual,
    no_unresolved = !any(feature_audit$resolution_class == "unresolved"),
    mutually_exclusive = all(vapply(split(feature_audit$resolution_class,
      feature_audit$original_mirna_id), function(x) length(unique(x)) == 1L, logical(1))))
  list(ok = all(checks), checks = checks,
       reason = if (all(checks)) "" else paste(names(checks)[!checks], collapse = ", "))
}

evomics_default_db_path <- function(script_dir, filename) {
  candidates <- c(
    file.path(script_dir, "../../../db/rda", filename),
    file.path(script_dir, "../../db/rda", filename),
    file.path(script_dir, "../db/rda", filename)
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) normalizePath(hit[1]) else NULL
}

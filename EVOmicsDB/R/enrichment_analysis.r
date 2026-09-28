#!/usr/bin/env Rscript

# EVOmicsDB enrichment-input preparation
# Maps gene/protein identifiers, or projects miRNAs to targets using the local
# SQLite target database. The actual GO/KEGG over-representation analysis is
# performed by the downstream plotting scripts.

required_packages <- c("argparse", "AnnotationDbi", "org.Hs.eg.db")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them in the R environment used by the EVOmicsDB backend."
  )
}

## ------------------------------------------------------------------ helpers
ensure_parent_directory <- function(path) {
  parent <- dirname(normalizePath(path, mustWork = FALSE))
  if (!dir.exists(parent)) {
    ok <- dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(parent)) {
      stop("Unable to create output directory: ", parent)
    }
  }
  invisible(NULL)
}

get_script_directory <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command, value = TRUE)
  if (length(file_argument) == 0L) {
    return(getwd())
  }
  dirname(normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_argument[[1L]]))))
}

detect_gene_id_type <- function(ids) {
  ids <- trimws(as.character(ids))
  ids <- ids[nzchar(ids)]
  if (length(ids) == 0L) {
    stop("Feature identifiers are empty.")
  }

  ensembl_ratio <- mean(grepl(
    "^ENSG[0-9]+(?:\\.[0-9]+)?$",
    ids,
    perl = TRUE
  ))
  entrez_ratio <- mean(grepl("^[0-9]+$", ids))
  uniprot_ratio <- mean(grepl(
    "^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])(?:-[0-9]+)?$",
    ids,
    perl = TRUE
  ))

  detected <- c(
    ENSEMBL = ensembl_ratio,
    ENTREZID = entrez_ratio,
    UNIPROT = uniprot_ratio
  )
  if (max(detected) >= 0.80) {
    return(names(which.max(detected)))
  }
  "SYMBOL"
}

normalize_mapping_key <- function(ids, id_type) {
  ids <- trimws(as.character(ids))
  if (id_type == "ENSEMBL") {
    ids <- sub("\\.[0-9]+$", "", ids)
  }
  if (id_type == "UNIPROT") {
    ids <- sub("-[0-9]+$", "", ids)
  }
  ids
}

map_gene_identifiers <- function(ids, id_type) {
  mapping_key <- normalize_mapping_key(ids, id_type)
  valid_keys <- unique(mapping_key[nzchar(mapping_key) & !is.na(mapping_key)])
  if (length(valid_keys) == 0L) {
    return(data.frame(
      mapping_key = character(),
      SYMBOL = character(),
      ENTREZID = character(),
      stringsAsFactors = FALSE
    ))
  }

  columns <- unique(c("SYMBOL", "ENTREZID"))
  mapped <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = valid_keys,
      columns = columns,
      keytype = id_type
    )
  )
  colnames(mapped)[colnames(mapped) == id_type] <- "mapping_key"
  if (!"SYMBOL" %in% colnames(mapped) && id_type == "SYMBOL") {
    mapped$SYMBOL <- mapped$mapping_key
  }
  if (!"ENTREZID" %in% colnames(mapped) && id_type == "ENTREZID") {
    mapped$ENTREZID <- mapped$mapping_key
  }
  mapped$mapping_key <- as.character(mapped$mapping_key)
  mapped$SYMBOL <- as.character(mapped$SYMBOL)
  mapped$ENTREZID <- as.character(mapped$ENTREZID)
  mapped <- mapped[
    !is.na(mapped$ENTREZID) & nzchar(mapped$ENTREZID),
    c("mapping_key", "SYMBOL", "ENTREZID"),
    drop = FALSE
  ]
  unique(mapped)
}

resolve_target_database <- function(explicit_path, script_directory) {
  candidates <- c(
    explicit_path,
    Sys.getenv("MIRNA_TARGET_DB", unset = ""),
    Sys.getenv("EVOMICS_TARGET_DB", unset = ""),
    file.path(script_directory, "../../../db/rda/mirna_targets.sqlite"),
    file.path(script_directory, "../../local_mirna_target_db/mirna_targets.sqlite")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop(
      "Local miRNA target database not found. Supply --target_db or set ",
      "EVOMICS_TARGET_DB. Checked: ",
      paste(candidates, collapse = ", ")
    )
  }
  normalizePath(existing[[1L]])
}

query_mirna_targets <- function(
  database_path,
  mirnas,
  resolution,
  evidence,
  min_mirdb_score,
  min_predicted_databases
) {
  database_packages <- c("DBI", "RSQLite")
  missing <- database_packages[
    !vapply(database_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0L) {
    stop(
      "miRNA target projection requires R package(s): ",
      paste(missing, collapse = ", ")
    )
  }

  mirnas <- unique(tolower(trimws(as.character(mirnas))))
  mirnas <- mirnas[nzchar(mirnas) & !is.na(mirnas)]
  if (length(mirnas) == 0L) {
    return(data.frame())
  }

  connection <- DBI::dbConnect(
    RSQLite::SQLite(),
    database_path,
    flags = RSQLite::SQLITE_RO
  )
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  DBI::dbExecute(connection, "PRAGMA query_only=ON")

  # Retrieve annotations at the configured measurement resolution. A precursor
  # stays one measured feature; mature-arm annotations supply a target union.
  annotated_ids <- DBI::dbGetQuery(connection,
    "SELECT DISTINCT mirna_id FROM web_target_map WHERE species = 'hsa'")$mirna_id
  feature_keys <- resolve_mirna_feature_key(mirnas, resolution)
  mirnas <- annotated_ids[
    resolve_mirna_feature_key(annotated_ids, resolution) %in% feature_keys
  ]
  if (!length(mirnas)) return(data.frame())

  required_columns <- c(
    "species", "mirna_id", "target_symbol", "target_entrez",
    "is_validated", "predicted_database_count", "mirdb_score", "sources"
  )
  table_info <- DBI::dbGetQuery(
    connection,
    "PRAGMA table_info(web_target_map)"
  )
  missing_columns <- setdiff(required_columns, table_info$name)
  if (length(missing_columns) > 0L) {
    stop(
      "The local target database has an incompatible web_target_map schema. ",
      "Missing: ", paste(missing_columns, collapse = ", ")
    )
  }

  evidence_sql <- switch(
    evidence,
    validated = "is_validated = 1",
    predicted = "predicted_database_count >= ?",
    recommended = paste0(
      "(is_validated = 1 OR mirdb_score >= ? ",
      "OR predicted_database_count >= ?)"
    )
  )
  evidence_parameters <- switch(
    evidence,
    validated = list(),
    predicted = list(as.integer(min_predicted_databases)),
    recommended = list(
      as.numeric(min_mirdb_score),
      as.integer(min_predicted_databases)
    )
  )

  chunks <- split(mirnas, ceiling(seq_along(mirnas) / 500L))
  results <- lapply(chunks, function(chunk) {
    placeholders <- paste(rep("?", length(chunk)), collapse = ",")
    sql <- paste0(
      "SELECT species, mirna_id, target_symbol, target_entrez, ",
      "is_validated, predicted_database_count, mirdb_score, sources ",
      "FROM web_target_map WHERE species = 'hsa' AND mirna_id IN (",
      placeholders, ") AND ", evidence_sql
    )
    DBI::dbGetQuery(
      connection,
      sql,
      params = c(as.list(chunk), evidence_parameters)
    )
  })

  result <- do.call(rbind, results)
  if (is.null(result) || nrow(result) == 0L) {
    return(data.frame())
  }
  result$mirna_id <- tolower(trimws(as.character(result$mirna_id)))
  result$target_symbol <- trimws(as.character(result$target_symbol))
  result$target_entrez <- trimws(as.character(result$target_entrez))
  result <- result[
    nzchar(result$mirna_id) &
      nzchar(result$target_symbol) &
      nzchar(result$target_entrez) &
      !is.na(result$target_entrez),
    ,
    drop = FALSE
  ]
  unique(result)
}

## ------------------------------------------------------------------ arguments
parser <- argparse::ArgumentParser(
  description = "Prepare mapped EVOmicsDB differential results for enrichment analysis."
)
parser$add_argument("--input_rda", type = "character", required = TRUE)
parser$add_argument("--output_rda", type = "character", required = TRUE)

# Backward-compatible web parameter: TRUE means a gene/protein-like feature ID.
parser$add_argument("--is_mrna", type = "logical", default = TRUE)

# Optional production parameters. Existing backend calls may omit all of them.
parser$add_argument("--input_type", type = "character", default = NULL)
parser$add_argument("--id_type", type = "character", default = "auto")
parser$add_argument("--species", type = "character", default = "hsa")
parser$add_argument("--min_mapping_rate", type = "numeric", default = 0.10)
parser$add_argument("--target_db", type = "character", default = NULL)
parser$add_argument("--target_evidence", type = "character", default = "recommended")
parser$add_argument("--min_mirdb_score", type = "numeric", default = 80)
parser$add_argument("--min_predicted_databases", type = "integer", default = 2)
parser$add_argument("--mirna_resolution", type = "character", default = NULL)

args <- parser$parse_args()

## ------------------------------------------------------------------ validation
if (!file.exists(args$input_rda)) {
  stop("Input RDA does not exist: ", args$input_rda)
}
if (!nzchar(trimws(args$output_rda))) {
  stop("--output_rda cannot be empty.")
}

input_type <- if (is.null(args$input_type)) {
  if (isTRUE(args$is_mrna)) "gene" else "mirna"
} else {
  tolower(trimws(args$input_type))
}
if (!input_type %in% c("gene", "mirna")) {
  stop("--input_type must be gene or mirna.")
}

id_type <- toupper(trimws(args$id_type))
supported_id_types <- c("AUTO", "SYMBOL", "ENSEMBL", "ENTREZID", "UNIPROT")
if (!id_type %in% supported_id_types) {
  stop(
    "--id_type must be one of: ",
    paste(supported_id_types, collapse = ", "),
    "."
  )
}
if (tolower(trimws(args$species)) != "hsa") {
  stop("This production script currently supports human data only (--species hsa).")
}
if (length(args$min_mapping_rate) != 1L ||
    is.na(args$min_mapping_rate) ||
    !is.finite(args$min_mapping_rate) ||
    args$min_mapping_rate < 0 ||
    args$min_mapping_rate > 1) {
  stop("--min_mapping_rate must be between 0 and 1.")
}
target_evidence <- tolower(trimws(args$target_evidence))
if (!target_evidence %in% c("recommended", "validated", "predicted")) {
  stop("--target_evidence must be recommended, validated, or predicted.")
}
if (length(args$min_mirdb_score) != 1L ||
    is.na(args$min_mirdb_score) ||
    !is.finite(args$min_mirdb_score) ||
    args$min_mirdb_score < 0) {
  stop("--min_mirdb_score must be a non-negative finite value.")
}
if (length(args$min_predicted_databases) != 1L ||
    is.na(args$min_predicted_databases) ||
    args$min_predicted_databases < 1L) {
  stop("--min_predicted_databases must be a positive integer.")
}

## ------------------------------------------------------------------ load input
input_environment <- new.env(parent = emptyenv())
load(args$input_rda, envir = input_environment)
if (!exists("deg", envir = input_environment, inherits = FALSE)) {
  stop("The input RDA is missing the required `deg` table.")
}

deg_input <- as.data.frame(input_environment$deg)
if (nrow(deg_input) == 0L ||
    is.null(rownames(deg_input)) ||
    anyNA(rownames(deg_input)) ||
    any(rownames(deg_input) == "")) {
  stop("The DEG table is empty or lacks valid feature row names.")
}
if (anyDuplicated(rownames(deg_input))) {
  stop("The DEG table contains duplicated feature identifiers.")
}

feature_ids <- if ("feature_id" %in% colnames(deg_input)) {
  trimws(as.character(deg_input$feature_id))
} else {
  rownames(deg_input)
}
if (anyNA(feature_ids) || any(feature_ids == "") || anyDuplicated(feature_ids)) {
  stop("The DEG table contains empty or duplicated feature_id values.")
}
if ("feature_id" %in% colnames(deg_input)) {
  base_deg <- deg_input
  base_deg$feature_id <- feature_ids
} else {
  base_deg <- data.frame(
    feature_id = feature_ids,
    deg_input,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
base_deg$.feature_order <- seq_len(nrow(base_deg))

## ------------------------------------------------------------------ mapping
if (input_type == "gene") {
  if (id_type == "AUTO") {
    id_type <- detect_gene_id_type(feature_ids)
  }
  mapping_key <- normalize_mapping_key(feature_ids, id_type)
  base_deg$mapping_key <- mapping_key
  identifier_map <- map_gene_identifiers(feature_ids, id_type)

  deg <- merge(
    base_deg,
    identifier_map,
    by = "mapping_key",
    all.x = TRUE,
    sort = FALSE
  )
  deg <- deg[order(deg$.feature_order), , drop = FALSE]
  mapped_features <- unique(deg$feature_id[
    !is.na(deg$ENTREZID) & nzchar(deg$ENTREZID)
  ])
  mapping_rate <- length(mapped_features) / length(feature_ids)
  universe_entrez <- unique(deg$ENTREZID[
    !is.na(deg$ENTREZID) & nzchar(deg$ENTREZID)
  ])

  mapping_audit <- data.frame(
    metric = c(
      "input_features", "mapped_features", "mapping_rate",
      "mapped_entrez_ids", "expanded_mapping_rows"
    ),
    value = c(
      length(feature_ids),
      length(mapped_features),
      mapping_rate,
      length(universe_entrez),
      nrow(deg)
    ),
    stringsAsFactors = FALSE
  )
  enrichment_metadata <- list(
    input_type = "gene",
    species = "hsa",
    id_type = id_type,
    mapping_is_one_to_many = nrow(deg) > nrow(base_deg)
  )
} else {
  # Fail closed for miRNA: omission must never silently turn a precursor into
  # a mature measurement. The API supplies the audited per-entry config value.
  if (is.null(args$mirna_resolution) ||
      !args$mirna_resolution %in% c("mirbase_mature", "mirbase_precursor")) {
    stop("miRNA input requires --mirna_resolution mirbase_mature or mirbase_precursor.")
  }
  source(file.path(get_script_directory(), "mirna_resolution_utils.R"))
  normalized_mirna <- tolower(trimws(sub("/.*$", "", feature_ids)))
  base_deg$miRNA <- resolve_mirna_feature_key(normalized_mirna, args$mirna_resolution)
  target_database <- resolve_target_database(
    args$target_db,
    get_script_directory()
  )
  target_map <- query_mirna_targets(
    database_path = target_database,
    mirnas = normalized_mirna,
    resolution = args$mirna_resolution,
    evidence = target_evidence,
    min_mirdb_score = args$min_mirdb_score,
    min_predicted_databases = args$min_predicted_databases
  )
  if (nrow(target_map) == 0L) {
    stop("No miRNA-target interactions matched the supplied features and evidence rule.")
  }

  target_map <- target_map[
    ,
    c(
      "mirna_id", "target_symbol", "target_entrez", "is_validated",
      "predicted_database_count", "mirdb_score", "sources"
    ),
    drop = FALSE
  ]
  colnames(target_map)[colnames(target_map) == "mirna_id"] <- "miRNA"
  colnames(target_map)[colnames(target_map) == "target_symbol"] <- "SYMBOL"
  colnames(target_map)[colnames(target_map) == "target_entrez"] <- "ENTREZID"

  # Collapse by stable target identifier, preserving all supporting evidence.
  # Screening occurs before collapse, identically for foreground and background.
  target_map$target_gene <- target_map$ENTREZID
  feature_resolution_audit <- mirna_feature_matching_audit(
    normalized_mirna, target_map, args$mirna_resolution
  )
  # Keep every measured ID, including distinct probes that share a normalized
  # annotation key; the feature audit helper itself deduplicates input keys.
  feature_resolution_audit <- feature_resolution_audit[
    match(normalized_mirna, feature_resolution_audit$original_mirna_id), , drop = FALSE
  ]
  feature_resolution_audit$original_mirna_id <- feature_ids
  rownames(feature_resolution_audit) <- NULL
  annotation_rows <- target_map
  target_map <- resolve_mirna_target_annotations(
    target_map, normalized_mirna, args$mirna_resolution
  )
  annotation_keys <- paste(
    resolve_mirna_feature_key(annotation_rows$miRNA, args$mirna_resolution),
    annotation_rows$ENTREZID, sep = "\t"
  )
  grouped_annotations <- split(seq_len(nrow(annotation_rows)), annotation_keys)
  target_map$miRNA <- resolve_mirna_feature_key(target_map$miRNA, args$mirna_resolution)
  output_keys <- paste(target_map$miRNA, target_map$ENTREZID, sep = "\t")
  finite_max <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
  summarize_column <- function(column, fun) {
    values <- vapply(grouped_annotations, function(index) fun(annotation_rows[[column]][index]),
      FUN.VALUE = if (column %in% c("sources", "miRNA")) character(1) else numeric(1))
    unname(values[output_keys])
  }
  target_map$is_validated <- summarize_column("is_validated", function(x) as.numeric(any(x == 1, na.rm = TRUE)))
  target_map$predicted_database_count <- summarize_column("predicted_database_count", finite_max)
  target_map$mirdb_score <- summarize_column("mirdb_score", finite_max)
  combine_unique <- function(x) paste(sort(unique(x)), collapse = ";")
  target_map$sources <- summarize_column("sources", combine_unique)
  target_map$supporting_annotation_mirnas <- summarize_column("miRNA", combine_unique)
  target_map$target_gene <- NULL

  deg <- merge(
    base_deg,
    target_map,
    by = "miRNA",
    all.x = TRUE,
    sort = FALSE
  )
  deg <- deg[order(deg$.feature_order), , drop = FALSE]
  mapped_features <- unique(deg$feature_id[
    !is.na(deg$ENTREZID) & nzchar(deg$ENTREZID)
  ])
  mapping_rate <- length(mapped_features) / length(feature_ids)
  universe_entrez <- unique(deg$ENTREZID[
    !is.na(deg$ENTREZID) & nzchar(deg$ENTREZID)
  ])

  mapping_audit <- data.frame(
    metric = c(
      "input_mirnas", "mapped_mirnas", "mapping_rate",
      "unique_target_entrez_ids", "mirna_target_rows"
    ),
    value = c(
      length(feature_ids),
      length(mapped_features),
      mapping_rate,
      length(universe_entrez),
      sum(!is.na(deg$ENTREZID) & nzchar(deg$ENTREZID))
    ),
    stringsAsFactors = FALSE
  )
  enrichment_metadata <- list(
    input_type = "mirna",
    mirna_resolution = args$mirna_resolution,
    precursor_target_rule = if (args$mirna_resolution == "mirbase_precursor")
      "union of qualifying unarmed/3p/5p annotations; no mature-arm expression inferred" else NULL,
    feature_resolution_audit = feature_resolution_audit,
    species = "hsa",
    target_database = target_database,
    target_evidence = target_evidence,
    min_mirdb_score = args$min_mirdb_score,
    min_predicted_databases = args$min_predicted_databases,
    target_projection_does_not_imply_regulatory_direction = TRUE
  )
}

if (mapping_rate < args$min_mapping_rate) {
  stop(
    "Identifier/target mapping rate is too low: ",
    format(100 * mapping_rate, digits = 4),
    "%; required minimum is ",
    format(100 * args$min_mapping_rate, digits = 4),
    "%. Verify --input_type and --id_type."
  )
}

deg$.feature_order <- NULL
rownames(deg) <- make.unique(as.character(deg$feature_id))
unmapped_features <- setdiff(feature_ids, mapped_features)
enrichment_metadata$mapping_rate <- mapping_rate
enrichment_metadata$unmapped_features <- unmapped_features

ensure_parent_directory(args$output_rda)
save(
  deg,
  universe_entrez,
  mapping_audit,
  enrichment_metadata,
  file = args$output_rda
)

message(
  "Enrichment input prepared: ", length(feature_ids), " input features; ",
  length(mapped_features), " mapped (",
  format(100 * mapping_rate, digits = 4), "%); ",
  length(universe_entrez), " unique background ENTREZ IDs; mode = ",
  input_type,
  if (input_type == "gene") paste0("; ID type = ", id_type, ".") else
    paste0("; target evidence = ", target_evidence, ".")
)

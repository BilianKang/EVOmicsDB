get_script_directory <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command, value = TRUE)
  if (!length(file_argument)) return(getwd())
  dirname(normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_argument[[1L]]))))
}

canonicalize_evomics_group <- function(input_environment, sample_names) {
  if (exists("group_list", envir = input_environment, inherits = FALSE)) {
    raw_group <- input_environment$group_list
    if (!is.null(names(raw_group))) {
      index <- match(sample_names, names(raw_group))
      if (anyNA(index)) stop("group_list does not contain every expression-matrix sample.")
      raw_group <- raw_group[index]
    } else if (length(raw_group) != length(sample_names)) {
      stop("Unnamed group_list length does not match the expression matrix.")
    }
    key <- tolower(trimws(as.character(raw_group)))
    group <- ifelse(
      key %in% c("control", "normal", "non-cancer", "noncancer", "n"),
      "Control",
      ifelse(key %in% c("case", "cancer", "tumor", "tumour", "c"), "Tumor", NA_character_)
    )
    if (anyNA(group)) stop("group_list contains unrecognized case/control labels.")
  } else {
    warning("Legacy RDA detected: inferring groups from ^N/^C sample prefixes.")
    normal <- grepl("^N", sample_names, ignore.case = TRUE)
    tumor <- grepl("^C", sample_names, ignore.case = TRUE)
    if (any(normal & tumor) || any(!normal & !tumor)) {
      stop("Legacy prefix grouping failed; the RDA must provide group_list.")
    }
    group <- ifelse(normal, "Control", "Tumor")
  }
  group <- factor(group, levels = c("Control", "Tumor"))
  if (any(table(group) == 0L)) stop("Both Control and Tumor groups are required.")
  names(group) <- sample_names
  group
}

load_evomics_model_data <- function(input_rda) {
  input_environment <- new.env(parent = emptyenv())
  load(input_rda, envir = input_environment)
  if (!exists("deg", envir = input_environment, inherits = FALSE)) {
    stop("The input RDA is missing the required `deg` table.")
  }
  expression <- if (exists("exp_model", envir = input_environment, inherits = FALSE)) {
    input_environment$exp_model
  } else if (exists("exp", envir = input_environment, inherits = FALSE)) {
    warning("Legacy RDA detected: using `exp` because `exp_model` is absent.")
    input_environment$exp
  } else {
    stop("The input RDA contains neither `exp_model` nor legacy `exp`.")
  }
  expression <- as.matrix(expression)
  storage.mode(expression) <- "numeric"
  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("The model expression matrix must have feature and sample names.")
  }

  group <- canonicalize_evomics_group(input_environment, colnames(expression))

  list(
    exp_model = expression,
    deg = as.data.frame(input_environment$deg),
    group = group,
    schema_version = if (exists("preprocessing", envir = input_environment, inherits = FALSE)) {
      input_environment$preprocessing$schema_version
    } else {
      "legacy"
    }
  )
}

load_evomics_visual_data <- function(input_rda) {
  input_environment <- new.env(parent = emptyenv())
  load(input_rda, envir = input_environment)
  if (!exists("deg", envir = input_environment, inherits = FALSE)) {
    stop("The input RDA is missing the required `deg` table.")
  }
  expression <- if (exists("exp_visual", envir = input_environment, inherits = FALSE)) {
    input_environment$exp_visual
  } else if (exists("exp", envir = input_environment, inherits = FALSE)) {
    warning("Legacy RDA detected: using `exp` because `exp_visual` is absent.")
    input_environment$exp
  } else {
    stop("The input RDA contains neither `exp_visual` nor legacy `exp`.")
  }
  expression <- as.matrix(expression)
  storage.mode(expression) <- "numeric"
  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("The visualization expression matrix must have feature and sample names.")
  }
  group <- canonicalize_evomics_group(input_environment, colnames(expression))
  list(
    exp_visual = expression,
    deg = as.data.frame(input_environment$deg),
    group = group,
    schema_version = if (exists("preprocessing", envir = input_environment, inherits = FALSE)) {
      input_environment$preprocessing$schema_version
    } else {
      "legacy"
    }
  )
}

feature_complete_cases <- function(expression, group, feature, min_per_group = 2L) {
  values <- as.numeric(expression[feature, ])
  keep <- is.finite(values) & !is.na(group)
  retained_group <- droplevels(group[keep])
  valid <- sum(keep) >= 4L && nlevels(retained_group) == 2L &&
    all(table(retained_group) >= min_per_group) && length(unique(values[keep])) >= 2L
  list(valid = valid, values = values[keep], group = retained_group, keep = keep)
}

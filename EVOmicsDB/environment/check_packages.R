#!/usr/bin/env Rscript
# Inspect dependencies without changing the user's installation.
args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args)) args[[1L]] else "demo"
profiles <- list(
  demo = c("argparse", "digest", "edgeR", "limma", "dplyr", "ggplot2", "ggrepel", "pROC", "scales"),
  figures = c("jsonlite", "data.table", "dplyr", "ggplot2", "patchwork", "scales", "png", "ragg", "svglite", "enrichplot")
)
if (!profile %in% names(profiles)) stop("Choose demo or figures; other analyses have script-specific dependencies.")
packages <- profiles[[profile]]
present <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
versions <- vapply(seq_along(packages), function(i) {
  if (present[[i]]) as.character(utils::packageVersion(packages[[i]])) else "MISSING"
}, character(1))
cat(R.version.string, "\n")
print(data.frame(package = packages, version = versions), row.names = FALSE)
if (any(!present)) {
  cat("Missing:", paste(packages[!present], collapse = ", "), "\n")
  quit(status = 1L)
}
cat("PASS: dependencies present for", profile, "(not a version-lock or font check).\n")

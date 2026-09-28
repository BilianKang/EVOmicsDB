#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(patchwork)
  library(tidyverse)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(gsub("~\\+~", " ", sub("^--file=", "", script_arg[1])), mustWork = FALSE))
} else getwd()
evomics_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) return(getwd())
  dirname(normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg[1])), mustWork = FALSE))
}
source(file.path(script_dir, "evomics_explore_utils.R"))
source(file.path(script_dir, "mirna_resolution_utils.R"))

parser <- ArgumentParser(
  description = "Matched-feature cross-omics effect-size correlation analysis."
)

# Existing web parameters: keep these names stable for backward compatibility.
parser$add_argument("--input_file1", type = "character", required = TRUE)
parser$add_argument("--input_file2", type = "character", required = TRUE)
parser$add_argument("--omics1_name", type = "character", required = TRUE)
parser$add_argument("--omics2_name", type = "character", required = TRUE)
parser$add_argument("--output_png", type = "character", required = TRUE)
parser$add_argument("--output_csv", type = "character", default = NULL)
parser$add_argument("--output_pdf", type = "character", default = NULL)

parser$add_argument("--main_title", type = "character", default = "")
parser$add_argument("--x_axis_title", type = "character", default = "log2FC")
parser$add_argument("--y_axis_title", type = "character", default = "Omics2 log2FC")
parser$add_argument("--main_title_size", type = "numeric", default = 14)
parser$add_argument("--axis_title_font_size", type = "numeric", default = 12)
parser$add_argument("--axis_text_font_size", type = "numeric", default = 10)

parser$add_argument("--dot_size", type = "numeric", default = 2)
parser$add_argument("--dot_alpha", type = "numeric", default = 0.7)
parser$add_argument(
  "--fill_colors", type = "character",
  default = "#1f78b4,#33a02c,#e31a1c,#ff7f00"
)
parser$add_argument("--border_color", type = "character", default = "black")
parser$add_argument("--border_size", type = "numeric", default = 0.5)
parser$add_argument("--is_show_border", type = "logical", default = TRUE)

parser$add_argument("--density_line_size", type = "numeric", default = 1)
parser$add_argument("--density_color1", type = "character", default = "#1f78b4")
parser$add_argument("--density_color2", type = "character", default = "#e31a1c")
parser$add_argument("--is_show_density", type = "logical", default = TRUE)

parser$add_argument("--is_show_legend", type = "logical", default = TRUE)
parser$add_argument("--legend_title", type = "character", default = "Direction concordance")
parser$add_argument("--legend_title_size", type = "numeric", default = 10)
parser$add_argument("--legend_text_size", type = "numeric", default = 9)
parser$add_argument("--legend_position", type = "character", default = "right")
parser$add_argument("--width", type = "numeric", default = 8.5)
parser$add_argument("--height", type = "numeric", default = 5)

# New backend parameters. Defaults preserve simple gene/protein workflows while
# preventing invalid row-order correlation for mRNA-miRNA comparisons.
parser$add_argument(
  "--match_mode", type = "character", default = "auto",
  choices = c("auto", "shared_id", "mirna_target", "lncrna_target")
)
parser$add_argument(
  "--mirna_omics", type = "character", default = "auto",
  choices = c("auto", "omics1", "omics2"),
  help = paste(
    "Which input is miRNA in mirna_target mode.",
    "Set this explicitly for web jobs; auto uses names and then target-map overlap."
  )
)
parser$add_argument("--target_map", type = "character", default = NULL)
parser$add_argument("--lncrna_target_db", type = "character", default = NULL)
parser$add_argument("--species", type = "character", default = "hsa")
parser$add_argument("--include_derived_cerna", type = "logical", default = FALSE)
parser$add_argument("--lncrna_omics", type = "character", default = "auto",
                    choices = c("auto", "omics1", "omics2"))
parser$add_argument("--mirna_column", type = "character", default = "miRNA")
parser$add_argument("--target_column", type = "character", default = "target_gene")
parser$add_argument("--mirna_resolution", type = "character", default = NULL,
                    choices = c("mirbase_mature", "mirbase_precursor"))
parser$add_argument(
  "--cor_method", type = "character", default = "spearman",
  choices = c("spearman", "pearson")
)
parser$add_argument("--min_pairs", type = "integer", default = 10)
parser$add_argument(
  "--contrast1_multiplier", type = "numeric", default = 1,
  help = "Use -1 if dataset 1 uses the opposite contrast direction."
)
parser$add_argument(
  "--contrast2_multiplier", type = "numeric", default = 1,
  help = "Use -1 if dataset 2 uses the opposite contrast direction."
)
parser$add_argument("--output_pairs_csv", type = "character", default = NULL)
parser$add_argument("--status_json", type = "character", default = NULL,
                    help = "Structured Explore execution status JSON")

args <- parser$parse_args()
if (!is.null(args$status_json) && file.exists(args$status_json)) unlink(args$status_json)
write_cross_status <- function(status, code, message_text, counts = list(), artifacts = list()) {
  evomics_status_write(args$status_json, "cross_omics_trend_analysis", status, code,
                       message_text, counts, artifacts)
}
cross_empty <- function(code, message_text, counts = list(), artifacts = list()) {
  # Preserve the response download contract for early empty exits.  The pair
  # table is still written separately below; this typed main CSV is what the
  # HTTP response exposes as ``data``.
  if (!is.null(args$output_csv) && nzchar(trimws(args$output_csv)) &&
      !file.exists(args$output_csv)) {
    empty_summary <- data.frame(
      pair_id = character(), x = numeric(), y = numeric(),
      stringsAsFactors = FALSE
    )
    utils::write.table(empty_summary, file = args$output_csv, sep = ",",
                       quote = FALSE, row.names = FALSE)
    artifacts <- c(list(csv = args$output_csv), artifacts)
  }
  evomics_status_empty(args$status_json, "cross_omics_trend_analysis", code,
                       message_text, counts, artifacts)
}
if (!is.finite(args$contrast1_multiplier) || !is.finite(args$contrast2_multiplier) ||
    !args$contrast1_multiplier %in% c(-1, 1) || !args$contrast2_multiplier %in% c(-1, 1)) {
  stop("Contrast multipliers must be either +1 or -1.")
}
if (args$min_pairs < 3) stop("--min_pairs must be at least 3.")

read_omics_rda <- function(path, multiplier = 1) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists("deg", envir = env, inherits = FALSE)) {
    stop("Object 'deg' was not found in: ", path)
  }
  df <- as.data.frame(env$deg)
  if (!("logFC" %in% names(df))) {
    stop("Column 'logFC' was not found in deg from: ", path)
  }
  data <- df %>%
    mutate(
      ID = trimws(as.character(rownames(df))),
      logFC = suppressWarnings(as.numeric(logFC)) * multiplier
    ) %>%
    filter(!is.na(ID), nzchar(ID), is.finite(logFC)) %>%
    arrange(desc(abs(logFC))) %>%
    distinct(ID, .keep_all = TRUE) %>%
    select(ID, logFC, everything())
  preprocessing <- if (exists("preprocessing", envir = env, inherits = FALSE)) env$preprocessing else list()
  list(
    data = data,
    organism = as.character(preprocessing$organism %||% ""),
    taxid = as.character(preprocessing$organism_taxid %||% "")
  )
}

normalise_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  sub("\\.[0-9]+$", "", x)
}
normalise_lncrna <- function(x) toupper(sub("\\.[0-9]+$", "", trimws(as.character(x))))

omics1_input <- read_omics_rda(args$input_file1, args$contrast1_multiplier)
omics2_input <- read_omics_rda(args$input_file2, args$contrast2_multiplier)
if (nzchar(omics1_input$taxid) && nzchar(omics2_input$taxid) &&
    omics1_input$taxid != omics2_input$taxid) {
  stop("Cross-omics inputs have different organism TaxIDs: ",
       omics1_input$taxid, " versus ", omics2_input$taxid, ".")
}
omics1 <- omics1_input$data
omics2 <- omics2_input$data

name_is_mirna1 <- grepl("mirna|micro.?rna", args$omics1_name, ignore.case = TRUE)
name_is_mirna2 <- grepl("mirna|micro.?rna", args$omics2_name, ignore.case = TRUE)
name_is_lncrna1 <- grepl("lncrna|long.?non.?coding", args$omics1_name, ignore.case = TRUE)
name_is_lncrna2 <- grepl("lncrna|long.?non.?coding", args$omics2_name, ignore.case = TRUE)
mode <- args$match_mode
mirna_side_used <- NA_character_
if (mode == "auto") {
  # A supplied target map is an explicit signal that projection was requested.
  # This avoids silently falling back to shared-ID matching when display names
  # are accession-only labels such as "GSE220445".
  has_target_map <- !is.null(args$target_map) && nzchar(trimws(args$target_map))
  has_lnc_map <- !is.null(args$lncrna_target_db) && nzchar(trimws(args$lncrna_target_db))
  mode <- if (has_lnc_map && xor(name_is_lncrna1, name_is_lncrna2)) {
    "lncrna_target"
  } else if (has_target_map || xor(name_is_mirna1, name_is_mirna2)) {
    "mirna_target"
  } else {
    "shared_id"
  }
}

if (mode == "shared_id") {
  paired <- inner_join(
    omics1 %>% transmute(pair_id = ID, x = logFC),
    omics2 %>% transmute(pair_id = ID, y = logFC),
    by = "pair_id"
  ) %>%
    mutate(n_targets = NA_integer_)
  unit_label <- "matched molecular features"
} else if (mode == "lncrna_target") {
  taxid <- if (nzchar(omics1_input$taxid)) omics1_input$taxid else omics2_input$taxid
  if (nzchar(taxid) && taxid != "9606") {
    stop("The deployed lncRNA target index is human-only; received TaxID ", taxid, ".")
  }
  lnc_side <- args$lncrna_omics
  if (lnc_side == "auto") {
    if (xor(name_is_lncrna1, name_is_lncrna2)) {
      lnc_side <- if (name_is_lncrna1) "omics1" else "omics2"
    } else {
      stop("Cannot infer which input is lncRNA. Set --lncrna_omics omics1 or omics2 explicitly.")
    }
  }
  lnc_map <- evomics_read_lnc_target_map(
    args$lncrna_target_db, evidence = "all", species = args$species,
    include_derived_cerna = args$include_derived_cerna
  ) %>%
    transmute(lncRNA = normalise_lncrna(lncRNA), target_gene = normalise_gene(target)) %>%
    filter(nzchar(lncRNA), nzchar(target_gene)) %>% distinct()
  if (lnc_side == "omics1") {
    lnc_effect <- omics1 %>% transmute(lncRNA = normalise_lncrna(ID), lnc_fc = logFC)
    target_effect <- omics2 %>% transmute(target_gene = normalise_gene(ID), target_fc = logFC)
    paired <- lnc_map %>% inner_join(lnc_effect, by = "lncRNA") %>%
      inner_join(target_effect, by = "target_gene") %>% group_by(lncRNA) %>%
      summarise(x = first(lnc_fc), y = median(target_fc, na.rm = TRUE),
                n_targets = n_distinct(target_gene), .groups = "drop") %>% rename(pair_id = lncRNA)
  } else {
    target_effect <- omics1 %>% transmute(target_gene = normalise_gene(ID), target_fc = logFC)
    lnc_effect <- omics2 %>% transmute(lncRNA = normalise_lncrna(ID), lnc_fc = logFC)
    paired <- lnc_map %>% inner_join(target_effect, by = "target_gene") %>%
      inner_join(lnc_effect, by = "lncRNA") %>% group_by(lncRNA) %>%
      summarise(x = median(target_fc, na.rm = TRUE), y = first(lnc_fc),
                n_targets = n_distinct(target_gene), .groups = "drop") %>% rename(pair_id = lncRNA)
  }
  mirna_side_used <- lnc_side
  unit_label <- "lncRNAs with detected target features"
} else {
  if (is.null(args$mirna_resolution) || !nzchar(args$mirna_resolution)) {
    stop("miRNA target projection requires an explicit --mirna_resolution from dataset configuration")
  }
  taxid <- if (nzchar(omics1_input$taxid)) omics1_input$taxid else omics2_input$taxid
  if (nzchar(taxid) && taxid != "9606") {
    stop("The deployed miRNA target index is human-only; received TaxID ", taxid, ".")
  }
  target_map <- evomics_read_target_map(args$target_map, args$mirna_column, args$target_column) %>%
    transmute(miRNA = as.character(miRNA), target_gene = normalise_gene(target_gene)) %>%
    filter(nzchar(miRNA), nzchar(target_gene)) %>% distinct()

  mirna_side <- args$mirna_omics
  if (mirna_side == "auto") {
    if (xor(name_is_mirna1, name_is_mirna2)) {
      mirna_side <- if (name_is_mirna1) "omics1" else "omics2"
    } else {
      map_mirnas <- unique(target_map$miRNA)
      overlap1 <- sum(unique(resolve_mirna_feature_key(omics1$ID, args$mirna_resolution)) %in% unique(resolve_mirna_feature_key(target_map$miRNA, args$mirna_resolution)))
      overlap2 <- sum(unique(resolve_mirna_feature_key(omics2$ID, args$mirna_resolution)) %in% unique(resolve_mirna_feature_key(target_map$miRNA, args$mirna_resolution)))
      if (max(overlap1, overlap2) == 0 || overlap1 == overlap2) {
        stop(
          "Cannot infer which input is miRNA from names or target-map overlap. ",
          "Set --mirna_omics omics1 or --mirna_omics omics2 explicitly."
        )
      }
      mirna_side <- if (overlap1 > overlap2) "omics1" else "omics2"
      message(
        "Auto-inferred miRNA input from target-map overlap: ", mirna_side,
        " (omics1=", overlap1, ", omics2=", overlap2, ")."
      )
    }
  }

  # Align annotations to the configured measurement resolution.  In precursor
  # mode this collapses 3p/5p records to canonical hairpins and retains the
  # measured precursor as the expression key.
  mirna_ids <- if (mirna_side == "omics1") omics1$ID else omics2$ID
  target_map <- resolve_mirna_target_annotations(
    target_map %>% transmute(miRNA = miRNA, target_gene = target_gene),
    mirna_ids, args$mirna_resolution
  ) %>% mutate(miRNA = resolve_mirna_feature_key(miRNA, args$mirna_resolution)) %>% distinct()

  if (mirna_side == "omics1") {
    mirna_effect <- omics1 %>%
      transmute(miRNA = resolve_mirna_feature_key(ID, args$mirna_resolution), mirna_fc = logFC)
    target_effect <- omics2 %>%
      transmute(target_gene = normalise_gene(ID), target_fc = logFC)
    paired <- target_map %>%
      inner_join(mirna_effect, by = "miRNA") %>%
      inner_join(target_effect, by = "target_gene") %>%
      group_by(miRNA) %>%
      summarise(
        x = first(mirna_fc),
        y = median(target_fc, na.rm = TRUE),
        n_targets = n_distinct(target_gene),
        .groups = "drop"
      ) %>%
      rename(pair_id = miRNA)
  } else {
    target_effect <- omics1 %>%
      transmute(target_gene = normalise_gene(ID), target_fc = logFC)
    mirna_effect <- omics2 %>%
      transmute(miRNA = resolve_mirna_feature_key(ID, args$mirna_resolution), mirna_fc = logFC)
    paired <- target_map %>%
      inner_join(target_effect, by = "target_gene") %>%
      inner_join(mirna_effect, by = "miRNA") %>%
      group_by(miRNA) %>%
      summarise(
        x = median(target_fc, na.rm = TRUE),
        y = first(mirna_fc),
        n_targets = n_distinct(target_gene),
        .groups = "drop"
      ) %>%
      rename(pair_id = miRNA)
  }
  mirna_side_used <- mirna_side
  unit_label <- "miRNAs with detected target genes"
}

paired <- paired %>% filter(is.finite(x), is.finite(y))
if (nrow(paired) < args$min_pairs) {
  if (!is.null(args$output_pairs_csv) && nzchar(trimws(args$output_pairs_csv))) {
    write.csv(paired, args$output_pairs_csv, row.names = FALSE, quote = TRUE)
  }
  cross_empty(
    "insufficient_matched_units",
    paste0("Only ", nrow(paired), " valid matched units were found; at least ",
           args$min_pairs, " are required. Check ID types, target mapping, and contrasts."),
    counts = list(matched_units = nrow(paired), min_pairs = args$min_pairs),
    artifacts = list(pairs_csv = args$output_pairs_csv)
  )
}
if (sd(paired$x) == 0 || sd(paired$y) == 0) {
  if (!is.null(args$output_pairs_csv) && nzchar(trimws(args$output_pairs_csv))) {
    write.csv(paired, args$output_pairs_csv, row.names = FALSE, quote = TRUE)
  }
  cross_empty(
    "zero_variance",
    "Correlation is undefined because one matched effect-size vector has zero variance. Check the target map, ID mapping, and differential-analysis input.",
    counts = list(matched_units = nrow(paired)),
    artifacts = list(pairs_csv = args$output_pairs_csv)
  )
}

cor_res <- if (args$cor_method == "spearman") {
  suppressWarnings(cor.test(paired$x, paired$y, method = "spearman", exact = FALSE))
} else {
  cor.test(paired$x, paired$y, method = "pearson")
}

group_levels <- c(
  "Both up", "Both down",
  paste0(args$omics1_name, " up / ", args$omics2_name, " down"),
  paste0(args$omics1_name, " down / ", args$omics2_name, " up")
)
paired <- paired %>%
  mutate(
    Direction = case_when(
      x >= 0 & y >= 0 ~ group_levels[1],
      x < 0 & y < 0 ~ group_levels[2],
      x >= 0 & y < 0 ~ group_levels[3],
      TRUE ~ group_levels[4]
    ),
    Direction = factor(Direction, levels = group_levels)
  )

fill_values <- trimws(strsplit(args$fill_colors, ",", fixed = TRUE)[[1]])
if (length(fill_values) != 4) {
  stop("--fill_colors must contain exactly four comma-separated colors.")
}
names(fill_values) <- group_levels

legend_pos <- if (!args$is_show_legend || tolower(args$legend_position) == "none") {
  "none"
} else {
  args$legend_position
}
main_title <- if (nzchar(args$main_title)) {
  args$main_title
} else {
  paste0("Cross-omics effect-size concordance: ", args$omics1_name, " vs ", args$omics2_name)
}
x_title <- if (args$x_axis_title == "log2FC") paste0(args$omics1_name, " log2FC") else args$x_axis_title
y_title <- if (grepl("Random Index", args$y_axis_title, ignore.case = TRUE) ||
               args$y_axis_title == "Omics2 log2FC") {
  paste0(args$omics2_name, " log2FC")
} else {
  args$y_axis_title
}
if (mode == "mirna_target") {
  if (mirna_side_used == "omics1") {
    if (args$x_axis_title == "log2FC") x_title <- paste0(args$omics1_name, " log2FC")
    if (args$y_axis_title == "Omics2 log2FC") y_title <- "Median target-gene log2FC"
  } else {
    if (args$x_axis_title == "log2FC") x_title <- "Median target-gene log2FC"
    if (args$y_axis_title == "Omics2 log2FC") y_title <- paste0(args$omics2_name, " log2FC")
  }
}
if (mode == "lncrna_target") {
  if (mirna_side_used == "omics1") {
    if (args$x_axis_title == "log2FC") x_title <- paste0(args$omics1_name, " log2FC")
    if (args$y_axis_title == "Omics2 log2FC") y_title <- "Median target-feature log2FC"
  } else {
    if (args$x_axis_title == "log2FC") x_title <- "Median target-feature log2FC"
    if (args$y_axis_title == "Omics2 log2FC") y_title <- paste0(args$omics2_name, " log2FC")
  }
}

p <- ggplot(paired, aes(x = x, y = y, fill = Direction)) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed", colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", colour = "grey55")

if (args$is_show_border) {
  p <- p + geom_point(
    shape = 21, size = args$dot_size, alpha = args$dot_alpha,
    colour = args$border_color, stroke = args$border_size
  )
} else {
  p <- p + geom_point(
    shape = 21, size = args$dot_size, alpha = args$dot_alpha,
    colour = NA, stroke = 0
  )
}

stat_symbol <- if (args$cor_method == "spearman") "rho" else "r"
stat_label <- paste0(
  stat_symbol, " = ", formatC(unname(cor_res$estimate), digits = 3, format = "f"),
  "\nP = ", format.pval(cor_res$p.value, digits = 3, eps = 1e-300),
  "\nn = ", nrow(paired)
)

p <- p +
  annotate(
    "text", x = -Inf, y = Inf, label = stat_label,
    hjust = -0.08, vjust = 1.15, size = 3.8
  ) +
  scale_fill_manual(
    values = fill_values,
    labels = function(z) paste0(z, " (n=", table(paired$Direction)[z], ")"),
    drop = FALSE
  ) +
  labs(
    x = x_title, y = y_title, fill = args$legend_title,
    title = main_title,
    caption = paste0(
      "Matched units: ", unit_label,
      "; contrasts were harmonized before correlation.",
      if (mode == "mirna_target") " miRNA effects were projected through curated target relationships; this is descriptive effect-size concordance, not a sample-level paired association test." else ""
    )
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(size = args$main_title_size, hjust = 0.5),
    axis.title = element_text(size = args$axis_title_font_size),
    axis.text = element_text(size = args$axis_text_font_size),
    legend.title = element_text(size = args$legend_title_size),
    legend.text = element_text(size = args$legend_text_size),
    legend.position = legend_pos,
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    plot.caption = element_text(hjust = 0, margin = margin(t = 3)),
    plot.margin = margin(6, 8, 6, 12)
  )

final_plot <- p
if (args$is_show_density && nrow(paired) >= 5) {
  density_df <- bind_rows(
    data.frame(logFC = paired$x, Omics = args$omics1_name),
    data.frame(logFC = paired$y, Omics = args$omics2_name)
  )
  density_values <- setNames(
    c(args$density_color1, args$density_color2),
    c(args$omics1_name, args$omics2_name)
  )
  p_density <- ggplot(density_df, aes(x = logFC, colour = Omics)) +
    geom_density(linewidth = args$density_line_size, adjust = 1) +
    scale_colour_manual(values = density_values) +
    labs(x = "log2FC distribution", y = "Relative density", colour = NULL) +
    theme_classic(base_family = "Arial") +
    theme(
      axis.title = element_text(size = max(7, args$axis_title_font_size - 2)),
      axis.text = element_text(size = max(6, args$axis_text_font_size - 2)),
      legend.position = "top",
      legend.text = element_text(size = args$legend_text_size),
      plot.margin = margin(2, 3, 2, 3)
    )
  final_plot <- p / p_density + plot_layout(heights = c(4, 1.15))
}

ggsave(
  args$output_png, plot = final_plot, width = args$width, height = args$height,
  units = "in", dpi = 300, bg = "white"
)
if (!is.null(args$output_pdf) && nzchar(trimws(args$output_pdf))) {
  ggsave(
    args$output_pdf, plot = final_plot, width = args$width, height = args$height,
    units = "in", device = cairo_pdf, bg = "white"
  )
}

if (!is.null(args$output_csv) && nzchar(trimws(args$output_csv))) {
  summary_df <- data.frame(
    omics1 = args$omics1_name,
    omics2 = args$omics2_name,
    match_mode = mode,
    correlation_method = args$cor_method,
    correlation = unname(cor_res$estimate),
    R_squared = unname(cor_res$estimate)^2,
    R_squared_note = if (args$cor_method == "spearman") {
      "squared Spearman rho, not coefficient of determination"
    } else {
      "squared Pearson correlation coefficient; descriptive, not model R-squared"
    },
    p_value = cor_res$p.value,
    n_matched_units = nrow(paired),
    omics1_features = nrow(omics1),
    omics2_features = nrow(omics2),
    contrast1_multiplier = args$contrast1_multiplier,
    contrast2_multiplier = args$contrast2_multiplier
  )
  summary_df$target_map <- if (mode == "mirna_target") {
    normalizePath(args$target_map, mustWork = FALSE)
  } else if (mode == "lncrna_target") {
    normalizePath(args$lncrna_target_db, mustWork = FALSE)
  } else NA_character_
  summary_df$min_pairs <- args$min_pairs
  summary_df$mirna_omics <- if (mode %in% c("mirna_target", "lncrna_target")) {
    mirna_side_used
  } else {
    NA_character_
  }
  summary_df$mirna_resolution <- if (mode == "mirna_target") args$mirna_resolution else NA_character_
  summary_df$resolution_source <- if (mode == "mirna_target") "dataset_config" else NA_character_
  summary_df$correlation_p_value_note <- "Descriptive effect-size concordance; not a sample-level paired multi-omics association test."
  write.csv(summary_df, args$output_csv, row.names = FALSE, quote = TRUE)
}
if (!is.null(args$output_pairs_csv) && nzchar(trimws(args$output_pairs_csv))) {
  write.csv(paired, args$output_pairs_csv, row.names = FALSE, quote = TRUE)
}

write_cross_status(
  "success", "completed", "Cross-omics trend analysis generated successfully",
  counts = list(matched_units = nrow(paired), mode = mode),
  artifacts = list(png = args$output_png, csv = args$output_csv, pdf = args$output_pdf,
                   pairs_csv = args$output_pairs_csv)
)

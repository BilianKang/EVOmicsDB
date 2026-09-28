# Reproducing EVOmicsDB analyses

## Environment

The Python tools use Python 3.10+ and its standard library. The statistical revision was run with R 4.5.2. [environment/R-packages.tsv](../environment/R-packages.tsv) records the R package versions. It is not a complete environment lock, and it is not a list of packages that every user must install.

The demo needs `argparse`, `digest`, `edgeR`, `limma`, `dplyr`, `ggplot2`, `ggrepel`, `pROC` and `scales`. In the R installation intended for the analysis, a starting installation is:

```r
install.packages(c("BiocManager", "argparse", "digest", "dplyr", "ggplot2", "ggrepel", "pROC", "scales"))
BiocManager::install(c("edgeR", "limma"))
```

These commands install versions available to that R/Bioconductor environment; they do not pin the historical versions. For numerical reproduction, match the recorded inventory and annotation snapshots. Consult the inventory and script requirements for analyses beyond the demo.

Figure redraws additionally use `jsonlite`, `data.table`, `patchwork`, `scales`, `png`, `ragg`, `svglite`, and `enrichplot` (Figure 4g). Some scripts use Arial and Cairo graphics; operating-system libraries and fonts may also be required. A dependency check reads the environment and does not install packages:

```sh
Rscript environment/check_packages.R demo
Rscript environment/check_packages.R figures
```

## A registered entry

`data/entries.tsv` is the single human-readable 110-entry registry. `config/differential_analysis.json` supplies the execution settings. Use `scripts/run_entry.py` to pass these settings to the original differential-analysis script.

```sh
python3 scripts/run_entry.py --id 19 --data-root /absolute/local_inputs --out outputs/ID19 --dry-run
python3 scripts/run_entry.py --id 19 --data-root /absolute/local_inputs --out outputs/ID19
```

The dry run verifies required input identities and prints the command. It does not create the output directory. A real run writes its command, effective configuration, input hashes and output hashes to `run.json`, with R output in `run.log`. Both commands require the input files to be present.

[data/input_access.tsv](../data/input_access.tsv) records endpoint/request recipes and the expected local paths for matrix and metadata retrieval. The current expected hashes are in [data/input_assets.json](../data/input_assets.json). Endpoint availability and response identity were not rechecked here; an original source-study download is not necessarily the processed input used by these analyses. See the [input guide](../data/README.md).

## Analysis modules

The 27 files in `R/` retain the original analytical implementations; resource-path and version-label corrections are recorded in provenance. Shared `*_utils.R` files are libraries. `search_canonical.r` is an adapter that reads JSON from standard input and does not expose `--help`; use `search_export.r` for command-line exports from a canonical RDA. Other analysis entry points expose `--help` when their dependencies are installed.

| Task | Main entry points |
|---|---|
| Differential analysis | `differential_analysis.r` via `scripts/run_entry.py` |
| Volcano, heatmap and PCA | `volcano_plot.r`, `heatmap.r`, `pca.r` |
| GO/KEGG ORA and GSEA | `enrichment_analysis.r`, `go_bar_plot.r`, `go_bubble_plot.r`, `kegg_bar_plot.r`, `kegg_bubble_plot.r`, `gsea.r` |
| Marker ROC, ranking and panel fitting | `roc.r`, `roc_rf.r`, `roc_svm.r`, `roc_logistic.r` |
| PPI and cross-omics | `ppi.r`, `cross_omics_trend_analysis.r`, `regulatory_network.r`, `lncrna_target_network.r`, `go_network.r`, `kegg_network.r` |
| Canonical search values | `search_export.r`, `search_canonical.r` |

For ROC/ranking commands supply `--utils_r R/roc_utils.R`, the required output arguments and explicit case/control mapping when labels are ambiguous. Preserve the manuscript's selection rules. Figure 4h names WNT2B, TPSAB1 and TPSB2 explicitly; it is not an automatically selected Top-3 panel. The SVM-ranked panel is SERPINH1, ITGA3, FCGBP, NDUFA13 and ACADVL; the RF-ranked panel is CKB, SERPINH1, FCGBP, UQCRC1 and CEMIP2. Both are fitted by logistic regression and evaluated on the same 70 samples. This statement is not a replacement for missing sample-level model-reproduction records.

Web orchestration, cache/index services, upload processing and production deployment are outside this standalone analysis distribution and are not supplied by the companion analysis archive. Original R resource fallback paths refer to the backend layout; use explicit arguments below in this reorganized distribution. Standalone R defaults may differ from webpage settings, especially lncRNA thresholds; pass manuscript settings explicitly.

## External annotations

Full matrices, canonical expression objects and complete miRNA/lncRNA/KEGG databases are not bundled. The expected external identities are in `data/input_assets.json`; source accessions and publication identifiers are in `data/entries.tsv`.

| Script or workflow | Explicit resource argument |
|---|---|
| `gsea.r`, `kegg_bar_plot.r`, `kegg_bubble_plot.r` | `--kegg_mapping_rda /absolute/kegg_mapping.rda` |
| miRNA `enrichment_analysis.r` | `--target_db /absolute/mirna_targets.sqlite` and `--mirna_resolution` from configuration |
| `regulatory_network.r`, `cross_omics_trend_analysis.r` | `--target_map /absolute/mirna_targets.sqlite` and `--mirna_resolution` from `config/mirna_resolution.json` for miRNA; for the applicable trend branch, `--lncrna_target_db /absolute/external_targets.sqlite` |
| `lncrna_target_network.r` | `--master_lnc_target_db /absolute/external_targets.sqlite` |
| mRNA/protein `kegg_network.r` | `--kegg_mapping_rda /absolute/kegg_mapping.rda` |
| Metabolite-containing `kegg_network.r` | `--kegg_joint_snapshot_dir resources/kegg_joint_snapshot`, `--metabolite_mapping_file resources/kegg_joint_snapshot/release_metabolite_mapping.tsv`, and applicable `--lncrna_target_db` / `--mirna_target_db` / `--mirna_resolution` |
| `ppi.r` | `--string_cache_dir /absolute/string_cache` as needed |

Check each script's help for the required dataset, type, threshold and output arguments. Explicit paths locate inputs; they do not download or reconstruct them. A hash identifies a file but is not sufficient to rebuild it. Exact recovery of the authors' processed target databases still requires the matching inputs or a verified construction procedure; one is not invented by this packaging change.

`resources/kegg_joint_snapshot/` contains the manifest, author metabolite mapping and org.Hs.eg.db version. `gene_symbols.tsv`, four reference tables (`gene_pathway.tsv`, `compound_pathway.tsv`, `compound_synonyms.tsv`, `pathways.tsv`) and six raw captures are external. Applicable upstream access/reuse conditions apply. To restore a lawfully held, version-matched local copy:

```sh
python3 scripts/restore_kegg_snapshot.py --source /absolute/authorized_snapshot
python3 scripts/restore_kegg_snapshot.py --check
```

The tool makes no network requests and refuses mismatched hashes. Current KEGG REST responses may differ from the historical snapshot; do not replace the expected hashes to make a new snapshot appear identical. Keep restored external tables out of public ZIPs. See `THIRD_PARTY_NOTICES.md` and the snapshot manifest for provenance and upstream links.

## Recount and redraw

Recount Figure 1d from all version-matched registered matrices; use a new output directory:

```sh
python3 scripts/recount_figure1d.py --data-root /absolute/local_inputs --out outputs/Figure1d_recount
```

Redraw selected figures from the included frozen inputs:

```sh
Rscript scripts/draw_Figure1d.R data/figures/Figure1d/counts.tsv outputs/Figure1d
Rscript scripts/draw_Figure4g.R data/figures/Figure4g/GSEA.rds outputs/Figure4g/Figure4g
Rscript scripts/draw_Figure5.R data/figures/Figure5 outputs/Figure5
Rscript scripts/draw_Figure6d.R data/figures/Figure6d/ID9_ID19_GO_full.tsv outputs/Figure6d/Figure6d
```

Figure 1d and Figure 5 take output directories; Figure 4g and Figure 6d take output filename prefixes. Figure 5 refuses an existing output directory. Run these commands once with new destinations. Original plot scripts can overwrite files at existing destinations; select new paths when preserving previous work.

| Figure | Included support and boundary |
|---|---|
| 1a–c | Registry available; supporting count/ledger tables in the companion archive; native artwork outside this repository |
| 1d | Counts and counting rule, recount and standalone redraw; full identifier union in full archive |
| 1e | Frozen and recovered source values in the companion archive |
| 2–3 | Framework/interface schematics; Figure 2 counts have registry support, while Figure 3 thumbnail source mappings remain incomplete |
| 4a–j | Core analysis code; archived summaries and differential tables in full archive. Exact Figure 4b display matrix and complete 4i–j model-reproduction records are not supplied |
| 4g | Corrected GSEA object and five-pathway redraw supplied |
| 5 | Frozen HSP90AA1 Search values, ROC coordinates, summaries and redraw supplied |
| 6a–c | Core code; trend/network source records in full archive |
| 6d | Full GO export and all-131-term redraw supplied |

Redrawing frozen values and recomputing statistical analyses are different operations. Original publication artwork need not be pixel-identical across operating systems or fonts. This package does not claim complete end-to-end reconstruction of every panel.

## Companion archive

The companion record is [Zenodo DOI 10.5281/zenodo.22918243](https://doi.org/10.5281/zenodo.22918243), analytical version 1.1.0. Paths refer to the extracted `EVOmicsDB_v1.1.0/` archive root:

| Archive location | Contents |
|---|---|
| `01_figures/FIGURE_MAP.tsv` | Complete panel-to-evidence map and remaining gaps |
| `02_datasets/ENTRY_INDEX.tsv` | Entry-to-result and execution-record index |
| `03_results/` | Complete differential, GSEA, joint GO, metabolite KEGG and miRNA ORA results |
| `05_methods/` | Reproduction/validation notes and consolidated execution records |

These archive paths do not imply the same directory structure in this GitHub repository. In the local registry, `zenodo_differential_results` uses this archive-relative convention. Complete negative and empty results remain in the full archive.

## Verification and provenance

```sh
python3 scripts/verify_archive.py
python3 tests/validate_package.py
```

[provenance.json](provenance.json) records the supplied GitHub archive and the sources of corrected content. Its retained-file mappings distinguish unchanged copies from adapted files; `SHA256SUMS.txt` records all current release files. `docs/VERIFICATION.md` records what was actually checked for this layout. Historical analysis checks remain in the paired full archive. Neither syntax checks nor the artificial demo establish full clinical-result replication, a clean-environment installation or a website deployment.

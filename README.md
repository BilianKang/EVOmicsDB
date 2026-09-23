# EVOmicsDB

Analytical code for **EVOmicsDB: a clinically annotated computational framework for pan-cancer extracellular vesicle multi-omics analysis**.

**Analytical version:** 1.1.0 · **Registry:** 110 matrix entries · **Analysis revision:** 22 September 2026 · **Figure 1d recount:** 23 September 2026.

This repository contains the R analyses, entry-specific configurations, a small artificial example and selected figure inputs. A matrix entry is not an independent study or a unique donor. RF/SVM outputs rank exploratory candidates; panel ROC estimates describe apparent within-entry performance.

## Start here

Use Python 3.10+ and R 4.5.2 with the required packages. Run from this directory:

```sh
python3 scripts/verify_archive.py
Rscript environment/check_packages.R demo
python3 scripts/run_demo.py --out outputs/demo
```

The example runs differential analysis, a volcano plot and ROC analysis on **artificial log2 values**. Use a new output directory for each run. Its outputs have no biological interpretation. See [reproduction instructions](docs/REPRODUCING.md) for installation, real inputs and figure commands.

## Repository guide

| Folder | Contents |
|---|---|
| `R/` | 27 original analysis scripts and shared helpers, kept together |
| `scripts/` | Entry/demo runners, figure redraws and integrity tools |
| `config/` | Differential-analysis and miRNA-resolution settings |
| `data/` | 110-entry registry, input identities, grouping labels, artificial example and selected figure inputs |
| `resources/` | KEGG snapshot manifest and retained mapping tables; external reference tables are required |
| `environment/` | Recorded R versions and package checks |
| `tests/` | Configuration, syntax and bounded regression checks |
| `docs/` | Reproduction guide, file provenance and current verification scope |

For an input-ready registered entry:

```sh
python3 scripts/run_entry.py --id 19 --data-root /absolute/local_inputs --out outputs/ID19 --dry-run
python3 scripts/run_entry.py --id 19 --data-root /absolute/local_inputs --out outputs/ID19
```

The data root must contain the matching `db/...` files listed in `data/input_assets.json`. The runner validates input hashes before analysis. Full expression matrices and complete third-party target databases are external. Some workflows require explicit annotation paths; see the guide.

## Manuscript figures and full archive

Small frozen inputs are included to redraw Figure 1d, Figure 4g, Figure 5 and Figure 6d. Redrawing these inputs does not recompute their underlying statistical analyses. Complete reconstruction of all Figure 4 panels is not claimed; the original Figure 4b display matrix and the complete Figure 4i–j model-reproduction records remain unavailable in this package.

The paired full analysis archive is `EVOmicsDB_v1.1.0_20260923.zip`; it contains the result tree, supplementary tables, additional figure records, historical validation and backend integration source. The current organized archive is identified by `current_paired_archive` in [provenance](docs/provenance.json); `source_archives` preserves the earlier input archive identities. The organized archive uses `code/`, `data/`, `results/`, `figures/` and `docs/`, with original Python integration under `code/reference/`. This GitHub package does not deploy the production website.

## Version and citation

Project: <https://github.com/BilianKang/EVOmicsDB>. Website: <https://evomicsdb.com>.

This folder is a locally prepared distribution of analytical version 1.1.0. The new version-specific Zenodo DOI has not been supplied in these files. Add the published version DOI to `CITATION.cff` and the manuscript before claiming public availability. The historical DOI `10.5281/zenodo.19534750` identifies the older 102-entry release and must not be cited as this revision.

The authors' code is covered by [MIT](LICENSE). Experimental data and third-party annotations retain their source terms; see [third-party notices](THIRD_PARTY_NOTICES.md).

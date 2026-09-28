# EVOmicsDB

Analytical code for **EVOmicsDB: a clinically annotated computational framework for pan-cancer extracellular vesicle multi-omics analysis**.

**Software version: 1.1.0 · Registry: 110 matrix entries · [Zenodo record](https://doi.org/10.5281/zenodo.22918243)**

This repository contains the R analysis modules, entry-specific settings, launchers, artificial example data and selected figure redraw inputs. A matrix entry is not an independent study or a unique donor. RF/SVM outputs rank exploratory candidates; retained panel ROC estimates describe apparent within-entry performance.

## Quick start

Use Python 3.10+ and R with the required packages. From the repository root:

```sh
python3 scripts/verify_archive.py
python3 tests/validate_package.py
Rscript environment/check_packages.R demo
python3 scripts/run_demo.py --out outputs/demo
```

The demo runs differential analysis, a volcano plot and marker ROC on artificial log2 values. It requires a new output directory and has no biological interpretation. For package installation and registered datasets, start with [reproduction instructions](docs/REPRODUCING.md).

## Repository guide

| Folder | Role |
|---|---|
| `R/` | 27 statistical modules and their shared helpers |
| `scripts/` | Analysis launchers, four figure redraws, recounting and integrity/resource tools |
| `config/` | Differential-analysis and miRNA-resolution settings |
| `data/` | [Entry registry and input guide](data/README.md), artificial example and selected figure inputs |
| `resources/` | KEGG snapshot identities and author metabolite mappings. The gene-symbol table is external |
| `environment/` | Package inventory and dependency checks |
| `tests/` | Configuration, syntax and focused helper checks |
| `docs/` | Reproduction instructions, validation scope and file provenance |

## Run a registered entry

```sh
python3 scripts/run_entry.py --id 19 --data-root /absolute/local_inputs --out outputs/ID19 --dry-run
python3 scripts/run_entry.py --id 19 --data-root /absolute/local_inputs --out outputs/ID19
```

Both commands require the matching external inputs. The runner validates file hashes before execution. See the [input guide](data/README.md).

## Selected figures and complete evidence

The [figure input guide](data/figures/README.md) maps Figure 1d, Figure 4g, Figure 5 and Figure 6d to their input files and scripts. Saved-value redraws do not recompute the underlying analyses. Exact reconstruction of all manuscript panels is not claimed.

Complete statistical outputs, additional panel source tables and consolidated execution records are in the companion [Zenodo record](https://doi.org/10.5281/zenodo.22918243). Resolve the registry's `zenodo_differential_results` paths from the extracted archive root, `EVOmicsDB_v1.1.0/`. The archive guide in [REPRODUCING.md](docs/REPRODUCING.md#companion-archive) explains the other directories. Native artwork and the production website application are outside this code repository.

## Citation and reuse

Use [CITATION.cff](CITATION.cff) and analytical version **1.1.0**. Author software and documentation use [MIT](LICENSE); source data and annotations retain the [applicable third-party terms](THIRD_PARTY_NOTICES.md).

[Project repository](https://github.com/BilianKang/EVOmicsDB) · [Database website](https://evomicsdb.com) · [Checks and limitations](docs/VERIFICATION.md)

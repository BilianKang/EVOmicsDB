# Data used by this code repository

| Path | Purpose | Used by |
|---|---|---|
| `entries.tsv` | 110-entry registry, study accessions and analysis settings | Entry runner, configuration checks and Figure 1d recount |
| `input_assets.json` | Expected identities of external matrices, metadata and annotation inputs | Entry runner and recount |
| `input_access.tsv` | Recorded website endpoint/request patterns, paired with current expected local input paths | Manual retrieval; verify against the input catalogue |
| `example/` | Artificial expression values and group labels | `scripts/run_demo.py` |
| `figures/` | Frozen inputs for four supplied redraw scripts | [Figure input guide](figures/README.md) |

Paths in `input_assets.json` and `input_access.tsv` resolve below the runner's `--data-root`. The endpoint recipes were retained from the supplied package; they were not contacted or checked for current availability during this revision. HTTP success alone does not establish a match to the expected processed matrix. Compare retrieved bytes with the hashes in `input_assets.json`. Earlier HTTP status records have not been presented as current evidence.

For example, the recorded matrix route uses POST to `https://evomicsdb.com/api/analysis/browse/download_datafile` with a JSON body such as `{"dataset_id":"19"}`. Its expected local filename is recorded in `input_access.tsv`. No automatic data download is performed by the launchers.

Full expression matrices and full clinical workbooks are external. Group labels come from those external inputs.

The `zenodo_differential_results` column points inside the companion archive, not this repository. The registry contains 7,120 matrix-level sample records; assays can share cohorts, so this is not a unique-donor count.

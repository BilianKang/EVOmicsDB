# Selected figure inputs

Drawing entry points remain in `scripts/`; their input files remain here. Run the commands below from the repository root, using fresh output destinations.

| Figure | Input directory | Drawing entry point | Scope |
|---|---|---|---|
| 1d | `data/figures/Figure1d/` | `scripts/draw_Figure1d.R` | Six molecular-category counts and the counting rule |
| 4g | `data/figures/Figure4g/` | `scripts/draw_Figure4g.R` | Five selected pathways from the saved ID19 GSEA object |
| 5 | `data/figures/Figure5/` | `scripts/draw_Figure5.R` | HSP90AA1 values, entry-level labels, ROC coordinates and summaries |
| 6d | `data/figures/Figure6d/` | `scripts/draw_Figure6d.R` | Full GO input table and its 131 displayed terms |

```sh
Rscript scripts/draw_Figure1d.R data/figures/Figure1d/counts.tsv outputs/Figure1d
Rscript scripts/draw_Figure4g.R data/figures/Figure4g/GSEA.rds outputs/Figure4g/Figure4g
Rscript scripts/draw_Figure5.R data/figures/Figure5 outputs/Figure5
Rscript scripts/draw_Figure6d.R data/figures/Figure6d/ID9_ID19_GO_full.tsv outputs/Figure6d/Figure6d
```

Figures 1d and 5 take output directories; 4g and 6d take filename prefixes. Figure 5 rejects an existing output directory. These commands redraw saved values and require the R packages/fonts described in the [reproduction guide](../../docs/REPRODUCING.md).

The complete panel-evidence map and other source tables are in the companion archive's `01_figures/`. Figure 2's framework counts can be checked against the registry. Figure 3's illustrative statistical thumbnails lack an exact input-to-panel mapping in the supplied material. The missing Figure 4b display input and complete Figure 4i–j model-reproduction records are not supplied here. These limits are not resolved by the presence of generic analysis modules.

# Verification and scope

This document describes checks of the corrected GitHub distribution. It does not claim full real-data reruns or complete manuscript reconstruction.

| Check | Scope |
|---|---|
| Source preservation | Analysis modules, four figure entry points and metabolite mappings retained. The gene-symbol table and minimal grouping exports are external |
| Scientific consistency | Registry and configuration use the final 110-entry specification; Figure 1d uses the corrected six-category counts |
| Input contracts | 237 external-input identity records; all configured input paths and 220 retrieval-recipe paths have identities |
| Group labels | Case and control labels remain in the external inputs. This repository does not include separate grouping exports |
| Syntax and helper regression | Full original validation entry point checks Python/R syntax, three R helpers, external KEGG restoration, registry/configuration, retained Figure 1d counts and provenance |
| Integrity | SHA256SUMS.txt covers every release file except itself; the verification tool checks listed file identities |
| Artificial example | Differential analysis, volcano plot and marker ROC complete using data/example/ inputs |
| Repository navigation | Local documentation links, script paths and selected figure inputs resolve |

The bounded validation command and artificial example passed in the available local environment. Demo dependencies were present. The figure dependency check reported svglite missing, so complete execution of all four figure exporters was not assessed. Recorded R/package versions remain references, not a tested clean-environment lock.

## Repeat the checks

From the repository root:

```sh
python3 scripts/verify_archive.py
python3 tests/validate_package.py
Rscript environment/check_packages.R demo
python3 scripts/run_demo.py --out outputs/new_demo
```

The package validator accepts `--rscript /absolute/path/to/Rscript`. It does not download data or install packages. Use fresh output directories.

Full expression matrices, complete clinical workbooks and several annotations remain external. Recorded endpoint recipes do not guarantee availability or byte identity of current server responses.

Statistical outputs retained for the figures describe the original analyses. Apparent within-entry ROC estimates are not independent validation. Missing exact Figure 4 inputs/model records and other panel reconstruction gaps remain documented in the reproduction guide and companion archive. The production website is outside this repository's scope.

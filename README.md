# CRC Manuscript Code Repository

This working copy combines the original manuscript code with a framework for future processing and figure-generation workflows.

All 15 files tracked in the original repository at commit `08b22acb28622d4470b5567101c9f2f4d3187aa6` are preserved unchanged under [additional_analyses/](additional_analyses/). This includes all 13 notebooks, the ZL R script, and the [original README](additional_analyses/README.md).

Original code by Hongcheng Yao, Alex Tang, Yiwen He, and Zhengyue Li. The original figure descriptions and contributor information remain in the original README.

The first plotting notebook, [plot/scrna_global_umap.ipynb](plot/scrna_global_umap.ipynb), is an unchanged copy of [HY_fig1c_globalUmap.ipynb](additional_analyses/HY_fig1c_globalUmap.ipynb) with a filename based on its analysis topic. Its source cells, saved outputs, and metadata are preserved; processing steps have not yet been separated from plotting, and execution from the new location has not been validated.

The `process/`, `R/`, `envs/`, and `config/` locations remain scaffolding only. No analysis code has been extracted from the original notebooks, no datasets have been added, and no analysis has been executed as part of this reorganization.

## Directory layout

```text
crcmetmanuscript2026/
|-- README.md
|-- LICENSE                         # Empty placeholder
|-- CITATION.cff                    # Empty placeholder
|-- .gitignore
|-- additional_analyses/
|   |-- README.md                   # Original GitHub README, unchanged
|   |-- HY_*.ipynb                  # All seven original HY notebooks
|   |-- AT-*.ipynb                  # All five original AT notebooks
|   |-- YH-*.ipynb                  # Original YH notebook
|   `-- ZL-*.R                     # Original ZL script
|-- envs/                           # Reproducible environment specifications
|-- config/                         # Data paths and analysis configuration
|-- process/
|   |-- README.md
|   |-- scrna/                      # Single-cell RNA analysis
|   |-- tcr/                        # T-cell receptor analysis
|   `-- xenium/                     # Spatial transcriptomics analysis
|-- plot/
|   |-- README.md                   # Plotting notebook guidance and status
|   `-- scrna_global_umap.ipynb     # Unchanged copy; execution not yet validated
|-- R/                             # Shared R functions
|-- data/
|   `-- README.md                   # Data access and provenance guidance
|-- docs/
|   |-- architecture_guide.md       # Organizational recommendations
|   |-- additional_analyses_guide.md
|   |-- migration_notes.md         # Checks and remaining path dependencies
|   |-- migration_manifest.tsv     # Original paths, new paths, and checksums
|   |-- migration_provenance.json  # Source commit and migration metadata
|   `-- figure_map.tsv             # Header only; dependencies not yet mapped
`-- results/                       # Generated tables and figures
```

Empty directories contain `.gitkeep` files so that Git can retain their structure. These files have no computational role. Git metadata remains in `.git/` at the repository root.

## Start here

Read [the migration notes](docs/migration_notes.md) for the completed checks and outstanding execution dependencies, and [the architecture guide](docs/architecture_guide.md) before adding new workflows. Establish the relationship between each figure input and its actual producing script, then record it in [the figure map](docs/figure_map.tsv).

Use `process/` for upstream computation, `plot/` for figure-oriented entry points, and `R/` for shared functions. Keep all documentation, comments, and explanatory notebook text in English.

The [preservation guide](docs/additional_analyses_guide.md) explains the role of `additional_analyses/`. The directory contains the entire original repository snapshot, including the HY notebooks. Keep this preserved baseline identifiable as new workflows are developed.

`LICENSE` and `CITATION.cff` are intentionally empty. They do not yet specify a license or valid citation metadata and must be completed before a public release. The empty `envs/` and `config/` directories likewise do not define a runnable environment.

# CRC Manuscript Code Repository

Code for the analyses and figures in the CRC manuscript.

Code by Hongcheng Yao, Alex Tang, Yiwen He, and Zhengyue Li.

## Find the code

- [plot/](plot/): plotting notebooks, starting with [global UMAP](plot/scrna_global_umap.ipynb).
- [process/](process/): upstream processing workflows; no runnable scripts are available yet.
- [additional_analyses/](additional_analyses/): original analysis notebooks, scripts, and figure descriptions.

## Getting started

1. Find your notebook's figure panels and required inputs in the [figure map](docs/figure_map.tsv).
2. Check the [data manifest](docs/data_manifest.tsv) for data versions, download links, and checksums.
3. Follow the [data access guide](docs/data_access.md) to obtain and verify inputs and review software requirements before running the notebook.

Use `data/` for inputs and `results/` for generated outputs, relative to the repository root.

## Current availability

Data are excluded from GitHub. The Zenodo data release is pending, so public download links are not yet available.

The input mapping currently covers the global UMAP notebook only. Its analysis has not been rerun after the input path update, and a complete environment specification is pending.

License and citation information are pending.

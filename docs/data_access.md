# Data Access

## Availability

Data are stored under the local `data/` directory, which is excluded from Git. The two files below have been verified locally; a fresh clone contains no study data. Zenodo publication is pending, so no public download links are available yet.

After publication, use the specific version DOI and download links recorded in the [data manifest](data_manifest.tsv).

## Current files

| Path relative to the repository root | Use |
| --- | --- |
| `data/scrna/all_seurat_obj.RDS` | Required input to the [global UMAP notebook](../plot/scrna_global_umap.ipynb). |
| `data/scrna/metadata/18_144_metadata.txt.gz` | Fixed collaborator-supplied metadata for upstream processing; not required by the global UMAP notebook. |

The current inventory covers this notebook and its upstream metadata only. File sizes, checksums, versions, and source information are recorded in the data manifest.

## Prepare and check inputs

1. Find your notebook's `input_files` in the [figure map](figure_map.tsv), then match each path to `relative_path` in the [data manifest](data_manifest.tsv).
2. For files you already have, create the parent directories and place each file at its listed path, relative to the repository root. Preserve filename case.

From the repository root, check the global UMAP input:

```bash
wc -c < data/scrna/all_seurat_obj.RDS
sha256sum data/scrna/all_seurat_obj.RDS
```

Compare the results with `bytes` and `sha256` in the matching manifest row; both must match. To check metadata, substitute its path. On macOS, use `shasum -a 256` instead of `sha256sum`. These checks verify file integrity only.

## Global UMAP software

Open the notebook from the repository root or `plot/`. It uses `here` to locate its input; the inspected version is 1.0.1. It also requires Seurat, singlecellmethods, uwot, and the other packages loaded in the notebook.

Historical notebook metadata records R 4.1.2. Available evidence supports `immunogenomics/singlecellmethods` at commit `58ca323b336a24c8b13b8dfbbbd3516fd89bc845`, using `weighted_pca(do_corr=TRUE)`. This package version is a historical inference, not an execution-time environment lock.

## Reproduction scope

Start from the saved Seurat object. The original count matrices are unavailable, and the historical object-generation notebook has not been imported into `process/`.

The notebook retains its original analysis cells and historical outputs and has not been rerun after the input path update. A complete environment specification and end-to-end validation remain pending.

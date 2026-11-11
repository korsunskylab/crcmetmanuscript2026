# Local Repository Migration Notes

## Source and scope

- Source repository: https://github.com/korsunskylab/crcmetmanuscript2026
- Source branch: `main`
- Source commit: `08b22acb28622d4470b5567101c9f2f4d3187aa6`
- Full history at clone: 14 commits; the clone is not shallow.
- Local working branch: `local/restructure-code`
- Framework source: the neighboring `ccr_ver/` directory.

All 15 original tracked files were moved under `additional_analyses/` with their original names and contents. The original README remains `additional_analyses/README.md`. The framework's preservation guide was placed at `docs/additional_analyses_guide.md` to avoid that filename collision.

The new root-level framework contains documentation, placeholders, and a header-only figure map. It does not yet contain a newly implemented plotting or processing workflow. No notebooks or R analysis scripts were executed, and no datasets were downloaded. No changes were pushed to GitHub.

After the initial relocation, [plot/scrna_global_umap.ipynb](../plot/scrna_global_umap.ipynb) was copied from [additional_analyses/HY_fig1c_globalUmap.ipynb](../additional_analyses/HY_fig1c_globalUmap.ipynb). The two files have identical contents, including source cells, saved outputs, and metadata. The new filename describes the analysis independently of figure numbering. Processing has not been separated from plotting, and execution from the new location has not been validated. The migration manifest continues to describe the 15 preserved original files.

## File-integrity checks

The [migration manifest](migration_manifest.tsv) records all original paths, new paths, source Git blobs, sizes, and SHA-256 checksums. Original checked-out files were verified against the source commit before relocation, then checked again after relocation. The neighboring framework source was copied without modification.

The [provenance record](migration_provenance.json) records the source version and initial migration scope. The existing Git history and root `.git/` directory were retained. The local branch for this reorganization is `local/restructure-code`.

## Static path and dependency observations

These observations were obtained by reading notebook source cells and the R script. They are not execution results or a complete dependency audit. Existing relative relationships among the original tracked files have been retained, but the required data files are not included in this clone.

| Original analysis | Observed dependency or path assumption | Follow-up before execution |
| --- | --- | --- |
| `HY_fig1c_globalUmap.ipynb` | Reads `all_seurat_obj.RDS` relative to the working directory. | Identify and provide the correct input object and document the working directory. |
| HY clonotype and diversity notebooks | Read files under `tcr/`; the diversity notebook uses both `obj.RDS` and `obj.rds`. | Preserve the distinct input identities and filename case, and trace their producing workflows. |
| `HY_fig5d_figs7ab_xeUmap.ipynb` | Reads `xenium/cell_meta.rds` and separate `obj.RDS` files under `all/`, `myeloid/`, and `t/`. | Retain these distinct object locations when documenting data layout. |
| Other HY spatial notebooks | Read metadata and analysis objects from `.`. | Document the actual data locations and intended working directory for each notebook. |
| `HY_fig6acde_figs9abcedfghi.ipynb` | Source cell 3 loads and documents `singlecellmethods` from a personal absolute filesystem path. | Provide a usable, versioned dependency before execution or refactoring. |
| AT notebooks | Include `gs://` input locations and local CSV or H5AD inputs. | Establish data availability, access requirements, and expected local inputs. |
| `YH-2025-Treg-analysis_SuppFig6A-D.ipynb` | Source cell 6 reads `crc_analysis_v2.h5ad` and `foxp3cd8.h5ad`. | Obtain the exact input versions and define their location. |
| `ZL-2025-Treg-analysis_SuppFig6E-N.R` | Lines 5, 279, and 387 set the working directory to a fixed `/projects/...` location. | Review the original working-directory and input assumptions before a portability update. |

For an unchanged notebook that originally resolved paths from the original repository root, `additional_analyses/` now occupies that original relative position. This does not supply missing data or guarantee that the notebook will run. Avoid treating the new repository root and `additional_analyses/` as interchangeable working directories.

No path strings were edited inside the preserved code. Future functional changes should be tracked separately from the relocation. The new `docs/figure_map.tsv` remains header-only until actual input-to-producer relationships have been established.

## Review the local changes

After committing the reorganization, inspect the working tree and most recent local commit with:

```bash
git status --short --branch
git show --stat --oneline HEAD
git show --summary --find-renames HEAD
```

The empty `LICENSE` and `CITATION.cff` files are framework placeholders. They have not been populated with a license choice or publication metadata.

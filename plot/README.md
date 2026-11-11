# Figure Generation

This directory contains notebooks or scripts named by analysis topic. Record current manuscript figure and panel assignments in `../docs/figure_map.tsv` so that filenames can remain stable when figure numbering changes.

All original GitHub notebooks, including the HY notebooks, are preserved under [additional_analyses/](../additional_analyses/).

[scrna_global_umap.ipynb](scrna_global_umap.ipynb) is a byte-for-byte copy of [HY_fig1c_globalUmap.ipynb](../additional_analyses/HY_fig1c_globalUmap.ipynb). Only its filename and location differ. It retains the original QC, normalization, PCA, UMAP, and plotting code, together with saved outputs and metadata. Upstream processing has not yet been separated, input paths have not been adjusted, and the copied notebook has not been executed or validated from this location.

Each figure entry point should document its input files, upstream producers, environment, configuration, expected outputs, and execution instructions. Use a fresh session and a documented working directory when validating a notebook.

Prefer stable upstream results for shared or expensive computations. Figure-specific filtering, formatting, and short transparent calculations may remain in the notebook when their role is documented. Keep reusable functions in `../R/`.

The figure map in `../docs/figure_map.tsv` contains a header only; populate it as the actual input-to-producer relationships are verified.

See [the architecture guide](../docs/architecture_guide.md) for output handling and reproduction guidance.

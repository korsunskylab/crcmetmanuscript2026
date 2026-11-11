# Processing and Analysis

This directory is reserved for code that produces analysis results and figure inputs. It includes both preprocessing and statistical analysis.

The original GitHub code is preserved under [additional_analyses/](../additional_analyses/). No upstream script has been extracted from that code or added to this directory during the reorganization.

Organize scripts by data modality:

- `scrna/`: quality control, normalization, integration, annotation, and dimensionality reduction.
- `tcr/`: clonotype preparation, sampling, diversity summaries, and related analysis.
- `xenium/`: spatial data preparation, coordinate handling, distances, scores, and statistical models.

These subdirectories currently contain placeholders only. Add scripts after tracing the inputs and outputs used for the manuscript. Document prerequisites and execution order; numeric filename prefixes may help but do not replace dependency documentation.

Keep shared R functions in `../R/`, configuration in `../config/`, and environment specifications in `../envs/`. Record the actual producing script for every figure input in `../docs/figure_map.tsv`.

See [the architecture guide](../docs/architecture_guide.md) for the intended processing and plotting boundary.

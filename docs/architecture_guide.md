# Repository Architecture and Completion Guide

## Scope and current status

This document records recommendations for organizing a manuscript code repository. All original versioned files from GitHub commit `08b22acb28622d4470b5567101c9f2f4d3187aa6` have been preserved unchanged under `additional_analyses/`, including every HY notebook and the original README. The first notebook in `plot/`, [scrna_global_umap.ipynb](../plot/scrna_global_umap.ipynb), is an unchanged copy of the original global UMAP notebook. The processing, shared-function, environment, and configuration locations remain scaffolding. Upstream dependencies have not yet been fully mapped, and no analysis workflow has been executed or validated as part of this reorganization.

Keep all repository documentation, explanatory notebook text, code comments, and user-facing messages in English. Preserve scientific identifiers and names where needed.

The proposed directory names are organizational choices, not a directory standard imposed by Nature journals. Journal guidance emphasizes accessible code, identifiable versions, installation and execution documentation, and reproducibility. Requirements for a target journal should be checked when preparing the submission.

## Directory responsibilities

| Location | Intended responsibility |
| --- | --- |
| `additional_analyses/` | All original versioned repository files, including HY, AT, YH, ZL, and the original README, preserved unchanged. |
| `process/scrna/` | Single-cell RNA preprocessing and analysis that produce reusable results or figure inputs. |
| `process/tcr/` | T-cell receptor preprocessing, clonotype summaries, sampling, and statistical analysis. |
| `process/xenium/` | Spatial preprocessing, coordinate handling, annotations, distances, expression summaries, and models. |
| `plot/` | Figure-oriented entry points that consume documented inputs and produce manuscript panels. |
| `R/` | Shared R functions used by processing scripts or plotting notebooks. |
| `envs/` | Environment specifications and software version information corresponding to the manuscript analysis. |
| `config/` | Portable data-path configuration and explicit scientific parameters. |
| `data/` | Data-access documentation and, when intentionally added, documented local data locations. |
| `docs/` | Figure dependencies, workflow explanations, and reproduction instructions. |
| `results/` | Generated result tables, figures, and other run outputs. |

Organize upstream work by modality because a processing result may support several figures. Name plotting entry points by analysis topic and record current figure and panel assignments in `figure_map.tsv` so readers can locate a manuscript result directly even when figure numbering changes. Put shared functions in one maintained location instead of copying them across notebooks.

## Preserving the original repository contents

`additional_analyses/` holds the complete original set of 15 tracked files: 13 notebooks, one R script, and the original README. This directory includes all seven HY notebooks. Its name does not imply that its contents are obsolete.

Original filenames, file contents, attribution, and figure descriptions were preserved. The source commit and per-file checksums are recorded in [the migration manifest](migration_manifest.tsv) and [provenance record](migration_provenance.json). Git metadata remains at the repository root.

The original README is retained at `additional_analyses/README.md`. The framework's separate preservation guidance is now [additional_analyses_guide.md](additional_analyses_guide.md), so it does not replace the original README. The repository-level README links to both the original material and the new framework.

Moving files can affect relative paths and working-directory assumptions. Static observations are documented in [the migration notes](migration_notes.md); no analysis was executed. Document any future path or analytical changes separately from this initial relocation.

When synchronizing framework updates into this checkout, explicitly preserve `additional_analyses/` and `.git/`. Directory names and `.gitignore` rules are not automatic protection against deletion by external synchronization tools. Inspect planned additions, replacements, and deletions before applying a synchronization operation.

The `plot/` directory currently contains its README and `scrna_global_umap.ipynb`, copied from `additional_analyses/HY_fig1c_globalUmap.ipynb` without content changes. Processing remains embedded in that notebook, and its execution from the new location has not been validated. The `process/` directory contains guidance and placeholders only; upstream scripts have not been inferred or generated. Empty license and citation placeholders are not publication metadata.

## The processing and plotting boundary

`process/` includes both preprocessing and analysis. It is not limited to file conversion or quality control.

Move shared, expensive, or reusable computations upstream when an explicit saved result provides a useful boundary. Examples include dimensionality reduction, repeated sampling, nearest-neighbor distances, reusable expression scores, and model fitting. These are examples of responsibilities, not scripts already supplied by this scaffold.

`plot/` should provide an understandable route from defined figure inputs to manuscript panels. Figure-specific formatting, filtering, and brief transparent calculations may remain in a notebook. Do not split every short transformation into a separate script merely to make the plotting directory appear computationally empty.

When separating an existing notebook, preserve its analytical behavior and compare its outputs with the established results. Record the exact input and configuration version used on each side of the boundary. Do not silently substitute a different processed dataset while reorganizing paths.

## Establish dependencies before moving code

For each manuscript panel, trace the following chain:

```text
Documented starting data
    -> Actual upstream scripts and configuration
    -> Named, versioned analysis results or figure inputs
    -> Figure notebook or plotting script
    -> Expected tables and figure panels
```

Populate `figure_map.tsv` only after inspecting the actual files and scripts. Do not infer the producing script solely from a matching filename. Record all prerequisites where a figure has multiple inputs or a processing step has multiple dependencies.

The figure map is initially header-only. Its columns are:

| Column | Information to record later |
| --- | --- |
| `figure_panel` | Manuscript figure or panel identifiers. |
| `plot_notebook` | Repository-relative path to the figure entry point. |
| `input_files` | Required input files, distinguished by their full relative paths. |
| `producer_scripts` | Actual scripts that produce those inputs. |
| `prerequisites` | Earlier steps, required annotations, and other execution dependencies. |
| `output_files` | Expected generated tables or figure files. |
| `data_reference` | Accession, DOI, or a link to documented data access. |
| `environment_reference` | The relevant environment specification or version record. |
| `config_reference` | The configuration and parameters used for the result. |
| `notes` | Additional context, including manual steps, resource requirements, or limitations. |

Use a consistent delimiter, such as a semicolon, when a field contains multiple paths. More detailed relationships can be documented separately if a single row becomes difficult to read.

## Paths, parameters, and shared dependencies

Define a consistent project root and document the expected working directory for each entry point. Put machine-specific paths in local configuration rather than embedding personal filesystem paths in analysis code.

Record scientific settings explicitly, including relevant thresholds, sample exclusions, identifier mappings, random seeds, and thread settings. Keep these settings separate from machine-specific locations so that changing a data directory does not change the analysis.

For custom packages, provide a usable installation source and a fixed version or commit, or include the required reusable functions in the repository. An instruction that loads an unpublished package from a personal development directory is not a portable dependency specification.

If an upstream workflow uses a scheduler, include the actual computational scripts and their arguments as well as the submission wrappers. Document scheduler-specific settings separately from scientific parameters.

If a workflow includes manual annotations or region selection, retain the corresponding mapping tables or region definitions and document how they were produced. The presence of a script alone does not capture those decisions.

## Two reproduction entry points

Plan two documented starting points where the available data support them:

1. **Figure reproduction:** obtain fixed, versioned figure inputs and regenerate the published panels.
2. **Full analysis reproduction:** begin from a clearly specified earlier data stage, run the upstream processing and analysis, and then regenerate the panels.

These entry points are recommendations for later implementation. No runner is implemented in the scaffold. State the actual starting point, expected outputs, runtime, and resource requirements for each route. If only processed inputs are available, describe that scope explicitly rather than claiming reproduction from raw data.

Small representative examples can help users verify installation. They should be distinguished from the full analysis and should have documented expected results.

## Environments and validation

Record the language, package, and external-tool versions used to produce the manuscript results. Add appropriate environment or dependency files under `envs/` once those versions have been established. Do not fill this directory with guessed versions or assume that the latest release reproduces the original analysis.

Validate notebooks from a fresh session in the documented order. Check that required objects are created or loaded explicitly and that prior interactive state is unnecessary. Verify relevant sample counts, intermediate summaries, statistical results, and figure outputs rather than relying only on successful execution.

For stochastic or platform-dependent steps, document seeds and settings and describe the expected level of numerical or visual agreement. A fixed figure-input archive can preserve the exact inputs used for the publication while the upstream workflow documents how those inputs were generated.

A full software package, a workflow engine, containers, and continuous integration are optional implementation choices. Add them when they solve a demonstrated dependency, execution, or maintenance problem; the scaffold does not require them.

## Generated files and notebook outputs

Keep generated result tables and figures under `results/` or another explicitly documented output location. The scaffold excludes generated results from Git while retaining the directory placeholder.

Keep notebooks readable and easy to review. Selected lightweight outputs may be useful; large embedded images can instead be distributed as separate rendered notebooks, HTML/PDF reports, or release artifacts. Document which outputs were generated by the publication version of the code.

The default data ignore rules can be narrowed later for deliberately included small examples. Large datasets and intermediate objects should have documented external storage, versions, and retrieval instructions.

## Publication metadata and versioning

Before a public release, replace the empty `LICENSE` placeholder with the selected license text and complete `CITATION.cff` with valid project and citation metadata. Empty placeholders are not a license declaration or valid citation information.

After validating the manuscript workflow, identify the corresponding code version with a release or tag and archive that version in a service that provides a persistent identifier, such as Zenodo. Link the repository and archived version in the manuscript's Code Availability statement. Maintain access to that version when later development continues.

Use version control to preserve historical work. Keep the primary reproduction path focused on the scripts and settings that actually support the manuscript results, and clearly distinguish any retained exploratory material.

## Suggested completion sequence

1. Inventory the final figure entry points and their actual inputs.
2. Trace each input back to its generating scripts, configuration, and starting data.
3. Populate the figure map and resolve missing or ambiguous dependencies.
4. Place upstream scripts by modality and centralize shared functions without changing the analysis.
5. Add portable configuration, environment specifications, and data-access documentation.
6. Validate a representative workflow, then validate the complete supported reproduction route.
7. Complete publication metadata and archive the validated manuscript version.

## References and examples

- [Nature: Guidelines for authors submitting code and software](https://media.nature.com/full/nature-cms/documents/GuidelinesCodePublication.pdf).
- [Nature: Reporting standards and availability of data, materials, code and protocols](https://www.nature.com/nature/editorial-policies/reporting-standards).
- [SpaceBar manuscript repository](https://github.com/grantkinsler/SpatialBarcodes): separate extraction, plotting, helper, and intermediate-data locations.
- [scMultiSim manuscript repository](https://github.com/ZhangLabGT/scMultiSim_manuscript): data generation and benchmarking code with dependency and resource information.
- [AUTOENCODIX reproducibility repository](https://github.com/jan-forest/autoencodix-reproducibility): configuration, environments, and execution entry points for reproducing results.

These repositories illustrate possible organizational approaches. They are not a substitute for the target journal's current requirements, and their organization should not be copied without considering this project's actual dependencies.

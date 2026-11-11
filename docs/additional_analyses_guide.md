# Guide to the Preserved Original Analyses

[additional_analyses/](../additional_analyses/) contains the complete set of versioned files from the original GitHub repository at commit `08b22acb28622d4470b5567101c9f2f4d3187aa6`. This includes HY code as well as code from the other contributors.

## Preserved contents

| Original content | File count |
| --- | ---: |
| HY notebooks | 7 |
| AT notebooks | 5 |
| YH notebook | 1 |
| ZL R script | 1 |
| Original README | 1 |
| Total | 15 |

The [original README](../additional_analyses/README.md), filenames, notebook source cells, saved outputs, metadata, and R script contents were retained without modification. Git metadata remains at the repository root and was not moved into this directory.

This location does not classify the retained analyses as obsolete or exploratory. It preserves the original version while `plot/` and `process/` remain available for future organized workflows.

## Verification and execution scope

The [migration manifest](migration_manifest.tsv) records each original path, its new path, the source commit, Git blob identifier, file size, and SHA-256 checksum. Each relocated file was checked against its original committed content.

No notebook or R analysis was executed. Relative data paths, external storage references, and machine-specific package or working-directory paths still require attention before execution. See [the migration notes](migration_notes.md) for static observations and the limits of the checks performed.

## Future updates

Keep the preserved files identifiable when developing revised workflows. Document how new scripts and notebooks relate to their original versions, and record any analytical or path changes separately from this initial relocation.

When synchronizing framework updates into this Git checkout, explicitly preserve `additional_analyses/` and `.git/`. Directory names and Git ignore rules do not automatically protect files from deletion by external synchronization tools.

This guide was placed under `docs/` to avoid overwriting the original `additional_analyses/README.md`. See [the architecture guide](architecture_guide.md) for the wider organization and reproduction recommendations.

# Data Access and Provenance

This scaffold contains no datasets, accession numbers, download instructions for actual study files, or data-access assertions.

When the data dependencies have been established, document:

- The starting point for each workflow: raw data, a defined preprocessing output, or a fixed figure-input dataset.
- The repository, accession, DOI, or documented access procedure for each dataset.
- The expected local destination, file format, version, size, and checksum.
- The script and configuration that produce each intermediate file.
- The meaning of required metadata fields, identifier mappings, annotations, and coordinate units.
- Any access restrictions and the reproduction steps that depend on restricted inputs.

Use descriptive filenames or explicit subdirectories to distinguish different analysis objects. Avoid ambiguous dependencies on several unrelated files with the same generic object name.

The scaffold's `.gitignore` excludes data files by default. Store large datasets and intermediate objects in an appropriate data repository or archive, and document their retrieval here. If a small example dataset is later included in Git, narrow the ignore rules explicitly and explain its scope. Do not present example data as the complete manuscript dataset.

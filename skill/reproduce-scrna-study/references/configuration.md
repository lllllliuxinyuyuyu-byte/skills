# Configuration reference

## Metadata CSV

Provide one row per deposited matrix with these columns:

| Column | Meaning |
|---|---|
| `sample_id` | Stable accession or sample identifier |
| `file_name` | H5 filename or 10X matrix directory |
| `sample_group` | Plotting label |
| `cohort` | Experimental cohort; do not merge unrelated cohorts |
| `day` | Numeric time point or blank for references |
| `cell_line` | Cell line or tissue source |
| `protocol` | Differentiation or reference protocol |
| `stage` | Biological stage label |
| `sample_kind` | `derived` or `reference` |

Additional columns are retained as Seurat metadata.

## YAML sections

- `project`: name, data directory, metadata path, and output directory.
- `qc`: minimum detected genes, maximum UMIs, and maximum mitochondrial percentage.
- `analysis`: seed, variable-feature count, PCs, graph dimensions, resolution, reductions, and optional covariates to regress.
- `figures`: grouping variables, marker genes, color palettes, and paper-comparison labels.
- `gene_sets`: named gene vectors used for relative module scores.
- `assertions`: expected cell totals for declared sample cohorts.
- `narrative_claims`: claims whose stated sample IDs and time points are audited against metadata.

Resolve relative paths against the directory containing the YAML file. Environment variables written as `${NAME}` are expanded at runtime.

## Recommended cohort design

Represent a primary longitudinal series and an additional late-stage cohort as separate assertions even if both contain a Day 8 sample. This prevents a Day 12 conclusion from being attributed to a Day 0–8 core analysis.

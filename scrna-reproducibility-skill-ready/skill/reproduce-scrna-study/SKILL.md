---
name: reproduce-scrna-study
description: Reproduce and audit single-cell RNA-seq studies from 10X matrices with an R/Seurat workflow. Use when Codex needs to inspect GEO files, build sample metadata, run QC, PCA, clustering, tSNE or UMAP, score gene modules, compare generated plots with a paper, write figure captions, or verify that sample IDs, time points, cell totals, methods, and reported conclusions are internally consistent.
---

# Reproduce an scRNA-seq study

## Workflow

1. Inspect the supplied paper, repository metadata, and local matrix filenames before running analysis.
2. Create a metadata CSV and YAML configuration from the templates described in `references/configuration.md`.
3. Separate exact paper reproduction from exploratory additions. Reproduce the paper's dimensionality reduction first when its method is known.
4. Run the deterministic R workflow:

```text
Rscript scripts/run_scrna_workflow.R --config <config.yaml> --stage all
```

5. Run the independent metadata and narrative audit:

```text
Rscript scripts/audit_consistency.R --config <config.yaml>
```

6. Inspect every output figure. Do not infer cell types solely from clustering geometry.
7. Report deviations from the paper next to the affected result, including software-version, filtering, regression, and parameter differences.
8. Label every figure as either a paper comparison or a newly added analysis. Name the original figure when it is a comparison.

## Guardrails

- Never treat cells as independent biological replicates for inferential differential-expression tests.
- Never silently substitute UMAP for a paper's tSNE analysis.
- Never merge distinct cohorts merely because they share a time point.
- Keep raw matrices, large RDS objects, copyrighted PDFs, credentials, and personal application materials outside a public repository.
- Treat module scores as relative summaries, not proof of maturity, quiescence, or causal regulation.
- Define endothelial quiescence as a viable and reversible physiological state, not simply low proliferation.
- Preserve a machine-readable table linking samples, cohorts, time points, figures, and expected cell totals.

## Resources

- Read `references/configuration.md` when preparing or modifying a project configuration.
- Read `references/interpretation.md` before writing biological conclusions or figure captions.
- Use `scripts/run_scrna_workflow.R` for analysis and figure generation.
- Use `scripts/audit_consistency.R` to detect contradictions such as mentioning Day 12 when only Day 0–8 matrices were analyzed.

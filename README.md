# scRNA-seq Reproducibility Skill

[![R](https://img.shields.io/badge/R-%E2%89%A54.3-276DC3?logo=r)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-5.x-2C7FB8)](https://satijalab.org/seurat/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A configuration-driven **Codex Skill + R/Seurat workflow** for reproducing and auditing single-cell RNA-seq studies from 10X matrices. It was developed from a reconstruction of endothelial differentiation dataset **GSE131736**, while the workflow itself is dataset-independent.

中文简介：这是一个面向论文复现的单细胞 RNA 测序工作流，可完成10X矩阵质控、Seurat聚类、tSNE/UMAP、标志基因展示、状态评分、原文图形对照，以及样本—时间点—细胞数—文字结论的一致性检查。

## Why this repository is different

- Reproduces the reduction reported by the paper before adding exploratory views.
- Keeps biological cohorts explicit instead of silently merging samples with the same time point.
- Audits cell totals and narrative claims against the actual matrix metadata.
- Records deviations from published filters or incompletely disclosed parameters.
- Avoids cell-level pseudo-replication when biological replicates are insufficient.
- Packages the workflow as a reusable Codex Skill rather than a one-off notebook.

## Example result

GSE131736 contains a four-matrix H9 Day 0/4/6/8 core series with 21,369 deposited cells and a separate Day 8 replicate 3–Day 12 late cohort with 13,657 deposited cells. Keeping these cohorts separate prevents Day 12 conclusions from being attributed to the core trajectory.

![Core tSNE reproduction](examples/GSE131736/results/figures/core_tsne.png)

The example also adds a clearly labelled UMAP and relative module-score summary for hypothesis generation.

| Output | Purpose |
|---|---|
| QC retention | Show the effect of filtering by sample |
| tSNE | Compare with the dimensionality reduction used in the paper |
| UMAP | Explore broader state relationships; not labelled as an original-figure reproduction |
| Marker plots | Inspect endothelial, mesenchymal, pluripotency, and proliferation programs |
| Module heatmap | Summarize relative state programs across samples |
| Consistency audit | Verify matrix IDs, time points, cohorts, and reported cell totals |

## Repository layout

```text
skill/reproduce-scrna-study/   installable Codex Skill
examples/GSE131736/            public configuration and example outputs
tests/                         lightweight validation tests
.github/workflows/             GitHub Actions checks
```

Raw matrices, RDS objects, and the article PDF are intentionally excluded.

The example analysis environment is recorded in [`sessionInfo_example.txt`](examples/GSE131736/results/sessionInfo_example.txt).

## Quick start

Requirements: R 4.3 or newer and enough memory for the selected dataset. Large scRNA-seq datasets commonly require 16–32 GB RAM.

```bash
Rscript skill/reproduce-scrna-study/scripts/install_dependencies.R
```

Download the required 10X H5 matrices into a local directory, then set an environment variable:

```bash
export SCRNA_DATA_DIR=/path/to/GSE131736_RAW
```

On PowerShell:

```powershell
$env:SCRNA_DATA_DIR = "C:\path\to\GSE131736_RAW"
```

Validate without loading matrices:

```bash
Rscript skill/reproduce-scrna-study/scripts/run_scrna_workflow.R \
  --config examples/GSE131736/config.yaml --stage validate --dry-run
```

Run the complete workflow:

```bash
Rscript skill/reproduce-scrna-study/scripts/run_scrna_workflow.R \
  --config examples/GSE131736/config.yaml --stage all

Rscript skill/reproduce-scrna-study/scripts/audit_consistency.R \
  --config examples/GSE131736/config.yaml
```

Outputs are written to the configured `results/` directory.

## Install the Codex Skill

Copy `skill/reproduce-scrna-study` into your Codex skills directory, then invoke it with prompts such as:

```text
Use $reproduce-scrna-study to inspect these GEO 10X matrices, reproduce the paper's tSNE, add a UMAP, and audit every reported sample and cell total.
```

## Reproducibility boundaries

The example uses current R/Seurat versions rather than the historical Seurat release used by the article. Exact visual coordinates can therefore differ even when the major structure is reproduced. Module scores are descriptive; functional maturity and causal ECM effects require independent experiments.

## License

Code and Skill instructions are released under the MIT License. Dataset and article rights remain with their original providers.

## Data and primary publication

- GEO accession: **GSE131736**
- McCracken IR et al. *Transcriptional dynamics of pluripotent stem cell-derived endothelial cell differentiation revealed by single-cell RNA-seq*. European Heart Journal. 2020;41:1024–1036. DOI: [10.1093/eurheartj/ehz351](https://doi.org/10.1093/eurheartj/ehz351)

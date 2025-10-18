# scDECAF: Single-Cell Disentanglement by Canonical Factors

<img width="1627" height="780" alt="scDECAF overview" src="https://github.com/user-attachments/assets/396833cf-9716-4860-ac09-84d6fa37f6b3" />


---

## Overview

**scDECAF** is a statistical learning algorithm for single-cell RNA-seq analysis.
It identifies **gene signatures, states, and transcriptional programs** by learning vector representations of gene sets.

Through sparse selection, scDECAF highlights the most **biologically relevant programs**, improving interpretability of single-cell data.

---

## Installation

Requires **R ≥ 4.0.0**.

```r
install.packages("devtools")
devtools::install_github("DavisLaboratory/scDECAF")
```

---

## Quick Start

Here is a **minimal runnable example** with toy data:

```r
library(scDECAF)

# Simulated expression matrix (100 genes x 200 cells)
set.seed(123)
x <- matrix(rpois(100*200, lambda = 5), nrow = 100, ncol = 200)
rownames(x) <- paste0("Gene", 1:100)
colnames(x) <- paste0("Cell", 1:200)

# Define toy gene sets
genesetlist <- list(
  Pathway_A = rownames(x)[1:15],
  Pathway_B = rownames(x)[16:30],
  Pathway_C = rownames(x)[31:45]
)

# Highly variable genes (here, all genes for simplicity)
hvg <- rownames(x)

# Dummy 2D embedding (e.g., from PCA/UMAP)
cell_embedding <- matrix(rnorm(200*2), ncol = 2)
rownames(cell_embedding) <- colnames(x)

# Sparse selection of relevant gene sets
selected_gs <- pruneGenesets(
  data = x,
  genesetlist = genesetlist,
  hvg = hvg,
  embedding = cell_embedding,
  min_gs_size = 5,
  lambda = exp(-3)
)

# Build gene–set assignment matrix
target <- genesets2ids(
  x[match(hvg, rownames(x)), ],
  genesetlist[selected_gs]
)

# Compute gene-set scores
ann_res <- scDECAF(
  data = x,
  gs = target,
  hvg = hvg,
  k = 5,
  embedding = cell_embedding,
  n_components = min(2, ncol(target) - 1),
  max_iter = 2,
  thresh = 0.5
)

# Extract per-cell scores
scores <- attributes(ann_res)$raw_scores
head(scores[, 1:3])  # preview first few components
```

You can now add `scores` to your single-cell object (`SingleCellExperiment`, `Seurat`, or `AnnData`) and visualize them per cell.

---

## Input Requirements

* **data**: log-normalised single-cell expression matrix
* **genesetlist**: named list of gene sets
* **hvg**: highly variable genes
* **embedding**: reduced dimension embedding (UMAP, PCA, PHATE, etc.)
* **min_gs_size**: minimum gene set size
* **lambda**: shrinkage penalty
* **n_components**: number of CCA components
* **k**: nearest neighbors for refinement

---

## Demo 

| Workflow | Notebook |
|-----------|------|
| Pathway/gene signature screening and scoring | [Kang et al. (2018): 25K PBMC single cells](https://github.com/DavisLaboratory/scDECAF-reproducibility/blob/master/kang_pbmc/kang_pbmc.ipynb) |
| Optimization of sparsity operator | [Experimentation with the shrinkage penalty and gene set screening results in Kang et al. (2018)](https://github.com/DavisLaboratory/scDECAF-reproducibility/blob/master/kang_pbmc/sparse_mode_effect.ipynb)|
|PMBC COVID-19 analysis| [Combining reference atlas mapping and Milo analysis with scDECAF gene set screening](https://github.com/DavisLaboratory/scDECAF-reproducibility/blob/master/cite_pbmc/TotalVI_scDECAF_analysis-addMilo.ipynb)|
|Drug2cell analysis with scDECAF| [Running scDECAF with pre-computed Drug2cell scores in HECOA Organoid Atlas]|

---

## Reproducibility

Full analysis notebooks reproducing the manuscript are available in the [reproducibility repository](https://github.com/DavisLaboratory/scDECAF-reproducibility).

---

## Citation

If you use scDECAF, please cite:

> Hediyehzadeh, Whitfield, et al., *Identification of cell types, states and programs by learning gene set representations*, bioRxiv (2023).
> [https://doi.org/10.1101/2023.09.08.556842](https://doi.org/10.1101/2023.09.08.556842)

---

## License

This project is released under the same license as the Davis Laboratory repositories.

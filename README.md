![GitHub all releases](https://img.shields.io/github/downloads/DavisLaboratory/scDECAF/total)

scDECAF - single cell disentanglement by canonical factors
=================================================================

![Fig1_scDECAF_github](https://github.com/DavisLaboratory/scDECAF/assets/7257233/61723efa-2d7c-47d9-af0b-c7205b8b5644)



scDECAF is a statistical learning algorithm to identify cell types, states and programs in single-cell gene expression data using vector representation of gene sets. scDECAF improves biological interpretation by selecting a subset of most biologically relevant programs.



## Installation
(Requires R >= 4.0.0)

```
install.packages("devtools")
devtools::install_github("DavisLaboratory/scDECAF")
```
## Examples
See notebooks in the [reproducibility repository](https://github.com/DavisLaboratory/scDECAF-reproducibility)

## Quick start
```{r}
# sparse selection of most relevant genesets
# also plots number of genesets surviving the sparsity threshold
selected_gs <- pruneGenesets(data = x, genesetlist = my_genesets, hvg = hvg,
                            embedding = cell_embedding, min_gs_size = 3, lambda = exp(-3))
                            


# print selected genesets
as.character(selected_gs)



# print ranking/importance of geneset
head(attributes(selected_gs)$"glmnet_coef")



# subset on selected genesets from the full genesets list and prepare gene-geneset assignment binary matrix
rownames(cell_embedding) = cell_names
target <- genesets2ids(x[match(hvg, rownames(x)),], my_genesets[selected_gs])



# compute geneset scores per cell for the sparse set of selected genesets 
ann_res <- scDECAF(data = x, gs = target, standardize = FALSE, 
                   hvg = hvg, k = 10, embedding = cell_embedding,
                   n_components = ncol(target) - 1, max_iter = 2, thresh = 0.5)
                   


# get geneset scores per cell for the sparse set of genesets
scores_constrained = attributes(ann_res)$raw_scores

# plot scores
```

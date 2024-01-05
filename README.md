
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

scDECAF takes the followings as input:

**data** : A numeric matrix of log-normalised single cell gene expression (SCT normalisation from seurat, scran- or scanpy- normalised data). Rows are genes, columns are the cells.

**genesetlist**: A list of lists. Each element of the list is a list of gene IDs or symbols (depending on `rownames(data)`) in a gene set. The outer list has to be named.

**hvg**: Character vector of highly variable genes in `data`. If the data is already subsetted on HVGs, then set this to `rownames(data)`

**embedding**: A numeric matrix 2-D or higher dimensional embedding of the cells, e.g. UMAP, PCA, PHATE, Diffusion components etc. Rows are cells, columns are the dimension of the data in the reduced dimension space.

**min_gs_size** : Scalar. Minimum number of genes in a gene set (after considering hvgs).

**lambda**: Shrinkage regulariser penalty

**K** Scalar. This iss number of components in the CCA model. Has to be smaller than the number of gene sets in `genesetlist` or the prunned `genesetlist`.




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



# compute geneset scores per cell for the sparse set of selected genesets. `K` is number of components in the CCA model. 
ann_res <- scDECAF(data = x, gs = target, standardize = FALSE, 
                   hvg = hvg, k = 10, embedding = cell_embedding,
                   n_components = ncol(target) - 1, max_iter = 2, thresh = 0.5)
                   


# get geneset scores per cell for the sparse set of genesets
scores_constrained = attributes(ann_res)$raw_scores


```
You can now store the scores to your data container (`sce`, `seuratObj`, `anndata` etc) and visualise the scores per cell.

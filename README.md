scDECAF - Celltype and phenotype annotation in single cell RNAseq 
=================================================================


<img src="https://user-images.githubusercontent.com/7257233/107848582-ad5a2980-6e48-11eb-8590-ddd00223e9c5.png" width="700px" align="center">



What can you do with scDECAF?
---------------------
*scDECAF* is a tool for *mapping phenotype and celltype similarities* in single cell RNAseq from a collection of genesets or makers such as those available from [cellMarker](http://biocc.hrbmu.edu.cn/CellMarker/), [PanglaoDB](https://panglaodb.se/) and [MSigDB](http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp). Hence, *scDECAF* does not require reference datasets or any quantitative information about desired phenotypes or celltypes. *scDECAF* learns a vector representation for cells and genesets and measures their similarity in the vector space. *scDECAF* is not an integration nor a classification method, and is not designed for geneset enrichment analysis. It is rather a method for learning similarities between gene expression profiles of the cells and a geneset collection.



What are the other similar methods?
------------------------
**Celltype annotation** is typically done by examining the expression of canonical markers, or using reference or atlas datasets such as the
[Human Lung Cell Atlas](https://hlca.ds.czbiohub.org/). The reference-based methods may involve integration of reference with the query data. 
Examples include:

* [scVI](https://www.nature.com/articles/s41592-018-0229-2) which uses Deep Generative Models
* [scArches](https://www.biorxiv.org/content/10.1101/2020.07.16.205997v1) also uses Generative Models with Transfer Learning

Some popular reference-based methods do not involve direct integration of reference and query datasets. They annotate cells based on the correlation between reference and query:

* [SingleR](https://bioconductor.org/packages/release/bioc/html/SingleR.html)
* [scmap](https://www.nature.com/articles/nmeth.4644)


The well known annotation tool, [Seurat](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8), offers reference-based annotation strategies both with and without implicit integration of datasets. Another category of methods use pre-trained classifiers for celltype annotation. [scPred](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1862-5), for example, uses SVM classifier to predict the identity of the cells in a query dataset. 

Amongs the mentioned methods, *scDECAF* is only comparable to [garnett](https://www.nature.com/articles/s41592-019-0535-3); they both take genesets/marker lists as input and use no quantitative information other than the expression profiles in the dataset of interest.


#### Install from Github
Requires R >= 4.0.0

```
devtools::install_github("DavisLaboratory/scDECAF")
```



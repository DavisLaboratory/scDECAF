scDECAF - Celltype and phenotype annotation in single cell RNAseq 
=================================================================


<img src="https://user-images.githubusercontent.com/7257233/107848582-ad5a2980-6e48-11eb-8590-ddd00223e9c5.png" width="700px" align="center">



What can you do with scDECAF?
---------------------
scDECAF is a tool for mapping phenotype and celltype similarities in single cell RNAseq from a collection of genesets or makers such as those available from [cellMarker](http://biocc.hrbmu.edu.cn/CellMarker/), [PanglaoDB](https://panglaodb.se/) and [MSigDB](http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp). Hence, scDECAF does not require reference datasets or any quantitative information about desired phenotypes or celltypes. scDECAF learns a vector representation for cells and genesets and measures their similarity in the vector space. scDECAF is not an integration nor a classification method, and is not designed for geneset enrichment analysis. It is rather a method for learning similarities between gene expression profiles of the cells and a geneset collection.



What are the other similar methods?
------------------------



### Install from Github
Requires R >= 4.0.0

```
devtools::install_github("DavisLaboratory/scDECAF")
```



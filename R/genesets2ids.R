#' Creates a binary matrix from a list of genes
#'
#' @param y numeric matrix. The gene expression matrix from Highly variable genes
#' @param genesetlist list of genesets. Each element of the list is a geneset.
#'     Each geneset can contain genes or gene ids. However, the genes in the genesets must match the rownames in \code{y}.
#'
#' @return numeric matrix. A binary matrix is returned where rows are genes, columns are genesets.
#' @author Soroor Hediyeh-zadeh
#' @export
genesets2ids <- function(y, genesetlist){
  INDX <- matrix(0, nrow = nrow(y), ncol= length(genesetlist),
                 dimnames = list(rownames(y), names(genesetlist)))
  gs_idx <- limma::ids2indices(genesetlist, rownames(y))

  for (i in names(gs_idx)){
    INDX[gs_idx[[i]],i] <- 1
  }

  INDX[rowSums(INDX)>0,]
}

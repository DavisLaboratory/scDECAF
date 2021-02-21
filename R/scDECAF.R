#' Reference-free cell type and phenotype annotation by learning geneset representations
#'
#'
#' @param data numeric matrix of log-normalised values. Rows are genes, columns are cells or samples
#' @param gs binary matrix created by \code{genesets2ids}. Rows are genes, columns are genesets. Each column represents a geneset.
#' @param hvg character vector of Highly Variable Genes.
#' @param cca.k numeric. Number of canonical factors used to learn geneset representations.
#' @param standardize logical. Should input expression data be standardized? default to TRUE.
#' @param k numeric. Number of nearest neighbors used to refine cell-label assignment. Default to 30.
#' @param thresh numeric. Assignment confidence probability threshold. Default to 0.5.
#' @param embedding numeric matrix. UMAP, tSNE or other cell embedding that captures transcriptional similarity of the cells.
#' @param n_components numeric. Number of canonical factors to use. Default to all.
#' @param max_iter numeric. currently not used.
#'
#' @return Data frame of annotation results
#'
#' @details
#' The cell-geneset assignments are determined by cosine similarity of vector representations
#' of cells and genesets by CCA. There is the a “re-assignment” step where the umap (embedding) coordinates of the cells are
#' used to refine labels based on a local (\code{k}=3-15) or global (\code{k}>15) neighborhood of
#' the cell, as the lower-dimensional embedding of the cells captures the transcriptional heterogeneity of
#' cells, and cells with similar transcriptional profiles should ideally have same labels/phenotypes.
#' Users have their own choice of embedding, as there are various methods to obtain low-dimensional embedding
#' for the cells, and each offer unique properties in terms of preservation of local or global structures, or transcriptional
#' heterogeneity. The input gene expression matrix is preferred to contain filtered, log-normalised values.
#' However, the algorithm is expected to work with RPKM and/or log-RPKM values as well.
#'
#' @author Soroor Hediyeh-zadeh
#'
#' @importFrom stats aggregate median
#' @export
scDECAF <- function(data, gs, hvg, cca.k= NULL, standardize=TRUE,
                    k = 30, thresh = 0.5, embedding,
                    n_components = ncol(gs) - 1, max_iter = 1){
  # To Do:
  # add nperms and trace to main function argument
  # change default for standardize - should there be a difference for numeric vs binary genesets?
  # if coeffs are zero, return NA instead of removing the cells -- done
  # NAs and NaNs in input
  # evalution of the warning statement

  if(class(data)[1] %in% c("dgCMatrix", "dgTMatrix")) data <- as.matrix(data)
  if(any(colSums(data) == 0)) stop("Empty cells detected. Please remove cells with zero counts.")
  if(any(colSums(abs(gs)) == 0)) stop("Empty genesets detected. Please remove genesets with no genes.")
  if(any(!hvg %in% rownames(data))) stop("Some or all of HVG features are not present in data matrix.")
  if(all(is.null(rownames(embedding)))) stop("Row names are required for embedding.")
  if(ncol(data) != nrow(embedding)) stop("Number of cells differ between the input and embedding.")
  if (is.null(cca.k)) cca.k <- min(dim(gs)) - 1



  xt <- data[match(rownames(gs), rownames(data)),]
  #xt <- xt[, colSums(xt) > 0]


  cell_names <- colnames(xt)
  gs_names <- colnames(gs)


  message(paste("Learning geneset representations by CCA using ", cca.k, " dimensions...", sep=""))
  perm.out <- PMA::CCA.permute(x=xt, z=gs, typex = "standard", typez = "standard",
                               nperms = 10, standardize = standardize, trace = TRUE)
  CCA.out <- PMA::CCA(x=xt, z=gs, typex = "standard", typez = "standard",
                      penaltyx = 0.7, penaltyz = 0.7, K= cca.k, v=perm.out$v.init[,1], standardize = standardize,
                      trace = TRUE)


  cell.coef <- CCA.out$u
  gs.coef.1 <- CCA.out$v



  invalidCell <- (rowSums(cell.coef) == 0)
  n_invalid <-  sum(invalidCell)

  if (any(invalidCell)) {
    warning(paste("Detected", n_invalid , "cells with 0 canonical coefficients")) # sum(rowSums(cell.coef) == 0)

    unassigned_cells <- cell_names[invalidCell]
    cell.coef <- cell.coef[!invalidCell,]
    cell_names <- cell_names[!invalidCell]

  }



  P_t <- gs.coef.1
  P_s_1 <- cell.coef


  P_t <- P_t/sqrt(rowSums(P_t^2))
  P_s_1 <- P_s_1/sqrt(rowSums(P_s_1^2))


  message("Assigning cells to genesets in latent space...")


  # note indexing fails if number of genesets < 10
  nn <- FNN::get.knnx(data =  P_t[,1:n_components],
                      query =  P_s_1[,1:n_components], k = 1)



  transported_labels <- colnames(gs)[nn$nn.index[,1]]
  names(transported_labels) <- cell_names




  message("Verifying assignments by weighted k-nn...")
  message(paste("using", k, "nearest neighbors and confidence threshold", thresh))

  #xhvg <- data[match(hvg, rownames(data)), match(cell_names, colnames(data))]
  # cell.nn.hvg <- FNN::get.knn(t(xhvg), k = k)


  embedding <- embedding[match(cell_names, rownames(embedding)),]
  cell.nn.latent <- FNN::get.knn(embedding, k = k)
  #cell.nn.latent <- FNN::get.knn(t(xhvg), k = k)




  Sigma <- array(matrixStats::rowSds(cell.nn.latent$nn.dist), dim(cell.nn.latent$nn.dist))



  cell.nn.weights <- exp(-0.5*((cell.nn.latent$nn.dist)^2)/Sigma)

  reassigned_celltype <- transported_labels
  uncertainty <- rep(0, length(transported_labels))

  while(max_iter > 0){
    max_iter = max_iter - 1
    cell.nn.latent <- FNN::get.knn(embedding, k = k) # note there can be a difference in #cells retained by CCA and total num cells
    n_unique_nn_labels_per_cell <- apply(cell.nn.latent$nn.index, 1, FUN = function(x) length(unique(reassigned_celltype[x])))

    #I <- array(0,dim(cell.nn.weights))
    for ( i in which(n_unique_nn_labels_per_cell > 1)){
      #query_label <- transported_labels[i]
      # compute probabilities for each label in the neighborhood
      nn_labels <- reassigned_celltype[cell.nn.latent$nn.index[i,]]
      query_label <- unique(nn_labels)
      nn_weight <- cell.nn.weights[i,]

      label_score <- vector("numeric")

      if(!all(nn_weight == 0)){
        for (j in query_label) {

          #score <- sum(ifelse(nn_labels==query_label,1,0)*nn_weight)/sum(nn_weight)
          label_score[j] <- sum(ifelse(nn_labels==j,1,0)*nn_weight)/sum(nn_weight)

        }

        #message(paste("cell", i))
        #print(label_score)
        reassigned_celltype[i] <- names(label_score)[which.max(label_score)]
        uncertainty[i] <- 1 - label_score[which.max(label_score)]

        #I[i,] <- ifelse(nn_labels==query_label,1,0)

      }else{
        uncertainty[i] <- 1
      }
    }
  }
  # D <- rowSums(indics*cell.nn.weights)/rowSums(cell.nn.weights)
  # p_uncertain <- 1-D



  reassigned_celltype[uncertainty > thresh] <- "Unknown"

  centroids <- aggregate(.~label, FUN= median, data = data.frame(embedding,
                                                                 label = reassigned_celltype))

  nn_centroids <- FNN::get.knnx(centroids[, -1], embedding, k = nrow(centroids))
  labels_by_nearest_centroid <- centroids[nn_centroids$nn.index[,1],1] # first nearest centroid

  # uncertainty is proportional to the distance of the cell from label centroid
  dist_to_centroid <- rep(0, length(transported_labels))
  s <- which(reassigned_celltype!= labels_by_nearest_centroid)
  for(cell in s){
    idx <- match(reassigned_celltype[cell], centroids$label)
    dist_to_centroid[cell] <- 1 - exp(-0.5*(nn_centroids$nn.dist[cell, idx]^2))
  }



  message("cell-geneset assignment completed.")





  ## percent neighbors preserved in CCA and hvg subspace/latent spaces -----
  N = 200
  cell.nn.latent <- FNN::get.knn(embedding, k = N)
  cell.nn.cca.latent <- FNN::get.knn(P_s_1, k = N)
  nn_indices <- cbind(cell.nn.latent$nn.index, cell.nn.cca.latent$nn.index)
  prop_preserved_nns <- apply(nn_indices, 1, FUN=function(x) {
    length(intersect(x[1:N], x[(N+1):length(x)]))/N  # at the moment this is the intersection of neighbors, should we compute knn purity?
    #length(intersect(x[1:N], x[(N+1):length(x)]))/length(union(x[1:N], x[(N+1):length(x)]))
  })


  if(n_invalid == 0){
    unassigned_cells <- NULL
    n_unassigned <- 0
  } else {
    n_unassigned <- n_invalid
  }



  out <- data.frame(cell = c(names(transported_labels), unassigned_cells),
                      pred_celltype = c(transported_labels, rep("unassigned", n_unassigned)),
                      score = c(0.5*(2-nn$nn.dist^2), rep(NA, n_unassigned)),
                      reassigned_celltype = c(reassigned_celltype,  rep("unassigned", n_unassigned)),
                      uncertainty = c(uncertainty, rep(NA, n_unassigned)),
                      prop_neighbourhood_preserved = c(prop_preserved_nns, rep(NA, n_unassigned)),
                      dist_to_centroid = c(dist_to_centroid, rep(NA, n_unassigned)))



  out <- out[match(colnames(data), out$cell),]

  # score for all genesets
  nngs <- FNN::get.knnx(data =  P_t[,1:n_components],
                        query =  P_s_1[,1:n_components],
                        k = nrow(P_t))



  scores <- 0.5*(2-(nngs$nn.dist^2))

  indicies <- nngs$nn.index

  for (i in seq_len(nrow(scores))) {
    scores[i,] <- scores[i,match(gs_names, gs_names[indicies[i,]])]
    }
  rownames(scores) <- names(transported_labels)
  colnames(scores) <- gs_names
  attr(out,"raw_scores") <- scores

  return(out)
}






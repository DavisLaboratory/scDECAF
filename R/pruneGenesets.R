#' Penalized generalized linear models for model-based selection of relevant genesets
#'
#'
#'
#'
#' @param data numeric matrix of log-normalised values. Rows are genes, columns are cells or samples
#' @param genesetlist list of genesets. Each element of the list is a geneset.
#'     Each geneset can contain genes or gene ids. However,
#'     the genes in the genesets must match the rownames in \code{data}.
#'
#' @param embedding numeric matrix. UMAP, PHATE or other cell embedding that captures transcriptional similarity of the cells.
#' @param hvg character vector of Highly Variable Genes.
#' @param min_gs_size numeric. Minimum number of genes per geneset.
#' @param suppress_plot logical. Should the plot of glm cross-validation results be suppressed? default to FALSE.
#' @param lambda numeric or character. The lasso penalty. This is set to optimal value selected
#' by cross-validation by default. Should be set to exp(log(lambda)), where
#' log(lambda) is chosen based on the diagnostic plot returned by the function
#' (requires \code{suppress_plot = FALSE}).
#' @param gamma numeric. Controls the behavior of shrinkage operator e.g. \code{gamma=0.5}
#'     is equivalent to the Elastic Net model
#' @param nfolds numeric. Number of folds for cross-validation.
#'
#' @return a character vector of selected genesets. Regression coefficients from the glmnet model
#' can be accessed via the \code{glmnet_coef} attribute.
#' @author Soroor Hediyeh-zadeh
#' @importFrom graphics plot
#' @export
pruneGenesets <- function(data, genesetlist,
                    embedding,
                    hvg = NULL,
                    min_gs_size = 3,
                    suppress_plot=FALSE,
                    lambda = 'lambda.1se',
                    gamma=0, nfolds = 10){

  if(!is.null(hvg) & any(!hvg %in% rownames(data))) stop("Some or all of HVG features are not present in data matrix.")
  if(class(data)[1] %in% c("dgCMatrix", "dgTMatrix")) data <- as.matrix(data)

  if(!is.null(hvg)){
    target <- genesets2ids(data[match(hvg, rownames(data)),], genesetlist)
  } else{
    target <- genesets2ids(data, genesetlist)
  }


  target <- target[,colSums(target) > min_gs_size]

  mean_expr_per_gs <-  t(data[match(rownames(target), rownames(data)),]) %*% target
  mean_expr_per_gs <- apply(mean_expr_per_gs, 1, FUN=function(x) x/colSums(target))

  embedding <- data.matrix(embedding)

  message("Computing optimal shrinkage value by cross-validation")
  cvfit <- glmnet::cv.glmnet(x = t(mean_expr_per_gs), y=embedding ,
                             family = "mgaussian", nfolds = nfolds,
                             gamma = gamma, standardize = FALSE) # gamma = 0.5 elastic net

  if(!suppress_plot) plot(cvfit)
  if(lambda == 'lambda.1se') lambda <- cvfit$lambda.1se

  message(paste("Fitting penalized multi-response gaussian GLM with alpha", round(lambda,3)))
  mfit <- glmnet::glmnet(x = t(mean_expr_per_gs), y=embedding , family = "mgaussian",
                 gamma = gamma, lambda = lambda, standardize = FALSE)

  message("Returning selected genesets with non-zero regression coefficients")
  coefs <- glmnet::coef.glmnet(mfit)[[1]][-1,"s0"] # -1 to drop coef for intercept
  selected_gs <- names(coefs)[coefs!=0]

  attr(selected_gs, "glmnet_coef") <- glmnet::coef.glmnet(mfit)


  return(selected_gs)




}


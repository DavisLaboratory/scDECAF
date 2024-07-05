#' Fit multiple shrinkage values for Penalized generalized linear modelling for selection of relevant pre-computed Drug2Cell scores
#'
#'
#'
#'
#' @param data numeric matrix of Drug2Cell scores. columns are Drug2Cell scores, rows are cells or samples.
#' @param embedding numeric matrix. UMAP, PHATE or other cell embedding that captures transcriptional similarity of the cells.
#' @param suppress_plot logical. Should the plot of glm cross-validation results be suppressed? default to FALSE.
#' @param lambda numeric or character. The lasso penalty. This is set to optimal value selected
#' by cross-validation by default. Should be set to exp(log(lambda)), where
#' log(lambda) is chosen based on the diagnostic plot returned by the function
#' (requires \code{suppress_plot = FALSE}).
#' @param gamma numeric. Controls the behavior of shrinkage operator e.g. \code{gamma=0.5}
#'     is equivalent to the Elastic Net model
#' @param nfolds numeric. Number of folds for cross-validation.
#' @details
#' This function should be used to determine the optimal shrinkage value \code{lambda} for geneset selection. 
#' This function is applicable to any pre-computed gene set and signature scores.
#' 
#' @return a character vector of selected genesets. Regression coefficients from the glmnet model
#' can be accessed via the \code{glmnet_coef} attribute.
#' @author Soroor Hediyeh-zadeh
#' @importFrom graphics plot
#' @export
pruneDrug2CellScores.fitLambda <- function(data, 
                    embedding,
                    suppress_plot=FALSE,
                    lambda = 'lambda.1se',
                    gamma=0, nfolds = 10, trace.it = 0){

  
  if(class(data)[1] %in% c("dgCMatrix", "dgTMatrix")) data <- as.matrix(data)

  # embedding <- data.matrix(embedding)

  message("Computing optimal shrinkage value by cross-validation")
  # x_t <- t(data) 
  cvfit <- glmnet::cv.glmnet(x = data, y=embedding ,
                             family = "mgaussian", nfolds = nfolds,
                             gamma = gamma, standardize = FALSE, trace.it = trace.it) # gamma = 0.5 elastic net

  if(!suppress_plot) plot(cvfit)
  if(lambda == 'lambda.1se') lambda <- cvfit$lambda.1se

  message(paste("Fitting penalized multi-response gaussian GLM with alpha", round(lambda,3)))
  mfit <- glmnet::glmnet(x = data, y=embedding , family = "mgaussian",
                 gamma = gamma, lambda = lambda, standardize = FALSE)

  message("Returning selected genesets with non-zero regression coefficients")
  coefs <- glmnet::coef.glmnet(mfit)[[1]][-1,"s0"] # -1 to drop coef for intercept
  selected_gs <- names(coefs)[coefs!=0]

  attr(selected_gs, "glmnet_coef") <- glmnet::coef.glmnet(mfit)


  return(selected_gs)




}





#' Penalized generalized linear models for model-based selection of relevant pre-computed Drug2Cell scores
#'
#'
#'
#'
#' @param data numeric matrix of Drug2Cell scores. columns are Drug2Cell scores, rows are cells or samples.
#' @param embedding numeric matrix. UMAP, PHATE or other cell embedding that captures transcriptional similarity of the cells.
#' @param lambda numeric or character. The lasso penalty. This is set to optimal value selected
#' by cross-validation by default. Should be set to exp(log(lambda)), where
#' log(lambda) is chosen based on the diagnostic plot returned by the function
#' (requires \code{suppress_plot = FALSE}), for example using \code{pruneDrug2CellScores.fitLambda}
#' @param gamma numeric. Controls the behavior of shrinkage operator e.g. \code{gamma=0.5}
#'     is equivalent to the Elastic Net model
#' @details
#' This function is applicable to any pre-computed gene set and signature scores.
#' 
#' @return a character vector of selected genesets. Regression coefficients from the glmnet model
#' can be accessed via the \code{glmnet_coef} attribute.
#' @author Soroor Hediyeh-zadeh
#' @importFrom graphics plot
#' @export
pruneDrug2CellScores <- function(data, 
                    embedding,
                    lambda = 'lambda.1se',
                    gamma=0, trace.it = 0){

  
  if(class(data)[1] %in% c("dgCMatrix", "dgTMatrix")) data <- as.matrix(data)

  # embedding <- data.matrix(embedding)

  if(lambda == 'lambda.1se') lambda <- cvfit$lambda.1se

  message(paste("Fitting penalized multi-response gaussian GLM with alpha", round(lambda,3)))
  mfit <- glmnet::glmnet(x = data, y=embedding , family = "mgaussian",
                 gamma = gamma, lambda = lambda, standardize = FALSE)

  message("Returning selected genesets with non-zero regression coefficients")
  coefs <- glmnet::coef.glmnet(mfit)[[1]][-1,"s0"] # -1 to drop coef for intercept
  selected_gs <- names(coefs)[coefs!=0]

  attr(selected_gs, "glmnet_coef") <- glmnet::coef.glmnet(mfit)


  return(selected_gs)

}
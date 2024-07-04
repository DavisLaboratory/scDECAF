#' Penalized generalized linear models for model-based selection of relevant pre-computed Drug2Cell scores
#'
#'
#'
#'
#' @param data numeric matrix of Drug2Cell scores. Rows are Drug2Cell scores, columns are cells or samples
#' @param embedding numeric matrix. UMAP, PHATE or other cell embedding that captures transcriptional similarity of the cells.
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
#'
#' @importFrom graphics plot
#' @export
pruneDrug2CellScores <- function(data, 
                    embedding,
                    suppress_plot=FALSE,
                    lambda = 'lambda.1se',
                    gamma=0, nfolds = 10){

  
  if(class(data)[1] %in% c("dgCMatrix", "dgTMatrix")) data <- as.matrix(data)

  embedding <- data.matrix(embedding)

  message("Computing optimal shrinkage value by cross-validation")
  cvfit <- glmnet::cv.glmnet(x = t(data), y=embedding ,
                             family = "mgaussian", nfolds = nfolds,
                             gamma = gamma, standardize = FALSE) # gamma = 0.5 elastic net

  if(!suppress_plot) plot(cvfit)
  if(lambda == 'lambda.1se') lambda <- cvfit$lambda.1se

  message(paste("Fitting penalized multi-response gaussian GLM with alpha", round(lambda,3)))
  mfit <- glmnet::glmnet(x = t(data), y=embedding , family = "mgaussian",
                 gamma = gamma, lambda = lambda, standardize = FALSE)

  message("Returning selected genesets with non-zero regression coefficients")
  coefs <- glmnet::coef.glmnet(mfit)[[1]][-1,"s0"] # -1 to drop coef for intercept
  selected_gs <- names(coefs)[coefs!=0]

  attr(selected_gs, "glmnet_coef") <- glmnet::coef.glmnet(mfit)


  return(selected_gs)




}

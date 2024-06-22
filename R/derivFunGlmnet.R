derivFunGlmnet <- function(u, parms) {
  d <- length(u)
  du <- vapply(
    seq_len(d),
    \(j) {
      degVecs <- parms$glmnetFits[[j]]$degVec
      coefs <- parms$glmnetFits[[j]]$coef
      features <- DEEButil::evaluateMonomials(matrix(u, nrow=1), degVecs)
      if (length(coefs) == length(features) + 1) {
        sum(c(1, features) * coefs) # add intercept
      } else {
        sum(features * coefs)
      }
    },
    numeric(1))
  return(du)
}


prepareParmsGlmnet <- function(parms, opts) {

  opts <- asOpts(opts, c("Glmnet", "DerivFun"))

  degree <- opts$polyDeg
  d <- getDim(parms$trajs)
  degVecs <- DEEButil::getMonomialExponents(d, degree)
  featureValues <- DEEButil::evaluateMonomials(parms$trajs$state, degVecs)

  fits <- lapply(
    seq_len(d),
    \(j) {
      fit <- glmnet::cv.glmnet(
        featureValues,
        parms$trajs$deriv[,j],
        lambda = opts$lambda,
        typeMeasure = opts$typeMeasure,
        nFolds = opts$nFolds,
        relax = opts$relax,
        gamma = opts$gamma,
        alpha = opts$alpha,
        thresh = opts$thresh,
        maxit = opts$maxit,
        standardize = opts$standardize,
        intercept = opts$intercept)
      coefs <-
        if (opts$relax) {
          glmnet::coef.relaxed(fit, s = opts$lambdaChoice, gamma = opts$gammaChoice)
        } else {
          glmnet::coef.glmnet(fit, s = opts$lambdaChoice, gamma = opts$gammaChoice)
        }
      list(degVec = degVecs[coefs@i, ], coef = coefs@x)
    }
  )

  parms$glmnetFits <- fits
  return(parms)
}

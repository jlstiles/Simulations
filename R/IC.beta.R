#' @title IC.beta
#' @description computes IC for beta coefficients of logistic regression
#' mainly a helper function.  
#' @param W, matrix or data.frame of covariates
#' @param A, a binary vector of treatment assignments
#' @param Y, a binary vector of outcomes
#' @param Qform, a formula for Y in terms of the covariates as input in glm
#' 
#' @return  a list with elements IC_beta and fit, the glm fit object.  
#' @export
IC.beta = function(data,OC=NULL, Ynode, Qform, verbose = FALSE, parallelize = FALSE) {
  n = nrow(data)
  # This option is for when feeding in sequential regression
  if (!is.null(OC)) data[,Ynode] = OC
  # we only fit on non deaths or uncensored
  cens = is.na(data[,Ynode])
  data = data[!cens,]
  n1 = nrow(data)
  # form the design matrix based on the formula
  X = model.matrix(Qform,data)
  X = as.data.frame(X[,-1])
  # fit the regression
  Y = data[,Ynode]
  if (!verbose) {
    fit = suppressWarnings(stats::glm(Y~.,data=X,
                                      family='binomial'))
  } else {
    fit = stats::glm(Y~.,data=X,family='binomial')
  }
  goods = 1:(ncol(X)+1)
  if (any(is.na(coef(fit)))) {
    print(paste0("you have a singular covariance matrix so we will refit without these variables",
                 paste(names(coef(fit))[is.na(coef(fit))], collapse = " ")))
    goods = which(!is.na(coef(fit)))
    X = X[,(goods-1)]
    if (!verbose) {
      fit = suppressWarnings(stats::glm(Y~.,data=X,
                                        family='binomial'))
    } else {
      fit = stats::glm(Y~.,data=X,family='binomial')
    }
  }
  # predictions over data
  Qk = predict(fit,type='response')
  
  X$Y = NULL
  X=cbind(int = rep(1,n1),X)
  
  # calculate the score
  # score_beta = sapply(1:n1,FUN = function(x) {
  #   X[x,]*(Y[x]-Qk[x])
  # })
  
  if (parallelize) cores = getOption("mc.cores",parallel::detectCores()) else cores = 1L
  score_beta = mclapply(1:n1,FUN = function(x) {
    X[x,]*(Y[x]-Qk[x])
  }, mc.cores = getOption("mc.cores", cores))
  
  LL = length(score_beta[[1]])
  score_beta = vapply(1:length(score_beta), FUN = function(x){
    return(unlist(score_beta[[x]]))
  }, FUN.VALUE = rep(1,LL))
  
  # averaging hessians to approx the true average then invert as per IC
  hessian = mclapply(1:n1,FUN = function(x) {
    mat = -(1-Qk[x])*Qk[x]*as.numeric(X[x,])%*%t(as.numeric(X[x,]))
    return(mat)
  }, mc.cores = getOption("mc.cores", cores))
  
  # M1 = summary(fit)$cov.unscaled*n1
  fisher = -Reduce('+', hessian)/n1
  M = solve(fisher)
  
  Xfull = matrix(rep(NA,(n*ncol(X))),nrow = n)
  Xfull = as.data.frame(Xfull)
  Xfull[!cens,] = X
  colnames(Xfull) = colnames(X)
  # calculate the IC for beta
  IC_beta = matrix(rep(0, nrow(M)*n), nrow = nrow(M))
  IC_beta[,!cens] = apply(score_beta,2,FUN = function(x) M%*%as.numeric(x))
  IC_beta = IC_beta
  return(list(IC_beta = IC_beta, fit = fit, X = Xfull, hessian = M, goods = goods, cens = cens))
  
}

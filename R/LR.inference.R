
#' @title LR.inference
#' @description Function that gives inference for logistic regression plug-in
#' estimators of ATE and treatment effect function variance (VTE)  
#' @param W, matrix or data.frame of covariates
#' @param A, a binary vector of treatment assignments
#' @param Y, a binary vector of outcomes
#' @param Qform, a formula for Y in terms of the covariates as input in glm
#' @param alpha, significance level for the (1-alpha)100 percent CI's. 0.05 is default
#' @param simultaneous.inference, TRUE if user wants simultaneous confidence
#' bounds for both ATE and blip variance at level alpha. default is FALSE
#' 
#' @return  if simultaneous.inference is specified as TRUE then will return a vector giving
#' pt estimate, left and right bound for ATE, simultaneous ATE CI, blip variance, 
#' and simultaneous blip variance.  Otherwise gives pt estimate, left and right bound for ATE
#' and blip variance.  
#' @export
#' @example /inst/examples/example_LR_inference.R
LR.inference = function(W, A, Y, Qform, alpha = .05, simultaneous.inference = FALSE) {
  
  n = length(Y)
  X = as.data.frame(cbind(A,W,Y))
  X0 = X1 = X
  X0$A = 0
  X1$A = 1
  
  newdata = rbind(X, X1, X0)
  newdata = model.matrix(Qform,newdata)
  newdata = as.data.frame(newdata[,-1])
  colnames(newdata)[2:ncol(newdata)] = paste0("X",2:ncol(newdata))
  
  # fit the regression
  Qfit = stats::glm(Y~.,data=newdata[1:n,],
                    family='binomial')
  # predictions over data, A=1 and A=0
  Qk = predict(Qfit,type='response')
  Q1k = predict(Qfit,newdata=newdata[(n+1):(2*n),],type='response')
  Q0k = predict(Qfit,newdata=newdata[(2*n+1):(3*n),],type='response')
  
  # covariates and treatment for convenient use
  X = newdata
  X$Y = NULL
  X=cbind(int = rep(1,n),X)
  head(X)
  # calculate the score
  score_beta = sapply(1:n,FUN = function(x) {
    X[x,]*(Y[x]-Qk[x])
  })
  
  # averaging hessians to approx the deriv of hessian and mean then inverse
  hessian = lapply(1:n,FUN = function(x) {
    mat = -(1-Qk[x])*Qk[x]*as.numeric(X[x,])%*%t(as.numeric(X[x,]))
    return(mat)
  })
  fisher = -Reduce('+', hessian)/n
  M = solve(fisher)
  
  # calculate the IC for beta
  IC_beta = apply(score_beta,2,FUN = function(x) M%*%as.numeric(x))
  
  # SE_test = apply(IC_beta,1,sd)*sqrt(n-1)/n
  # SE_test
  
  blip = Q1k-Q0k
  ate = mean(Q1k-Q0k)
  # calculate the deriv to mult by IC_beta
  deriv1 = rowMeans(vapply(1:n, FUN = function(x) {
    return((1-Q1k[x])*Q1k[x]*as.numeric(X[(n+x),])-(1-Q0k[x])*Q0k[x]*
             as.numeric(X[(2*n+x),]))
  }, FUN.VALUE=rep(1,ncol(X))))
  
  deriv = rowMeans(vapply(1:n, FUN = function(x) {
    return(2*(blip[x]-ate)*((1-Q1k[x])*Q1k[x]*as.numeric(X[(n+x),])-(1-Q0k[x])*Q0k[x]*
                              as.numeric(X[(2*n+x),])))
  }, FUN.VALUE=rep(1,ncol(X))))
  
  psi = var(blip)
  # connect both parts of IC to form the full one
  
  IC = apply(IC_beta,2,FUN = function(x) t(deriv)%*%x) + (blip - ate)^2 - psi
  IC1 = apply(IC_beta, 2, FUN = function(x) t(deriv1)%*%x) + blip -ate
  # standard error
  SE = sd(IC)*sqrt((n-1))/n
  SE1 = sd(IC1)*sqrt((n-1))/n
  
  qq = qnorm(1-alpha/2)
  CI = c(bv_delta = psi, left = psi - qq*SE, right = psi + qq*SE)
  
  CI_ate = c(ate_delta = ate, left = ate - qq*SE1, right = ate + qq*SE1)
  
  if (simultaneous.inference) {
    corM = stats::cor(data.frame(IC=IC, IC1=IC1))
    Z = rmvnorm(1000000,c(0,0),corM)
    zabs = apply(Z,1,FUN = function(x) max(abs(x)))
    zscore = quantile(zabs, 1-alpha)
    CI_simul_ate = c(ate_deltasimul = ate, left = ate - zscore*SE1, right = ate + zscore*SE1)  
    CI_simul_bv = c(bv_deltasimul = psi, left = psi - zscore*SE, right = psi + zscore*SE) 
    
    return(c(CI_ate, CI_simul_ate, CI, CI_simul_bv))
  } else {
    return(c(CI_ate, CI))
  }
}

#' @export
sim.longTSM = function(n, dag, gform, Qform, formulas, formulags, setA, T_end, 
                       Lnodes, Anodes, Ynodes, tmle = TRUE, gcomp)
{
  
  OdatL = sim(Ddyn, n = n)
  OdatL$ID = NULL
  data = OdatL
  data_ltmle = data
  if (T_end >1){
    for (t in 1:(T_end-1)) {
      data_ltmle[data_ltmle[,Ynodes[t]]==1,Ynodes[t+1]] = 1 
    }
  }
  
  nombre = Ynodes[T_end]
  Yend = grep(nombre, colnames(data))
  
  res = ltmle(data=data_ltmle[,1:Yend], Anodes=Anodes, Lnodes = Lnodes, 
              Ynodes=Ynodes,survivalOutcome = TRUE, abar = setA, 
              Qform = Qform, gform = gform, 
              gbounds = c(0.000001,1),deterministic.g.function = NULL,  
              estimate.time = TRUE, gcomp = gcomp, iptw.only = FALSE, stratify = FALSE,
              deterministic.Q.function = NULL,variance.method = "ic", 
              observation.weights = NULL, id = NULL)
  
  TSMinfo = long.TSM(data = data, Ynodes = Ynodes, Anodes = Anodes, 
                     formulas = formulas, formulas_g = formulags, setA = setA, tmle = tmle)
  
  TSMinfo$CI
  sd(TSMinfo$IC)*sqrt(n-1)/n
  summary(res)[[1]]$std.dev
  sd(res$IC$iptw)*sqrt(n-1)/n
  
  c(summary(res)[[1]]$estimate, summary(res)[[1]]$CI)
  
  CIs = c(c(summary(res)[[1]]$estimate, summary(res)[[1]]$CI),summary(res)[[1]]$std.dev,
          TSMinfo$CI, sd(TSMinfo$IC)*sqrt(n-1)/n,sd(res$IC$iptw)*sqrt(n-1)/n)
  names(CIs)[c(2:3,6:7)] = c("left", "right")
  if (tmle & !gcomp) {
    names(CIs)[c(1,4,5,8,9)] = c("tmle", "SE tmle","LRdelta_tmle","SE LR_tmle", "SE iptw")
  } else if (tmle & gcomp){
    names(CIs)[c(1,4,5,8,9)] = c("gcomp", "SE gcomp","LRdelta_tmle","SE LR_tmle", "SE iptw")
  } else if (!tmle & gcomp){
    names(CIs)[c(1,4,5,8,9)] = c("gcomp", "SE gcomp","LRdelta","SE LR", "SE iptw")  
  } else {
    names(CIs)[c(1,4,5,8,9)] = c("tmle", "SE tmle","LRdelta","SE LR", "SE iptw") 
    }
  return(CIs)
}


#' @export
sim_single = function(n, g0, Q0, form) {
  X = gendata(n, g0, Q0)
  X$Y = NULL
  X1 = X0 = X
  X1$A = 1
  X0$A = 0
  newX = rbind(X,X1,X0)
  W = X
  W$A = NULL
  
  results_glm <- glm(form,data = simdata,family = binomial())
  Qk = predict(results_glm, newdata = X, type = 'response')
  Q1k = predict(results_glm, newdata = X1,type = 'response')
  Q0k = predict(results_glm, newdata = X0,type = 'response')
  initest = var(Q1k - Q0k)
  
  gfit = glm(A~.,data = X, family = binomial())
  gk = predict(gfit, type = 'response')
  
  initdata = data.frame(A = X$A, Y = simdata$Y, gk = gk, Qk = Qk, Q1k = Q1k, Q0k = Q0k)
  
  sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
  
  ci_tmle = ci_gentmle(sigma_info)[c(2,4,5)]
  return(c(ci_tmle, initest))
}

#' @export
sim_lr = function(n, g0, Q0, formQ, formG) {
  
  data = gendata(n, g0, Q0)
  
  X=data
  X$Y = NULL
  X1 = X0 = X
  X1$A = 1
  X0$A = 0
  newX = rbind(X,X1,X0)
  W = X
  W$A = NULL
  
  results <- glm(formQ,data = data,family=binomial())
  Qk = predict(results, newdata = X, type = 'response')
  Q1k = predict(results, newdata = X1,type = 'response')
  Q0k = predict(results, newdata = X0,type = 'response')
  
  gfit = glm(formG,data = X, family = 'binomial')
  gk = predict(gfit, type = 'response')
  initest = var(Q1k - Q0k)
  initest_ATE = mean(Q1k - Q0k)
  
  initdata = data.frame(A = X$A, Y = data$Y, gk = gk, Qk = Qk, 
                        Q1k = Q1k, Q0k = Q0k)
  
  sigmait_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100,g.trunc = 1e-2)
  sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
  
  simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE,param_sigmaATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000, g.trunc = 1e-2,
                                 simultaneous.inference = TRUE)
  
  ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                               submodel = submodel_logit, loss = loss_loglik,
                               approach = "full", max_iter = 100, g.trunc = 1e-2)
  
  steps = c(sigma_info$steps,sigmait_info$steps, simul_info$steps, ATE_info$steps)
  converge = c(sigma_info$converge, sigmait_info$converge,simul_info$converge, ATE_info$converge)
  
  ci_sig = ci_gentmle(sigma_info)[c(2,4,5)]
  ci_sigit = ci_gentmle(sigmait_info)[c(2,4,5)]
  ci_simul = ci_gentmle(simul_info)[2,c(2,4,5)]
  ci_simulATE = ci_gentmle(simul_info)[1,c(2,4,5)]
  ci_ATE = ci_gentmle(ATE_info)[c(2,4,5)]
  
  cis = c(ci_sig, ci_sigit, ci_simul, ci_simulATE, 
          ci_ATE)
  names(cis)[c(1,4,7,10,13)] =
    c("sig", "sigit","simul","simulATE","ATE")
  names(converge) = names(steps) = names(cis)[c(1:3,5)]
  results = c(cis, initest = initest, initest_ATE = initest_ATE,steps = steps, 
              converge = converge)
  
  return(results)
}


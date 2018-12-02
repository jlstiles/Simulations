library(boot)

#' @export
sim_hal = function(data, gform = NULL, Qform = NULL, V = 10, single = FALSE, estimator, method, gn = NULL,
                   cvhal = TRUE, parallel = FALSE) {
  # n=100
  # single = TRUE
  # V = 10
  # g0 = g0_1
  # Q0 = Q0_2
  # Qform = formula("Y~A*(W1+W2+W3+W4)")
  # gform = formula("A~.")
  # data = gendata(n, g0, Q0)
  # head(simdata)
  n = nrow(data)
  Y=data$Y
  A=data$A
  X=data
  X$Y = NULL
  X1 = X0 = X
  X1$A = 1
  X0$A = 0
  newdata = rbind(X,X1,X0)
  W = X
  W$A = NULL
  
  if (all(A==1 | A == 0)) familyG = binomial() else familyG = gaussian()
  if (all(Y==1 | Y == 0)) familyQ = binomial() else familyQ = binomial()
  
  folds = origami::make_folds(n, V=V)
  stack = lapply(folds, FUN = function(x) {
    # x=folds[[5]]
    if (V==1) {tr = val = 1:n} else {
      tr = x$training_set
      val = x$validation_set
    }
    nt=length(tr)
    nv = length(val)
    
    data = data[tr,]
    Y = Y[tr]
    X = X[tr,]
    newtr = c(val, (n+val),(2*n+val))
    newdata = newdata[newtr,]
    
    if (cvhal){
      halresults <- hal::hal(Y = Y,newX = newdata,
                             X = X, family = familyQ,
                             verbose = FALSE, parallel = parallel)
      
      Qk = halresults$pred[1:nv]
      Q1k = halresults$pred[nv+1:nv]
      Q0k = halresults$pred[2*nv+1:nv]} else {
        halresults = stats::glm(Qform, data = data, family = familyQ)
        Qk = predict(halresults, newdata = newdata[1:nv,], type = 'response')
        Q1k = predict(halresults, newdata = newdata[nv + 1:nv,], type = 'response')
        Q0k = predict(halresults, newdata = newdata[2*nv + 1:nv,], type = 'response')
      }
    
    
    A = A[tr]
    W1 = W[tr,]
    newW = W[val,]
    
    if (is.null(gn)) {
      if (cvhal){
        halresultsG <- hal::hal(Y = A,newX = newW,
                                X = W1, family = familyG,
                                verbose = FALSE, parallel = parallel)
        gk = halresultsG$pred[1:nv]} else {
          halresultsG = stats::glm(gform, data = X, family = familyG)
          gk = predict(halresultsG, newdata = newW, type = 'response')
        }
    } else {
      gk = gn[val]
    }
    
    return(list(Qk = Qk, Q0k = Q0k, Q1k = Q1k, gk = gk,inds = x$validation_set))
  })
  
  Qk = unlist(lapply(stack, FUN = function(x) x$Qk))
  Q1k = unlist(lapply(stack, FUN = function(x) x$Q1k))
  Q0k = unlist(lapply(stack, FUN = function(x) x$Q0k))
  gk = unlist(lapply(stack, FUN = function(x) x$gk))
  inds = unlist(lapply(stack, FUN = function(x) x$inds))
  
  initest = var(Q1k - Q0k)
  initest_ATE = mean(Q1k - Q0k)
  
  if (any(estimator %in% c("simul 1step", "simul line", "simul full"))) {
    simul_lr = TRUE
  } else {simul_lr = FALSE}
  if (!is.null(Qform)) {
    ci_lr = LR.inference(W = W[inds,], A = A[inds], Y = Y[inds], Qform = Qform,
                         simultaneous.inference = simul_lr)
  }
  initdata = data.frame(A = A[inds], Y = Y[inds], gk = gk, Qk = Qk, Q1k = Q1k, Q0k = Q0k)
  
  if ("single 1step" %in% estimator) {
    sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
    
    ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000,g.trunc = 1e-2)
    steps = c(steps_bv1step = sigma_info$steps, steps_ate1step = ATE_info$steps)
    converge = c(sigma_info$converge, ATE_info$converge)
    cis = c(ci_gentmle(sigma_info)[c(2,4,5)], ci_gentmle(ATE_info)[c(2,4,5)])  
    
    names(cis)[c(1,4)] = names(steps) = names(converge) = c("bv1step", "ate1step")
  } else {
    cis = c()
    converge = c()
    steps = c()}
  
  if ("single iterative" %in% estimator) {
    sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 10000, g.trunc = 1e-2)
    
    ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "full", max_iter = 10000,g.trunc = 1e-2)
    steps1 = c(sigma_info$steps, ATE_info$steps)
    converge1 = c(sigma_info$converge, ATE_info$converge)
    cis1 = c(ci_gentmle(sigma_info)[c(2,4,5)], ci_gentmle(ATE_info)[c(2,4,5)])  
    names(cis1)[c(1,4)] = names(steps1) = names(converge1) = c("bv", "ate")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if ("simul 1step" %in% estimator) {
    simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "recursive", max_iter = 10000, g.trunc = 1e-2,
                                   simultaneous.inference = TRUE)
    
    steps1 = c(simul = simul_info$steps)
    converge1 = c(simul = simul_info$converge)
    cis1 = c(ci_gentmle(simul_info)[2,c(2,4,5)], ci_gentmle(simul_info)[1,c(2,4,5)])
    names(cis1)[c(1,4)] = c("bv_simul", "ate_simul")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if ("simul line" %in% estimator) {
    simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "line", max_iter = 100, g.trunc = 1e-2,
                                   simultaneous.inference = TRUE)
    
    steps1 = c(steps_line = simul_info$steps)
    converge1 = c(con_line = simul_info$converge)
    cis1 = c(ci_gentmle(simul_info)[2,c(2,4,5)], ci_gentmle(simul_info)[1,c(2,4,5)])
    names(cis1)[c(1,4)] = c("bv_line", "ate_line")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if ("simul line" %in% estimator) {
    simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100, g.trunc = 1e-2,
                                   simultaneous.inference = TRUE)
    
    steps1 = c(steps_full = simul_info$steps)
    converge1 = c(con_full = simul_info$converge)
    cis1 = c(ci_gentmle(simul_info)[2,c(2,4,5)], ci_gentmle(simul_info)[1,c(2,4,5)])
    names(cis1)[c(1,4)] = c("bv_full", "ate_full")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if (method == "method.NNloglik"){
    Qrisk = with(initdata, -mean(Y*log(Qk) + (1 - Y)*log(1 - Qk)))
    Grisk = with(initdata, -mean(A*log(gk) + (1 - A)*log(1 - gk)))
  } else {
    Qrisk = with(initdata, mean((Y - Qk)^2))
    Grisk = with(initdata, -mean(A*log(gk) + (1 - A)*log(1 - gk)))
  }
  if (!is.null(Qform)) {
    cis = c(cis, ci_lr)
  }
  
  results = c(cis, initest = initest, initest_ATE = initest_ATE, steps = steps, converge = converge, 
              Qrisk = Qrisk, Grisk = Grisk)
}

#' @title sim_cv
#' @description Function that simulates data and performs TMLE estimates
#' data is 4 covariates and binary treatment and outcome.  The covariates
#' are generated according to the gendata function.
#' @param n, sample size
#' @param g0, treatment mechanism function, call g0_linear to see the format
#' @param Q0, outcome model function, call Q0 linear to see the format--WILL
#' be more generalized but for now rather limited
#' @param SL.library, SuperLearner library for outcome predictions
#' @param SL.libraryG, SuperLearner Library for treatment mechanism
#' @param method, SuperLearner meta fitting method
#' @param cv, set to TRUE for CV-TMLE
#' @param V, number of folds for the CV tmle
#' @param SL, number of folds for each superlearner
#' @param gform, a linear form to specify for estimating pscore
#' @param Qform, a linear form to specify for estimating outcome prediction
#' @param estimator, a character vector containing any set of "single 1step" for
#' one step tmle single param estimates for ATE and blip variance, "single iterative"
#' for the same with iterative tmle, or "simul 1 step", "simul line", "simul full"
#' to compute simultaneous estimates and CI's for ATE and blip variance.  line, full
#' and 1 step are just different targeting methods for tmle
#' @param dgp, a list containing an element named DF for the data.frame with A, Y and 
#' covariates which are named whatever, BV0 and ATE0 for true blip variance and 
#' average treatment effect respectively.  
#' @return  a vector with the following elements in this order:
#' TMLE pt estimates and confidence intervals each with estimate, 
#' left and right bounds and initial estimates for BV and ATE
#' Superlearner coefficients and risks for both pscore and outcome
#' estimation, tmle pt estimates, CI's for BV and ATE for Qform model 
#' assumed plus initial estimates for these as well as pt estimates and
#' CI's for BV and ATE using the delta method, ie sandwich estimator
#' under non-parametric model
#' @export
#' @example /inst/examples/example_sim_cv.R
sim_cv = function(n, g0, Q0, SL.library, SL.libraryG, method = "method.NNLS", 
                  cv = TRUE, V = 10, SL = 10L, gform, Qform, estimator, dgp = NULL, gendata.fcn) {
  
  if (!is.null(dgp)) {
    data = dgp$DF
    BV0 = dgp$BV0
    ATE0 = dgp$ATE0
    blip_n = dgp$blip_n
  } else {
    data = gendata.fcn(n, g0, Q0)
  }
  
  if (is.vector(g0)) gn = g0 else gn = NULL
  
  Y = data$Y
  A = data$A
  
  X = data
  X$Y = NULL
  X1 = X0 = X
  X0$A = 0
  X1$A = 1
  newdata = rbind(X,X1,X0)
  
  W = X
  W$A = NULL
  
  single_info = sim_hal(data = data, gform = gform, Qform = Qform, V = 10, estimator = estimator,
                        method = method, gn = gn, cvhal = FALSE)
  
  stack = SL.stack1(Y=Y, X=X, A=A, W=W, newdata=newdata, method=method, 
                    SL.library=SL.library, SL.libraryG=SL.libraryG,cv = cv, V = V, SL = SL, gn = gn)
  # proc.time() - time
  
  initdata = stack$initdata
  
  # stack
  initest = with(initdata,var(Q1k - Q0k))
  initest_ATE = with(initdata, mean(Q1k - Q0k))
  
  if ("single 1step" %in% estimator) {
    sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
    
    ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000,g.trunc = 1e-2)
    steps = c(steps_bv1step = sigma_info$steps, steps_ate1step = ATE_info$steps)
    converge = c(sigma_info$converge, ATE_info$converge)
    cis = c(ci_gentmle(sigma_info)[c(2,4,5)], ci_gentmle(ATE_info)[c(2,4,5)])  
    
    names(cis)[c(1,4)] = names(steps) = names(converge) = c("bv1step", "ate1step")
  } else {
    cis = c()
    converge = c()
    steps = c()
  }
  
  if ("single iterative" %in% estimator) {
    sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100, g.trunc = 1e-2)
    
    ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "full", max_iter = 100,g.trunc = 1e-2)
    steps1 = c(sigma_info$steps, ATE_info$steps)
    converge1 = c(sigma_info$converge, ATE_info$converge)
    cis1 = c(ci_gentmle(sigma_info)[c(2,4,5)], ci_gentmle(ATE_info)[c(2,4,5)])  
    names(cis1)[c(1,4)] = names(steps1) = names(converge1) = c("bv", "ate")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if ("simul 1step" %in% estimator) {
    simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "recursive", max_iter = 10000, g.trunc = 1e-2,
                                   simultaneous.inference = TRUE)
    
    steps1 = c(simul = simul_info$steps)
    converge1 = c(simul = simul_info$converge)
    cis1 = c(ci_gentmle(simul_info)[2,c(2,4,5)], ci_gentmle(simul_info)[1,c(2,4,5)])
    names(cis1)[c(1,4)] = c("bv_simul", "ate_simul")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if ("simul line" %in% estimator) {
    simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "line", max_iter = 100, g.trunc = 1e-2,
                                   simultaneous.inference = TRUE)
    
    steps1 = c(steps_line = simul_info$steps)
    converg1 = c(con_line = simul_info$converge)
    cis1 = c(ci_gentmle(simul_info)[2,c(2,4,5)], ci_gentmle(simul_info)[1,c(2,4,5)])
    names(cis1)[c(1,4)] = c("bv_line", "ate_line")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if ("simul full" %in% estimator) {
    simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100, g.trunc = 1e-2,
                                   simultaneous.inference = TRUE)
    
    steps1 = c(steps_full = simul_info$steps)
    converg1 = c(con_full = simul_info$converge)
    cis1 = c(ci_gentmle(simul_info)[2,c(2,4,5)], ci_gentmle(simul_info)[1,c(2,4,5)])
    names(cis1)[c(1,4)] = c("bv_full", "ate_full")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if (!is.null(dgp)) {
    results = list(res = c(cis, initest = initest, initest_ATE = initest_ATE, steps = steps, converge = converge, 
                           Qcoef = stack$Qcoef, Gcoef = stack$Gcoef, Qrisk = stack$Qrisk, 
                           Grisk = stack$Grisk, single = single_info, BV0 = BV0, ATE0 = ATE0), blip_n = blip_n)
  } else {
    results = c(cis, initest = initest, initest_ATE = initest_ATE, steps = steps, converge = converge, 
                Qcoef = stack$Qcoef, Gcoef = stack$Gcoef, Qrisk = stack$Qrisk, 
                Grisk = stack$Grisk, single = single_info)
  }
  return(results)
  
}

# input data.frame with A, Y and covariates spit out lr CI based on delta method
# remember the order of vars is A, mainterms, then interactions
#' @export
sim_cv4 = function(n, g0, Q0, SL.library, SL.libraryG, method = "method.NNLS", 
                   cv = TRUE, V = 10, SL = 10L, gform, Qform, estimator, dgp = NULL) {
  
  # g0 = dgps[[3]]$PGn 
  # Q0 = Q0_trig1
  # SL.library = SL.library
  # SL.libraryG = SL.libraryG 
  # method = "method.NNLS"
  # cv = TRUE
  # V = 10
  # SL = 10L
  # gform = gform 
  # Qform = Qform 
  # estimator = c("single 1step")
  # dgp = dgps[[3]]
  
  if (!is.null(dgp)) {
    data = dgp$DF
    BV0 = dgp$BV0
    ATE0 = dgp$ATE0
    blip_n = dgp$blip_n
  } else {
    data = gendata(n, g0, Q0)
  }
  
  if (is.vector(g0)) gn = g0 else gn = NULL
  
  X = data
  X1 = X0 = X
  X0$A = 0
  X1$A = 1
  Y = data$Y
  A = data$A
  W = X
  W$A = NULL
  W$Y = NULL
  newdata = rbind(X,X1,X0)
  
  if (any(c("simul 1step", "simul line", "simul full") %in% estimator)) {
    lrsingle = FALSE
  } else {
    lrsingle = TRUE
  }
  
  single_info = sim_hal(data = data, gform = gform, Qform = Qform, V = 10, estimator = estimator,
                        method = method, gn = gn, cvhal = FALSE)
  
  mainform = paste0(paste(colnames(data)[2:4],"+",collapse=""),colnames(data)[5])
  mainform
  squares = paste0(paste0("I(",colnames(data)[2:5]),"^2)")
  squares
  squares = paste0(paste(squares[1:3],"+",collapse=""),squares[4])
  squares
  mainsq = paste0(mainform,"+",squares)
  mainsq
  mainsq.int = paste0("Y~A*(",mainsq,")")
  mainsq.int = formula(mainsq.int) 
  mainsq.int
  
  newdata = model.matrix(mainsq.int,newdata)
  newdata = as.data.frame(newdata[,-1])
  colnames(newdata)[2:ncol(newdata)] = paste0("X",2:ncol(newdata))
  head(newdata)
  X = newdata[1:n,]
  X$Y = NULL
  
  # time = proc.time()
  
  stack = SL.stack1(Y=Y, X=X, A=A, W=W, newdata=newdata, method=method, 
                    SL.library=SL.library, SL.libraryG=SL.libraryG,cv = cv, V = V, SL = SL, gn = gn)
  # proc.time() - time
  
  initdata = stack$initdata
  
  # stack
  initest = with(initdata,var(Q1k - Q0k))
  initest_ATE = with(initdata, mean(Q1k - Q0k))
  
  if ("single 1step" %in% estimator) {
    sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
    
    ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000,g.trunc = 1e-2)
    steps = c(steps_bv1step = sigma_info$steps, steps_ate1step = ATE_info$steps)
    converge = c(sigma_info$converge, ATE_info$converge)
    cis = c(ci_gentmle(sigma_info)[c(2,4,5)], ci_gentmle(ATE_info)[c(2,4,5)])  
    
    names(cis)[c(1,4)] = names(steps) = names(converge) = c("bv1step", "ate1step")
  } else {
    cis = c()
    converge = c()
    steps = c()
  }
  
  if ("single iterative" %in% estimator) {
    sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100, g.trunc = 1e-2)
    
    ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "full", max_iter = 100,g.trunc = 1e-2)
    steps1 = c(sigma_info$steps, ATE_info$steps)
    converge1 = c(sigma_info$converge, ATE_info$converge)
    cis1 = c(ci_gentmle(sigma_info)[c(2,4,5)], ci_gentmle(ATE_info)[c(2,4,5)])  
    names(cis1)[c(1,4)] = names(steps1) = names(converge1) = c("bv", "ate")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if ("simul 1step" %in% estimator) {
    simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "recursive", max_iter = 10000, g.trunc = 1e-2,
                                   simultaneous.inference = TRUE)
    
    steps1 = c(simul = simul_info$steps)
    converge1 = c(simul = simul_info$converge)
    cis1 = c(ci_gentmle(simul_info)[2,c(2,4,5)], ci_gentmle(simul_info)[1,c(2,4,5)])
    names(cis1)[c(1,4)] = c("bv_simul", "ate_simul")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if ("simul line" %in% estimator) {
    simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "line", max_iter = 100, g.trunc = 1e-2,
                                   simultaneous.inference = TRUE)
    
    steps1 = c(steps_line = simul_info$steps)
    converg1 = c(con_line = simul_info$converge)
    cis1 = c(ci_gentmle(simul_info)[2,c(2,4,5)], ci_gentmle(simul_info)[1,c(2,4,5)])
    names(cis1)[c(1,4)] = c("bv_line", "ate_line")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if ("simul full" %in% estimator) {
    simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100, g.trunc = 1e-2,
                                   simultaneous.inference = TRUE)
    
    steps1 = c(steps_full = simul_info$steps)
    converg1 = c(con_full = simul_info$converge)
    cis1 = c(ci_gentmle(simul_info)[2,c(2,4,5)], ci_gentmle(simul_info)[1,c(2,4,5)])
    names(cis1)[c(1,4)] = c("bv_full", "ate_full")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if (!is.null(dgp)) {
    results = list(res = c(cis, initest = initest, initest_ATE = initest_ATE, steps = steps, converge = converge, 
                           Qcoef = stack$Qcoef, Gcoef = stack$Gcoef, Qrisk = stack$Qrisk, 
                           Grisk = stack$Grisk, single = single_info, BV0 = BV0, ATE0 = ATE0), blip_n = blip_n)
  } else {
    results = c(cis, initest = initest, initest_ATE = initest_ATE, steps = steps, converge = converge, 
                Qcoef = stack$Qcoef, Gcoef = stack$Gcoef, Qrisk = stack$Qrisk, 
                Grisk = stack$Grisk, single = single_info)
  }
  return(results)
  
}

####
####
####
####

#' @export
sim_cvnew = function(n, g0, Q0, SL.library, SL.libraryG, method = "method.NNLS",
                     V = 10, SL = 10L, gform, Qform, estimator, dgp = NULL, gendata.fcn) {
  
  if (!is.null(dgp)) {
    data = dgp$DF
    BV0 = dgp$BV0
    ATE0 = dgp$ATE0
    blip_n = dgp$blip_n
  } else {
    data = gendata.fcn(n, g0, Q0)
  }
  
  if (is.vector(g0)) gn = g0 else gn = NULL
  
  Y = data$Y
  A = data$A
  
  X = data
  X$Y = NULL
  X1 = X0 = X
  X0$A = 0
  X1$A = 1
  newdata = rbind(X,X1,X0)
  
  W = X
  W$A = NULL
  
  single_info = sim_hal(data = data, gform = gform, Qform = Qform, V = 10, estimator = estimator,
                        method = method, gn = gn, cvhal = FALSE)
  
  stack = SL.stack1(Y=Y, X=X, A=A, W=W, newdata=newdata, method=method, 
                    SL.library=SL.library, SL.libraryG=SL.libraryG,cv = TRUE, V = V, SL = SL, gn = gn)
  # proc.time() - time
  
  
  initdata = stack$initdata
  
  # stack
  initest = with(initdata,var(Q1k - Q0k))
  initest_ATE = with(initdata, mean(Q1k - Q0k))
  
  if ("single 1step" %in% estimator) {
    sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_EB2), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
    
    ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000,g.trunc = 1e-2)
    steps = c(steps_EB21step = sigma_info$steps, steps_ate1step = ATE_info$steps)
    converge = c(sigma_info$converge, ATE_info$converge)
    D_blipvar = sigma_info$Dstar[[1]] - 2*ATE_info$tmleests*ATE_info$Dstar[[1]]
    psi_blipvar = sigma_info$tmleests - ATE_info$tmleests^2
    SE_blipvar = sd(D_blipvar)*sqrt(n-1)/n
    ci_blipvar = c(-1.96*SE_blipvar, 1.96*SE_blipvar) + psi_blipvar
    
    cis = c(psi_blipvar, ci_blipvar , ci_gentmle(ATE_info)[c(2,4,5)])  
    
    names(cis)[c(1,4)] = names(steps) = names(converge) = c("bv1step", "ate1step")
  } else {
    cis = c()
    converge = c()
    steps = c()
  }
  
  if ("single iterative" %in% estimator) {
    sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_EB2), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100, g.trunc = 1e-2)
    
    ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "full", max_iter = 100,g.trunc = 1e-2)
    steps1 = c(sigma_info$steps, ATE_info$steps)
    converge1 = c(sigma_info$converge, ATE_info$converge)
    D_blipvar = sigma_info$Dstar[[1]] - 2*ATE_info$tmleests*ATE_info$Dstar[[1]]
    psi_blipvar = sigma_info$tmleests - ATE_info$tmleests^2
    SE_blipvar = sd(D_blipvar)*sqrt(n-1)/n
    ci_blipvar = c(-1.96*SE_blipvar, 1.96*SE_blipvar) + psi_blipvar
    
    cis1 = c(psi_blipvar, ci_blipvar , ci_gentmle(ATE_info)[c(2,4,5)])  
    names(cis1)[c(1,4)] = names(steps1) = names(converge1) = c("bv", "ate")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if ("simul 1step" %in% estimator) {
    simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_EB2), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "recursive", max_iter = 10000, g.trunc = 1e-2,
                                   simultaneous.inference = FALSE)
    
    steps1 = c(simul = simul_info$steps)
    converge1 = c(simul = simul_info$converge)
    
    D_blipvar = simul_info$Dstar[[2]] - 2*simul_info$tmleests[1]*simul_info$Dstar[[1]]
    psi_blipvar = simul_info$tmleests[2] - simul_info$tmleests[1]^2
    SE_blipvar = sd(D_blipvar)*sqrt(n-1)/n
    ci_blipvar = c(-1.96*SE_blipvar, 1.96*SE_blipvar) + psi_blipvar
    
    cis1 = c(psi_blipvar, ci_blipvar, ci_gentmle(simul_info)[1,c(2,4,5)])
    names(cis1)[c(1,4)] = c("bv_simul", "ate_simul")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if ("simul line" %in% estimator) {
    simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_EB2), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "line", max_iter = 100, g.trunc = 1e-2,
                                   simultaneous.inference = FALSE)
    
    steps1 = c(simul = simul_info$steps)
    converge1 = c(simul = simul_info$converge)
    
    D_blipvar = simul_info$Dstar[[2]] - 2*simul_info$tmleests[1]*simul_info$Dstar[[1]]
    psi_blipvar = simul_info$tmleests[2] - simul_info$tmleests[1]^2
    SE_blipvar = sd(D_blipvar)*sqrt(n-1)/n
    ci_blipvar = c(-1.96*SE_blipvar, 1.96*SE_blipvar) + psi_blipvar
    
    cis1 = c(psi_blipvar, ci_blipvar, ci_gentmle(simul_info)[1,c(2,4,5)])
    names(cis1)[c(1,4)] = c("bv_simul", "ate_simul")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if ("simul full" %in% estimator) {
    simul_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE, param_EB2), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100, g.trunc = 1e-2,
                                   simultaneous.inference = FALSE)
    
    steps1 = c(simul = simul_info$steps)
    converge1 = c(simul = simul_info$converge)
    
    D_blipvar = simul_info$Dstar[[2]] - 2*simul_info$tmleests[1]*simul_info$Dstar[[1]]
    psi_blipvar = simul_info$tmleests[2] - simul_info$tmleests[1]^2
    SE_blipvar = sd(D_blipvar)*sqrt(n-1)/n
    ci_blipvar = c(-1.96*SE_blipvar, 1.96*SE_blipvar) + psi_blipvar
    
    cis1 = c(psi_blipvar, ci_blipvar, ci_gentmle(simul_info)[1,c(2,4,5)])
    names(cis1)[c(1,4)] = c("bv_simul", "ate_simul")
    
    cis = c(cis, cis1)
    converge = c(converge, converge1)
    steps = c(steps, steps1)
  }
  
  if (!is.null(dgp)) {
    results = list(res = c(cis, initest = initest, initest_ATE = initest_ATE, steps = steps, converge = converge, 
                           Qcoef = stack$Qcoef, Gcoef = stack$Gcoef, Qrisk = stack$Qrisk, 
                           Grisk = stack$Grisk, single = single_info, BV0 = BV0, ATE0 = ATE0), blip_n = blip_n)
  } else {
    results = c(cis, initest = initest, initest_ATE = initest_ATE, steps = steps, converge = converge, 
                Qcoef = stack$Qcoef, Gcoef = stack$Gcoef, Qrisk = stack$Qrisk, 
                Grisk = stack$Grisk, single = single_info)
  }
  return(results)
  
}

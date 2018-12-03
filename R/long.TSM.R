#' @title long.TSM
#' @description computes delta method inferencec for logistic regression plug-in 
#' estimator of treatment specific mean for survival data in wide form, 
#' including pt treatment. Also does logistic regression plug-in using clever
#' covariate in the regression, inference computed using delta method as well.
#' Note, in performing tmle, pscores are obtained via user supplied regression
#' formulas for each time point.  Formulas for the conditional means are also
#' user supplied.  
#' @param data, data.frame of variables in time ordering from left to right
#' @param Ynodes, character vector of time-ordered Ynodes
#' @param Anodes, character vector of time-ordered Anodes
#' @param formulas, list of formulas for the conditional means
#' @param setA the value to which you intervene on A,vector of length that of Anodes
#' @param alpha significance level for two-sided CI.
#' @param formulas_g NULL by default but must be filled in if tmle = TRUE
#' @param tmle option to put a clever covariate in the regression.  If the propensity
#' score is well-specified, will give good coverage asymptotically despite mispecifying
#' the outcome model.  Otherwise gives delta method inference for the "closest" linear
#' model fit using clever covariate in your regression.  Set to TRUE for RCT's.
#' @param parallelize FALSE by default (still in development)
#' @return  a list with elements, CI for the confidence interval and  IC for the 
#' influence curve. 
#' @export
#' @example /inst/examples/example_longTSM.R
long.TSM = function(data, Ynodes, Anodes, formulas, setA, alpha = .05,
                    formulas_g = NULL, tmle = FALSE, parallelize = FALSE)
{
  
  n = nrow(data)
  endtime = length(setA)
  times = 1:endtime
  # we will work backwards in time
  times = times[order(times,decreasing = TRUE)]
  
  # if tmle is TRUE we need to get all the IC's for the pscore estimation and the new design
  #  which includes the clever covariates as well as adds them to the formula for the regression' 
  
  if (tmle) {
    # ic list for each pscore regression
    ICg = list()
    # set up data frames for both H and H when intervened upon
    H = data.frame(matrix(rep(0,(length(Anodes)*n)), nrow = n, ncol = length(Anodes)))
    Ha = data.frame(matrix(rep(0,(length(Anodes)*n)), nrow = n, ncol = length(Anodes)))
    # Habar is the intervened upon clever covariate for that conditional mean
    Habar = list()
    # prob of observed treatment at time t
    gaw = list()
    # prob of receiving intervened-upon treatment up to time t
    gabar = list()
    # observed clever covariates at time t
    HAW = list()
    
    for (t in 1:endtime) {
      # grab the IC information for time t for pscore regression at time t
      ICg[[t]] = IC.beta(data,OC=NULL, Anodes[t], Qform = formulas_g[[t]], 
                         verbose = FALSE, parallelize = FALSE)
      
      # get values where the outcome existed for treatment
      reals = !ICg[[t]]$cens
      # set up prediction to be for A=1 or A=0 depending on setA
      gaw[[t]] = rep(NA,n)
      if (setA[t]==1){
        gaw[[t]][reals] = predict(ICg[[t]]$fit, type = 'response') 
      } else {
        gaw[[t]][reals] = (1- predict(ICg[[t]]$fit, type = 'response')) 
      }
      
      A = data[reals,Anodes[t]]
      A_index = grep(Anodes[t], colnames(data))
      Lta = data[,1:A_index]
      for (a in 1:t) Lta[,Anodes[a]] = setA[a]
      Lta = model.matrix(formulas_g[[t]], Lta)
      Lta = Lta[,ICg[[t]]$goods]
      gabar[[t]] = rep(NA,n)
      if (setA[t]==1) {
        gabar[[t]][reals] = plogis(Lta %*% ICg[[t]]$fit$coef)
      } else {
        gabar[[t]][reals] = 1 - plogis(Lta %*% ICg[[t]]$fit$coef)
      }
      
      # the clever covariate piece at time t
      H[reals,t] = I(A==setA[t])/gaw[[t]][reals]
      HAW[[t]] = rowProds(as.matrix(H[,1:t]))
      
      # the intervened upon clever cov
      Ha[reals,t] = 1/gabar[[t]][reals]
      Habar[[t]] = rowProds(as.matrix(Ha[,1:t]))
      
      # add the clever cov to the formula
      temp = as.character(formulas[[t]])
      temp = paste(temp[c(2, 1, 3:length(temp))],"",collapse = "")
      Hname = paste0("H",t)
      temp = paste0(temp, " + ", Hname)
      formulas[[t]] = formula(temp)
      # slide it in the data.frame just before the outcome
      slide = grep(Anodes[t], colnames(data))
      dd = ncol(data)
      # take the rowProds to get the pscore and indicator (the clever cov)
      labels = colnames(data)
      data = cbind(data[,1:slide], HAW[[t]], data[,(slide+1):dd])
      colnames(data)[c(1:slide, (slide+2):ncol(data))] = labels
      colnames(data)[(slide+1)] = Hname
    }
  }
  # getting outcome indicies and then proceed as if estimated g is the truth
  Yinds = vapply(Ynodes, FUN = function(x) grep(x,colnames(data)), FUN.VALUE = 1)
  for (t in times) {
    # If t is the end time we just do a regression on the outcome
    if (t == max(times)) {
      IC_tplus1 = 0
      OC = NULL
      # if (t == 1) design = data else design = data[,-Yinds[1:(t-1)]]
      ICinfo_t = IC.beta(data = data, OC = OC, Ynode = Ynodes[t],
                         Qform = formulas[[t]], parallelize = parallelize)
      IC_t = ICinfo_t$IC_beta
      
      if (tmle) {
        QAW = predict(ICinfo_t$fit, type = 'response')
        # keep real outcome indices
        keeps = !ICinfo_t$cens
        
        X_t = ICinfo_t$X[keeps,]
        hess_qinv = ICinfo_t$hessian
        # take gradient of beta X_g(t) wrt betas for g
        ###
        ###
        # get the beta_g coefficient
        H_ind = grep(paste0("H",t), names(ICinfo_t$fit$coef))
        beta_g = ICinfo_t$fit$coef[H_ind]
        # Get design matrices for all previous g fits--doesn't include 
        
        for (a in 1:t) {
          if (a==1) {
            L = ICg[[1]]$X[keeps,]
          } else {
            M = ICg[[a]]$X[keeps,]
            L = cbind(L,M)
          }
        }
        
        Y_t = data[keeps,Ynodes[t]]
        # add to this the little piece for the clever beta coefficient
        gpred_mat = do.call(cbind,gabar[1:t])[keeps,]
        if (t==1) {
          gpred_mat = matrix(gpred_mat, ncol = 1)
        } else {
          gpred_mat = lapply(1:t, FUN = function(col){
            no.covs = length(ICg[[t]]$IC_beta[,1])
            new = matrix(rep(NA,length(Y_t)*no.covs), ncol = no.covs)
            if (setA[col]==1) {
              new = rep(gpred_mat[,col] - 1, no.covs)
              new = matrix(new, ncol = no.covs)
            } else {
              new = rep(1 - gpred_mat[,col], no.covs)
              new = matrix(new, ncol = no.covs)
            }
            return(new)
          })
          
          gpred_mat = do.call(cbind, gpred_mat)
        }
        gpred_dotLg = vapply(1:nrow(gpred_mat), FUN = function(x){
          return(as.numeric(gpred_mat[x,]*L[x,]))
        }, FUN.VALUE = rep(1,ncol(L)))
        # gpred_dotLg
        
        # average the covariances for part 2 of the IC
        hess_g = mclapply(1:nrow(L),FUN = function(x) {
          mat = -beta_g*HAW[[t]][keeps][x]*(1-QAW[x])*
            QAW[x]*as.numeric(X_t[x,])%*%t(gpred_dotLg[,x])
          return(mat)
        })
        
        hess_g = Reduce('+', hess_g)/nrow(L)
        # taking average
        piece_g = gpred_dotLg %*% ((Y_t-QAW)*HAW[[t]][keeps])/nrow(L)
        # update hess_g to add this to the row for coefficient of H_t
        hess_gplus = hess_g
        hess_gplus[H_ind,] = hess_g[H_ind,] + piece_g
        # finally, we can multiply by the stacked IC
        if (t == 1) ICg_t = ICg[[1]]$IC_beta else {
          ICg_t = ICg[[1]]$IC_beta
          for (i in 2:t) ICg_t = rbind(ICg_t, ICg[[i]]$IC_beta)
        }
        IC_temp = (hess_qinv %*% hess_gplus) %*% ICg_t
        IC_t = IC_t + IC_temp
      }  
    } else {
      IC_tplus1 = IC_t
      ICinfo_tplus1 = ICinfo_t
      # get the indice to make the design matrix based on the formula
      Yind = grep(Ynodes[t+1], colnames(data))
      # for the outcomes at t which are 0 we use prev regression o.w use Y = 1
      # goods are those that did not die
      Y_t = data[,Ynodes[t]]
      goods = vapply(Y_t, FUN = function(x) {
        t = ifelse(!is.na(x), x==0, FALSE) 
      }, FUN.VALUE = TRUE)
      reals = vapply(Y_t, FUN = function(x) {
        ifelse(!is.na(x), x==1, FALSE) 
      }, FUN.VALUE = TRUE)
      # make the outcome predictions based on previous beta by grabbing previous 
      # design, intervening then setting up the design--should prob 
      # switch to datatable
      Xa_tplus1 = data[goods,-Yinds[1:t]]
      for (i in 1:(t+1)) {
        col = grep(Anodes[i], colnames(Xa_tplus1))
        Xa_tplus1[,col] = setA[i]
        if (tmle) Xa_tplus1[,paste0("H",(t+1))]=Habar[[(t+1)]][goods]
      }
      Xa_tplus1 = model.matrix(formulas[[t+1]],Xa_tplus1)
      Xa_tplus1 = Xa_tplus1[,(ICinfo_tplus1$goods)]
      OC = rep(NA,n)
      OC[goods] = plogis(Xa_tplus1 %*% ICinfo_tplus1$fit$coef)
      OC[reals] = 1
      # get the new beta
      # if (t == 1) design = data else design = data[,-Yinds[1:(t-1)]]
      ICinfo_t = IC.beta(data = data, OC = OC, Ynode = Ynodes[t], 
                         Qform = formulas[[t]],parallelize = parallelize)
      X_t = ICinfo_t$X[goods,]
      OCgoods = OC[goods]
      # This is to create the M matrix from the paper
      hess = mclapply(1:length(OCgoods),FUN = function(x) {
        mat = (1-OCgoods[x])*OCgoods[x]*as.numeric(X_t[x,])%*%t(as.numeric(Xa_tplus1[x,]))
        return(mat)
      }, mc.cores = getOption("mc.cores", parallel::detectCores()))
      M = Reduce('+', hess)/length(OCgoods)
      M = ICinfo_t$hessian %*% M 
      # form the IC
      IC_temp = apply(IC_tplus1,2,FUN = function(x) M%*%as.numeric(x))
      
      IC_t = ICinfo_t$IC_beta + IC_temp
      # We need to add extra pieces if we do tmle
      if (tmle) {
        QAW = predict(ICinfo_t$fit, type = 'response')
        keeps = !ICinfo_t$cens
        X_t = ICinfo_t$X[keeps,]
        hess_qinv = ICinfo_t$hessian
        
        # take gradient of beta X_g(t)^abar wrt betas for g
        ###
        ###
        # get the beta_g coefficient
        H_ind = grep(paste0("H",t), names(ICinfo_t$fit$coef))
        beta_g = ICinfo_t$fit$coef[H_ind]
        # Get intervened upon matrices for all previous g fits--doesn't include H
        
        for (a in 1:t) {
          if (a==1) {
            L = ICg[[1]]$X[keeps,]
          } else {
            M = ICg[[a]]$X[keeps,]
            L = cbind(L, M)
          }
        }
        
        Y_t = data[keeps,Ynodes[t]]
        # add to this the little piece for the clever beta coefficient\
        gpred_mat = do.call(cbind,gabar[1:t])[keeps,]
        if (t==1) gpred_mat = matrix(gpred_mat, ncol = 1)
        gpred_mat = lapply(1:t, FUN = function(col){
          no.covs = length(ICg[[t]]$IC_beta[,1])
          new = matrix(rep(NA,length(Y_t)*no.covs), ncol = no.covs)
          if (setA[col]==1) {
            new = rep(gpred_mat[,col] - 1, no.covs)
            new = matrix(new, ncol = no.covs)
          } else {
            new = rep(1 - gpred_mat[,col], no.covs)
            new = matrix(new, ncol = no.covs)
          }
          return(new)
        })
        
        gpred_mat = do.call(cbind, gpred_mat)
        gpred_dotLg = vapply(1:nrow(gpred_mat), FUN = function(x){
          return(as.numeric(gpred_mat[x,]*L[x,]))
        }, FUN.VALUE = rep(1,ncol(L)))
        # gpred_dotLg
        
        # average the covariances for part 2 of the IC
        hess_g = mclapply(1:nrow(L),FUN = function(x) {
          mat = -beta_g*HAW[[t]][keeps][x]*(1-QAW[x])*
            QAW[x]*as.numeric(X_t[x,])%*%t(gpred_dotLg[,x])
          return(mat)
        }, mc.cores = getOption("mc.cores", detectCores()))
        
        hess_g = Reduce('+', hess_g)/nrow(L)
        
        #####
        #####
        #####
        # not Y_t in there next line
        #####
        #####
        
        piece_g = gpred_dotLg %*% ((OC[reals|goods]-QAW)*HAW[[t]][keeps])/nrow(L)
        # update hess_g to add this to the row for coefficient of H_t
        hess_gplus = hess_g
        hess_gplus[H_ind,] = hess_g[H_ind,] + piece_g
        # finally, we can multiply by the stacked IC
        if (t == 1) ICg_t = ICg[[1]]$IC_beta else {
          ICg_t = ICg[[1]]$IC_beta
          for (i in 2:t) ICg_t = rbind(ICg_t, ICg[[i]]$IC_beta)
        }
        # ICg_t = ICg[[1]]$IC_beta
        # 
        # if (t > 1) {
        #   for (a in 2:t) ICg_t = rbind(ICg_t, ICg[[a]]$IC_beta)
        # }
        IC_temp = (hess_qinv %*% hess_gplus) %*% ICg_t
        
        # We have to go up one time point, grab the Lg and gpreds
        H_indtplus1 = grep(paste0("H",t+1), names(ICinfo_tplus1$fit$coef))
        beta_gtplus1 = ICinfo_tplus1$fit$coef[H_indtplus1]
        keeps1 = !is.na(data[,Yinds[t+1]])
        L1 = ICg[[t+1]]$X[keeps1,]
        L0 = as.data.frame(matrix(rep(NA, ncol(L)*n), ncol = ncol(L)))
        L0[keeps,] = L
        L1 = cbind(L0[keeps1,],L1)
        
        gpred1 = gabar[[t+1]][keeps1]
        no.covs = length(ICg[[t+1]]$IC_beta[,1])
        new = matrix(rep(NA,length(OCgoods)*no.covs), ncol = no.covs)
        if (setA[t+1]==1) {
          new = rep(gpred1 - 1, no.covs)
          new = matrix(new, ncol = no.covs)
        } else {
          new = rep(1 - gpred1, no.covs)
          new = matrix(new, ncol = no.covs)
        }
        
        gtemp = matrix(rep(NA,ncol(gpred_mat)*n), ncol = ncol(gpred_mat))
        gtemp[keeps] = gpred_mat
        gpred_mat1 = cbind(gtemp[keeps1,], new)
        
        gpred_dotLg1 = vapply(1:nrow(gpred_mat1), FUN = function(x){
          return(as.numeric(gpred_mat1[x,]*L1[x,]))
        }, FUN.VALUE = rep(1,ncol(L1)))
        # gpred_dotLg
        
        # average the covariances for part 2 of the IC
        hess_g1 = mclapply(1:nrow(L1),FUN = function(x) {
          mat = beta_gtplus1*HAW[[t+1]][keeps1][x]*(1-OCgoods[x])*
            OCgoods[x]*as.numeric(X_t[x,])%*%t(gpred_dotLg1[,x])
          return(mat)
        }, mc.cores = getOption("mc.cores", detectCores()))
        
        hess_g1 = Reduce('+', hess_g1)/nrow(L1)
        
        ICg1_t = rbind(ICg_t, ICg[[t+1]]$IC_beta)
        
        IC_temp1 = (hess_qinv %*% hess_g1) %*% ICg1_t
        
        IC_t = IC_t + IC_temp + IC_temp1
      }
    }
  }
  # if tmle we need to add the g pieces of the IC
  
  n = nrow(data)
  # if (t == 1) XA = data else XA = data[,-Yinds[1:(t-1)]]
  XA = data
  XA[,Anodes[1]] = setA[1]
  XA = model.matrix(formulas[[1]],XA)
  XA = XA[,(ICinfo_t$goods)]
  if (tmle) {
    if ("H1" %in% names(ICinfo_t$fit$coef)) XA[,"H1"] = Habar[[1]]
  }
  QAk = plogis(XA %*% ICinfo_t$fit$coef)
  
  # score = rowMeans(sapply(1:n,FUN = function(x) {
  #   QAk[x]*(1 - QAk[x])*XA[x,]
  # })) 
  score = t(XA) %*% (QAk*(1-QAk))/nrow(XA)
  psi = mean(QAk)
  if (tmle) {
    score_g = gpred_dotLg %*% (QAk*(1 - QAk)*Habar[[t]][keeps])*beta_g/length(Y_t)
    ICextra = t(ICg_t) %*% score_g
    IC = apply(IC_t, 2, FUN = function(x) sum(score*x)) + QAk - psi + ICextra
  } else {
    IC = apply(IC_t, 2, FUN = function(x) sum(score*x)) + QAk - psi
  }
  SE = sd(IC)*sqrt(n-1)/n
  
  qq = qnorm(1-alpha/2)
  CI = c(psi = psi, left = psi - qq*SE, right = psi + qq*SE)
  return(list(CI = CI, IC = IC))
}


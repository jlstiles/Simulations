
#' @title get.dgp
#' @description randomly creates a dgp and attempts to satisfy user specs. Number of covariates
#' is not limited but may take a while beyond d = 20 with too many terms. Limit time with depth
#' and maxterms parameters.
#' @param n, sample size
#' @param d, dimension of potential confounders
#' @param pos, a small value to make sure prop scores are in (pos, 1 - pos)
#' @param minATE, minimum causal risk difference for the population.  
#' @param minBV, minimum blip variance for population
#' @param depth, specify depth of interaction--must be less than or equal d.  
#' @param maxterms, maximum terms per interaction.  For example, this would limit 
#' two way interactions to maximally 10 terms as well as three way or main terms.
#' With high dimension it is wise to set this low because it might take a while
#' otherwise.  Still in development--perhaps future will set this for each depth
#' @param minterms sets a minimum number of total covariate terms, including 
#' interactions with eachother--do not set lower than 1.
#' @param mininters sets the minimum number of interactions with treatment to include
#' This must be bigger or equal to minterms
#' @param num.binaries specifies number of main terms you want as binaries, must be 
#' less than d.
#' @param force.confounding forces variables used for p-score to overlap with those
#' used for outcome regression. 
#' @param skewing randomly skews an otherwise centered dgp for generating binary treatment
#' default is c(-1, 1).  Set to c(-5,-1) to deliberately skew more regularly or widen to 
#' c(-3, 3) to skew more randomly.
#' @return  a sample DF, the true average treatment effect, ATE0 and blip variance
#' BV0, the sample pscores, PGn, the sample true blips, blip_n, the sample 
#' true prob of death under treatment, PQ1n, and prob of death under control
#' PQ0n
#' @export
#' @example /inst/examples/example_get.dgp.R
remakeDGP = function(n, object, limit_inter = NULL) 
{
  f_Aforms = object$f_Aforms
  f_Yforms = object$f_Yforms
  terms = object$terms
  skewage = object$skewage
  termsQW = object$termsQW
  terms_inter = object$terms_inter
  coef_Q = object$coef_Q
  U_Wdist = object$U_Wdist
  d = object$d
  num.binaries = object$num.binaries
  depth = object$depth
  coef_G = object$coef_G
  coef_Q = object$coef_Q
  pos = object$pos
  # sample size of population

  # the population matrix of potential confounders, consisting of normals and binaries
  types = list(function(x) sin(x), function(x) cos(x), 
               function(x) x^2, function(x) x, function(x) x^3, function(x) exp(x))
  
  no.types = length(types)
  
   # Create the W matrix
  Wmat = lapply(1:d, FUN = function(i) {
      W = effect(n, dist = U_Wdist[[i]]$dist, params = U_Wdist[[i]]$params)
      return(W)})
  
  U_W = do.call(cbind, Wmat)
  
  f_A = lapply(1:length(f_Aforms), FUN = function(x) {
    if (x<=num.binaries) return(Wmat[[x]]) else {
      return(f_Aforms[[x]](Wmat[[x]]))
    }})
                    
  f_A = do.call(cbind, f_A)
  
  f_Y = lapply(1:length(f_Aforms), FUN = function(x) {
    if (x<=num.binaries) return(Wmat[[x]]) else {
      return(f_Yforms[[x]](Wmat[[x]]))
    }})
  
  f_Y = do.call(cbind, f_Y)
  
  # All of the interaction combos of columns possible up to the depth of interaction
  # user specifies
  choos = lapply(1:depth, FUN = function(x) {
    c = combn(1:d, x)
    if (!is.matrix(c)) c = as.matrix(c)
    return(c)
  })
  # combine specified columns as to randomly chosen interactions
  col.comb = lapply(1:length(terms), FUN = function(a) {
    col.choos = terms[[a]]
    if (length(col.choos) == 0) {
      return(integer(0))
    } else {
      df = vapply(col.choos, FUN = function(x) {
        col.inds = choos[[a]][,x]
        v = rep(1, n)
        for (c in col.inds) v = v*f_A[,c]
        return(v)
      }, FUN.VALUE = rep(1, n))
      return(df)
    }
  })
  
  # put the cols in one matrix used for p-score
  dfG = do.call(cbind, col.comb)
  
  # transform the columns by plugging into randomly drawn functions and standardize
  # so no variable dominates unnecessarily
  dfG = apply(dfG, 2, FUN = function(col) {
    if (all(col == 1 | col ==0)) {
      v = col
      return(v)
    } else {
      v = (col - mean(col))/sd(col)
      return(v)
    }
  })
  
  # create an intercept for skewing deliberately
  dfG = cbind(dfG, rep(1, n))
  PG = plogis(dfG %*% c(coef_G, skewage))
  
  # Creating A  based on p-scores for whole population of 1e6
  PG = pmin(pmax(PG, pos), 1-pos)
  # hist(PG0, breaks = 100)
  A = rbinom(n, 1, PG)

  # combine W interactions and mains for OC
  col.combQ = lapply(1:length(termsQW), FUN = function(a) {
    col.choos = termsQW[[a]]
    if (length(col.choos) == 0) {
      return(integer(0))
    } else {
      df = vapply(col.choos, FUN = function(x) {
        col.inds = choos[[a]][,x]
        v = rep(1, n)
        for (c in col.inds) v = v*f_Y[,c]
        return(v)
      }, FUN.VALUE = rep(1, n))
      return(df)
    }
  })

  
  # combine cols used for interaction with A
  col.comb_inter = lapply(1:length(terms_inter), FUN = function(a) {
    col.choos = terms_inter[[a]]
    if (length(col.choos) == 0) {
      return(integer(0))
    } else {
      df = vapply(col.choos, FUN = function(x) {
        col.inds = choos[[a]][,x]
        v = rep(1, n)
        for (c in col.inds) v = v*f_Y[,c]
        return(v)
      }, FUN.VALUE = rep(1, n))
      return(df)
    }
  })
  
  # put the cols in one matrix for W interactions and mains
  dfQWA = do.call(cbind, col.combQ)
  dfQWA = cbind(dfQWA, A)
  
  # put the cols in one matrix for interactions with A = 1 
  dfQ_inter = do.call(cbind, col.comb_inter)
  # and for population
  dfQ_interA = apply(dfQ_inter, 2, FUN = function(col) A*col)
  
  # OC df cols for W interactions and A plugged into randomly drawn functions (types)
  dfQWA = apply(dfQWA, 2, FUN = function(col) {
    if (all(col == 1 | col ==0)) {
      return(col)
    } else {
      # v = types[[sample(1:no.types, 1)]](col)
      v = (col - mean(col))/sd(col)
      return(v)
    }
  })
  
  # This skips interactions appendages to dfQWA if no interactions are chosen
  no.inters = sum(unlist(lapply(terms_inter, sum)))
  
  if (no.inters != 0) {
    # apply the fcns as per setting A to its draw, A =1 and A=0
  
    dfQ_inter0 = vapply(1:ncol(dfQ_interA), FUN = function(col) rep(0,n), FUN.VALUE = rep(1,n))
    # We standardize these columns too for the population as is
    means = apply(dfQ_interA, 2, FUN = function(col) mean(col))
    sds = apply(dfQ_interA, 2, FUN = function(col) sd(col))
    dfQ_interA = apply(dfQ_interA, 2, FUN = function(col) (col - mean(col))/sd(col))
    
    # We apply this to the pop under A =1 and A = 0 so we apply the same fcn for these as for 
    # the true observed population as is
    dfQ_inter = vapply(1:ncol(dfQ_inter), FUN = function(col) {
      (dfQ_inter[,col] - means[col])/sds[col]
    }, FUN.VALUE = rep(1,n))
    dfQ_inter0 = vapply(1:ncol(dfQ_inter0), FUN = function(col) {
      (dfQ_inter0[,col] - means[col])/sds[col]
    }, FUN.VALUE = rep(1,n))
    
    # standardize the treatment column too, to be fair!!
    dfQ = cbind(dfQWA, dfQ_interA)
    dfQW1 = dfQWA
    dfQW1[, ncol(dfQW1)] = (1 - mean(A))/sd(A)
    dfQ1 = cbind(dfQW1, dfQ_inter)
    
    dfQW0 = dfQWA
    dfQW0[, ncol(dfQW0)] = - mean(A)/sd(A)
    dfQ0 = cbind(dfQW0, dfQ_inter0)
    
    # else we have no interactions with A shenanigans and everything is very simple
  } else {
    dfQ = cbind(dfQWA)
    dfQW1 = dfQWA
    dfQW1[, ncol(dfQW1)] = (1 - mean(A))/sd(A)
    dfQ1 = cbind(dfQW1)
    dfQW0 = dfQWA
    dfQW0[, ncol(dfQW0)] = - mean(A)/sd(A)
    dfQ0 = cbind(dfQW0)
  }
  
  # compute true probs under A = 1 and A = 0 and related truths
  PQ1 = plogis(dfQ1 %*% coef_Q)
  PQ0 = plogis(dfQ0 %*% coef_Q)
  blip_true = PQ1 - PQ0
  ATE0 = mean(blip_true)
  BV0 = var(blip_true)
  
  # finally we create the population probs of death
  PQ = plogis(dfQ %*% coef_Q)

  PQ = pmin(pmax(PQ, .00001), 1-.00001)
  # hist(PQ)
  # take the draw for the population
  Y = rbinom(n, 1, PQ)
  # make sure our loglikelihood loss is bounded reasonably, no one gets super lucky or unlucky!
  # mean(Y*A/PG0-Y*(1-A)/(1-PG0))
  # ATE0
  # take a sample of size n and return sample probs blips, the dataframe 
  # with covariates but the user never sees the formula. Now they can use DF
  # to try and recover the truth
  blip_n = PQ1 - PQ0
  An = A
  Yn = Y
  Wn = U_W
  DF = cbind(Wn, An, Yn)
  colnames(DF)[c((d + 1), (d + 2))] = c("A", "Y")
  colnames(DF)[1:d] = paste0("W",1:d)
  return(list(DF = DF, blip_n = blip_n, Wn = Wn,
              PQ1n = PQ1, PQ0n = PQ0, PQn = PQ, PGn = PG))
}


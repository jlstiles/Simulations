
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
#' @param forcedist Can be used for specifying the prob of drawing a dist in Wdist
#' @param Wdist These are dist for generating W's  
#' @return  a sample DF, the true average treatment effect, ATE0 and blip variance
#' BV0, the sample pscores, PGn, the sample true blips, blip_n, the sample 
#' true prob of death under treatment, PQ1n, and prob of death under control
#' PQ0n
#' @export
#' @example /inst/examples/example_get.dgp.R
get.dgp = function(n, d, pos = 0.01, minATE = -2, minBV = 0, depth, maxterms, minterms, 
                   mininters, num.binaries = floor(d/4), skewing = c(-1,1), 
                   force.confounding = TRUE, limit_inter = NULL, forcedist = c(.2,.2,.2,.2,.2),
                   Wdist = list(normal = list(dist = rnorm, params = list(mean = 0, sd = 1)),
                                beta = list(dist = rbeta, params = list(shape1= 2, shape2 = 1)),
                                beta = list(dist = rbeta, params = list(shape1= 10, shape2 = 11)),
                                chisquared = list(dist = rchisq, params = list(df = 1, ncp = 5)),
                                uniform = list(dist = runif, params = list(min = -1, max = 1)))
) 
{
  # n = 1000; d = 4; pos = .01; minATE = -2; minBV = .03; depth = 2; maxterms = 2; minterms = 1; mininters = 1
  # num.binaries = 2; force.confounding = TRUE
  if (minterms == 0) 
    stop("minterms must be atleast 1")
  if (mininters > minterms) 
    stop("minimum number of interactions cannot exceed number of covariate terms, obviously, hello!!!")
  
  # sample size of population
  N = 1e+06
  
  # hist(rbeta(1000, 10,11))
  # mean(rbeta(1000,10,11))
  # 
  # hist(rchisq(1000, 1,5), 100)
  # mean(rchisq(1000,1,5))
  
  # randomly drawn binary distributions on random columns, cols don't need to be random here
  r = runif(num.binaries, 0.3, 0.7)
  Wdist_bin = lapply(r, FUN = function(x) list(dist = rbinom, params = list(size = 1, prob = x)))
  # the population matrix of potential confounders, consisting of normals and binaries
  types = list(function(x) sin(x), function(x) cos(x), 
               function(x) x^2, function(x) x, function(x) x^3, function(x) exp(x))
  
  no.types = length(types)
  
   # Create the W matrix
  Wmat = lapply(1:d, FUN = function(col) {
    if (col <= num.binaries) {
      W = rbinom(N, 1, r[col])
      i = col
      return(list(W, dist=i))
    } else {
      i = which(rmultinom(1,1,forcedist)==1)
      W = effect(N, dist = Wdist[[i]]$dist, params = Wdist[[i]]$params)
      return(list(W, dist=i))
    }
  })
  
  U_W = do.call(cbind, lapply(Wmat, FUN = function(x) x[[1]]))
  U_Wdist = append(Wdist_bin, lapply(Wmat[(num.binaries+1):d], FUN = function(x) Wdist[x[[2]]][[1]]))
  
  # select the fcn of W for g
  Wmat_g = lapply(1:d, FUN = function(col) {
    if (col <= num.binaries) {
      return(list(U_W[,col], "none-binary"))
    } else {
      fcn = types[[sample(1:no.types, 1)]]
      v = fcn(U_W[,col])
      return(list(v, fcn))
    }
  })

  f_A = do.call(cbind, lapply(Wmat_g, FUN = function(x) x[[1]]))
  f_Aforms = lapply(Wmat_g, FUN = function(x) x[[2]])
  
  Wmat_Q = lapply(1:d, FUN = function(col) {
    if (col <= num.binaries) {
      return(list(U_W[,col], "none-binary"))
    } else {
      if (force.confounding) {
        fcn = f_Aforms[[col]]
        v = fcn(U_W[,col])
        return(list(v, fcn))} else {
        fcn = types[[sample(1:no.types, 1)]]
        v = fcn(U_W[,col])
        return(list(v, fcn))
      }
    }
  })
  
  f_Y = do.call(cbind, lapply(Wmat_Q, FUN = function(x) x[[1]]))
  f_Yforms = lapply(Wmat_Q, FUN = function(x) x[[2]])
  
  # All of the interaction combos of columns possible up to the depth of interaction
  # user specifies
  choos = lapply(1:depth, FUN = function(x) {
    c = combn(1:d, x)
    if (!is.matrix(c)) c = as.matrix(c)
    return(c)
  })
  
  # types of transformations to apply, can add many more

  ##
  # begin p-score construction for population
  ##
  
  # select the interaction terms to be included according to maxterms and minterms
  # maxterms is maxterms for any depth, so for if maxterms is 5, you cannot have more 
  # than 5 2 way interactions or more than 5 8 way interactions, etc. and minterms makes 
  # sure we have a model of certain complexity minterms must be greater than equal 1
  s = -1
  
  while (s < minterms) {
    no.terms = sample(0:min(maxterms, ncol(choos[[1]])),1)
    select.cols = sample(1:ncol(choos[[1]]), no.terms)
    terms = list()
    terms[[1]] = select.cols
    terms = append(terms, lapply(choos[2:length(choos)], FUN = function(x) {
      if (!is.null(limit_inter)) L = min(limit_inter, ncol(x)) else L = ncol(x)
      no.terms = sample(0:min(maxterms, L),1)
      select.cols = sample(1:ncol(x), no.terms)
      return(select.cols)
    }))
    s = sum(unlist(lapply(terms, FUN = function(x) length(x))))
  }
  # combine specified columns as to randomly chosen interactions
  col.comb = lapply(1:length(terms), FUN = function(a) {
    col.choos = terms[[a]]
    if (length(col.choos) == 0) {
      return(integer(0))
    } else {
      df = vapply(col.choos, FUN = function(x) {
        col.inds = choos[[a]][,x]
        v = rep(1, N)
        for (c in col.inds) v = v*f_A[,c]
        return(v)
      }, FUN.VALUE = rep(1, N))
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
  
  skewage = runif(1, skewing[1], skewing[2])  
  dfG = cbind(dfG, rep(1, N))
  coef_G = c(runif(ncol(dfG)-1, -1, 1))
  
  # satisfying positivity constraints, we don't want too high a percentage beyond
  # the user specified positivity probs of pos and 1-pos. Since .8^20 is small we
  # can usually satisfy this constraint but remains to be seen for ridiculous scenarios
  tol = TRUE
  its = 0
  while (tol & its < 20) {
    PG0 = plogis(dfG %*% c(coef_G, skewage))
    coef_G = 0.8 * coef_G
    its = its + 1
    tol = mean(PG0 < pos) > 0.01 | mean(PG0 > (1 - pos)) > 
      0.01
  }
  
  # Creating A  based on p-scores for whole population of 1e6
  PG0 = pmin(pmax(PG0, pos), 1-pos)
  # hist(PG0, breaks = 100)
  A = rbinom(N, 1, PG0)
  
  ###
  # NOW WE TAKE CARE OF TRUE OC Probs
  ###
  
  if (force.confounding) {
    termsQW = terms
    s = -1 
    while (s < mininters) {
      terms_inter = lapply(terms, FUN = function(x) {
        no.terms = sample(0:length(x),1)
        select.cols = sample(x, no.terms)
        return(select.cols)
      })
      if (mininters == 0) s = Inf else s = sum(unlist(lapply(terms_inter, sum)))
    }
  } else {
    # for barQ first select interaction terms for just the W's
    s = -1
    while (s < minterms) {
      termsQW = lapply(choos, FUN = function(x) {
        no.terms = sample(0:min(maxterms, ncol(x)),1)
        select.cols = sample(1:ncol(x), no.terms)
        return(select.cols)
      })
      s = sum(unlist(lapply(termsQW, sum)))
    }
    s = -1 
    while (s < mininters) {
      terms_inter = lapply(choos, FUN = function(x) {
        no.terms = sample(0:min(maxterms, ncol(x)),1)
        select.cols = sample(1:ncol(x), no.terms)
        return(select.cols)
      })
      if (mininters == 0) s = Inf else s = sum(unlist(lapply(terms_inter, sum)))
    }
  }
  
  # if (force.confounding) use same terms in OC interactions
  # as for p-scores or a high percentage there of
  # for barQ select interaction terms of W's that will interact with A

  
  # combine W interactions and mains for OC
  col.combQ = lapply(1:length(termsQW), FUN = function(a) {
    col.choos = termsQW[[a]]
    if (length(col.choos) == 0) {
      return(integer(0))
    } else {
      df = vapply(col.choos, FUN = function(x) {
        col.inds = choos[[a]][,x]
        v = rep(1, N)
        for (c in col.inds) v = v*f_Y[,c]
        return(v)
      }, FUN.VALUE = rep(1, N))
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
        v = rep(1, N)
        for (c in col.inds) v = v*f_Y[,c]
        return(v)
      }, FUN.VALUE = rep(1, N))
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

    dfQ_inter0 = vapply(1:ncol(dfQ_interA), FUN = function(col) rep(0,N), FUN.VALUE = rep(1,N))
    
    # We standardize these columns too for the population as is
    means = apply(dfQ_interA, 2, FUN = function(col) mean(col))
    sds = apply(dfQ_interA, 2, FUN = function(col) sd(col))
    dfQ_interA = apply(dfQ_interA, 2, FUN = function(col) (col - mean(col))/sd(col))
    
    # We apply this to the pop under A =1 and A = 0 so we apply the same fcn for these as for 
    # the true observed population as is
    dfQ_inter = vapply(1:ncol(dfQ_inter), FUN = function(col) {
      (dfQ_inter[,col] - means[col])/sds[col]
    }, FUN.VALUE = rep(1,N))
    dfQ_inter0 = vapply(1:ncol(dfQ_inter0), FUN = function(col) {
      (dfQ_inter0[,col] - means[col])/sds[col]
    }, FUN.VALUE = rep(1,N))
    
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
  
  # treatment node position for convenience
  TXpos = ncol(dfQWA)  
  # draw an extent of coefficient variation
  a = runif(1, 0, 1)
  # use that extent to select coefs for all the OC columns
  coef_Q = runif(ncol(dfQ), -a, a)
  # compute true probs under A = 1 and A = 0 and related truths
  PQ1 = plogis(dfQ1 %*% coef_Q)
  PQ0 = plogis(dfQ0 %*% coef_Q)
  blip_true = PQ1 - PQ0
  ATE0 = mean(blip_true)
  BV0 = var(blip_true)
  
  # we then tweak coefs to satisfy the user specs on true ATE
  jj = 1
  while (abs(ATE0) <= minATE & jj <= 20) {
    coef_Q[TXpos] = 1.2*coef_Q[TXpos]
    PQ1 = plogis(dfQ1 %*% coef_Q)
    PQ0 = plogis(dfQ0 %*% coef_Q)
    blip_true = PQ1 - PQ0
    ATE0 = mean(blip_true)
    jj = jj + 1
  } 
  
  # we then tweak coefs to satisfy the user specs on true blip var, BV0
  jj = 1 
  if (no.inters != 0) {
    while (BV0 <= minBV & jj <= 20) {
      coef_Q[(ncol(dfQWA) + 1):ncol(dfQ)] = 1.2 * coef_Q[(ncol(dfQWA) + 1):ncol(dfQ)]
      PQ1 = plogis(dfQ1 %*% coef_Q)
      PQ0 = plogis(dfQ0 %*% coef_Q)
      blip_true = PQ1 - PQ0
      BV0 = var(blip_true)
      jj = jj + 1
    } 
  } 
  
  ATE0 = mean(blip_true)
  # finally we create the population probs of death
  PQ = plogis(dfQ %*% coef_Q)
  # max(PQ)
  # min(PQ)
  # 
  PQ = pmin(pmax(PQ, .00001), 1-.00001)
  
  # hist(PQ)
  # take the draw for the population
  Y = rbinom(N, 1, PQ)
  # make sure our loglikelihood loss is bounded reasonably, no one gets super lucky or unlucky!
  # mean(Y*A/PG0-Y*(1-A)/(1-PG0))
  # ATE0
  # take a sample of size n and return sample probs blips, the dataframe 
  # with covariates but the user never sees the formula. Now they can use DF
  # to try and recover the truth
  S = sample(1:N, n)
  PQ1n = PQ1[S]
  PQ0n = PQ0[S]
  PQn = PQ[S]
  PGn = PG0[S]
  blip_n = PQ1n - PQ0n
  An = A[S]
  Yn = Y[S]
  Wn = as.data.frame(U_W[S, ])
  DF = cbind(Wn, An, Yn)
  colnames(DF)[c((d + 1), (d + 2))] = c("A", "Y")
  colnames(DF)[1:d] = paste0("W",1:d)
  return(list(BV0 = BV0, ATE0 = ATE0, DF = DF, blip_n = blip_n, 
              PQ1n = PQ1n, PQ0n = PQ0n, PQn = PQn, PGn = PGn, U_Wdist = U_Wdist,
              f_Aforms = f_Aforms, f_Yforms = f_Yforms, terms = terms, skewage = skewage, its = its,
            termsQW = termsQW, terms_inter = terms_inter, coef_Q = coef_Q, coef_G = coef_G,
            d = d, num.binaries = num.binaries, depth = depth, pos = pos))
}


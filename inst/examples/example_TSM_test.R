
# using built-in package functions, g0_linear and define Q0_linear to specify
# pscore and outcome model probabilities
g0_linear
Q0_linear = function(A,W1,W2,W3,W4) plogis(A + W1 + W2 + A*(W3 + W4)+W3+W4)
big = gendata(1e6,g0_linear, Q0_linear)
setA = 1
truth1 = mean(with(big, Q0_linear(setA, W1,W2,W3,W4)))
truth1
truth0 = mean(with(big, Q0_linear(0, W1,W2,W3,W4)))
truth0
# get a randomly drawn dataframe under the specified model
simmie = function(n, truth) {
  # n=1000
  data = gendata(n, g0_linear, Q0_linear)
  # well-specified model
  
  formulas = list(formula("Y ~ W1 + W2 + A*(W3 + W4)"))
  Ynodes = c("Y")
  Anodes = c("A")
  setA = 1
  
  TSMinfo = long.TSM(data = data, Ynodes = Ynodes, Anodes = Anodes, 
                     formulas = formulas, setA = setA, alpha = .05)
  
  cover = (TSMinfo$CI[2] < truth) & (TSMinfo$CI[3] > truth)
  # cover
  # info$CI
  
  return(cover)
}

detectCores()
cl = makeCluster(detectCores(), type = "SOCK")
registerDoSNOW(cl)

B = 100
ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
            .errorhandling = "remove")%dopar%
            {simmie(1000, truth1)}
results =  mean(unlist(ALL))
results

# should cover each truth 95 percent of the time and both truths
# simultaneously 95 percent of the time for the simultaneous CI's


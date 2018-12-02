
# using built-in package functions, g0_linear and define Q0_linear to specify
# pscore and outcome model probabilities
g0_linear
Q0_linear = function(A,W1,W2,W3,W4) plogis(A + W1 + W2 + A*(W3 + W4)+W3+W4)
big = gendata(1e6,g0_linear, Q0_linear)
setA = 1
truth = mean(with(big, Q0_linear(setA, W1,W2,W3,W4)))
truth = get.truth(g0_linear, Q0_linear)
truth
# get a randomly drawn dataframe under the specified model
simmie = function(n, truth) {
  n=1000
  data = gendata(n, g0_linear, Q0_linear)
  # well-specified model
  Qform = formula("Y ~ W1 + W2 + A*(W3 + W4)")
  
  # specifying the covariates, treatment and outcome
  W = data[,2:5]
  A = data$A
  Y = data$Y
  
  # should cover each truth 95 percent of the time.
  info = LR.inference(W=W,A=A,Y=Y,Qform=Qform, alpha = .05)
  
  cover = (info[2] < truth[1]) & (info[3] > truth[1])
  cover1 = (info[5] < truth[2]) & (info[6] > truth[2])
  return(c(cover, cover1))
}

detectCores()
cl = makeCluster(detectCores(), type = "SOCK")
registerDoSNOW(cl)

B = 100
ALL=foreach(i=1:B,.packages=c("gentmle2","mvtnorm","hal","Simulations","SuperLearner"),
            .errorhandling = "remove")%dopar%
            {simmie(1000, truth)}
results =  do.call(rbind, ALL)
results[1:10,]

mean(results[,1])
mean(results[,2])

# should cover each truth 95 percent of the time and both truths
# simultaneously 95 percent of the time for the simultaneous CI's


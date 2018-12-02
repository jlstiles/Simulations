
# using built-in package functions, g0_linear and define Q0_linear to specify
# pscore and outcome model probabilities
g0_linear
Q0_linear = function(A,W1,W2,W3,W4) plogis(A + W1 + W2 + A*(W3 + W4) + W3 + W4)

# get a randomly drawn dataframe under the specified model
data = gendata(1000, g0_linear, Q0_linear)

# get the truth
truth = get.truth(g0_linear, Q0_linear)
truth

# well-specified model
Qform = formula("Y ~ W1 + W2 + A*(W3 + W4)")

# specifying the covariates, treatment and outcome
W = data[,2:5]
A = data$A
Y = data$Y

# should cover each truth 95 percent of the time.
info = LR.inference(W=W,A=A,Y=Y,Qform=Qform, alpha = .05)
info

# should cover each truth 95 percent of the time and both truths
# simultaneously 95 percent of the time for the simultaneous CI's
info1 = LR.inference(W=W,A=A,Y=Y,Qform=Qform, alpha = .05, 
                    simultaneous.inference = TRUE)
info1


# basic example with very simple SuperLearner library
# generate the data according to a couple of built-in functions here

data = gendata(1000, g0 = g0_linear, Q0 = Q0_trig1)

#drop Y
Y = data$Y
A = data$A
#drop Y to form X and form newdata
X = data[,-6]
X0 = X1 = X
X0$A = 0
X1$A = 1
newdata = rbind(X,X1,X0)
#form W
W = X[,-1]

#declare SL library
SL.library = SL.libraryG = c("SL.glm","SL.mean")

stack = SL.stack(Y, X, A, W, newdata, method = "method.NNloglik",
                 SL.library, SL.libraryG, V=10, mc.cores = 4)

# simultaneously run one-step tmle for ATE and blip variance with 
# simultaneous CI
tmle.info = gentmle(initdata=stack$initdata, 
                    params=list(param_ATE,param_sigmaATE), 
                    submodel = submodel_logit, loss = loss_loglik,
                    approach = "recursive", max_iter = 10000, 
                    g.trunc = 1e-2,
                    simultaneous.inference = TRUE)
tmle.info$steps
# get simultaneous CIs
ci_gentmle(tmle.info)



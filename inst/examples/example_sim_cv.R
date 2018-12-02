# basic example with very simple SuperLearner library
SL.library = c("SL.glm", "SL.mean")
SL.libraryG = c("SL.glm", "SL.mean")
n=1000
g0 = g0_linear
Q0 = Q0_trig
gform = formula("A ~.")
Qform = formula("Y ~ A*(W1 + W2 + W3 + W4)")
estimator = c("single iterative")
result = sim_cv(n, g0, Q0, SL.library, SL.libraryG, method = "method.NNLS", 
                  cv = TRUE, V = 2, SL = 2L, gform, Qform, estimator, dgp = NULL)
          


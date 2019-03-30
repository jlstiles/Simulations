# getting a random sample from a randomly drawn dgp.  We specify 3 covariates, a minimum population ATE of .1
# a minimum variance of blip of .03, up to 2 way interactions, limit positivity violations to less than 1% of the
# population having propensity scores below .05 or above .95, 1 binary, up to 1 interaction terms and at least 1
# term included as a covariate in outcome regression and in treatment mechanism
undebug(get.dgp)
dgp = get.dgp1(n = 1000, d = 5, pos = 0.05, minATE = 0, minBV = 0, depth = 4, maxterms = 10, minterms = 1, 
              mininters = 0, num.binaries = 2, force.confounding = TRUE, limit_inter = 5) 

# population blip variance (VTE)
dgp$BV0
# population average treatment effect
dgp$ATE0
# sample proportion with treatment
mean(dgp$DF$A)
# sample proportion who died
mean(dgp$DF$Y)
# min sample pscore
min(dgp$PGn)
# max sample pscore
max(dgp$PGn)
# sample dataframe
head(dgp$DF)
# histogram of blips
hist(dgp$blip_n)
# histogram of propensity scores
hist(dgp$PGn,50)

dgp$
###
# The following example is for pt treatment.  
###
# using built-in package functions, g0_linear and define Q0_linear to specify
# pscore and outcome model probabilities
data(longdata)
head(data_pt)
g0_linear
Q0_linear = function(A,W1,W2,W3,W4) plogis(A + W1 + W2 + A*(W3 + W4) + W3 +W4)

# get the truth setting A to 1
setA = 1
truth = mean(with(gendata(1e6, g0_linear, Q0_linear), Q0_linear(A=setA,W1,W2,W3,W4)))
truth

# well-specified model
n=1000

formulas = list(formula("Y ~ A + W1 + W2 + A*(W3 + W4) + W3 +W4"))
Ynodes = c("Y")
Anodes = c("A")
setA = 1
formulas_g = list(formula("A ~ W1 + W2 + W3 + W4"))

TSMinfo = long.TSM(data = data_pt, Ynodes = Ynodes, Anodes = Anodes, 
                   formulas = formulas, formulas_g = formulas_g, tmle = TRUE, setA = setA, alpha = .05)

TSMinfo1 = long.TSM(data = data_pt, Ynodes = Ynodes, Anodes = Anodes, 
                    formulas = formulas, formulas_g = formulas_g, tmle = FALSE, setA = setA, alpha = .05)
# get CI
TSMinfo$CI
TSMinfo1$CI
# get influence curve
TSMinfo$IC

# TMLE recovers truth from misspecified outcome model by getting pscore right
# non-tmle does not

# misspecified outcome regression
formulas = list(formula("Y ~ A"))
Ynodes = c("Y")
Anodes = c("A")
setA = 1
# correctly specified g
formulas_g = list(formula("A ~ W1 + W2 + W3 + W4"))

TSMinfo = long.TSM(data = data_pt, Ynodes = Ynodes, Anodes = Anodes, 
                   formulas = formulas, formulas_g = formulas_g, tmle = TRUE, setA = setA, alpha = .05)

# non-tmle
TSMinfo1 = long.TSM(data = data_pt, Ynodes = Ynodes, Anodes = Anodes, 
                    formulas = formulas, formulas_g = formulas_g, tmle = FALSE, setA = setA, alpha = .05)
# get CI
TSMinfo$CI
TSMinfo1$CI

####
# example with longitudinal intervention
####
# enter data with time-ordering correct
head(data_long)
Ynodes = c("Y_1", "Y_2","Y_3","Y_4")
Anodes = c("A1_0", "A1_1","A1_2","A1_3")

# specify the formulas
formula0 = formula("Y_1 ~ L2_0 + L1_0 + A1_0") 
formula1 = formula("Y_2 ~ L2_1 + A1_1")
formula2 = formula("Y_3 ~ L2_2 + A1_2")
formulas = list(formula0, formula1, formula2)
formula_g0 = formula("A1_0 ~ L2_0 + L1_0") 
formula_g1 = formula("A1_1 ~ L2_1 + A1_0")
formula_g2 = formula("A1_2 ~ L2_2 + A1_1")
formulas_g = list(formula_g0, formula_g1, formula_g2)
setA = c(0,1,1)

# tmle (takes about 6 seconds on a mac with 4 cores)
time = proc.time()
TSMinfo = long.TSM(data = data_long, Ynodes = Ynodes, Anodes = Anodes, 
                   formulas = formulas, formulas_g = formulas_g, tmle = TRUE, setA = setA, alpha = .05,
                   parallel = TRUE)
proc.time() - time
# non-tmle
TSMinfo1 = long.TSM(data = data_long, Ynodes = Ynodes, Anodes = Anodes, 
                    formulas = formulas, formulas_g = formulas_g, tmle = FALSE, setA = setA, alpha = .05)
# get CI
TSMinfo$CI

# ignoring propensity score gives much smaller variance but gives small variance
# and high bias if confounding is present
TSMinfo1$CI


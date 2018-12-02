# getting a random sample from a randomly drawn dgp.  We specify 3 covariates, a minimum population ATE of .1
# a minimum variance of blip of .03, up to 2 way interactions, limit positivity violations to less than 1% of the
# population having propensity scores below .05 or above .95, 1 binary, up to 1 interaction terms and at least 1
# term included as a covariate in outcome regression and in treatment mechanism
dgp = get.dgp(n = 1000, d = 3, pos = 0.05, minATE = .1, minBV = .03, depth = 2, maxterms = 3, minterms = 1, 
              mininters = 1, num.binaries = 1) 

# population blip variance
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
hist(dgp$PGn)

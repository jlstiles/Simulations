% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/long.TSM.R
\name{long.TSM}
\alias{long.TSM}
\title{long.TSM}
\usage{
long.TSM(data, Ynodes, Anodes, formulas, setA, alpha = 0.05,
  formulas_g = NULL, tmle = FALSE, parallelize = FALSE)
}
\arguments{
\item{data, }{data.frame of variables in time ordering from left to right}

\item{Ynodes, }{character vector of time-ordered Ynodes}

\item{Anodes, }{character vector of time-ordered Anodes}

\item{formulas, }{list of formulas for the conditional means}

\item{setA}{the value to which you intervene on A,vector of length that of Anodes}

\item{alpha}{significance level for two-sided CI.}

\item{formulas_g}{NULL by default but must be filled in if tmle = TRUE}

\item{tmle}{option to put a clever covariate in the regression.  If the propensity
score is well-specified, will give good coverage asymptotically despite mispecifying
the outcome model.  Otherwise gives delta method inference for the "closest" linear
model fit using clever covariate in your regression.  Set to TRUE for RCT's.}

\item{parallelize}{FALSE by default (still in development)}
}
\value{
a list with elements, CI for the confidence interval and  IC for the
influence curve.
}
\description{
computes delta method inferencec for logistic regression plug-in
estimator of treatment specific mean for survival data in wide form,
including pt treatment. Also does logistic regression plug-in using clever
covariate in the regression, inference computed using delta method as well.
Note, in performing tmle, pscores are obtained via user supplied regression
formulas for each time point.  Formulas for the conditional means are also
user supplied.
}
\examples{
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

}

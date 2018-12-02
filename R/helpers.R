#' @export
get.truth = function(g0, Q0, N=1e6) {
  testdata=gendata(N, g0=g0, Q0 = Q0)
  blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
  propensity = with(testdata, g0(W1,W2,W3,W4))
  ATE0 = mean(blip_true)
  var0 = var(blip_true)
  return(c(ATE0=ATE0, var0=var0))
}

#' @export
cov.check = function(data, truth, ind) {
  ans = vapply(ind,FUN = function(x){
    covs = data[,x+1]<=truth&data[,x+2]>=truth
    mean(covs)
  }, FUN.VALUE = 1, USE.NAMES = FALSE)
  names(ans) = colnames(data)[ind]
  return(ans)
}

#' @export
cov.simul = function(data, truth, ind) {
  covs = data[,ind[1]+1]<=truth[1]&data[,ind[1]+2]>=truth[1]
  covs1 = data[,ind[2]+1]<=truth[2]&data[,ind[2]+2]>=truth[2]
  mean(covs*covs1)
}

# function to check MSE
#' @export
perf=function(ests,truth){
  n=length(ests)
  var=((n-1)/n)*var(ests)
  bias=mean(ests)-truth
  mse=mean((ests-truth)^2)
  c(var=var,bias=bias,mse=mse)
}

#' @export
gendata=function(n,g0, Q0){
  W1 = runif(n,-3,3)
  W2=rnorm(n)
  W3=runif(n)
  W4=rnorm(n)
  A=rbinom(n,1,g0(W1,W2,W3,W4))
  Y=rbinom(n,1,Q0(A,W1,W2,W3,W4))
  data.frame(A,W1,W2,W3,W4,Y)
}

#' @export
Q0_trig1= function (A, W1, W2, W3, W4) 
{
  plogis(.14*(2* A + 2*A * W1 + 20*cos(W1) * A - 3*W1 * sin(2*W2)+ cos(W1)
              -3*W2+4*A*(W2^2) +3*cos(W4)*A +A*W1^2- 2 * sin(W2)*W4 - 6*A* W3 * W4-3))
}

#' @export
Q0_trig =function (A, W1, W2, W3, W4)
{
  plogis(.14*(2* A  + 20*cos(W1) * A +cos(W1)-4*A*(W2^2) +3*cos(W4)*A +A*W1^2))
}

#' @export
Q0_1 = function (A, W1, W2, W3, W4) 
{
  plogis(.14*(2* A + 2*A * W1 + 4*A*W3*W4+W2*W1+W3*W4+10*A*cos(W4)))
}

#' @export
Q0_2 = function (A, W1, W2, W3, W4) 
{
  plogis(.14*(2* A + 5*A * W1 + 4*A*W3*W4+W2*W1+W3*W4+10*A*cos(W4)))
}

#' @export
g0_linear= function (W1, W2, W3, W4) 
{
  plogis(.5*(-0.8 * W1 + 0.39 * W2 + 0.08 * W3 - 0.12 * W4 - 0.15))
}

#' @export
g0_1 = function (W1, W2, W3, W4) 
{
  plogis(.5*(-0.08 * W1^2*W2+.5*W1 + 0.49 * cos(W2)*W3 + 0.18 * W3^2 - 0.12 * sin(W4) - 0.15))
}

#' @export
Q0_noise =function (A, W1, W2, W3, W4)
{
  plogis(.2*(.1*A+2*A*W1-10*A*W2+3*A*W3^2+W1+W2+.4*W3+.3*W4))
}

#' @export
gendata_noise=function(n, g0, Q0){
  W1 = runif(n,-3,3)
  # W1= rnorm(n)
  # W1=rnorm(n)
  W2=rbinom(n,1,.5)
  W3=rnorm(n)
  W4=rnorm(n)
  A=rbinom(n,1,g0(W1,W2,W3,W4))
  Y=rbinom(n,1,Q0(A,W1,W2,W3,W4))
  data.frame(A,W1,W2,W3,W4,Y)
}

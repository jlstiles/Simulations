# Script to generate data.  Low dimension and not may covs that matter
devtools::install_github("jlstiles/Simulations")
library(Simulations)

dgps1 = lapply(1:20, FUN = function(x) {
  dgp = get.dgp(n = 500, d = 22, pos = 0.05, minATE0 = 0, minVTE0 = 0, depth = 3, maxterms = 5, minterms = 3, 
                mininters = 3, num.binaries = 5, force.confounding = TRUE, N = 1e5)
  temp = remakeDGP(n=1e6, dgp)
  ATE0 = mean(temp$TE_n)
  VTE0 = var(temp$TE_n)
  return(return(list(dgp = dgp, ATE0 = ATE0, VTE0 = VTE0)))
})

save(dgps1, file = "dgps1.RDATA")
rm(dgps1)
# Low dimension and many covs that matter and more depth

dgps2 = lapply(1:20, FUN = function(x) {
  dgp = get.dgp(n = 500, d = 22, pos = 0.05, minATE0 = 0, minVTE0 = 0, depth = 4, maxterms = 10, minterms = 8, 
                mininters = 8, num.binaries = 5, force.confounding = TRUE, N = 1e5) 
  temp = remakeDGP(n=1e6, dgp)
  ATE0 = mean(temp$TE_n)
  VTE0 = var(temp$TE_n)
  return(return(list(dgp = dgp, ATE0 = ATE0, VTE0 = VTE0)))
})
save(dgps2, file = "dgps2.RDATA")
rm(dgps2)
#Low dimension main terms with some interactions

dgps3 = lapply(1:20, FUN = function(x) {
  dgp = get.dgp(n = 500, d = 22, pos = 0.05, minATE0 = 0, minVTE0 = 0, depth = 2, maxterms = 22, minterms = 15, 
                mininters = 12, num.binaries = 5, force.confounding = TRUE, N = 1e5)   
  temp = remakeDGP(n=1e6, dgp)
  ATE0 = mean(temp$blip_n)
  VTE0 = var(temp$blip_n)
  return(return(list(dgp = dgp, ATE0 = ATE0, VTE0 = VTE0)))
})

save(dgps3, file = "dgps3.RDATA")
rm(dgps3)
# High dimension with not many terms that matter and 3-way interaction

dgps4 = lapply(1:20, FUN = function(x) {
  dgp = get.dgp(n = 2000, d = 180, pos = 0.05, minATE0 = 0, minVTE0 = 0, depth = 2, maxterms = 10, minterms = 5, 
                mininters = 5, num.binaries = 20, force.confounding = TRUE, N = 1e5)  
  temp = remakeDGP(n=1e6, dgp)
  ATE0 = mean(temp$blip_n)
  VTE0 = var(temp$blip_n)
  return(return(list(dgp = dgp, ATE0 = ATE0, VTE0 = VTE0)))
})

save(dgps4, file = "dgps4.RDATA")
rm(dgps4)
# High dimension with more terms that matter only 2-way interaction

dgps5 = lapply(1:20, FUN = function(x) {
  dgp = get.dgp(n = 2000, d = 180, pos = 0.05, minATE0 = 0, minVTE0 = 0, depth = 2, maxterms = 40, minterms = 15, 
                mininters = 15, num.binaries = 20, force.confounding = TRUE, N = 1e5) 
  temp = remakeDGP(n=1e6, dgp)
  ATE0 = mean(temp$blip_n)
  VTE0 = var(temp$blip_n)
  return(return(list(dgp = dgp, ATE0 = ATE0, VTE0 = VTE0)))
})

save(dgps5, file = "dgps5.RDATA")
rm(dgps5)

# # check confounding:
# mean(dgp$blip_n)
# mean(dgp$DF$Y[dgp$DF$A==1]) - mean(dgp$DF$Y[dgp$DF$A==0])
# 
# dgp$coef_Q
# dgp$coef_G
# dgp$ATE0
# 
# big  = remakeDGP(1e6, dgp)
# head(big$DF)
# mean(big$blip_n)
# mean(big$DF$Y[big$DF$A==1]) - mean(big$DF$Y[big$DF$A==0])
# 
# # draw from same pop as dgp:
# dgp1 = remakeDGP(n=1000, object = dgp)
# mean(dgp1$blip_n)
# mean(dgp1$DF$Y[dgp1$DF$A==1]) - mean(dgp1$DF$Y[dgp1$DF$A==0])
# 
# dgps = lapply(1:20, FUN = function(x) {
#   dgp = get.dgp(n = 500, d = 22, pos = 0.05, minATE = 0, minBV = 0, depth = 3, maxterms = 5, minterms = 3, 
#                 mininters = 3, num.binaries = 5, force.confounding = TRUE, N = 1e5)
#   # temp = remakeDGP(n=1e6, dgp)
#   # ATE0 = mean(temp$blip_n)
#   return(dgp)
# })
# 
# dgps_stats = lapply(dgps, FUN = function(x) {
#   temp = remakeDGP(n=1e6, x)
#   ATE0 = mean(temp$blip_n)
#   VTE0 = var(temp$blip_n)
#   return(c(ATE0 = ATE0, VTE0 = VTE0))
# })
# 
# dgps_stats
# 
# B = 20
# library(foreach)
# detectCores()
# cl = makeCluster(detectCores(), type = "SOCK")
# 
# ALL=foreach(i=1:B,.packages=c("Simulations"), 
#             .errorhandling = "remove")%dopar%
#             {dgp = get.dgp(n = 500, d = 22, pos = 0.05, minATE = 0, minBV = 0.01, depth = 3, maxterms = 5, minterms = 3, 
#                            mininters = 3, num.binaries = 5, force.confounding = TRUE, N = 1e5)
#              temp = remakeDGP(n=1e6, dgp)
#              ATE0 = mean(temp$blip_n)
#              return(list(dgp, ATE0))
#             }
# 
# length(ALL)


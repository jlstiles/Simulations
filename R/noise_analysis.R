#' @export
noise.1 = function(data,n,rate, biasQ, sdQ) {
  rnorm(nrow(data), 
        with(data, biasQ(1,W1,W2,W3,W4,n=n,rate=rate)),
        with(data, sdQ(1,W1,W2,W3,W4,n=n,rate=rate)))
}

#' @export
# noise on barQ(0,W) is correlated per draw so not any more variant
noise.0 = function(data,noise,n,rate, biasQ, sdQ) {
  .5*noise + 
    sqrt(.75)*rnorm(nrow(data), 
                    with(data, biasQ(0,W1,W2,W3,W4,n=n,rate=rate)),
                    with(data, sdQ(0,W1,W2,W3,W4,n=n,rate=rate)))
}

#' @export
noise = function(data, noise_1, noise_0, biasQ, sdQ) {
  with(data, A*noise_1+(1-A)*noise_0)
}

#' @export
simBlipvar = function(n,rate, g0, Q0, biasQ, sdQ){
  # tack on the noise and use as an initial estimate
  data = gendata_noise(n, g0, Q0)
  noise_1 = noise.1(data, n, rate, biasQ, sdQ)
  noise_0 = noise.0(data, noise_1, n, rate, biasQ, sdQ)
  noise_A = noise(data, noise_1, noise_0 ,biasQ, sdQ)
  # noise_G = noiseG(data,V,n,rate)
  
  Q1k = plogis(qlogis(with(data,Q0(1,W1,W2,W3,W4)))+noise_1)
  Q0k = plogis(qlogis(with(data,Q0(0,W1,W2,W3,W4)))+noise_0)
  Qk = plogis(qlogis(with(data,Q0(A,W1,W2,W3,W4)))+noise_A)
  # gk = plogis(qlogis(with(data,g0(W1,W2,W3,W4)))+noise_G)
  gk = with(data,g0(W1,W2,W3,W4))
  
  initdata = data.frame(Qk=Qk,Q1k=Q1k,Q0k=Q0k,gk=gk,A=data$A,Y=data$Y)
  
  ATE_info = gentmle2::gentmle(initdata=initdata, params=list(param_ATE), 
                               submodel = submodel_logit, loss = loss_loglik,
                               approach = "line", max_iter = 100,g.trunc = 1e-2)
  sigmait_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                   submodel = submodel_logit, loss = loss_loglik,
                                   approach = "full", max_iter = 100,g.trunc = 1e-2)
  sigma_info = gentmle2::gentmle(initdata=initdata, params=list(param_sigmaATE), 
                                 submodel = submodel_logit, loss = loss_loglik,
                                 approach = "recursive", max_iter = 10000, g.trunc = 1e-2)
  
  steps = c(ATE_steps = ATE_info$steps, sigmait_steps = sigmait_info$steps,
            sigma_info = sigma_info$steps)
  
  ATE_ci = gentmle2::ci_gentmle(ATE_info)[c(2,4,5)]
  sigmait_ci = gentmle2::ci_gentmle(sigmait_info)[c(2,4,5)]
  sigma_ci = gentmle2::ci_gentmle(sigma_info)[c(2,4,5)]
  
  converges = c(sigmait = sigmait_info$converge,sigma = sigma_info$converge,
                ATE = ATE_info$converge)
  
  initest = var(Q1k-Q0k)
  return(c(sigmait_ci = sigmait_ci, ci_sigma = sigma_ci, ATE_ci = ATE_ci,
           sigma_init = initest, ATE_init = ATE_info$initests,
           steps = steps, converges = converges))
}

#' @export
getRes = function(allresults,B, ATE0, var0, varind = c(4,10), 
                  ateind = c(7,11)) {
  # allresults = L[[3]]
  # B=1000
  results = vapply(1:B,FUN = function(x) (allresults[[x]]),
                   FUN.VALUE=allresults[[1]])
  results1=as.matrix(results)
  results=apply(results1,2,as.numeric)
  row.names(results)=row.names(results1)
  results=t(results)
  # results = results[results[,"converges.sigma"]==0,]
  cov.sig.it = mean(results[,2]<=var0&var0<=results[,3])
  cov.sig.1step = mean(results[,5]<=var0&var0<=results[,6])
  cov.ate = mean(results[,8]<=ATE0&ATE0<=results[,9])
  # cov.onestepEst = mean(results[,36]<=var0&var0<=results[,37])
  coverage = c(cov.sig.it=cov.sig.it,
               cov.sig.1step=cov.sig.1step,
               cov.ate=cov.ate)
  
  ests=c(results[,7],results[,11])
  type = c(rep("TMLE",B),rep("init",B))
  
  ateests = data.frame(ests=ests,type=type)
  ggover1 = ggplot2::ggplot(ateests,aes(ests, fill=type)) + 
    geom_density(alpha=.3)+
    scale_fill_manual(values=c("red", "blue"))+
    theme(axis.title.x = element_blank())+ggtitle("ATE sampling distributions")+
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size=8))
  ggover1 = ggover1+geom_vline(xintercept = c(ATE0),color=c("green"))+
    # annotate("text",x=c(.2,.3),y=0,label=c("Event1","Event2"),hjust=1,angle=-40)
    # geom_text(x=.1,y=-5,label="asdfs",angle=-40,size = 3,hjust=1)
    geom_vline(xintercept=mean(results[,7]),color = "red")+
    geom_vline(xintercept=mean(results[,11]),color = "blue")+
    geom_vline(xintercept=ATE0,color = "green")
  # breaks = round(c(seq(min(ests),max(ests),length.out=4),ATE0,
  #            mean(results[,11]),mean(results[,7])),3)
  # ggover1 + scale_x_continuous(labels = c(round(seq(min(ests),max(ests),length.out=4),3),"ATE0",
  #                                         "mean init", "mean of TMLE")[order(breaks)],
  #                              breaks=breaks[order(breaks)])
  ggover1=cowplot::ggdraw(add_sub(ggover1,"Truth at green line", 
                         x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                         vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                         colour = "black", size = 8, angle = 0, lineheight = 0.9))
  
  ests = c(results[,4],results[,10])
  type= c(rep("one step",B),rep("init",B))
  # ests = c(results[,23],results[,7],results[,34])
  # type= c(rep("LR",B),rep("one-step multi",B),rep("init",B))
  varests = data.frame(ests=ests,type=type)[1:B,]
  
  ggover2 = ggplot2:: ggplot(varests,aes(ests)) + 
    geom_density(alpha=.3)+
    theme(axis.title.x = element_blank())+
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size=8,angle=315))+
    ggtitle("CATE variance sampling distributions")
  ggover2 = ggover2+geom_vline(xintercept = var0,color="green")+
    geom_vline(xintercept=mean(results[,4]),color = "black")+
    geom_vline(xintercept=mean(results[,10]),color = "red")
  # breaks = round(c(seq(min(varests[,1]),max(varests[,1]),length.out=4),mean(results[,10])),3)
  # ggover2 = ggover2 + scale_x_continuous(labels = c(round(seq(min(varests[,1]),max(varests[,1]),length.out=4),3),
  #                                         "mean init")[order(breaks)],
  #                              breaks=breaks[order(breaks)])
  ggover2=cowplot::ggdraw(add_sub(ggover2,"truth at green line. mean of initial est at red line",x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                         vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                         colour = "black", size = 8, angle = 0, lineheight = 0.9))
  
  performance.sig = t(apply(results[,varind], 2, perf,var0))
  performance.ate = t(apply(results[,ateind], 2, perf,ATE0))
  
  res = list(performance.ate=performance.ate, performance.sig=performance.sig, 
             coverage=coverage,ggover1=ggover1, ggover2=ggover2)
  return(res)
}

#' @export
noise_analysis = function(n, rate, truth, coverage, Q0, biasQ, sdQ) {
  # rate = .33
  # # V=1
  # # truth = gendata(10000)
  # n
  Q1_true = with(truth, Q0(1,W1,W2,W3,W4))
  Q0_true = with(truth, Q0(0,W1,W2,W3,W4))
  Q_true = with(truth,Q0(A,W1,W2,W3,W4))
  # g_true = with(truth, g0(W1,W2,W3,W4))
  
  noise.1 = function(data,n,rate) {
    rnorm(nrow(data), with(data, biasQ(1,W1,W2,W3,W4,n=n,rate=rate)), 
          with(data, sdQ(1,W1,W2,W3,W4,n=n,rate=rate)))}
  # noise on barQ(0,W) is correlated per draw so not any more variant
  noise.0 = function(data,noise,n,rate) {
    .5*noise +
      sqrt(.75)*rnorm(nrow(data), with(data, biasQ(0,W1,W2,W3,W4,n=n,rate=rate)),
                      with(data, sdQ(0,W1,W2,W3,W4,n=n,rate=rate)))}
  
  noise = function(data, noise_1, noise_0) with(data, A*noise_1+(1-A)*noise_0)
  
  noise_1 = noise.1(truth,n,rate)
  noise_0 = noise.0(truth,noise_1,n,rate)
  noise_A = noise(truth, noise_1, noise_0)
  # noise_G = noiseG(truth,V,n,rate)
  var(noise_1-noise_0)
  var(noise_A)
  Q1_test = plogis(qlogis(Q1_true) + noise_1)
  Q0_test = plogis(qlogis(Q0_true) + noise_0)
  Q_test = plogis(qlogis(Q_true) + noise_A)
  # g_test = plogis(qlogis(g_true) + noise_G)
  
  blip_test = Q1_test - Q0_test
  blip_true = Q1_true - Q0_true
  
  L2_blip = sqrt(mean((blip_test - blip_true)^2))
  L2_Q = sqrt(mean((Q_true - Q_test)^2))
  # L2_G = sqrt(mean((g_true - g_test)^2))
  
  # hist(Q_test)
  # hist(Q_true)
  # hist(blip_test)
  # hist(blip_true)
  # hist(g_test)
  ate_bias = mean(blip_test) - ATE0
  var_bias = var(blip_test) - var0
  
  title = paste("sample size ",n, "variance bias = ",round(var_bias,5))
  df = data.frame(CATE = c(blip_test,blip_true), type = c(rep("test",1e6),rep("true",1e6)))
  gg_testvstrue = ggplot2::ggplot(df,aes(x=CATE,color=type))+geom_density()+ggtitle(title)+
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(size=8,angle=315)) 
  caption = paste0("coverage = ",coverage)
  gg_testvstrue = cowplot::ggdraw(add_sub(gg_testvstrue,caption, x= .05, y = 0.5, hjust = 0, vjust = 0.5, vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                                 colour = "black", size = 9, angle = 0, lineheight = 0.9))
  # var_bias
  # ate_bias
  
  results = list(L2_blip = n^.25*L2_blip, L2_Q = n^.25*L2_Q,
                 ate_bias = ate_bias, var_bias = var_bias,plot = gg_testvstrue)
  return(results)
}

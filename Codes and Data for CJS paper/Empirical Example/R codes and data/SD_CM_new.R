### Jan 12, 2020

### This program tests for SD: H_0: F dominates G. In terms of CDF: G(z) >= F(z)
### General idea: reject if F(z) > G(z) at some z value.
###  Dilemma: F(z) = G(z) at two tails. one may not want to discount evidence there.
###  Hence, the KS actually makes good sense. 
###  By ignoring differential variability, it down plays the difference at two tails

### Test statistics:  we regard KS and T as two methods in two extremes
###  One could use  KS/\{ F(1-F) + G(1-G) \}^a as a test statistics.
###   KS: a=0; TT: a = 1/2;  We add a = 1/8, 1/4 as two intermediate choices.
###   We try DD distributions here

### It is felt that all previous methods do not have tight type I error for forefrong config.
###  It indicates that the centralization idea does not work well.

### This program tries our LLFC idea: creating new F and G from which we bootstrap.


SDCM.func = function(x0,y0,Boot.n) {
  nn1 = length(x0)       ## sample size of xx
  nn2 = length(y0)        ## sample size of y
  
  BB = matrix(0, Boot.n, 4)
  
  B.qq = c(.90, .95,.99)
  
  
  
  out0l=matrix(NA,8,4)  ## results for SD test
  colnames(out0l)=c("KS","T1/8","T1/4","T1/2")
  rownames(out0l)=c("Stats","Pvalue","C.10","C.05","C.01","m","n","B")
  
  out0l[6,1]=nn1
  out0l[7,1]=nn2
  out0l[8,1]=Boot.n
  
  
  
  epsilon = 1/(nn1+nn2)^2 
  ####  constant to prevent 0-denominator.
  
  
  ########
  tt = sort(c(x0, y0))[3:(nn1+nn2-3)] 
  FFm = ecdf(x0)(tt)
  GGn = ecdf(y0)(tt)
  DFFGG1 = FFm - GGn
  ## Put on weight on two ends than staight KS.
#  ks.pos = DFFGG1*(DFFGG1>0)
  ### used for centralization.
  temp = 1/(GGn*(1-GGn)/nn2 + FFm*(1-FFm)/nn1+epsilon)
  TT0 = max(DFFGG1)
  TT1 = max(DFFGG1*temp^.125)   ### t-statistic
  TT2 = max(DFFGG1*temp^.25)   ### t-statistic
  TT3 = max(DFFGG1*temp^.5)   ### t-statistic
 ##
  for(bb in 1:Boot.n) {
    xx.b = sample(x0,replace=T);
    yy.b = sample(y0,replace=T)
    tt = sort(c(xx.b, yy.b))[3:(nn1+nn2-3)]
    FFm = ecdf(xx.b)(tt)
    GGn = ecdf(yy.b)(tt)
    DFFGG1 = FFm - GGn
    temp = 1/(GGn*(1-GGn)/nn2 + FFm*(1-FFm)/nn1+epsilon)
    ###
    FFm0 = ecdf(x0)(tt)
    GGn0 = ecdf(y0)(tt)
    DFFGG10 = FFm0 - GGn0
    ## Put on weight on two ends than staight KS.
    ks.pos = DFFGG10*(DFFGG10>0)
    ###
    TT0.b = max(DFFGG1 - ks.pos)
    TT1.b = max((DFFGG1 - ks.pos)*temp^.125)
    TT2.b = max((DFFGG1 - ks.pos)*temp^.25)
    TT3.b = max((DFFGG1 - ks.pos)*temp^.5)
    BB[bb,] = c(TT0.b, TT1.b, TT2.b, TT3.b)
  }  ## end of bb.
  
  
  pvalue.l=c(length(BB[BB[,1]>TT0,1]),
             length(BB[BB[,2]>TT1,2]),
             length(BB[BB[,3]>TT2,3]),
             length(BB[BB[,4]>TT3,4]))/Boot.n
  Bqq.l = c(quantile(BB[,1], B.qq), quantile(BB[,2], B.qq), quantile(BB[,3], B.qq),quantile(BB[,4], B.qq))
  
  out0l[1,]=c(TT0,TT1,TT2,TT3)
  out0l[2,]=pvalue.l
  out0l[3,]=Bqq.l[c(1,4,7,10)]
  out0l[4,]=Bqq.l[c(2,5,8,11)]
  out0l[5,]=Bqq.l[c(3,6,9,12)]
  
  out0l=round(out0l,4)
  
  list(out_c=out0l,pvalue.c=pvalue.l,Bqq.c=Bqq.l)
  
   
}
  
   
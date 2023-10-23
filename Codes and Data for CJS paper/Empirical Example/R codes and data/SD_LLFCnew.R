### June 26, 2022: This one shows promise
## LLFC method
### This program tests for SD: H_0: F dominates G. In terms of CDF: G(z) >= F(z)
### General idea: reject if F(z) > G(z) at some z value.
###
###  function of new LLFC method used in the empirical example
###


SDLLFCnew.func = function(x0,y0,Boot.n) {
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
  
  cc0=-log(log(nn1+nn2)) ##-loglog(nn1+nn2)
  
  ###large negative value as penalty
  nv0=-10000.0
  
  
  ## SD value allows us to echange F and G.
  tt = sort(c(x0, y0))
  FFm = ecdf(x0)(tt); FFm = (FFm*(nn1+nn2)+1)/(nn1+nn2+2)
  GGn = ecdf(y0)(tt); GGn = (GGn*(nn1+nn2)+1)/(nn1+nn2+2)
 ########
  DFFGG1 = (FFm - GGn)[3:(nn1+nn2-2)]
  temp = (GGn*(1-GGn)/nn2 + FFm*(1-FFm)/nn1+epsilon)[3:(nn1+nn2-2)]
  TT0 = max(DFFGG1)            ### KS statistic
  TT1 = max(DFFGG1/temp^.125)   
  TT2 = max(DFFGG1/temp^.25)   ### t-statistic
  TT3 = max(DFFGG1/temp^.5)   ### t-statistic
  #################
  for(bb in 1:Boot.n) {
    xx.b = sample(x0,replace=T);
    yy.b = sample(y0,replace=T)
    ######
    tt = sort(c(xx.b, yy.b,x0,y0))
    FFm = ecdf(xx.b)(tt); FFm = (FFm*(nn1+nn2)+1)/(nn1+nn2+2)
    GGn = ecdf(yy.b)(tt); GGn = (GGn*(nn1+nn2)+1)/(nn1+nn2+2)
    DFFGG1 = (FFm - GGn)[6:((nn1+nn2)*2-5)]
    temp = (GGn*(1-GGn)/nn2 + FFm*(1-FFm)/nn1+epsilon)[6:((nn1+nn2)*2-5)]
    #############
    FFm0 = ecdf(x0)(tt); FFm0 = (FFm0*(nn1+nn2)+1)/(nn1+nn2+2)
    GGn0 = ecdf(y0)(tt); GGn0 = (GGn0*(nn1+nn2)+1)/(nn1+nn2+2)
    ##############
    DFFGG10 = (FFm0 - GGn0)[6:((nn1+nn2)*2-5)]
    #########
    ############
    cen.b = DFFGG10
    ### used for centralization in bootstrap.  	
    ###### set up 4 sets
    fi1=rep(1,length(DFFGG10))
    fi1=fi1*(DFFGG10>cc0)
    fi1nv=(1-fi1)*nv0
    fi2=rep(1,length(DFFGG10))
    fi2=fi2*((DFFGG10/temp^.125)>cc0)
    fi2nv=(1-fi2)*nv0
    fi3=rep(1,length(DFFGG10))
    fi3=fi3*((DFFGG10/temp^.25)>cc0)
    fi3nv=(1-fi3)*nv0
    fi4=rep(1,length(DFFGG10))
    fi4=fi4*((DFFGG10/temp^.5)>cc0)
    fi4nv=(1-fi4)*nv0
    #########
    tmp0b=(DFFGG1-cen.b)*fi1+fi1nv
    TT0.b = max(tmp0b)            ### KS statistic
    tmp1b=((DFFGG1-cen.b)/temp^.125)*fi2+fi2nv
    TT1.b = max(tmp1b)
    tmp2b=((DFFGG1-cen.b)/temp^.25)*fi3+fi3nv
    TT2.b = max(tmp2b)   ### t-statistic
    tmp3b=((DFFGG1-cen.b)/temp^.5)*fi4+fi4nv
    TT3.b = max(tmp3b)   ### t-statistic
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
  
  list(out_l=out0l,pvalue.l=pvalue.l,Bqq.l=Bqq.l)
  
   
}
  
   
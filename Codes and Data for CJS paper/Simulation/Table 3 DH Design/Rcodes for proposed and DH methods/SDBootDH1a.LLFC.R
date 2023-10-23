### 26 June 2022  Table 3 DH1 Design
#
###   The main features are: create LLFC via log-odds ratio. F dominates G
###   Weight KS statistics by [ F(1-F)/m + G(1-G)/n ]^a: a = 0, 1/8, 1/4, 1/2
###### DH1 design

#setwd("C:/simuLN1/SD_new26June22/version3/")

source("Generator.Rd")

Nrep = 20000     ## replications
ratio = .7          ## ratio of two sample sizes
nn1 = 200;       ## sample size of xx
nn2 = ratio*nn1    ## sample size of yy
Boot.n= 501          ## bootstrap number of repetitions

epsilon = 1/(nn1+nn2)^2
####  constant to prevent 0-denominator.

cc0=-log(log(nn1+nn2)) ##-loglog(nn1+nn2)

###large negative value as penalty
nv0=-10000.0

DH = T
ll = c(.25, .55, .8)
DH.pv = matrix(NA, length(ll), 12)

for(jj in 1:length(ll)) {
	All.pv = matrix(NA, Nrep, 4)
	BB = matrix(NA, Boot.n, 4)

    for(i in 1:Nrep) {
    if(DH) {
    	xxyy = generate.DH1(ll[jj], nn1, nn2, SD=T)      ## data generation
	} else {
    	xxyy = generate.DH2(ll[jj], nn1, nn2, SD=T)      ## data generation
	}

        tt = sort(c(xxyy[[1]], xxyy[[2]]))
        FFm = ecdf(xxyy[[1]])(tt); FFm = (FFm*(nn1+nn2)+1)/(nn1+nn2+2)
        GGn = ecdf(xxyy[[2]])(tt); GGn = (GGn*(nn1+nn2)+1)/(nn1+nn2+2)
##############
        DFFGG1 = (FFm - GGn)[3:(nn1+nn2-2)]
        temp = (GGn*(1-GGn)/nn2 + FFm*(1-FFm)/nn1+epsilon)[3:(nn1+nn2-2)]
        TT0 = max(DFFGG1)            ### KS statistic
        TT1 = max(DFFGG1/temp^.125)
        TT2 = max(DFFGG1/temp^.25)   ### t-statistic
        TT3 = max(DFFGG1/temp^.5)   ### t-statistic
#####################
        #####    	
    for(bb in 1:Boot.n) {
      xx.b = sample(xxyy[[1]],replace=T);
      yy.b = sample(xxyy[[2]],replace=T)
      ######
      tt = sort(c(xx.b, yy.b,xxyy[[1]],xxyy[[2]]))
      FFm = ecdf(xx.b)(tt); FFm = (FFm*(nn1+nn2)+1)/(nn1+nn2+2)
      GGn = ecdf(yy.b)(tt); GGn = (GGn*(nn1+nn2)+1)/(nn1+nn2+2)
      DFFGG1 = (FFm - GGn)[6:((nn1+nn2)*2-5)]
      temp = (GGn*(1-GGn)/nn2 + FFm*(1-FFm)/nn1+epsilon)[6:((nn1+nn2)*2-5)]
      #############
      FFm0 = ecdf(xxyy[[1]])(tt); FFm0 = (FFm0*(nn1+nn2)+1)/(nn1+nn2+2)
      GGn0 = ecdf(xxyy[[2]])(tt); GGn0 = (GGn0*(nn1+nn2)+1)/(nn1+nn2+2)
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
    All.pv[i,] = colMeans(cbind(BB[,1]>TT0,BB[,2]>TT1,BB[,3]>TT2,BB[,4]>TT3))
    }  ##end of Nrep
DH.pv[jj,] = c(colMeans(All.pv<0.01),colMeans(All.pv<0.05),colMeans(All.pv<0.1))
}
write(c("LLFC, TT0-3;nn1,nn2,Nrep, DH",nn1,nn2,Nrep),"DH1pva20k",ncolumns=4,append=T)
write(c("ll in DH=", ll, DH), "DH1pva20k", ncolumns=8, append=T)
write(round(t(DH.pv),4), "DH1pva20k", ncolumns=12, append=T)


write.csv(DH.pv,file="LLFC_DH1a20k.csv")



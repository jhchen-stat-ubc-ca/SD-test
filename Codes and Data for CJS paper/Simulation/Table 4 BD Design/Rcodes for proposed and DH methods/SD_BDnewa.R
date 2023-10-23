### June 26, 2022: Table 4  BD Design

### This program tests for SD: H_0: F dominates G. In terms of CDF: G(z) >= F(z)
### General idea: reject if F(z) > G(z) at some z value.
###  Dilemma: F(z) = G(z) at two tails. one may not want to discount evidence there.
###  Hence, the KS actually makes good sense. 
###  By ignoring differential variability, it down plays the difference at two tails
###
###   BD design


#setwd("C:/simuLN1/")

#setwd("C:/simuLN1/SD_new26June22/")

#setwd("C:/simuLN1/SD_new26June22/version3/")


source("Generator.Rd")

source("data_generator1.R")


Nrep = 20000      ## replications
#ratio = 0.7         ## ratio of two sample sizes
nn1 = 200;       ## sample size of xx
nn2 = 140        ## sample size of y
#nn2 = ratio*nn1    ## sample size of yy
Boot.n= 501          ## bootstrap number of repetitions

epsilon = 1/(nn1+nn2)^2 
####  constant to prevent 0-denominator.


cc0=-log(log(nn1+nn2)) ##-loglog(nn1+nn2)


###large negative value as penalty
nv0=-10000.0

 
#### alternatives

  DHpv=matrix(2, 4, 12)
  
  colnames(DHpv)=c("T0", "T1/8",  "T1/4", "T1/2",
                   "T0", "T1/8",  "T1/4", "T1/2",
                   "T0", "T1/8",  "T1/4", "T1/2") 
  row.names(DHpv)=c("k1","k2","k3","k4")
  

  ###
  #####
  dnm1=c(1:Nrep)
  dnm2=c("T0", "T1/8",  "T1/4", "T1/2") 
  dnm3=c("k1","k2","k3","k4")
  
  Allpv=array(rep(0.5,Nrep*4*4),dim=c(Nrep,4,4),dimnames = list(dnm1,dnm2,dnm3))
  ########
	for(i in 1:Nrep) {
     	xxyy=generate.LN1(nn1,nn2)  ##data generated for X and Y1-Y5 where Y_{kk} is the sample of G for case kk
    ########
     	for (jj in 1:4) {   ##jj different cases
   #  	  All.pv = matrix(0.5, Nrep, 4)
     	  BB = matrix(0, Boot.n, 4)
     	 x0=xxyy$xx
     	 y0=xxyy$yy[,jj+1]
   ## SD = T:   G is uniform distribution.
		## SD value allows us to echange F and G.
    	tt = sort(c(x0, y0))
    	FFm = ecdf(x0)(tt); FFm = (FFm*(nn1+nn2)+1)/(nn1+nn2+2)
    	GGn = ecdf(y0)(tt); GGn = (GGn*(nn1+nn2)+1)/(nn1+nn2+2)

    	##############
    	DFFGG1 = (FFm - GGn)[3:(nn1+nn2-2)]
    	temp = (GGn*(1-GGn)/nn2 + FFm*(1-FFm)/nn1+epsilon)[3:(nn1+nn2-2)]
    	TT0 = max(DFFGG1)            ### KS statistic
    	TT1 = max(DFFGG1/temp^.125)   
    	TT2 = max(DFFGG1/temp^.25)   ### t-statistic
    	TT3 = max(DFFGG1/temp^.5)   ### t-statistic
    	############
    	#####
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
      #############
      FFm0 = ecdf(x0)(tt); FFm0 = (FFm0*(nn1+nn2)+1)/(nn1+nn2+2)
      GGn0 = ecdf(y0)(tt); GGn0 = (GGn0*(nn1+nn2)+1)/(nn1+nn2+2)
      ##############
      DFFGG10 = (FFm0 - GGn0)[6:((nn1+nn2)*2-5)]
      #########
      ############
      cen.b = DFFGG10
      ### used for centralization in bootstrap.
      #########
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
  #  All.pv[i,] = colMeans(cbind(BB[,1]>TT0,BB[,2]>TT1,BB[,3]>TT2,BB[,4]>TT3))
    Allpv[i,,jj][] = colMeans(cbind(BB[,1]>TT0,BB[,2]>TT1,BB[,3]>TT2,BB[,4]>TT3))
    }  ## end of jj
  ##  DD.pv[jj,] = colMeans(cbind(All.pv<0.01,All.pv<0.05,All.pv<0.1))
} ##end of Nrep

##  DD.pv[jj,] = colMeans(cbind(All.pv<0.01,All.pv<0.05,All.pv<0.1))
for (jj in 1:4){
    a1=colMeans(Allpv[,,jj][]<0.01)
    a2=colMeans(Allpv[,,jj][]<0.05)
    a3=colMeans(Allpv[,,jj][]<0.1)
    DHpv[jj,]=c(a1,a2,a3)
}    
write.csv(DHpv, file="SD_BDa_LLFC20k.csv")

write(c("Nrep and nn1, nn2"), "SD_BDa_LLFC20k.txt",append=T)
write(c(Nrep, nn1, nn2), "SD_BDa_LLFC20k.txt",ncolumns=3,append=T)

#write(c("LLFC, TT0,1,2,3;nn1,nn2,Nrep,null=",nn1,nn2,Nrep,Null),"DHpv",ncolumns=3,append=T)
##write(round(t(DD),3), "DDpv", ncolumns=8, append=T)
##write(round(t(DD.pv),4), "DDpv", ncolumns=12, append=T)

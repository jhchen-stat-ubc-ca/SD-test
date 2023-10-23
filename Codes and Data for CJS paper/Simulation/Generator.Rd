## This file contains functions that generate two sample data.
### Plot cdf's and others. 

generate.DD = function (seg, nn1, nn2, SD) {
	### generate two samples: both taking values in the unit interval
	### Sample sizes are nn1 and nn2.
	### G is uniform and distribution of YY.
        		## It decides the percentage of observations
        		##     falls into 8 segments
        if(SD) {
        mm = rmultinom(1, nn1, seg)
        xx = NULL
        for(i in 1:8) xx = c(xx, (i-1+runif(mm[i]))/8)
        yy = runif(nn2)
	} else {
        mm = rmultinom(1, nn2, seg)
        yy = NULL
        for(i in 1:8) yy = c(yy, (i-1+runif(mm[i]))/8)
        xx = runif(nn1)
	}
        return(list(xx, yy))
}

generate.J1 = function (nn1, nn2, SD) {   
        mm = rbinom(1, nn1, 0.3)
        xx1 = rbeta(mm, 1, 1)^.5*0.3
        xx2 = 0.3+rbeta(nn1-mm, 1, 1)^.5*0.7
        xx = c(xx1, xx2)
        yy = runif(nn2)
        ### F dominates G, but CDF touch. cdf are round.
        if(SD) output = list(xx, yy)
        if(!SD) output = list(yy, xx)
        output
        }
 
generate.J2 = function (nn1, nn2, SD) {      
        mm = rbinom(1, nn1, 0.2)
        xx1 = rbeta(mm, 1, 1)^.5*0.3
        xx2 = 0.3+rbeta(nn1-mm, 2, 2)^.5*0.7
        xx = c(xx1, xx2)
        yy = runif(nn2)
        if(SD) output = list(xx, yy)
        if(!SD) output = list(yy, xx)
        output
	}
    
    
generate.DH1  = function (ll, nn1, nn2, SD) {  
       temp = runif(nn1)
	   yy = temp^2*(temp<=ll)/ll + temp*(temp>ll)
       xx = runif(nn2)
        if(SD) output = list(xx, yy)
        if(!SD) output = list(yy, xx)
        output
       }


generate.DH2  = function (ll, nn1, nn2, SD) {  
    temp = runif(nn1)
	yy = temp^2*(temp<=ll)/ll + temp*(temp>ll)
    temp = runif(nn2)
    xx = temp *(temp <= ll) + (ll + (temp - ll)^2/(1-ll))*(temp>ll)
        if(SD) output = list(xx, yy)
        if(!SD) output = list(yy, xx)
        output
}



generate.LN1  = function (nn1, nn2, mu1,mu2,sd1,sd2) {  
    temp = rnorm(nn1)  ##Z_{1i}
	xx = exp(sd1*tmp+mu1)
    temp = rnorm(nn2)  ##Z_{2i}
  yy = exp(sd2*tmp+mu2)
  output = list(xx, yy)
}



### Plot two empirical distributions 
#install.packages("Hmisc")
#library(Hmisc)
configuration.plot = function(xxyy) {
	gg <- c(rep('F',length(xxyy[[1]])), rep('G',length(xxyy[[2]])))
	Ecdf(c(xxyy[[1]],xxyy[[2]]), group=gg, xlab='', ylab='',label.curves=list(keys=1:2), subtitles=F)
}






##
##studies specified in Barrett and Donald (2003 Econometrica)
##
generate.LN1  = function (nn1, nn2) {  
## output: xx = one sample of nn1 for F, 
##         yy= nn2 by 5 columns with each for a sample of nn2 from G1-G5 seperately   
  z1 = rnorm(nn1)  ##Z_{1i}
  z2 = rnorm(nn2)  ##Z_{2i}
  #  Case 1
    mu1=0.85
    sd1=0.6
    mu2=mu1
    sd2=sd1
  xx1 = exp(sd1*z1+mu1)
  yy1 = exp(sd2*z2+mu2)
  #  Case 2  H0 is false
    mu1=0.85
    sd1=0.6
    mu2=0.6
    sd2=0.8
  # xx2 = xx1
  yy2 = exp(sd2*z2+mu2)
  ##
  # Case 3  H0 fails by a very small amount
  #
    mu1=0.85
    sd1=0.6
    mu2=1.2
    sd2=0.2
  # xx3 = xx1
  yy3 = exp(sd2*z2+mu2)
  ##
  # Case 4  ##two CDF have a single crossing
  ##
    mu1=0.85
    sd1=0.6
    mu2=0.8
    sd2=0.5
    mu3=0.9
    sd3=0.9
    z3 = rnorm(nn2)  ##Z_{3i}
    u_i= runif(nn2) ##U_i a uniform [0,1]
  #  xx4 = xx1
    yy4=rep(0,nn2)
    for (i in 1:nn2){
      if (u_i[i]>=0.1) {
          yy4[i]=exp(sd2*z2[i]+mu2)} else {
          yy4[i]=exp(sd3*z3[i]+mu3)
        }
      }
  ##
  # Case 5  ## two CDF have multiple crossings. 
  #
    mu1=0.85
    sd1=0.6
    mu2=0.85
    sd2=0.4
    mu3=0.4
    sd3=0.9
  #  xx5 = xx1
  #z3 = rnorm(nn2)  ##Z_{3i}
  #u_i= runif(nn2) ##U_i a uniform [0,1]
    yy5=rep(0,nn2)
    for (i in 1:nn2){
      if (u_i[i]>=0.1) {
        yy5[i]=exp(sd2*z2[i]+mu2)} else {
          yy5[i]=exp(sd3*z3[i]+mu3)
        }
    }
    
    yy=cbind.data.frame(yy1,yy2,yy3,yy4,yy5)
    output = list(xx=xx1, yy=yy)
}


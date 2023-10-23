


# load the original data 
load("realdata.Rda")

s1=log_UK_HHRelInc[[1]]  ##1968

s1=exp(s1)

s2=log_UK_HHRelInc[[2]] ##1969



s2=log_UK_HHRelInc[[11]] ##1978


s2=log_UK_HHRelInc[[5]] ##1968 +4


s2=log_UK_HHRelInc[[6]] ##1968 +5


s2=exp(s2)


s3=log_UK_HHRelInc[[11]] ##1968 +10


s3=exp(s3)



s4=log_UK_HHRelInc[[16]] ##1968 +15


s4=exp(s4)


s5=log_UK_HHRelInc[[21]] ##1968 +20


s5=exp(s5)




#write.csv(s1,'S1.csv')

#write.csv(s2,'S2.csv')

#write.csv(s3,'S3.csv')

#write.csv(s4,'S4.csv')

#write.csv(s5,'S5.csv')


##plot 1
####
x0=sort(c(s1,s2))

xmax0=max(x0)
xmin0=min(x0)


#xmax0=0.5
#xmin0=-0.5


n00=length(x0)

#x=c(1:10000)
#x=(x/10000)*xm0

#x=(x/20000)*(xm0+100)


#######S1/S2

F=ecdf(s1)(x0)
G=ecdf(s2)(x0)

######
plot(x0,F,type="l",xlab="relative income",ylab="CDF",ylim=c(0,1),xlim=c(xmin0,xmax0),pch=19, col="red")

lines(x0,G,type="l",pch=18, col="blue",lty=2,ylim=c(0,1),xlim=c(xmin0,xmax0))

legend(2, 0.5, legend=c("F (S1)", "G (S2)"),box.lty=1,
       col=c("red", "blue"), lty=1:2, cex=0.8)

title(main="S1/S2")
####
DGF=G-F
####
x=x0
##
xx=seq(1,length(x),by=11)
xx=c(xx,length(x))

######
plot(x0,DGF,type="l",xlab="relative income",ylab="G-F",pch=19)
abline(h=0)

title(main="S1/S2")




##plot 2
####
x0=sort(c(s2,s3))

xmax0=max(x0)
xmin0=min(x0)


#xmax0=0.5
#xmin0=-0.5


n00=length(x0)

#x=c(1:10000)
#x=(x/10000)*xm0

#x=(x/20000)*(xm0+100)


#######S2/S3

F=ecdf(s2)(x0)
G=ecdf(s3)(x0)

######
plot(x0,F,type="l",xlab="relative income",ylab="CDF",ylim=c(0,1),xlim=c(xmin0,xmax0),pch=19, col="red")

lines(x0,G,type="l",pch=18, col="blue",lty=2,ylim=c(0,1),xlim=c(xmin0,xmax0))

legend(2, 0.5, legend=c("F (S2)", "G (S3)"),box.lty=1,
       col=c("red", "blue"), lty=1:2, cex=0.8)

title(main="S2/S3")
####
DGF=G-F
####
x=x0
##
xx=seq(1,length(x),by=11)
xx=c(xx,length(x))

######
plot(x0,DGF,type="l",xlab="relative income",ylab="G-F",pch=19)
abline(h=0)

title(main="S2/S3")

##################################


##plot 3
####
x0=sort(c(s3,s4))

xmax0=max(x0)
xmin0=min(x0)


#xmax0=0.5
#xmin0=-0.5


n00=length(x0)

#x=c(1:10000)
#x=(x/10000)*xm0

#x=(x/20000)*(xm0+100)


#######S3/S4

F=ecdf(s3)(x0)
G=ecdf(s4)(x0)

######
plot(x0,F,type="l",xlab="relative income",ylab="CDF",ylim=c(0,1),xlim=c(xmin0,xmax0),pch=19, col="red")

lines(x0,G,type="l",pch=18, col="blue",lty=2,ylim=c(0,1),xlim=c(xmin0,xmax0))

legend(2, 0.5, legend=c("F (S3)", "G (S4)"),box.lty=1,
       col=c("red", "blue"), lty=1:2, cex=0.8)

title(main="S3/S4")
####
DGF=G-F
####
x=x0
##
xx=seq(1,length(x),by=11)
xx=c(xx,length(x))

######
plot(x0,DGF,type="l",xlab="relative income",ylab="G-F",pch=19)
abline(h=0)

title(main="S3/S4")
##########################



##plot 4
####
x0=sort(c(s4,s5))

xmax0=max(x0)
xmin0=min(x0)


#xmax0=0.5
#xmin0=-0.5


n00=length(x0)

#x=c(1:10000)
#x=(x/10000)*xm0

#x=(x/20000)*(xm0+100)


#######S4/S5

F=ecdf(s4)(x0)
G=ecdf(s5)(x0)

######
plot(x0,F,type="l",xlab="relative income",ylab="CDF",ylim=c(0,1),xlim=c(xmin0,xmax0),pch=19, col="red")

lines(x0,G,type="l",pch=18, col="blue",lty=2,ylim=c(0,1),xlim=c(xmin0,xmax0))

legend(2, 0.5, legend=c("F (S4)", "G (S5)"),box.lty=1,
       col=c("red", "blue"), lty=1:2, cex=0.8)

title(main="S4/S5")
####
DGF=G-F
####
x=x0
##
xx=seq(1,length(x),by=11)
xx=c(xx,length(x))

######
plot(x0,DGF,type="l",xlab="relative income",ylab="G-F",pch=19)
abline(h=0)

title(main="S4/S5")
##########################
##############

par(mfrow=c(2,2))

############################################



par(mfrow=c(1,1))


##plot 2
####
x0=sort(c(s2,s3))

xmax0=max(x0)
xmin0=min(x0)


#xmax0=0.5
#xmin0=-0.5


n00=length(x0)

#x=c(1:10000)
#x=(x/10000)*xm0

#x=(x/20000)*(xm0+100)


#######S2/S3

F=ecdf(s2)(x0)
G=ecdf(s3)(x0)

######
plot(x0,F,type="l",xlab="relative income",ylab="CDF",ylim=c(0,1),xlim=c(xmin0,xmax0),pch=19, col="red")

lines(x0,G,type="l",pch=18, col="blue",lty=2,ylim=c(0,1),xlim=c(xmin0,xmax0))

legend(0.2, 0.3, legend=c("F (S2)", "G (S3)"),box.lty=1,
       col=c("red", "blue"), lty=1:2, cex=0.8)

title(main="S2/S3")
####
DGF=G-F
####
x=x0
##
xx=seq(1,length(x),by=11)
xx=c(xx,length(x))

######
plot(x0,DGF,type="l",xlab="relative income",ylab="G-F",pch=19)
abline(h=0)

title(main="S2/S3")


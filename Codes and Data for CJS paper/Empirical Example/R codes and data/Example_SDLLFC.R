
rm(list=ls(all=TRUE))
######
##  Code for the analysis of the empirical example based on the proposed method
#
setwd("C:/Research/UBC/Codes Data/Empirical Example/R codes and data/")
#
source('SD_LLFCnew.R')
#
# load the original data 
load("realdata.Rda")
##
## 1968/1972
##
s1=log_UK_HHRelInc[[1]]  ##1968
s1=exp(s1)
s2=log_UK_HHRelInc[[6]] ##1972
s2=exp(s2)
##
Boot.n = 5001 ## number of replicates for bootstrap (odd number)
##
## S1 vs S2
outs1s2=SDLLFCnew.func(s1,s2,Boot.n)
outs2s1=SDLLFCnew.func(s2,s1,Boot.n)
#
output_llfc1=cbind(outs1s2$out_l,outs2s1$out_l)
round(output_llfc1,3)
write.csv(output_llfc1,file="S1_S2_SDllfc1.csv")
#
####
##
## 1972/1978
##
s2=log_UK_HHRelInc[[6]] ##1972
s2=exp(s2)
s3=log_UK_HHRelInc[[11]] ##1978
s3=exp(s3)
#
## S2 vs S3
outs2s3=SDLLFCnew.func(s2,s3,Boot.n)
outs3s2=SDLLFCnew.func(s3,s2,Boot.n)
#
output_llfc1=cbind(outs2s3$out_l,outs3s2$out_l)
round(output_llfc1,3)
write.csv(output_llfc1,file="S2_S3_SDllfc1.csv")
####
##
## 1978/1983
##
s3=log_UK_HHRelInc[[11]] ##1978
s3=exp(s3)
s4=log_UK_HHRelInc[[16]] ##1983
s4=exp(s4)
## S3 vs S4
outs3s4=SDLLFCnew.func(s3,s4,Boot.n)
outs4s3=SDLLFCnew.func(s4,s3,Boot.n)
#
output_llfc1=cbind(outs3s4$out_l,outs4s3$out_l)
#
round(output_llfc1,3)
#
write.csv(output_llfc1,file="S3_S4_SDllfc1.csv")
####
##
## 1983/1988
##
s3=log_UK_HHRelInc[[11]] ##1983
s3=exp(s3)
s4=log_UK_HHRelInc[[16]] ##1988
s4=exp(s4)
#
## S3 vs S4
outs3s4=SDLLFCnew.func(s3,s4,Boot.n)
outs4s3=SDLLFCnew.func(s4,s3,Boot.n)
#
output_llfc1=cbind(outs3s4$out_l,outs4s3$out_l)
#
round(output_llfc1,3)
#
write.csv(output_llfc1,file="S3_S4_SDllfc1.csv")










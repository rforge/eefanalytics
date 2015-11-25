data(iwq)

###############################################
## MLM analysis of cluster randomised trials ##
###############################################

output1 <- crtFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",data=iwq)

ES1 <- output1$ES$Intervention1

ICC <- output1$covParm[1]/sum(output1$covParm)
ICC

################################################################
## MLM analysis of cluster randomised trials with permutation ##
## under the null hypothesis                                  ##
################################################################

date()
output2 <- crtFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nPerm=1000,data=iwq)
ES2 <- output2$ES$Intervention1
ES2

ICC2 <- output2$covParm[1]/sum(output2$covParm)
ICC2

date()

perm <- output2$Perm


#### Within variance

obsg <- output2$ES$Intervention1[1,1]

tp <- ifelse(mean(perm[,1] > obsg[1])==0,"<0.001",mean(perm[,1] > obsg[1]) )

hist(perm[,1], breaks=40, col="white", border="blueviolet", 
		xlab="Distribution Under Null Hypothesis", main=paste("P(X|NULL)= ",tp,sep=""),
		xlim=range(c(perm[,1],obsg),na.rm=TRUE));
abline(v=obsg[1],lwd=2,col=4)


#### Total variance

obsg <- output2$ES$Intervention1[2,1]

tp <- ifelse(mean(perm[,2] > obsg[1])==0,"<0.001",mean(perm[,2] > obsg[1]) )

hist(perm[,1], breaks=40, col="white", border="blueviolet",
		xlab="Distribution Under Null Hypothesis", 
		main=paste("P(X|NULL)= ",tp,sep=""), xlim=range(c(perm[,1],obsg),
				na.rm=TRUE))
abline(v=obsg[1],lwd=2,col=4)

###############################################
## MLM analysis of cluster randomised trials ##	 
## with bootstrap confidence intervals       ##
###############################################

output3 <- crtFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nBoot=1000,data=iwq)
ES3 <- output3$ES
ES3


ICC3 <- output3$covParm[1]/sum(output3$covParm)
ICC3

#########################################################################
## Bayesian MLM analysis of cluster randomised trials with permutation ##
## under the null hypothesis                                           ##
#########################################################################


output4 <- crtBayes(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nSim=10000,data=iwq)
ES4 <- output4$ES
ES4




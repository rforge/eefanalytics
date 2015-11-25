data(catcht)

###################################################################
## Analysis of simple randomised trials using Hedges Effect Size ##
###################################################################

output1 <- srtFREQ(Posttest~ Intervention+Prettest,
		intervention="Intervention",data=catcht )
ES1 <- output1$ES

###################################################################
## Analysis of simple randomised trials using Hedges Effect Size ## 
## with Permutation p-value                                      ##
###################################################################

output2 <- srtFREQ(Posttest~ Intervention+Prettest,
		intervention="Intervention",nPerm=1000,data=catcht )
ES2 <- output2$ES

perm <- output2$Perm

obsg <- output2$ES[,1]

## calculating Permutation p-value

permPvalue <- ifelse(mean(perm[,1] > obsg[1])==0,"<0.001",
		mean(perm[,1] > obsg[1]) )

## Distribution of ES under H0
hist(perm[,1], breaks=40, col="white", border="blueviolet", 
		xlab="Distribution Under Null Hypothesis", 
		main=paste("P(X|NULL)= ",permPvalue,sep=""),
		xlim=range(c(perm[,1],obsg),na.rm=TRUE));
abline(v=obsg[1],lwd=2,col=4)



###################################################################
## Analysis of simple randomised trials using Hedges Effect Size ##
## with non-parametric bootstrap CI                              ##
###################################################################

output3 <- srtFREQ(Posttest~ Intervention+Prettest,
		intervention="Intervention",nBoot=1000,data=catcht)
ES3 <- output3$ES

## distribution of bootstraped ES
boot <- output3$Bootstrap

## ES given the data
obsg <- output2$ES[,1]

#Bootstrap distribution
hist(boot [,1], breaks=40, col="white", border="blueviolet", 
		xlab="Bootstrap Distribution", main="Bootstrap Distribution", 
		xlim=range(c(boot[,1],obsg),na.rm=TRUE));
abline(v=obsg[1],lwd=2,col=4)


####################################################################
## Analysis of simple randomised trials using Hedges Effect Size. ##
## Schools as fixed effects                                       ##
####################################################################

output4 <- srtFREQ(Posttest~ Intervention+Prettest+as.factor(School),
		intervention="Intervention",data=catcht )
ES4 <- output4$ES

###################################################################
## Analysis of simple randomised trials using Hedges Effect Size ##
## with Permutation p-value. Schools as fixed effects            ##
###################################################################

output5 <- srtFREQ(Posttest~ Intervention+Prettest,
		intervention="Intervention",nPerm=1000,data=catcht )
ES5 <- output5$ES

perm <- output5$Perm

obsg <- output5$ES[,1]

## calculating Permutation p-value

permPvalue <- ifelse(mean(perm[,1] > obsg[1])==0,"<0.001",
		mean(perm[,1] > obsg[1]) )

## Distribution of ES under H0
hist(perm[,1], breaks=40, col="white", border="blueviolet", 
		xlab="Distribution Under Null Hypothesis", main=paste("P(X|NULL)= ",
				permPvalue,sep=""), xlim=range(c(perm[,1],obsg),na.rm=TRUE));
abline(v=obsg[1],lwd=2,col=4)



####################################################################
## Analysis of simple randomised trials using Hedges Effect Size  ##
## with non-parametric bootstrap CI. Schools as fixed effects     ##
####################################################################

output6 <- srtFREQ(Posttest~ Intervention+Prettest,
		intervention="Intervention",nBoot=1000,data=catcht)
ES6 <- output6$ES

## distribution of bootstraped ES
boot <- output6$Bootstrap

## ES given the data
obsg <- output6$ES[,1]

#Bootstrap distribution
hist(boot [,1], breaks=40, col="white", border="blueviolet", 
		xlab="Bootstrap Distribution", main="Bootstrap Distribution", 
		xlim=range(c(boot[,1],obsg),na.rm=TRUE));
abline(v=obsg[1],lwd=2,col=4)

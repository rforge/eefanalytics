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




#### Distribution under the null

perm <- output2$Perm


#### Permutation P-value using total variance


obsg <- output2$ES[1]


permPvalue <- ifelse(mean(perm$permES> obsg[1])==0,"<0.001",
		mean(perm$permES > obsg[1]) )
permPvalue 

## Distribution of ES under H0
hist(perm$permES, breaks=40, col="white", border="blueviolet", 
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
ES5

####################################################################
## Analysis of simple randomised trials using Hedges Effect Size  ##
## with non-parametric bootstrap CI. Schools as fixed effects     ##
####################################################################

output6 <- srtFREQ(Posttest~ Intervention+Prettest,
		intervention="Intervention",nBoot=1000,data=catcht)
ES6 <- output6$ES
ES6 

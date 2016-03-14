data(catcht)

########################################################
## MLM analysis of multisite trials + 1.96SE ##
########################################################

output1 <- mstFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",data=catcht)


### Fixed effects
beta <- output1$Beta
beta

### Effect size
ES1 <- output1$ES$Intervention1
ES1

## Covariance matrix
covParm <- output1$covParm
covParm

### random effects for schools

randOut <- output1$"SchEffects"
randOut <- randOut[order(randOut$Estimate),]
barplot(randOut$Estimate,ylab="Deviations from Overall Average",
		names.arg=randOut$Schools,las=2)

###############################################
## MLM analysis of multisite trials          ##	 
## with bootstrap confidence intervals       ##
###############################################

output2 <- mstFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nBoot=1000,data=catcht)


### Effect size

ES2 <- output2$ES
ES2


#######################################################################
## MLM analysis of mutltisite trials with permutation p-value##
#######################################################################

output3 <- mstFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nPerm=1000,data=catcht)


#### Distribution under the null

perm <- output3$Perm
str(perm )

#### Permutation P-value using total variance


obsg <- output3$ES$Intervention1[2,1]
obsg


p_value <- ifelse(mean(perm$"Intervention1Total" > obsg)==0,"<0.001",
		mean(perm$"Intervention1Total" > obsg) )
p_value

hist(perm$"Intervention1Total", breaks=40, col="white", 
		border="blueviolet",
		xlab="Distribution Under Null Hypothesis", 
		main=paste("P(X|NULL)= ",p_value,sep=""), 
		xlim=range(c(perm$"Intervention1Total",obsg),
		na.rm=TRUE))
abline(v=obsg,lwd=2,col=4)






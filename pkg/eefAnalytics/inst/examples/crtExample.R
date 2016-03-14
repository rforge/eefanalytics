data(iwq)

########################################################
## MLM analysis of cluster randomised trials + 1.96SE ##
########################################################

output1 <- crtFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",data=iwq)


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
## MLM analysis of cluster randomised trials ##	 
## with bootstrap confidence intervals       ##
###############################################

output2 <- crtFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nBoot=1000,data=iwq)


### Effect size

ES2 <- output2$ES
ES2


#######################################################################
## MLM analysis of cluster randomised trials with permutation p-value##
#######################################################################

output3 <- crtFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nPerm=1000,data=iwq)


#### Distribution under the null

perm <- output3$Perm


#### Permutation P-value using total variance


obsg <- output3$ES$Intervention1[2,1]

p_value <- ifelse(mean(perm$"InterventionTotal" > obsg)==0,"<0.001",
		mean(perm$"InterventionTotal" > obsg) )
p_value

hist(perm$"InterventionTotal", breaks=40, col="white", border="blueviolet",
		xlab="Distribution Under Null Hypothesis", 
		main=paste("P(X|NULL)= ",p_value,sep=""), 
		xlim=range(c(perm$"InterventionTotal",obsg), na.rm=TRUE))
abline(v=obsg,lwd=2,col=4)






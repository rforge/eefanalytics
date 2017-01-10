if(interactive()){

data(crtData)

########################################################
## MLM analysis of cluster randomised trials + 1.96SE ##
########################################################

output1 <- crtFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",data=crtData)


### Fixed effects
beta <- output1$Beta
beta

### Effect size
ES1 <- output1$ES$Intervention1
ES1

## Covariance matrix
covParm <- output1$covParm
covParm

### plot random effects for schools

plotObject(analyticObject=output1)

###############################################
## MLM analysis of cluster randomised trials ##	 
## with bootstrap confidence intervals       ##
###############################################

output2 <- crtFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nBoot=1000,data=crtData)


### Effect size

ES2 <- output2$ES
ES2

### plot bootstrapped values

plotObject(analyticObject=output2, group=1)

#######################################################################
## MLM analysis of cluster randomised trials with permutation p-value##
#######################################################################

output3 <- crtFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nPerm=1000,data=crtData)


#### Distribution under the null

perm <- output3$Perm


### plot permutated values

plotObject(analyticObject=output3, group=1)


}



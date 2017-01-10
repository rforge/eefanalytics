if(interactive()){

data(mstData)

########################################################
## MLM analysis of multisite trials + 1.96SE ##
########################################################

output1 <- mstFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",data=mstData)


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
## MLM analysis of multisite trials          ##	 
## with bootstrap confidence intervals       ##
###############################################

output2 <- mstFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nBoot=1000,data=mstData)


### Effect size

ES2 <- output2$ES
ES2

### plot bootstrapped values 

plotObject(analyticObject=output2, group=1)

#######################################################################
## MLM analysis of mutltisite trials with permutation p-value##
#######################################################################

output3 <- mstFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nPerm=1000,data=mstData)


#### plot permutated values 

plotObject(analyticObject=output3, group=1)


}



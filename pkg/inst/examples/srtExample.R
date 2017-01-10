if(interactive()){

data(mstData)

###################################################################
## Analysis of simple randomised trials using Hedges Effect Size ##
###################################################################

output1 <- srtFREQ(Posttest~ Intervention+Prettest,
		intervention="Intervention",data=mstData )
ES1 <- output1$ES

###################################################################
## Analysis of simple randomised trials using Hedges Effect Size ## 
## with Permutation p-value                                      ##
###################################################################

output2 <- srtFREQ(Posttest~ Intervention+Prettest,
		intervention="Intervention",nPerm=1000,data=mstData )
ES2 <- output2$ES




#### plot permutated values

plotObject(analyticObject=output2, group=1)



###################################################################
## Analysis of simple randomised trials using Hedges Effect Size ##
## with non-parametric bootstrap CI                              ##
###################################################################

output3 <- srtFREQ(Posttest~ Intervention+Prettest,
		intervention="Intervention",nBoot=1000,data=mstData)
ES3 <- output3$ES

### plot bootstrapped values

plotObject(analyticObject=output3, group=1)

####################################################################
## Analysis of simple randomised trials using Hedges Effect Size. ##
## Schools as fixed effects                                       ##
####################################################################

output4 <- srtFREQ(Posttest~ Intervention+Prettest+as.factor(School),
		intervention="Intervention",data=mstData )
ES4 <- output4$ES

###################################################################
## Analysis of simple randomised trials using Hedges Effect Size ##
## with Permutation p-value. Schools as fixed effects            ##
###################################################################

output5 <- srtFREQ(Posttest~ Intervention+Prettest,
		intervention="Intervention",nPerm=1000,data=mstData )

ES5 <- output5$ES
ES5

#### plot permutated values

plotObject(analyticObject=output5, group=1)

####################################################################
## Analysis of simple randomised trials using Hedges Effect Size  ##
## with non-parametric bootstrap CI. Schools as fixed effects     ##
####################################################################

output6 <- srtFREQ(Posttest~ Intervention+Prettest,
		intervention="Intervention",nBoot=1000,data=mstData)
ES6 <- output6$ES
ES6 

### plot bootstrapped values

plotObject(analyticObject=output6, group=1)

}

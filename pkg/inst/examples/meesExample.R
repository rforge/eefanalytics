data(iwq)
#########################################################################
## Bayesian MLM analysis of cluster randomised trials with permutation ##
## under the null hypothesis                                           ##
#########################################################################


output4 <- crtBayes(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nSim=10000,data=iwq)

mees.plot(output=output4,x.cord=0.0,y.cord=0.2)


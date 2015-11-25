data(iwq)

#########################################################################
## Bayesian MLM analysis of cluster randomised trials with permutation ##
## under the null hypothesis                                           ##
#########################################################################

date()
output4 <- crtBayes(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nSim=10000,data=iwq)
ES4 <- output4$ES
ES4

date()



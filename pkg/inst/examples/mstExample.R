data(catcht)

#######################################################################
### Analysis of multisite randomised trials using Hedges Effect size ##
#######################################################################

output1 <- mstFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",data=catcht)

ES1 <- output1$ES$Intervention1
ES1

ICC1 <- sum(output1$covParm[-3])/sum(output1$covParm)
ICC1

##############################################################
## Analysis of multisite randomised trials with permutation ##
## under the null hypothesis                                ##
##############################################################


output2 <- mstFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nPerm=1000,data=catcht)

ES2 <- output2$ES$Intervention1
ES2





ICC2 <- sum(output2$covParm[-3])/sum(output2$covParm)
ICC2

##############################################
### Analysis of multisite randomised trials ##
##############################################


output3 <- mstFREQ(Posttest~ Intervention+Prettest,random="School",
		intervention="Intervention",nBoot=1000,data=catcht)

ES3 <- output3$ES$Intervention1
ES3


ICC3 <- sum(output3$covParm[-3])/sum(output3$covParm)
ICC3


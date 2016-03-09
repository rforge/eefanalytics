data(iwq)

########################################################
## Bayesian analysis of cluster randomised trials     ##
########################################################

output <- crtBayes(Posttest~ Intervention+Prettest,
		random="School",intervention="Intervention",
		nSim=10000,data=iwq)

### Fixed effects
beta <- output$Beta
beta

### Effect size
ES1 <- output$ES
ES1

## Covariance matrix
covParm <- output$covParm
covParm

### random effects for schools

randOut <- output$"SchEffects"
randOut <- randOut[order(randOut$Estimate),]
barplot(randOut$Estimate,ylab="Deviations from Overall Average",
		names.arg=randOut$Schools,las=2)



### Posterior probability given a fixed threshold
probES <- output$ProbES
str(probES )

			
plot(probES[,1] ,probES[,2],ylim=c(0,max(probES)),
		ylab="Probability",cex.lab=1,cex.axis=1,
		type="n", xlab=expression("Effect size" >= "x"),
		cex=1)
lines(probES[,1],probES[,2],col="chartreuse3",cex=1.5,
		lwd=1.5,lty=2)
lines(probES[,1],probES[,3],col="violetred",cex=1.5,
		lwd=1.5,lty=3)
lines(probES[,1],probES[,4],col="cornflowerblue",cex=1.5,
		lwd=1.5,lty=1)
points(probES[,1],probES[,2],col="chartreuse3",cex=1.5,
		lwd=1.5,pch=7)
points(probES[,1],probES[,3],col="violetred",cex=1.5,
		lwd=1.5,pch=1)
points(probES[,1],probES[,4],col="cornflowerblue",
		cex=1.5,lwd=1.5,pch=12)
legend(0,0.4,legend=c("Within ","Between ","Total "),
		lty=c(2,3,1),cex=1.5, pch=c(7,1,12),
		col=c("chartreuse3","violetred","cornflowerblue"),
		title="Variance Type")





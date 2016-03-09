
data(catcht)

############# SRT
caceOutput3<- caceSRTBoot(Posttest~ Prettest+ Intervention,
			intervention="Intervention",
			compliance = "Percentage_Attendance",
			nBoot=1000,data=catcht)

### visualising CACE ES

require(metafor)
forest(x=caceOutput3$CACE$ES, ci.lb=caceOutput3$CACE$LB, 
		ci.ub=caceOutput3$CACE$UB, xlab="CACE ES",
		slab=as.character(caceOutput3$CACE$Compliance), 
		xlim=c(-1,2), alim=c(-1,7), cex=1,lwd=1.5)


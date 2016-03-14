data(iwq)
############# CRT
caceOutput<- caceCRTBoot(Posttest~ Prettest+ Intervention,
			random="School",intervention="Intervention",
			compliance = "Percentage_Attendance",
			nBoot=1000,data=iwq)

### visualising CACE ES

require(metafor)
forest(x=caceOutput$CACE$ES, ci.lb=caceOutput$CACE$LB, 
		ci.ub=caceOutput$CACE$UB, xlab="CACE ES",
		slab=as.character(caceOutput$CACE$Compliance), 
		xlim=c(-1,2), alim=c(-1,3.75), cex=1,lwd=1.5)



############# MST
caceOutput2<- caceMSTBoot(Posttest~ Prettest+ Intervention,
			random="School",intervention="Intervention",
			compliance = "Percentage_Attendance",
			nBoot=1000,data=catcht)

### visualising CACE ES

require(metafor)
forest(x=caceOutput2$CACE$ES, ci.lb=caceOutput2$CACE$LB, 
		ci.ub=caceOutput2$CACE$UB, xlab="CACE ES",
		slab=as.character(caceOutput2$CACE$Compliance),
		xlim=c(-1,2), alim=c(-1,7), cex=1,lwd=1.5)



# write.csv(caceOutput2$CACE,file="D:\\eefAnalytics\\wd\\example.csv",row.names=FALSE)

if(interactive()){

data(mstData)
############# SRT
caceOutput3<- caceSRTBoot(Posttest~ Prettest+ Intervention,
			intervention="Intervention",
			compliance = "Percentage_Attendance",nBoot=1000,data=mstData)

cace <- caceOutput3$CACE
cace

Complier <- caceOutput3$Compliers
Complier 

### visualising CACE effect size

plotObject(analyticObject=ccaceOutput3)

}

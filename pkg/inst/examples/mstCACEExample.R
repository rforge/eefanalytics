if(interactive()){

data(mstData)
############# MST
caceOutput2<- caceMSTBoot(Posttest~ Prettest+ Intervention,
			random="School",intervention="Intervention",compliance = "Percentage_Attendance",nBoot=1000,data=mstData)

cace <- caceOutput2$CACE
cace

Complier <- caceOutput2$Compliers
Complier 

### visualising CACE effect size

plotObject(analyticObject=caceOutput2)


}


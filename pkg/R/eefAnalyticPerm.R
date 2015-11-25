
#' Multisite trial data.
#'
#' A multisite education trial data containing 54 schools.
#'
#' \itemize{
#'   \item Posttest. posttest scores
#'   \item Prettest. prettests scores
#'   \item Intervention. indicator for intervention groups 
#'coded as 1 for the intervention and 0 for the control group.
#'   \item School. numeric school identifier
#' }
#'
#' @format A data frame with 210 rows and 4 variables
#' @name catcht
NULL

###iwq

#' Cluster randomised trial data.
#'
#' A cluster randomised trial data containing 22 schools.
#'
#' \itemize{
#'   \item Posttest. posttest scores
#'   \item Prettest. prettests scores
#'   \item Intervention. indicator for intervention groups 
#'coded as 1 for the intervention and 0 for the control group.
#'   \item School. numeric school identifier
#' }
#'
#' @format A data frame with 265 rows and 4 variables
#' @name iwq
NULL


## IMPORTS ##
#' @importFrom lme4 lmer ranef
#' @importFrom geoR rinvchisq
#' @importFrom mvtnorm rmvnorm
NULL

#############################################################################
############# SRT main functions ################################################


#' Analysis of simple randomised trials using linear model.
#' 
#' \code{srtFREQ} performs analysis of education trials assuming simple randomisation across schools.
#' 
#' @export
#' @param formula specifies in the model using form posttest ~ pretests+Intervention+....
#' @param intervention specifies the intervention variable in the formula.
#' @param nBoot specifies number of bootstrap samples for non-parametric bootstrapping. Default is NULL. 
#' @param nPerm specifies number of perumation under the null hypothesis.  Default is NULL. 
#' @param data specifies data frame containing the data to be analysed. 
#' @return S3 \code{mcpi} object; a list consisting of
#' \itemize{
#' \item \code{Beta}. The parameter estimates from linear model with the associatie confidence intervals.
#' \item \code{ES}. Estimate Hedges effect size. It is based on standard error (SE) if \code{nBoot=NULL}. Otherwise, it is based on non-parametric bootstrap confidence intervals.
#' \item \code{sigma2}. The residual variance.
#' \item \code{Perm}. A vector containing the distribution of effect size under the null hypothesis. It is only produced if \code{nPerm} is specified.
#' \item \code{Bootstrap}. A vector containing the bootstrap effect sizes. It is only prduced is \code{nBoot} is specified.
#' }
#' @example inst/examples/srtExample.R
srtFREQ <- function(formula,intervention,nBoot=NULL,nPerm=NULL,data)UseMethod("srtFREQ")

#' @export
srtFREQ.formula <- function(formula,intervention,nBoot=NULL,nPerm=NULL,data=data){
	
	if(!is.null(nPerm) & !is.null(nBoot)){stop("Either nPerm or nBoot must be unspecified")}
	if(is.null(nPerm)){nPerm <-0}
	if(is.null(nBoot)){nBoot <-0}
	tmp3 <- which(colnames(data)==intervention)
	data[,tmp3] <- as.factor(data[,tmp3])
	mf <- model.frame(formula=formula, data=data)
	fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
	tmp <- colnames(fixedDesignMatrix )
	tmp[1]  <- "Intercept"
	colnames(fixedDesignMatrix)<- tmp
	posttest <- model.response(mf)
      intervention <- intervention
	trt <- data[,which(colnames(data)==intervention)]

	if(length(tmp3)!= 1){stop("Intervention variable misspecified")}

	output <- srt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt)	


	if(nPerm>0){

		if(nPerm<1000){stop("nPerm must be greater than 1000")}

		g <- matrix(NA,nPerm,(length(unique(trt))-1))
		
		set.seed(1020252)
		for (i in 1:nPerm){
				
			data[,which(colnames(data)==intervention)]<- sample(trt)
			fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))

			p2CRTFREQ <-srt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt)	

			g[i,]  <-  p2CRTFREQ$ES[,1]
		}

		output$Perm <-g
	}

	if(nBoot > 0){

		if(nBoot<1000){stop("nBoot must be greater than 1000")}


		tid <- c(1:nrow(fixedDesignMatrix))

		set.seed(1020252)
		bootSamples <- sapply(c(1:nBoot),function(x)sample(tid,replace=TRUE))

		bootResults <- data.frame(apply(bootSamples ,2,function(bt)srt.srt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,bt=bt)))
		bootResults2 <- t(matrix(unlist(bootResults ),(length(unique(trt))-1),nBoot))
		bootES <- apply(bootResults2,2,function(x)quantile(x,prob=c(0.025,0.975)))
		bootES2 <- t(rbind(output$ES[,1],bootES ))
		row.names(bootES2)<- row.names(output$ES)
		colnames(bootES2)<- colnames(output$ES)
		colnames(bootResults2)<- row.names(output$ES)
		output$ES <-round(bootES2,2)
		output$Bootstrap <- bootResults2  

	}

	return(output) 
}



#############################################################################
############# CRT main functions ################################################
#' Analysis of Cluster Randomised Trials using MLM.
#' 
#' \code{crtFREQ} performs analysis of cluster randomised trials.
#' 
#' @export
#' @param formula specifies in the model using form posttest ~ pretests+Intervention +....
#' @param random a string variable specifying the "clustering"  as contained in the data. 
#' @param intervention specifies the intervention variable in the formula.
#' @param nBoot specifies number of bootstrap samples for non-parametric bootstrapping. Default is NULL. 
#' @param nPerm specifies number of perumation under the null hypothesis.  Default is NULL. 
#' @param data specifies data frame containing the data to be analysed. 
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}. The parameter estimates from linear model with the associatie confidence intervals.
#' \item \code{ES}. A matric of  Hedges effect size based on within variance and total variance. It is base is based on standard error (SE) if \code{nBoot=NULL}. Otherwise, it is based on non-parametric bootstrap confidence intervals.
#' \item \code{covParm}. Vector of variance decomposition into between-variance (Schools) and within-variance (Pupils).
#' \item \code{SchEffects}. Estimated random effects for each school.
#' \item \code{Perm}. A matrix containing the distribution of effect size under the null hypothesis. It is only produced if \code{nPerm} is specified. The first columns contains effect size based on within variance and the second column contains effect size based on total variance.
#' \item \code{Bootstrap}. A matrix containing the bootstrap effect sizes. It is only prduced is \code{nBoot} is specified. The first columns contains effect size based on within variance and the second column contains effect size based on total variance.
#' }
#' @example inst/examples/crtExample.R
crtFREQ<- function(formula,random,intervention,nPerm=NULL,nBoot=NULL,data)UseMethod("crtFREQ")

#' @export
crtFREQ.formula <- function(formula,random,intervention,nPerm=NULL,nBoot=NULL,data){
  
	data <- data[order(data[,which(colnames(data)==random)]),]
    
	intervention <- intervention
	trt <- data[,which(colnames(data)==intervention)]
	tmp2 <- which(colnames(data)==random)
	cluster2 <-  data[,tmp2]

	chk <- sum(rowSums(table(cluster2,trt)!=0)>1)
	if(chk >0){stop("This not a CRT design")}
 	stp <- as.character(row.names(table(cluster2,trt)))
	stp2 <- as.numeric(apply(table(cluster2,trt),1,function(x)colnames(table(cluster2,trt))[x!=0]))

	if(!is.null(nPerm) & !is.null(nBoot)){stop("Either nPerm or nBoot must be unspecified")}
	if(is.null(nPerm)){nPerm <-0}
	if(is.null(nBoot)){nBoot <-0}
	tmp3 <- which(colnames(data)==intervention)
	data[,tmp3] <- as.factor(data[,tmp3])

	mf <- model.frame(formula=formula, data=data)
	mf <- mf[order(cluster2),]
	cluster <- cluster2[order(cluster2)]
	trt <- trt[order(cluster2)]
	fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
	tmp <- colnames(fixedDesignMatrix )
	tmp[1]  <- "Intercept"
	colnames(fixedDesignMatrix)<- tmp
	posttest <- model.response(mf)


	if(length(tmp2)!= 1){stop("Clustering variable misspecified")}
	if(length(tmp3)!= 1){stop("Intervention variable misspecified")}


	output <- crt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster)	


	if(nPerm>0){
		if(nPerm<999){stop("nPerm must be atleast 999")}
		output$Perm<- crt.perm(formula,data,stp,stp2,intervention,cluster,nPerm,random)
	}


	if(nBoot >0){


		tid <- c(1:nrow(fixedDesignMatrix))

		set.seed(1020252)
		bootSamples <- NULL

		for(ii in 1:length(unique(cluster))){
			selID <- tid[cluster==unique(cluster)[ii]]
			if(length(selID)>0){
				selID2<- sapply(c(1:nBoot),function(x)sample(selID,length(selID),replace=TRUE))
				bootSamples <- rbind(bootSamples ,selID2)
			}
				
		}

		
		bootResults <- apply(bootSamples ,2,function(bt)crt.crt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster,bt=bt))
		bootES <- bootCompile(output=output,trt=trt,bootResults=bootResults)
		output$ES <- bootES
		output$Bootstrap <- bootResults  
	}

	return(output) 
}




#############################################################################
############# MST main functions ################################################

#############################################################################
############# CRT main functions ################################################
#' Analysis of multisite andomised trials using multilevel model.
#' 
#' \code{crtFREQ} performs analysis of cluster randomised trials.
#' 
#' @export
#' @param formula specifies in the model using form posttest ~ pretests+Intervention +....
#' @param random a string variable specifying the "clustering"  as contained in the data. 
#' @param intervention specifies the intervention variable in the formula.
#' @param nBoot specifies number of bootstrap samples for non-parametric bootstrapping. Default is NULL. 
#' @param nPerm specifies number of perumation under the null hypothesis.  Default is NULL. 
#' @param data specifies data frame containing the data to be analysed. 
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}. The parameter estimates from linear model with the associatie confidence intervals.
#' \item \code{ES}. A matric of  Hedges effect size based on within variance and total variance. It is base is based on standard error (SE) if \code{nBoot=NULL}. Otherwise, it is based on non-parametric bootstrap confidence intervals.
#' \item \code{covParm}. Vector of variance decomposition into between-variance (Schools) and within-variance (Pupils).
#' \item \code{SchEffects}. Estimated random effects for each school.
#' \item \code{Perm}. A matrix containing the distribution of effect size under the null hypothesis. It is only produced if \code{nPerm} is specified. The first columns contains effect size based on within variance and the second column contains effect size based on total variance.
#' \item \code{Bootstrap}. A matrix containing the bootstrap effect sizes. It is only prduced is \code{nBoot} is specified. The first columns contains effect size based on within variance and the second column contains effect size based on total variance.
#' }
#' @example inst/examples/mstExample.R
mstFREQ<- function(formula,random,intervention,nPerm=NULL,data,nBoot=NULL)UseMethod("mstFREQ")

#' @export
mstFREQ.formula <- function(formula,random,intervention,nPerm=NULL,data,nBoot=NULL){

	data <- data[order(data[,which(colnames(data)==random)],data[,which(colnames(data)==intervention)]),]
	intervention <- intervention
	trt <- data[,which(colnames(data)==intervention)]
	tmp2 <- which(colnames(data)==random)
	cluster2 = data[,tmp2]

	chk <- sum(rowSums(table(cluster2,trt)!=0)>1)
	if(chk ==0){stop("This not an MST design")}
 	
	if(!is.null(nPerm) & !is.null(nBoot)){stop("Either nPerm or nBoot must be unspecified")}
	if(is.null(nPerm)){nPerm <-0}
	if(is.null(nBoot)){nBoot <-0}


	tmp3 <- which(colnames(data)==intervention)
	data[,tmp3] <- as.factor(data[,tmp3])
	mf <- model.frame(formula=formula, data=data)
	mf <- mf[order(cluster2),]
	cluster <- cluster2[order(cluster2)]
	trt <- trt[order(cluster2)]
	fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
	tmp <- colnames(fixedDesignMatrix )
	tmp[1]  <- "Intercept"
	colnames(fixedDesignMatrix)<- tmp
	posttest <- model.response(mf)
      intervention <- intervention
	trt <- data[,which(colnames(data)==intervention)]

	if(length(tmp2)!= 1){stop("Clustering variable misspecified")}
	if(length(tmp3)!= 1){stop("Intervention variable misspecified")}


	output <- rbd(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt,cluster=cluster)	

	if(nPerm > 0){

		if(nPerm<999){stop("nPerm must be atleast 999")}
		output$Perm <- mst.perm(formula,data,trt,intervention,nPerm,random,cluster)

	}

	if(nBoot>0){


		tid <- c(1:nrow(fixedDesignMatrix))

		set.seed(1020252)
		bootSamples <- NULL

		for(ii in 1:length(unique(cluster))){
			selID <- tid[cluster==unique(cluster)[ii]]
			if(length(selID)>0){
				selID2<- sapply(c(1:nBoot),function(x)sample(selID,length(selID),replace=TRUE))
				bootSamples <- rbind(bootSamples ,selID2)
			}
				
		}
	
		bootResults <- apply(bootSamples ,2,function(bt)rbd.rbd(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt,cluster=cluster,bt=bt))

		bootES <- bootCompile(output=output,trt=trt,bootResults=bootResults)
		output$ES <- bootES
		output$Bootstrap <- bootResults  

	}

	return(output) 
}

###################################################################################################
############# Bayesian analysis of cluster randomised trials
##################################################################################################

#' Bayesia analysis of cluster randomised trials Using vague priors.
#' 
#' \code{crtBayes} performs Bayesian analysis of clustered randomised trials using vague priors.
#' 
#' @export
#' @param formula specifies in the model using form posttest ~ pretests+Intervention+....
#' @param random a string variable specifying the "clustering"  as contained in the data.
#' @param intervention specifies the intervention variable in the formula.
#' @param nSim specifies number of MCMC iterations. A minimum of 10,000 is recommended.
#'
#' @param data specifies data frame containing the data to be analysed. 
#' @return S3 \code{mcpi} object; a list consisting of
#' \itemize{
#' \item \code{Beta}. The parameter estimates from linear model with the associatie confidence intervals.
#' \item \code{ES}. A matric of  Hedges effect size based on within variance and total variance. It is base is based on standard error (SE) if \code{nBoot=NULL}. Otherwise, it is based on non-parametric bootstrap confidence intervals.
#' \item \code{covParm}. Vector of variance decomposition into between-variance (Schools) and within-variance (Pupils).
#' \item \code{SchEffects}. Estimated random effects for each school.
#' \item \code{ProbES}. A maxtrix containing the probability of observing ES greater than a pre-specified value. First column is for within-variance, second column for between-variance and the third column for total-variance.
#' }
#' @example inst/examples/bayesExample.R
crtBayes <- function(formula,random,intervention,nSim,data)UseMethod("crtBayes")

#' @export
crtBayes.formula <- function(formula,random,intervention,nSim=nSim,data){
	tmp3 <- which(colnames(data)==intervention)
	data[,tmp3] <- as.factor(data[,tmp3])

	mf <- model.frame(formula=formula, data=data)
	fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
	tmp <- colnames(fixedDesignMatrix )
	tmp[1]  <- "Intercept"
	colnames(fixedDesignMatrix)<- tmp
	posttest <- model.response(mf)
      intervention <- intervention
	trt <- data[,which(colnames(data)==intervention)]
	tmp2 <- which(colnames(data)==random)
	tmp3 <- which(colnames(data)==intervention)
	cluster = as.factor(as.character(data[,tmp2]))
	nsim=nSim
	if(length(tmp2)!= 1){stop("Clustering variable misspecified")}
	if(length(tmp3)!= 1){stop("Intervention variable misspecified")}

	if(nsim < 10000){stop("nsim >= 10000 is recommended")}

	BayesOutput <- erantBAYES(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster, nsim=nsim)
	output  <- errantSummary(bayesObject=BayesOutput,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention)

	return(output) 
}




#############################################################################
############# Bayesian functions ################################################

## Gibbs sampling - internal
erantBAYES<- function(posttest,fixedDesignMatrix,intervention,cluster, nsim) {
      
	cluster <- as.factor(cluster)   
	randomDesignMatrix <- as.matrix(as.data.frame(lm(posttest~cluster-1,x=TRUE)$x))
	nschools <- ncol(randomDesignMatrix )
	npar <- ncol(fixedDesignMatrix)
	totalObservations <- nrow(randomDesignMatrix )
	
	lgrp <- nchar(paste(intervention))
	sgrp <- substring(colnames(fixedDesignMatrix),1,lgrp)
	ssgrp <- colnames(fixedDesignMatrix)[sgrp==intervention]
	iindex <- which(sgrp==intervention)	

	ncluster <- tapply(cluster,cluster,length)

	###########
  	# Initial values  #
  	###########
	beta <- summary(lm(posttest~fixedDesignMatrix))$coef[,1]
	sigma2.ui <- 1
	b <- rnorm(nschools,0,1)
	bmat<-matrix(b,ncol=1,byrow=T)

  	#################
  	# Output matrix #
  	#################

  	nused <- floor(nsim/2)
  	nburnin <-  nsim-nused
  	outputbik <- matrix(0,nused,nschools)
  	outputCov <- matrix(0,nused,4)
      outputES <- matrix(0,nused,3*length(ssgrp) )
      outputBeta <- matrix(0,nused,npar)

	for (i in 1:nsim){

    		set.seed(123*i+i)


		# Update sigma2.e
		
		Ve <- totalObservations
		s2e0 <- fixedDesignMatrix %*% t(t(beta))
            s2e1 <-  t(matrix(bmat,nrow=nschools,ncol=totalObservations))
   		s2e2 <- randomDesignMatrix*s2e1 
           	s2e3 <- t(t(posttest))-s2e0 -rowSums(s2e2) 
		s2e <- (t(s2e3)%*%s2e3)/Ve
   		sigma2.e<-rinvchisq(1, Ve, s2e)


 		# Update sigma2.ui
		Vui <- nschools 
		s2ui <- (t(bmat)%*%bmat)/Vui
            sigma2.ui<-rinvchisq(1,Vui,s2ui)


 		# Update b
		usim1 <- t(t(posttest))-s2e0
  		usim2<- tapply(usim1,cluster,mean)
		uV0 <- (ncluster +  as.numeric(sigma2.e/sigma2.ui))^-1
		uV <-  diag(as.numeric(sigma2.e)*uV0 )
            ui <- ncluster*uV0*usim2
   		b <-rmvnorm(1,ui,uV)
   		bmat<-matrix(b,ncol=1,byrow=T)

		# Update Beta 
            beta1 <-solve(t(fixedDesignMatrix) %*% fixedDesignMatrix )
		beta2 <-  t(matrix(bmat,nrow=nschools,ncol=totalObservations))
   		beta3 <- randomDesignMatrix*beta2
		beta4 <- t(t(posttest)) -rowSums(beta3)
            beta5 <- t(fixedDesignMatrix)%*% beta4
		beta6 <- beta1 %*% beta5 
   		vbeta<- beta1*as.numeric(sigma2.e )
   		beta<-c(rmvnorm(1,beta6,vbeta))
            
            treatment=beta[iindex]


   		if(i > nburnin ){ 
  			outputbik[(i-nburnin),] <- as.numeric(b)
  			outputCov[(i-nburnin),] <- c(sigma2.e,sigma2.ui,sigma2.e+sigma2.ui,(sigma2.ui)/(sigma2.e+sigma2.ui))
  			outputBeta[(i-nburnin),] <- as.numeric(beta)
  			outputES[(i-nburnin),] <- c(sapply(treatment,function(x)x/sqrt(c(sigma2.e,(sigma2.ui),sigma2.e+sigma2.ui))))
		}
	}


      esF <- c("Within","Between","Total")

      es.names <- c(sapply(ssgrp,function(x)paste(x,esF,sep="") ))
      colnames(outputBeta)<- colnames(fixedDesignMatrix)
	colnames(outputCov) <- c("Within","Between","Total","ICC")
	colnames(outputES) <- es.names
	colnames(outputbik) <- unique(cluster )

 	Output <- list(randomEffects=outputbik,CovParameters=outputCov,Beta=outputBeta,ES=outputES)
  	return(Output)
}


## summarise covariance parameters - internal

covSummary <- function(bayesObject){

	covParm <- bayesObject$CovParameters
      covParm2 <- colMeans(covParm) # remove all roundings that appear in the middle of functions 
      covParm3 <- t(apply(covParm,2,function(x)quantile(x,prob=c(0.025,0.975))))
	covParm4 <- data.frame(cbind(covParm2,covParm3))
	covParm5 <- sqrt(covParm4)
	colnames(covParm4) <- c("Variance","95% LB","95% UB")
	rownames(covParm4) <- c("Pupils","Schools","Total","ICC")
	covParm5 <- covParm4[,1]
	return(covParm5)
}


## summarise beta parameters - internal
betaSummary <- function(bayesObject){

	betaParm <- bayesObject$Beta
      betaParm2 <- colMeans(betaParm)
      betaParm3 <- t(apply(betaParm,2,function(x)quantile(x,prob=c(0.025,0.975))))
	betaParm4 <- data.frame(cbind(betaParm2,betaParm3))
	colnames(betaParm4) <- c("Estimate","95% LB","95% UB")
	rownames(betaParm4) <- colnames(betaParm)
	return(betaParm4)
}

## summarise Effect sizes - internal
esSummary <- function(bayesObject,fixedDesignMatrix,intervention){

	esParm <- bayesObject$ES
      esParm2 <- colMeans(esParm)
      esParm3 <- t(apply(esParm,2,function(x)quantile(x,prob=c(0.025,0.975))))
	esParm4 <- data.frame(cbind(esParm2,esParm3))
	colnames(esParm4) <- c("Estimate","95% LB","95% UB")
	btp <- nrow(esParm4)
      if(btp <=3){return(round(esParm4,2))}

	if(btp >3){
		btp2 <- seq(3,btp,3)
		btp3 <- seq(1,btp,3)
		ouptut <- list()
		for(i in 1:length(btp2)){
			ouptut[[i]] <- round(esParm4[btp3[i]:btp2[i],],2)
		}
	rname <- colnames(fixedDesignMatrix)[substring(colnames(fixedDesignMatrix),1,nchar(intervention))==intervention]
	names(ouptut)<- rname
	return(ouptut)
	}
	
}

## summarise random effects - internal

schSummary <- function(bayesObject){

	schParm <- bayesObject$randomEffects
      schParm2 <- colMeans(schParm)
      schParm3 <- t(apply(schParm,2,function(x)quantile(x,prob=c(0.025,0.975))))
	schParm4 <- data.frame(cbind(schParm2,schParm3))
	colnames(schParm4) <- c("Estimate","95% LB","95% UB")
	schParm5<- schParm4[order(abs(schParm4$Estimate),decreasing=TRUE),]
	return(schParm5)

}


##summarise minimum expected effect size - internal
esProb <- function(bayesObject,esOutput){
      es <- seq(0,1,0.1)
	esParm <- bayesObject$ES 
      esParm2 <- sapply(es,function(x)colMeans(esParm>=x))
	esParm4 <- data.frame(cbind(ES=es,t(esParm2)))
	rownames(esParm4) <- NULL
	return(esParm4)
	
}


## plots of minimum expected effect size - external
#############################################################################
############# CRT main functions ################################################
#' Minimum expected effect size plot.
#' 
#' \code{crtFREQ} performs analysis of cluster randomised trials.
#' 
#' @export
#' @param output \code{crtBayes} outout object
#' @param x.cord position of legend on x-axis.
#' @param y.cord position of legend on y-axis.
#' @return A plot of minimum expected effect size.
#' @example inst/examples/meesExample.R
mees.plot <- function(output,x.cord=0.8,y.cord=max(esParm2)-.05){
	esOutput<-output$ES
	tmp <- output$ProbES
      es <- tmp$ES
      esParm2 <- t(tmp[,-1])
	ngrp <- nrow(esParm2)/3
      grpname <- names(esOutput)
	if(ngrp > 1){
		par(mfrow=c(1,ngrp))
		for(k in 1:ngrp){
			plot(es,esParm2[(k*3-2),],ylim=c(0,max(esParm2)),ylab="Probability",main=paste(grpname[k]),cex.lab=1,cex.axis=1,type="n", xlab="ES > X",cex=1)
			lines(es,esParm2[(k*3-2),],col="chartreuse3",cex=1.5,lwd=1.5,lty=2)
			lines(es,esParm2[(k*3-1),],col="violetred",cex=1.5,lwd=1.5,lty=3)
			lines(es,esParm2[(k*3-0),],col="cornflowerblue",cex=1.5,lwd=1.5,lty=1)
			points(es,esParm2[(k*3-2),],col="chartreuse3",cex=1.5,lwd=1.5,pch=7)
			points(es,esParm2[(k*3-1),],col="violetred",cex=1.5,lwd=1.5,pch=1)
			points(es,esParm2[(k*3-0),],col="cornflowerblue",cex=1.5,lwd=1.5,pch=12)     		
 		}

	}
	if(ngrp == 1){
		
		plot(es,esParm2[1,],ylim=c(0,max(esParm2)),ylab="Probability",cex.lab=1,cex.axis=1,type="n", xlab=expression("Effect size" >= "x"),cex=1)
		lines(es,esParm2[1,],col="chartreuse3",cex=1.5,lwd=1.5,lty=2)
		lines(es,esParm2[2,],col="violetred",cex=1.5,lwd=1.5,lty=3)
		lines(es,esParm2[3,],col="cornflowerblue",cex=1.5,lwd=1.5,lty=1)
		points(es,esParm2[1,],col="chartreuse3",cex=1.5,lwd=1.5,pch=7)
		points(es,esParm2[2,],col="violetred",cex=1.5,lwd=1.5,pch=1)
		points(es,esParm2[3,],col="cornflowerblue",cex=1.5,lwd=1.5,pch=12)
      	legend(x.cord,y.cord,legend=c("Within ","Between ","Total "),lty=c(2,3,1), cex=.8,
		pch=c(7,1,12),col=c("chartreuse3","violetred","cornflowerblue"),title="Variance")

	}
	
}

##summarise all bayesian parameters - internal
errantSummary <- function(bayesObject,fixedDesignMatrix,intervention){
	covValues <- covSummary(bayesObject=bayesObject)
	covValues <- covValues[c(2,1,3,4)]
	betaValues <- betaSummary(bayesObject=bayesObject)
	esValues <- esSummary(bayesObject,fixedDesignMatrix,intervention)
	schValues <-schSummary(bayesObject=bayesObject)
	es.prob <- esProb(bayesObject=bayesObject,esOutput=esValues)
      output <- list(Beta=betaValues,covParm= covValues,ES=esValues,ProbES=es.prob,SchEffects=schValues)
}




#############################################################################
############# Frequentist functions ################################################

## lmer - internal

#erantFREQ <- function(posttest,fixedDesignMatrix,intervention,trt,cluster,hedges,multiSite,nboot){
#
#
#
#	if(multiSite==TRUE){
#
#		output <- rbd(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt,cluster=cluster)	
#	
#		bootResults <- apply(bootSamples ,2,function(bt)rbd.rbd(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt,cluster=cluster,bt=bt))
#
#		bootES <- bootCompile(output=output,trt=trt,bootResults=bootResults)
#		output$ES <- bootES
#	}
#
#	
#	if(multiSite !=TRUE&hedges==TRUE){
#
#			output <- crt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster)	
#	
#	}
#
#	if(multiSite !=TRUE&hedges!=TRUE){
#
#		output <- crt.boot(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster)	
#
#		bootResults <- apply(bootSamples ,2,function(bt)crt.crt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster,bt=bt))
#		bootES <- bootCompile(output=output,trt=trt,bootResults=bootResults)
#		output$ES <- bootES
#	}
#
#
#	return(output)
#}


### compile bootstrap results - internal
#bootCompile2 <- function(output,trt,bootResults){
#
#
#
#	withinCI <- apply(bootResults,2,function(x)quantile(x,prob=c(0.025,0.975)))
#	
#	tmpES <- list()
#	for(kk in 1:length(output$ES)){
#		tmp1 <- rbind(withinCI[,kk], TotalCI[,kk])
#		tmp2 <- cbind(output$ES[[kk]][,1],tmp1)
#		colnames(tmp2 )<- c("Estimate","95% LB","95% UB")
#		row.names(tmp2)	<- c("Within","Total")
#		tmpES[[kk]] <- round(tmp2,2)
#	}
#	names(tmpES) <- names(output$ES)
#	return(tmpES)
#
#}


## compile bootstrap results - internal
bootCompile <- function(output,trt,bootResults){

	withinBoot <- matrix(NA,nrow=length(bootResults),ncol=(length(unique(trt ))-1))
	totalBoot <- matrix(NA,nrow=length(bootResults),ncol=(length(unique(trt ))-1))
	for(k in 1:length(bootResults)){
		tmp <- bootResults[[k]]
		tmpR <- NULL
		for(j in 1:length(tmp)){
					
			tmpR  <- c(tmpR,tmp[[j]][,1])
		}
		withinBoot[k,] <- tmpR[seq(1,2*length(tmp),2)]
		totalBoot[k,] <- tmpR[seq(2,2*length(tmp),2)]
	}

	withinCI <- apply(withinBoot,2,function(x)quantile(x,prob=c(0.025,0.975)))
	TotalCI <- apply(totalBoot,2,function(x)quantile(x,prob=c(0.025,0.975)))
	tmpES <- list()
	for(kk in 1:length(output$ES)){
		tmp1 <- rbind(withinCI[,kk], TotalCI[,kk])
		tmp2 <- cbind(output$ES[[kk]][,1],tmp1)
		colnames(tmp2 )<- c("Estimate","95% LB","95% UB")
		row.names(tmp2)	<- c("Within","Total")
		tmpES[[kk]] <- round(tmp2,2)
	}
	names(tmpES) <- names(output$ES)
	return(tmpES)

}

## random intercept model - internal
crt <- function(posttest,fixedDesignMatrix,intervention,cluster){
	
	freqFit <- lmer(posttest~ fixedDesignMatrix-1+ (1|cluster))
		
      cit <- try(confint(freqFit),silent=TRUE)
	betaB <- data.frame(cbind(summary(freqFit)$coefficients[,1],cit[-c(1,2),]))
	row.names(betaB)<- colnames(fixedDesignMatrix)
	colnames(betaB) <- c("Estimate","95% LB ","95% UB")
	betaB <- betaB
	var.B<- as.numeric(summary(freqFit)$varcor)
	var.W<- summary(freqFit)$sigma^2
	var.tt <- var.W+var.B
	ICC <- var.B/var.tt
	sigmaBE <- c(var.B,var.W)
	sigmaBE <- sigmaBE
	names(sigmaBE)<- c("Schools","Pupils")
      schRand <- data.frame(ranef(freqFit)$cluster)
	names(schRand)<- "Estimate"
	btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)
	      
      output2 <- list()
	for( i in 1:length(btp )){
	
		beta <- betaB[btp[i],1]
      	group <- fixedDesignMatrix[,btp[i]]

  		esWithin <- g.within(var.w=var.W, beta=beta, icc=ICC , group=group, schoolID=cluster)
  		esTotal <- g.total(var.tt=var.tt, beta=beta, icc=ICC , group=group, schoolID=cluster)
   	
		output1 <- data.frame(rbind(esWithin,esTotal))
		colnames(output1) <- c("Estimate","95% LB","95% UB")
		rownames(output1) <- c("Within","Total")
		output2[[i]] <- round(output1,2)
	}
	names(output2) <- row.names(betaB)[btp]
      output <- list(Beta=round(betaB,2),covParm=round(sigmaBE,2),ES=output2,SchEffects=round(schRand,2))

	return(output)

}

crtP <- function(posttest,fixedDesignMatrix,intervention,cluster){
	
	freqFit <- lmer(posttest~ fixedDesignMatrix-1+ (1|cluster))
		
	betaB <- data.frame(summary(freqFit)$coefficients[,1])
	row.names(betaB)<- colnames(fixedDesignMatrix)
	#colnames(betaB) <- c("Estimate","95% LB ","95% UB")
	betaB <- betaB
	var.B<- as.numeric(summary(freqFit)$varcor)
	var.W<- summary(freqFit)$sigma^2
	var.tt <- var.W+var.B
	ICC <- var.B/var.tt
	sigmaBE <- c(var.B,var.W)
	sigmaBE <- sigmaBE
	btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)
	      
      output2 <- list()
	for( i in 1:length(btp )){
	
		beta <- betaB[btp[i],1]
      	group <- fixedDesignMatrix[,btp[i]]

  		esWithin <- g.within(var.w=var.W, beta=beta, icc=ICC , group=group, schoolID=cluster)
  		esTotal <- g.total(var.tt=var.tt, beta=beta, icc=ICC , group=group, schoolID=cluster)
   	
		output1 <- data.frame(rbind(esWithin,esTotal))
		colnames(output1) <- c("Estimate","95% LB","95% UB")
		rownames(output1) <- c("Within","Total")
		output2[[i]] <- round(output1,2)
	}
	names(output2) <- row.names(betaB)[btp]
      output <- list(ES=output2)

	return(output)

}


#### crt.perm

crt.perm <- function(formula,data,stp,stp2,intervention,cluster,nPerm,random){

		data2 <- data[,-which(colnames(data)==intervention)]
		
		g <- matrix(NA,nPerm,2*(length(unique(stp2))-1))

		for(i in 1:nPerm){
			set.seed(12890*i+1)
			tp3 <- data.frame(stp,sample(stp2))
			names(tp3) <- c(paste(random),paste(intervention))
			data.tp4 <- merge(data2,tp3,by=random)
			data.tp4 <- data.tp4[order(data.tp4[,which(colnames(data.tp4)==random)]),]
			cluster = data.tp4[,which(colnames(data.tp4)==random)]
			mf <- model.frame(formula=formula, data=data.tp4)
			fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data.tp4)))
			tmp <- colnames(fixedDesignMatrix )
			tmp[1]  <- "Intercept"
			colnames(fixedDesignMatrix)<- tmp
			posttest <- model.response(mf)
      		intervention <- intervention
				
			p2CRTFREQ <-crtP(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster)	
			chkppp <- unlist(p2CRTFREQ$ES)[c(1,2)]
			chkppp2 <- c(seq(1,6*(length(unique(tp3[,2]))-1),6),seq(2,6*(length(unique(tp3[,2]))-1),6))
			chkppp3 <- chkppp2[order(chkppp2)] 
			g[i,]  <-  chkppp[chkppp3]
		}
			ntpp <- rep(names(p2CRTFREQ$ES),2)
			ntpp <- ntpp[order(ntpp )]
			wt <- rep(c("Within","Total"),length(names(p2CRTFREQ$ES)))
			colnames(g) <- paste(ntpp ,wt,sep="")
			return(g)
}


###### MST perm


mst.perm <- function(formula,data,trt,intervention,nPerm,random,cluster){

		data2 <- data
		g <- matrix(NA,nPerm,2*(length(unique(trt))-1))
		for(i in 1:nPerm){
			set.seed(12890*i+1)
			data2[,which(colnames(data)==intervention)]<-unlist(tapply(trt,cluster,function(x)sample(x)))
			data3 <- data2[order(data2[,which(colnames(data2)==random)],data2[,which(colnames(data2)==intervention)]),]
			cluster = data3[,which(colnames(data3)==random)]
			 mf <- model.frame(formula=formula, data=data3)
			fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data3)))
			tmp <- colnames(fixedDesignMatrix )
			tmp[1]  <- "Intercept"
			colnames(fixedDesignMatrix)<- tmp
			posttest <- model.response(mf)
      		intervention <- intervention
			trt2 <- data3[,which(colnames(data3)==intervention)]	
			p2CRTFREQ <-rbdP(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt2,cluster=cluster)	

			chkppp <- unlist(p2CRTFREQ$ES)
			chkppp2 <- c(seq(1,6*(length(unique(trt))-1),6),seq(2,6*(length(unique(trt))-1),6))
			chkppp3 <- chkppp2[order(chkppp2)] 
			g[i,]  <-  chkppp[chkppp3]
		}
			ntpp <- rep(names(p2CRTFREQ$ES),2)
			ntpp <- ntpp[order(ntpp )]
			wt <- rep(c("Within","Total"),length(names(p2CRTFREQ$ES)))
			colnames(g) <- paste(ntpp ,wt,sep="")
			return(g)

}



## random intercept model with bootstrap CI - internal
crt.crt<- function(posttest,fixedDesignMatrix,intervention,cluster,bt){
	
	posttest2 <- posttest[bt]
	fixedDesignMatrix2 <- fixedDesignMatrix[bt,]
	cluster2 <- cluster[bt]


	freqFit <- try(lmer(posttest2~ fixedDesignMatrix2-1+(1|cluster2)),silent=TRUE)
	output2 <- NULL
     if(attr(freqFit,"class")!="try-error"){

	betaB <- data.frame(summary(freqFit)$coefficients[,1])
	row.names(betaB)<- colnames(fixedDesignMatrix)
	betaB <- betaB
	var.B2<- as.matrix(summary(freqFit)$varcor)
	var.B3 <- c(matrix(attr(var.B2[[1]],"stddev")))
	var.B <- var.B3^2
	var.W<- summary(freqFit)$sigma^2
	var.tt <- var.W+var.B
	ICC <- var.B/var.tt
	sigmaBE <- c(var.B,var.W)
	sigmaBE <- sigmaBE
	names(sigmaBE)<- c("Schools","Pupils")

	btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)
	      
      output2 <- list()
	for( i in 1:length(btp )){
	
		beta <- betaB[btp[i],1]
      	group <- fixedDesignMatrix[,btp[i]]

  		esWithin <- beta/sqrt(var.W)
  		esTotal <- beta/sqrt(var.tt) 
   	
		output1 <- data.frame(rbind(esWithin,esTotal))
		names(output1) <- c("Estimate")
		rownames(output1) <- c("Within","Total")
		output2[[i]] <- output1
	}
	names(output2) <- row.names(betaB)[btp]
	
	}
	return(output2)
}

## random intercept and treatmebt-by-school interaction model - internal

rbd <- function(posttest,fixedDesignMatrix,intervention,trt,cluster){

	freqFit <- lmer(posttest~ fixedDesignMatrix-1+(1|trt:cluster)+(1|cluster))

      cit <- confint(freqFit)
	betaB <- data.frame(cbind(summary(freqFit)$coefficients[,1],cit[-c(1:3),]))
	row.names(betaB)<- colnames(fixedDesignMatrix)
	colnames(betaB) <- c("Estimate","95% LB ","95% UB")
	betaB <- betaB
	var.B2<- as.matrix(summary(freqFit)$varcor)
	var.B3 <- c(c(matrix(attr(var.B2[[1]],"stddev"))),c(matrix(attr(var.B2[[2]],"stddev"))))
	var.B3 <- var.B3^2
	var.sch <- var.B3[1]
	var.schTrt <- var.B3[2]	
	var.B <- sum(var.B3)
	var.W<- summary(freqFit)$sigma^2
	var.tt <- var.W+var.B
	ICC <- var.B/var.tt
	sigmaBE <- c(var.sch,var.schTrt,var.W)
	sigmaBE <- sigmaBE
	names(sigmaBE)<- c("Schools","Intervention:School","Pupils")
      schRand <- data.frame(ranef(freqFit)$cluster)
	names(schRand)<- "Estimate"
	btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)
	      
      
      output2 <- list()
	for( i in 1:length(btp )){
	
		beta <- betaB[btp[i],1]
      	group <- fixedDesignMatrix[,btp[i]]

  		esWithin <- g.within(var.w=var.W, beta=beta, icc=ICC , group=group, schoolID=cluster)
  		esTotal <- g.total(var.tt=var.tt, beta=beta, icc=ICC , group=group, schoolID=cluster)
   	
		output1 <- data.frame(rbind(esWithin,esTotal))
		colnames(output1) <- c("Estimate","95% LB","95% UB")
		rownames(output1) <- c("Within","Total")
		output2[[i]] <- round(output1,2)
	}
	names(output2) <- row.names(betaB)[btp]
      output <- list(Beta=round(betaB,2),covParm=round(sigmaBE,2),ES=output2,SchEffects=round(schRand,2))


	return(output)
}


####

rbdP <- function(posttest,fixedDesignMatrix,intervention,trt,cluster){

	freqFit <- lmer(posttest~ fixedDesignMatrix-1+(1|trt:cluster)+(1|cluster))

	betaB <- data.frame(summary(freqFit)$coefficients)
	row.names(betaB)<- colnames(fixedDesignMatrix)
	#colnames(betaB) <- c("Estimate","95% LB ","95% UB")
	betaB <- betaB
	var.B2<- as.matrix(summary(freqFit)$varcor)
	var.B3 <- c(c(matrix(attr(var.B2[[1]],"stddev"))),c(matrix(attr(var.B2[[2]],"stddev"))))
	var.B3 <- var.B3^2
	var.sch <- var.B3[1]
	var.schTrt <- var.B3[2]	
	var.B <- sum(var.B3)
	var.W<- summary(freqFit)$sigma^2
	var.tt <- var.W+var.B
	ICC <- var.B/var.tt
	sigmaBE <- c(var.sch,var.schTrt,var.W)
	sigmaBE <- sigmaBE
	names(sigmaBE)<- c("Schools","Intervention:School","Pupils")
      schRand <- data.frame(ranef(freqFit)$cluster)
	names(schRand)<- "Estimate"
	btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)
	      
      
      output2 <- list()
	for( i in 1:length(btp )){
	
		beta <- betaB[btp[i],1]
      	group <- fixedDesignMatrix[,btp[i]]

  		esWithin <- g.within(var.w=var.W, beta=beta, icc=ICC , group=group, schoolID=cluster)
  		esTotal <- g.total(var.tt=var.tt, beta=beta, icc=ICC , group=group, schoolID=cluster)
   	
		output1 <- data.frame(rbind(esWithin,esTotal))
		colnames(output1) <- c("Estimate","95% LB","95% UB")
		rownames(output1) <- c("Within","Total")
		output2[[i]] <- round(output1,2)
	}
	names(output2) <- row.names(betaB)[btp]
      output <- list(ES=output2)


	return(output)
}


## random intercept and treatmebt-by-school interaction model with bootstrap CI - internal

rbd.rbd <- function(posttest,fixedDesignMatrix,intervention,trt,cluster,bt){
	
	posttest2 <- posttest[bt]
	fixedDesignMatrix2 <- fixedDesignMatrix[bt,]
	trt2 <- trt[bt]
	cluster2 <- cluster[bt]


	freqFit <- try(lmer(posttest2~ fixedDesignMatrix2-1+(1|trt2:cluster2)+(1|cluster2)),silent=TRUE)
	output2 <- NULL
     if(attr(freqFit,"class")!="try-error"){

	betaB <- data.frame(summary(freqFit)$coefficients[,1])
	row.names(betaB)<- colnames(fixedDesignMatrix)
	betaB <- betaB
	var.B2<- as.matrix(summary(freqFit)$varcor)
	var.B3 <- c(c(matrix(attr(var.B2[[1]],"stddev"))),c(matrix(attr(var.B2[[2]],"stddev"))))
	var.B3 <- var.B3^2
	var.sch <- var.B3[1]
	var.schTrt <- var.B3[2]	
	var.B <- sum(var.B3)
	var.W<- summary(freqFit)$sigma^2
	var.tt <- var.W+var.B
	ICC <- var.B/var.tt
	sigmaBE <- c(var.sch,var.schTrt,var.W)
	sigmaBE <- sigmaBE
	names(sigmaBE)<- c("Schools","Intervention:School","Pupils")

	btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)
	      
      output2 <- list()
	for( i in 1:length(btp )){
	
		beta <- betaB[btp[i],1]
      	group <- fixedDesignMatrix[,btp[i]]

  		esWithin <- beta/sqrt(var.W)
  		esTotal <- beta/sqrt(var.tt) 
   	
		output1 <- data.frame(rbind(esWithin,esTotal))
		names(output1) <- c("Estimate")
		rownames(output1) <- c("Within","Total")
		output2[[i]] <- output1
	}
	names(output2) <- row.names(betaB)[btp]
	
	}
	return(output2)
}


## Hedges's from  within variance

g.within <- function(var.w, beta, icc, group, schoolID){
   t <- group; id <- schoolID
   d.w <- (beta/sqrt(var.w))
   n.it <- table(id[t==1]); n.ic <- table(id[t==0])
   m.t <- length(unique(id[t==1])); m.c <- length(unique(id[t==0]))
   M <- (m.t + m.c)
   N.t <- sum(table(id[t==1])); N.c <- sum(table(id[t==0]))
   N <- (N.t + N.c)
   n.sim.1 <- ((N.c * sum(n.it^2))/(N.t*N))
   n.sim.2 <- ((N.t * sum(n.ic^2))/(N.c*N))
   n.sim <- (n.sim.1 + n.sim.2)
   vterm1 <- ((N.t+N.c)/(N.t*N.c))
   vterm2 <- (((1+(n.sim-1)*icc))/(1-icc))
   vterm3 <- ((d.w^2)/(2*(N-M)))
   se <- sqrt(vterm1*vterm2+vterm3)
   LB <- (d.w-1.96*se); UB <- (d.w+1.96*se)
   output <- data.frame(d.w, LB, UB)
   names(output) <- c("g", "LB", "UB")
   return(output)
}


## Hedges' ES from total variance

g.total <- function(var.tt, beta, icc, group, schoolID){
   t <- group; id <- schoolID
   n.it <- table(id[t==1]); n.ic <- table(id[t==0])
   m.t <- length(unique(id[t==1])); m.c <- length(unique(id[t==0]))
   M <- (m.t + m.c)
   N.t <- sum(table(id[t==1])); N.c <- sum(table(id[t==0]))
   N <- (N.t + N.c)
   n.ut <- ((N.t^2-sum(n.it^2))/(N.t*(m.t-1)))
   n.uc <- ((N.c^2-sum(n.ic^2))/(N.c*(m.c-1)))
   dt.1 <- (beta/sqrt(var.tt))
   dt.2 <- sqrt(1-icc*(((N-n.ut*m.t-n.uc*m.c)+n.ut+n.uc-2)/(N-2)))
   d.t <- (dt.1*dt.2)
   
   n.sim.1 <- ((N.c * sum(n.it^2))/(N.t*N))
   n.sim.2 <- ((N.t * sum(n.ic^2))/(N.c*N))
   n.sim <- (n.sim.1 + n.sim.2)
   B <- (n.ut*(m.t-1)+n.uc*(m.c-1))
   A.t <- ((N.t^2*sum(n.it^2)+(sum(n.it^2))^2-2*N.t*sum(n.it^3))/N.t^2)
   A.c <- ((N.c^2*sum(n.ic^2)+(sum(n.ic^2))^2-2*N.c*sum(n.ic^3))/N.c^2)
   A <- (A.t + A.c)
   
   vterm1 <- (((N.t+N.c)/(N.t*N.c))*(1+(n.sim-1)*icc))
   vterm2 <- (((N-2)*(1-icc)^2+A*icc^2+2*B*icc*(1-icc))*d.t^2)
   vterm3 <- (2*(N-2)*((N-2)-icc*(N-2-B)))
   se <- sqrt(vterm1+vterm2/vterm3)
   LB <- (d.t-1.96*se); UB <- (d.t+1.96*se)
   output <- data.frame(d.t, LB, UB)
   names(output)<- c("g", "LB", "UB")
   return(output)
}






## simple randomised analysis
srt <- function(posttest,fixedDesignMatrix,intervention,trt){
	
	freqFit <- lm(posttest~ fixedDesignMatrix-1)
      cit <- confint(freqFit)
	citt <- rowSums(is.na(cit))
	betaB <- data.frame(cbind(summary(freqFit)$coefficients[which(citt==0),1],cit[which(citt==0),]))
	row.names(betaB)<- colnames(fixedDesignMatrix)[which(citt==0)]
	colnames(betaB) <- c("Estimate","95% LB ","95% UB")
	betaB <- betaB

	btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)
	btp2 <- which(substring(colnames(fixedDesignMatrix),1,nchar(intervention))==intervention)

	tmpTRT <- table(trt)      
      output2 <- NULL
	for( i in 1:length(btp )){
	
  		beta <- betaB[btp[i],1]
   		sd.pool <- summary(freqFit)$sigma
   		cd <- (beta/sd.pool)
             trt2 <- unique(fixedDesignMatrix[,btp2[i]])
		n.c <- tmpTRT[names(tmpTRT)==trt2[2]] 
		n.t <- tmpTRT[names(tmpTRT)==trt2[1]] 
   		var.cd <- ((n.t+n.c)/(n.t*n.c)+cd^2/(2*(n.t+n.c)))
   		se.cd <- sqrt(var.cd)
   		cd.lb <- (cd - 1.96*se.cd)
   		cd.ub <- (cd + 1.96*se.cd)
   		j.df <- (1 - (3/(4*(n.t+n.c-2)-1)))
   		g <- (j.df*cd)
   		var.g <- (j.df^2 * var.cd)
   		se.g <- sqrt(var.g)
   		g.lb <- (g - 1.96*se.g)
   		g.ub <- (g + 1.96*se.g)
  		gtmp <- round(data.frame(g, g.lb, g.ub),2)
		output2 <- rbind(output2,gtmp)
	}
	names(output2)<- c("Estimate","95% LB","95% UB")
	row.names(output2) <- row.names(betaB)[btp]
      output <- list(Beta=round(betaB,2),ES=round(output2,2),sigma2=round(summary(freqFit)$sigma^2,2))
	return(output)

}


## SRT with bootstrap CI - internal
srt.srt<- function(posttest,fixedDesignMatrix,intervention,bt){
	
	posttest2 <- posttest[bt]
	fixedDesignMatrix2 <- fixedDesignMatrix[bt,]



	freqFit <- try(lm(posttest2~ fixedDesignMatrix2-1),silent=TRUE)
	output2 <- NULL
      if(attr(freqFit,"class")!="try-error"){

		betaB <- data.frame(summary(freqFit)$coefficients[,1])
		row.names(betaB)<- colnames(fixedDesignMatrix)
		betaB <- betaB

		sd.pool <- summary(freqFit)$sigma

		btp <- which(substring(row.names(betaB),1,nchar(intervention))==intervention)
	      
      		output1 <- NULL
		for( i in 1:length(btp )){
	
  			beta <- betaB[btp[i],1]	
     
  	
			output1[i] <-beta/sd.pool

		}
		names(output1) <- row.names(betaB)[btp]
	
	}
	output1<- matrix(output1,1,length(btp))

	return(output1)
}



#' Multisite trial data.
#'
#' A multisite education trial data containing 54 schools.
#'
#' \itemize{
#'   \item Posttest. posttest scores
#'   \item Prettest. prettests scores
#'   \item Intervention. indicator for intervention groups 
#'coded as 1 for the intervention group and 0 for the control group.
#'   \item Compliance. percentage of intervention sessions attended by pupils
#'   \item Compliance. percentage of sessions attended by pupils
#'   \item School. numeric school identifier
#' }
#'
#' @format A data frame with 210 rows and 5 variables
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
#'coded as 1 for the intervention group and 0 for the control group.
#'   \item Compliance. percentage of sessions attended by pupils
#'   \item School. numeric school identifier
#' }
#'
#' @format A data frame with 265 rows and 5 variables
#' @name iwq
NULL


## IMPORTS ##
#' @importFrom lme4 lmer ranef
#' @importFrom geoR rinvchisq
#' @importFrom mvtnorm rmvnorm
NULL

#############################################################################
############# SRT main functions ################################################


#' Analysis of Simple Randomised Trial (SRT).
#' 
#' \code{srtFREQ} perfoms analysis of education trials under the assumption of independent errors between pupils. 
#'This can also be used with schools as fixed effects.

#' @export
#' @param formula specifies the model to be analysed. It is of the form y~x1+x2+..., 
#' where y is the outcome variable and X's are the predictors.
#' @param intervention the name of the intervention variable as appeared in formula.  
#' This must be put in quotes.  For example "intervention" or "treatment" or "group".
#' @param nBoot number of bootstrap required to generate bootstrap confidence interval. Default is NULL.
#' @param nPerm number of permutations required to generate permutation p-value. Default  is NULL.
#' @param data data frame containing the data to be analysed. 
#' @return S3 \code{mcpi} object; a list consisting of
#' \itemize{
#' \item \code{Beta}. Estimates and confidence intervals for the predictors specified in the model. 
#' It will be the slope for a continuous predictor and the mean difference for a dummy variable or a categorical predictor.
#' \item \code{ES}. Hedges' effect size for the intervention effect. 
#' If \code{nBoot} is not specified, the confidence intervals are classical 95% CIs based on standard error. 
#' If \code{nBoot} is specified, they are non-parametric bootstrapped confidence intervals.
#' \item \code{sigma2}. Residual variance. Its square root will generate a pooled standard deviation. 
#' \item \code{Perm}. A vector containing the distribution of effect size under the null hypothesis. It is produced only if \code{nPerm} is specified.
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

		permES<- matrix(NA,nPerm,(length(unique(trt))-1))
		
		set.seed(1020252)
		for (i in 1:nPerm){
				
			data[,which(colnames(data)==intervention)]<- sample(trt)
			fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))

			p2CRTFREQ <-srt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt)	

			permES[i,]  <-  p2CRTFREQ$ES[,1]
		}
		output$Perm<- data.frame(permES)
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
		output$ES <-bootES2
		#output$Bootstrap <- bootResults2  

	}

	return(output) 
}



#############################################################################
############# CRT main functions ################################################
#' Analysis of Cluster Randomised Trials using MLM.
#' 
#' \code{crtFREQ} is a frequentist method that can be used to calculate effect size for cluster randomised trials based on residual variance or total variance.
#' 
#' @export
#' @param formula specifies the model to be analysed. It is of the form form y ~ x1+x2 +..., 
#' where y is the outcome variable and X's are the predictors.
#' @param random a string variable specifying the "clustering" variable as contained in the data.
#' This must be put between  quotes. For example, "school".
#' @param intervention specifies the name of the intervention variable as appeared in formula.
#' This must be put between quotes.  For example "intervention" or "treatment" or "group"...
#' @param nBoot number of bootstrap required to generate bootstrap confidence interval. Default is NULL. 
#' @param nPerm number of permutations required to generate permutation p-value. Default  is NULL.
#' @param data the data frame to be analysed. 
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}. Estimates and confidence intervals for the predictors specified in the model. 
#' It will be the slope for a continuous predictor and the mean difference for a dummy variable or a categorical predictor.
#' \item \code{ES}. Hedges' effect size for the intervention effect. If nBoot is not specified, 
#' the confidence intervals are 95% CIs based on standard error. If nBoot is specified, 
#' they are non-parametric bootstrapped confidence intervals.
#' \item \code{covParm}. A vector of variance decomposition into between-variance (Schools), within-variance (Pupils) and total variance.
#' It also contains the intra-cluster correlation (ICC).
#' \item \code{SchEffects}. Individual school effects at baseline.
#' \item \code{Perm}. A matrix of the distribution of ES under the null hypothesis. 
#' The two columns in the matrix represents ES based on within variance and total variance. 
#' Produced only if nPerm is specified.
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
		output$Perm <- round(data.frame(output$Perm),2)
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
		output$Bootstrap <- round(data.frame(bootResults),2)  
	}

	return(output) 
}




#############################################################################
############# MST main functions ################################################

#############################################################################
############# MST main functions ################################################
#' Analysis of Multisite Randomised Trials using MMLM.
#' 
#' \code{mstFREQ} performs multilevel model for analysing randomised multisite trials with schools and school by intervention interactions specified as random effects
#' 
#' @export
#' @param formula specifies the model to be analysed. It is of the form form y ~ x1+x2 +..., 
#' where y is the outcome variable and X's are the predictors.
#' @param random a string variable specifying the "clustering" variable as contained in the data.
#' This must be put between  quotes. For example, "school".
#' @param intervention specifies the name of the intervention variable as appeared in formula.
#' This must be put between quotes.  For example "intervention" or "treatment" or "group"...
#' @param nBoot number of bootstrap required to generate bootstrap confidence interval. Default is NULL. 
#' @param nPerm number of permutations required to generate permutation p-value. Default  is NULL.
#' @param data the data frame to be analysed. 
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{Beta}. Estimates and confidence intervals for the predictors specified in the model. 
#' It will be the slope for a continuous predictor and the mean difference for a dummy variable or a categorical predictor.
#' \item \code{ES}. Hedges' effect size for the intervention effect. If nBoot is not specified, 
#' the confidence intervals are 95% CIs based on standard error. If nBoot is specified, 
#' they are non-parametric bootstrapped confidence intervals.
#' \item \code{covParm}. A vector of variance decomposition into between-variance (Schools), within-variance (Pupils) and total variance.
#' It also contains the intra-cluster correlation (ICC).
#' \item \code{SchEffects}. Individual school effects at baseline.
#' \item \code{Perm}. A matrix of the distribution of ES under the null hypothesis. 
#' The two columns in the matrix represents ES based on within variance and total variance. 
#' Produced only if nPerm is specified.
#' }
#' @example inst/examples/mstExample.R
mstFREQ<- function(formula,random,intervention,nPerm=NULL,data,nBoot=NULL)UseMethod("mstFREQ")

#' @export
mstFREQ.formula <- function(formula,random,intervention,nPerm=NULL,data,nBoot=NULL){

	data <- data[order(data[,which(colnames(data)==random)],data[,which(colnames(data)==intervention)]),]
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
		output$Perm <- round(data.frame(output$Perm),2)
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
		output$Bootstrap <- round(data.frame(bootResults),2)  

	}

	return(output) 
}

###################################################################################################
############# Bayesian analysis of cluster randomised trials
##################################################################################################

#' Bayesia analysis of cluster randomised trials Using vague priors.
#' 
#' \code{crtBayes} performs analysis of cluster randomised trial using multilevel model within the Bayesian framework 
#' assuming vague priors.
#' 
#' @export
#' @param formula specifies the model to be analysed. It is of the form form y ~ x1+x2 +..., 
#' where y is the outcome variable and X's are the predictors.
#' @param random a string variable specifying the "clustering" variable as contained in the data.
#' This must be put between  quotes. For example, "school".
#' @param intervention specifies the name of the intervention variable as appeared in formula.
#' This must be put between quotes.  For example "intervention" or "treatment" or "group"...
#' @param nSim number of MCMC simulations to generate samples from full conditional posterior distributions. 
#' A minimum of 10,000 is recommended.
#' @param data specifies data frame containing the data to be analysed. 
#' @return S3 \code{mcpi} object; a list consisting of
#' \itemize{
#' \item \code{Beta}. Estimates and confidence intervals for the predictors specified in the model. 
#' It will be a slope for a continuous predictor and a mean difference for a dummy variable or a categorical predictor.
#' \item \code{ES}. Effect size for the intervention effect.
#' \item \code{covParm}. A vector of variance decomposition into between-variance (Schools), within-variance (Pupils) and total variance.
#' It also contains the intra-cluster correlation (ICC).
#' \item \code{ProbES}. A maxtrix containing the probability of observing ES greater than a pre-specified value. First column is for within-variance, second column for between-variance and the third column for total-variance.
#' \item \code{SchEffects}. Individual school effects at baseline.
#' }
#' @example inst/examples/bayesExample.R
crtBayes <- function(formula,random,intervention,nSim,data)UseMethod("crtBayes")

#' @export
crtBayes.formula <- function(formula,random,intervention,nSim=nSim,data){
	data <- data[order(data[,which(colnames(data)==random)],data[,which(colnames(data)==intervention)]),]
	tmp3 <- which(colnames(data)==intervention)
	data[,tmp3] <- as.factor(data[,tmp3])
	tmp2 <- which(colnames(data)==random)
	cluster2 = data[,tmp2]

	mf <- model.frame(formula=formula, data=data)
	mf <- mf[order(cluster2),]
	fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
	tmp <- colnames(fixedDesignMatrix )
	tmp[1]  <- "Intercept"
	colnames(fixedDesignMatrix)<- tmp
	posttest <- model.response(mf)
      intervention <- intervention
	trt <- data[,which(colnames(data)==intervention)]
	tmp2 <- which(colnames(data)==random)
	cluster <- cluster2[order(cluster2)]
	nsim=nSim
	if(length(tmp2)!= 1){stop("Clustering variable misspecified")}
	if(length(tmp3)!= 1){stop("Intervention variable misspecified")}

	if(nsim < 10000){stop("nsim >= 10000 is recommended")}

	BayesOutput <- erantBAYES(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster, nsim=nsim)
	output  <- errantSummary(bayesObject=BayesOutput,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention)
	output$SchEffects <- data.frame(Schools=unique(cluster),output$SchEffects)
	colnames(output$SchEffects ) <- c("Schools","Estimate","95% LB", "95% UB")
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
	schParm5<- schParm4
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

##summarise all bayesian parameters - internal
errantSummary <- function(bayesObject,fixedDesignMatrix,intervention){
	covValues <- covSummary(bayesObject=bayesObject)
	covValues <- data.frame(covValues[c(2,1,3,4)])
	row.names(covValues ) <- c("Schools","Pupils","Total","ICC")
	covValues <- t(covValues )
	row.names(covValues ) <- NULL
	betaValues <- betaSummary(bayesObject=bayesObject)
	esValues <- esSummary(bayesObject,fixedDesignMatrix,intervention)
	schValues <-schSummary(bayesObject=bayesObject)
	es.prob <- esProb(bayesObject=bayesObject,esOutput=esValues)
      output <- list(Beta=round(betaValues,2),covParm= round(covValues,2),ES=round(esValues,2),ProbES=round(data.frame(es.prob),2),SchEffects=round(schValues,2))
}




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
		
      np<- row.names(summary(freqFit)$coef)
      cit <- confint(freqFit,np)
	betaB <- data.frame(cbind(summary(freqFit)$coefficients[,1],cit))
	row.names(betaB)<- colnames(fixedDesignMatrix)
	colnames(betaB) <- c("Estimate","95% LB ","95% UB")
	betaB <- betaB
	var.B<- as.numeric(summary(freqFit)$varcor)
	var.W<- summary(freqFit)$sigma^2
	var.tt <- var.W+var.B
	ICC <- var.B/var.tt
	sigmaBE <- c(var.B,var.W,var.B+var.W,(var.B/(var.B+var.W)))
	sigmaBE <- sigmaBE
	names(sigmaBE)<- c("Schools","Pupils","Total","ICC")
      schRand <- data.frame(unique(cluster),ranef(freqFit)$cluster)
	names(schRand)<- c("Schools","Estimate")
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

## internal
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


## - internal

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


## - internal


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



## - internal
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

## - internal

rbd <- function(posttest,fixedDesignMatrix,intervention,trt,cluster){

	freqFit <- lmer(posttest~ fixedDesignMatrix-1+(1|trt:cluster)+(1|cluster))
      np<- row.names(summary(freqFit)$coef)
      cit <- confint(freqFit,np)
	betaB <- data.frame(cbind(summary(freqFit)$coefficients[,1],cit))
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
	sigmaBE <- c(var.sch,var.schTrt,var.W,var.tt,ICC )
	sigmaBE <- sigmaBE
	names(sigmaBE)<- c("Schools","Intervention:School","Pupils","Total","ICC")
      schRand <- data.frame(unique(cluster),ranef(freqFit)$cluster)
	names(schRand)<- c("Schools","Estimate")
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


## - internal

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


## - internal

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


## - internal

g.within <- function(var.w, beta, icc, group, schoolID){
   t <- group; id <- schoolID
   d.w <- (beta/sqrt(var.w))
   n.it <- table(id[t==1]); n.ic <- table(id[t==0])
   m.t <- length(unique(id[t==1])); m.c <- length(unique(id[t==0]))
   M <- (m.t + m.c)
   N.t <- sum(table(id[t==1])); N.c <- sum(table(id[t==0]))
   N <- (N.t + N.c)
   n.sim.1 <- ((N.c * sum(n.it^2))/(as.numeric(N.t)*as.numeric(N)))
   n.sim.2 <- ((N.t * sum(n.ic^2))/(as.numeric(N.c)*as.numeric(N)))
   n.sim <- (n.sim.1 + n.sim.2)
   vterm1 <- ((N.t+N.c)/(as.numeric(N.t)*as.numeric(N.c)))
   vterm2 <- (((1+(n.sim-1)*icc))/(1-icc))
   vterm3 <- ((d.w^2)/(2*(N-M)))
   se <- sqrt(vterm1*vterm2+vterm3)
   LB <- (d.w-1.96*se); UB <- (d.w+1.96*se)
   output <- data.frame(d.w, LB, UB)
   names(output) <- c("g", "LB", "UB")
   return(output)
}

## - internal

g.total <- function(var.tt, beta, icc, group, schoolID){
   t <- group; id <- schoolID
   n.it <- table(id[t==1]); n.ic <- table(id[t==0])
   m.t <- length(unique(id[t==1])); m.c <- length(unique(id[t==0]))
   M <- (m.t + m.c)
   N.t <- sum(table(id[t==1])); N.c <- sum(table(id[t==0]))
   N <- (N.t + N.c)
   n.ut <- ((N.t^2-sum(n.it^2))/(as.numeric(N.t)*as.numeric(m.t-1)))
   n.uc <- ((N.c^2-sum(n.ic^2))/(as.numeric(N.c)*as.numeric(m.c-1)))
   dt.1 <- (beta/sqrt(var.tt))
   dt.2 <- sqrt(1-icc*(((N-n.ut*m.t-n.uc*m.c)+n.ut+n.uc-2)/(N-2)))
   d.t <- (dt.1*dt.2)
   
   n.sim.1 <- ((as.numeric(N.c) * sum(n.it^2))/(as.numeric(N.t)*as.numeric(N)))
   n.sim.2 <- ((as.numeric(N.t) * sum(n.ic^2))/(as.numeric(N.c)*as.numeric(N)))
   n.sim <- (n.sim.1 + n.sim.2)
   B <- (n.ut*(m.t-1)+n.uc*(m.c-1))
   A.t <- ((as.numeric(N.t)^2*sum(n.it^2)+(sum(n.it^2))^2-2*as.numeric(N.t)*sum(n.it^3))/as.numeric(N.t)^2)
   A.c <- ((as.numeric(N.c)^2*sum(n.ic^2)+(sum(n.ic^2))^2-2*as.numeric(N.c)*sum(n.ic^3))/as.numeric(N.c)^2)
   A <- (A.t + A.c)
   
   vterm1 <- (((N.t+N.c)/(as.numeric(N.t)*as.numeric(N.c)))*(1+(n.sim-1)*icc))
   vterm2 <- (((N-2)*(1-icc)^2+A*icc^2+2*B*icc*(1-icc))*d.t^2)
   vterm3 <- (2*(N-2)*((N-2)-icc*(N-2-B)))
   se <- sqrt(vterm1+vterm2/vterm3)
   LB <- (d.t-1.96*se); UB <- (d.t+1.96*se)
   output <- data.frame(d.t, LB, UB)
   names(output)<- c("g", "LB", "UB")
   return(output)
}


## - internal
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


## - internal
srt.srt<- function(posttest,fixedDesignMatrix,intervention,bt){
	
	posttest2 <- posttest[bt]
	fixedDesignMatrix2 <- fixedDesignMatrix[bt,]



	freqFit <- try(lm(posttest2~ fixedDesignMatrix2-1),silent=TRUE)
	output2 <- NULL
      if(attr(freqFit,"class")!="try-error"){

		betaB <- data.frame(summary(freqFit)$coefficients[,1])
            ntpp <- as.character(sapply(as.character(rownames(betaB)),function(x)substring(x,(nchar("fixedDesignMatrix2")+1),nchar(x))))
		row.names(betaB)<- ntpp
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



#################
############# CRT main functions ################################################
#' CACE Analysis of Cluster Randomised Trials using MLM.
#' 
#' \code{caceCRTBoot} performs CACE analysis of cluster randomised trials. 
#' Intervention variable must be coded as dummy with multiple analysis for multi-arm trials. 
#' 
#' @export
#' @param formula specifies the model to be analysed.  
#' It is of the form y~x1+x2+...,  where y is the outcome variable and X's are the predictors.
#' @param intervention the name of the intervention variable as appeared in formula.  
#' This must be put in quotes.  For example "intervention" or "treatment" or "group".
#' @param random a string variable specifying the "clustering" variable as contained in the data.
#' This must be put between  quotes. For example, "school".
#' @param compliance percentages of sessions attended by pupils. 
#' @param nBoot number of bootstrap required to generate bootstrap confidence interval. This must be specified.
#' @param data data frame containing the data to be analysed.
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{CACE}. Estimated CACE effect size based on percentages of sessions attended by pupils. 
#' The percentage data is converted into the following grids (0, 10, 20, 30, 40, 50, 60, 70, 80, 90) 
#' and CACE effect size is calculated for each grid.
#' \item \code{Compliers}. A summary table of the percentage of pupils in the intervention and control groups that 
#' attended more than a pre-specified percentage of sessions. The values for the control group should be zeros if 
#' there is no dilution in which a pupil or school in the control group receives intervention.
#' }
#' @example inst/examples/crtCACEExample.R

caceCRTBoot <- function(formula,random,intervention,compliance,nBoot,data){
  
	data <- data[order(data[,which(colnames(data)==random)]),]
       
	intervention <- intervention
	trt <- data[,which(colnames(data)==intervention)]

	if(length(table(trt))!=2){stop("Applicable only to two-arms trials")}
	comp1 <- data[,which(colnames(data)==compliance)]
	tmp2 <- which(colnames(data)==random)
	cluster2 = data[,tmp2]
	tmp3 <- which(colnames(data)==intervention)
	data[,tmp3] <- as.factor(data[,tmp3])
	mf <- model.frame(formula=formula, data=data)
	mf <- mf[order(cluster2),]
	comp2 <- comp1[order(cluster2)]
	cluster <- cluster2[order(cluster2)]
	fixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(mf, "terms"), data=data)))
	tmp <- colnames(fixedDesignMatrix )
	tmp[1]  <- "Intercept"
	colnames(fixedDesignMatrix)<- tmp
	posttest <- model.response(mf)

	output <- crt.cace(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,cluster=cluster)	
	
	theta <- seq(0,ifelse(max(comp2)==100,90,max(comp2)),10)
	alpha_ITT_E <- output$ES[[1]][2,1]

	d <- sapply(theta,function(x)ifelse(comp2 > x,1,0))
	ibit <- data[,tmp3]
      tc.pp <- apply(d,2,function(x) table(x,ibit))
      c.tpp1 <- tc.pp[2,]/(tc.pp[1,]+tc.pp[2,])
	t.tpp1 <-tc.pp[4,]/(tc.pp[3,]+tc.pp[4,])
	p1 <- t.tpp1 - c.tpp1 
	caceES  <- alpha_ITT_E/(p1)
	p1Table <- rbind(t.tpp1,c.tpp1,p1)
	p1Table <- data.frame(round(p1Table,2))
	row.names(p1Table )<- c("pT","pC","P=PT-pC")
	colnames(p1Table ) <- paste("P > ", theta,sep="")

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
	row.names(bootSamples ) <- NULL
	bootSamples <- data.frame(bootSamples )

	bESOutput <- NULL

	for(i in 1:ncol(bootSamples )){

		bData <- data[bootSamples[,i],]

		bcluster <- bData[,tmp2]
		bmf <- model.frame(formula=formula, data=bData)
		bfixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(bmf, "terms"), data=bData)))
		bposttest <- model.response(bmf)
		boutput <- crt.cace(posttest=bposttest,fixedDesignMatrix=bfixedDesignMatrix,intervention=intervention,cluster=bcluster)	
	
		balpha_ITT_E <- boutput$ES[[1]][2,1]

		bcomp <- bData[,which(colnames(bData)==compliance)]
		bd <- sapply(theta,function(x)ifelse(bcomp > x,1,0))
		bibit <- bData[,tmp3]
      	btc.pp <- apply(bd,2,function(x) table(x,bibit))
      	bc.tpp1 <- btc.pp[2,]/(btc.pp[1,]+btc.pp[2,])
		bt.tpp1 <- btc.pp[4,]/(btc.pp[3,]+btc.pp[4,])
		bp1 <- bt.tpp1 - bc.tpp1 
		bcaceES  <- balpha_ITT_E/(bp1)
		bESOutput <- rbind(bESOutput ,bcaceES)  
	}

	bootES <- apply(bESOutput ,2,function(x)quantile(x,prob=c(0.025,0.975)))
	output <- data.frame(Compliance=paste("P>",as.character(theta)),ES=round(caceES,2),LB=round(bootES[1,],2),UB=round(bootES[2,],2))
	output2 <- list(CACE=output,Compliers=p1Table )
	return(output2 ) 
}

############


## - internal
crt.cace <- function(posttest,fixedDesignMatrix,intervention,cluster){
	
	freqFit <- lmer(posttest~ fixedDesignMatrix-1+ (1|cluster))
		
	betaB <- data.frame(cbind(summary(freqFit)$coefficients[,1]))
	row.names(betaB)<- colnames(fixedDesignMatrix)
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
      output <- list(ES=output2)

	return(output)
}

#####################
#################

#############################################################################
############# MST main functions ################################################
#' CACE Analysis of Multisite Randomised Trials.
#' 
#' \code{caceMSTBoot} performs CACE analysis of multisite randomised trials.
#' Intervention variable must be coded as dummy with multiple analysis for multi-arm trials. 
#' 
#' @export
#' @param formula specifies the model to be analysed.  
#' It is of the form y~x1+x2+...,  where y is the outcome variable and X's are the predictors.
#' @param intervention the name of the intervention variable as appeared in formula.  
#' This must be put in quotes.  For example "intervention" or "treatment" or "group".
#' @param random a string variable specifying the "clustering" variable as contained in the data.
#' This must be put between  quotes. For example, "school".
#' @param compliance percentages of sessions attended by pupils. 
#' @param nBoot number of bootstrap required to generate bootstrap confidence interval. This must be specified.
#' @param data data frame containing the data to be analysed. 
#' @return S3 object; a list consisting of
#' \itemize{
#' \item \code{CACE}. Estimated CACE effect size based on percentages of sessions attended by pupils. 
#' The percentage data is converted into the following grids (0, 10, 20, 30, 40, 50, 60, 70, 80, 90) 
#' and CACE effect size is calculated for each grid.
#' \item \code{Compliers}. A summary table of the percentage of pupils in the intervention and control groups that 
#' attended more than a pre-specified percentage of sessions. The values for the control group should be zeros if 
#' there is no dilution in which a pupil or school in the control group receives intervention.
#' }
#' @example inst/examples/mstCACEExample.R

caceMSTBoot <- function(formula,random,intervention,compliance,nBoot,data){

	data <- data[order(data[,which(colnames(data)==random)],data[,which(colnames(data)==intervention)]),]
	intervention <- intervention
	trt <- data[,which(colnames(data)==intervention)]
	tmp2 <- which(colnames(data)==random)
	cluster2 = data[,tmp2]

	chk <- sum(rowSums(table(cluster2,trt)!=0)>1)
	if(chk ==0){stop("This not an MST design")}
 	
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
	if(length(table(trt))!=2){stop("Applicable only to two-arms trials")}
	comp1 <- data[,which(colnames(data)==compliance)]
	comp2 <- comp1[order(cluster2)]

	if(length(tmp2)!= 1){stop("Clustering variable misspecified")}
	if(length(tmp3)!= 1){stop("Intervention variable misspecified")}


	output <- rbd.cace(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt,cluster=cluster)	


	theta <- seq(0,ifelse(max(comp2)==100,90,max(comp2)),10)
	alpha_ITT_E <- output$ES[[1]][2,1]

	d <- sapply(theta,function(x)ifelse(comp2 > x,1,0))
	ibit <- data[,tmp3]
      tc.pp <- apply(d,2,function(x) table(x,ibit))
      c.tpp1 <- tc.pp[2,]/(tc.pp[1,]+tc.pp[2,])
	t.tpp1 <-tc.pp[4,]/(tc.pp[3,]+tc.pp[4,])
	p1 <- t.tpp1 - c.tpp1 
	caceES  <- alpha_ITT_E/(p1)
	p1Table <- rbind(t.tpp1,c.tpp1,p1)
	p1Table <- data.frame(round(p1Table,2))
	row.names(p1Table )<- c("pT","pC","P=PT-pC")
	colnames(p1Table ) <- paste("P > ", theta,sep="")

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
	row.names(bootSamples ) <- NULL
	bootSamples <- data.frame(bootSamples )

	bESOutput <- NULL

	for(i in 1:ncol(bootSamples )){

		bData <- data[bootSamples[,i],]

		bcluster <- bData[,tmp2]
		bmf <- model.frame(formula=formula, data=bData)
		bfixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(bmf, "terms"), data=bData)))
		bposttest <- model.response(bmf)
		btrt <- bData[,which(colnames(bData)==intervention)]
		boutput <-rbd.cace(posttest=bposttest,fixedDesignMatrix=bfixedDesignMatrix,intervention=intervention,trt=btrt,cluster=bcluster)
	
		balpha_ITT_E <- boutput$ES[[1]][2,1]

		bcomp <- bData[,which(colnames(bData)==compliance)]
		bd <- sapply(theta,function(x)ifelse(bcomp > x,1,0))
		bibit <- bData[,tmp3]
      	btc.pp <- apply(bd,2,function(x) table(x,bibit))
      	bc.tpp1 <- btc.pp[2,]/(btc.pp[1,]+btc.pp[2,])
		bt.tpp1 <- btc.pp[4,]/(btc.pp[3,]+btc.pp[4,])
		bp1 <- bt.tpp1 - bc.tpp1 
		bcaceES  <- balpha_ITT_E/(bp1)
		bESOutput <- rbind(bESOutput ,bcaceES)  
	}

	bootES <- apply(bESOutput ,2,function(x)quantile(x,prob=c(0.025,0.975)))
	output <- data.frame(Compliance=paste("P>",as.character(theta)),ES=round(caceES,2),LB=round(bootES[1,],2),UB=round(bootES[2,],2))
	output2 <- list(CACE=output,Compliers=p1Table )
	return(output2 ) 

}


##########


##  - internal

rbd.cace <- function(posttest,fixedDesignMatrix,intervention,trt,cluster){

	freqFit <- lmer(posttest~ fixedDesignMatrix-1+(1|trt:cluster)+(1|cluster))

	betaB <- data.frame(cbind(summary(freqFit)$coefficients[,1]))
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



####

#' CACE Analysis of Simple Randomised Trials.
#' 
#' \code{caceSRTBoot} performs CACE analysis of simple randomised trials. 
#' Intervention variable must be coded as dummy with multiple analysis for multi-arm trials. 
#' 
#' @export
#' @param formula specifies the model to be analysed. 
#' It is of the form y~x1+x2+...,  where y is the outcome variable and X's are the predictors.
#' @param intervention the name of the intervention variable as appeared in formula.  
#' This must be put in quotes.  For example "intervention" or "treatment" or "group".
#' @param compliance percentages of sessions attended by pupils. 
#' @param nBoot number of bootstrap required to generate bootstrap confidence interval. This must be specified.
#' @param data data frame containing the data to be analysed. 
#' @return S3 \code{mcpi} object; a list consisting of
#' \itemize{
#' \item \code{CACE}. Estimated CACE effect size based on percentages of sessions attended by pupils. 
#' The percentage data is converted into the following grids (0, 10, 20, 30, 40, 50, 60, 70, 80, 90) 
#' and CACE effect size is calculated for each grid.
#' \item \code{Compliers}. A summary table of the percentage of pupils in the intervention and control groups that 
#' attended more than a pre-specified percentage of sessions. The values for the control group should be zeros if 
#' there is no dilution in which a pupil or school in the control group receives intervention.
#' }
#' @example inst/examples/srtCACEExample.R

caceSRTBoot <- function(formula,intervention,compliance,nBoot,data){

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

	comp2 <- data[,which(colnames(data)==compliance)]


	output <- srt(posttest=posttest,fixedDesignMatrix=fixedDesignMatrix,intervention=intervention,trt=trt)


	theta <- seq(0,ifelse(max(comp2)==100,90,max(comp2)),10)
	alpha_ITT_E <- output$ES[1,1]

	d <- sapply(theta,function(x)ifelse(comp2 > x,1,0))
	ibit <- data[,tmp3]
      tc.pp <- apply(d,2,function(x) table(x,ibit))
      c.tpp1 <- tc.pp[2,]/(tc.pp[1,]+tc.pp[2,])
	t.tpp1 <-tc.pp[4,]/(tc.pp[3,]+tc.pp[4,])
	p1 <- t.tpp1 - c.tpp1 
	caceES  <- alpha_ITT_E/(p1)
	p1Table <- rbind(t.tpp1,c.tpp1,p1)
	p1Table <- data.frame(round(p1Table,2))
	row.names(p1Table )<- c("pT","pC","P=PT-pC")
	colnames(p1Table ) <- paste("P > ", theta,sep="")

	tid <- c(1:nrow(fixedDesignMatrix))
	set.seed(1020252)
	bootSamples <- NULL


	tid <- c(1:nrow(fixedDesignMatrix))
	bootSamples <- sapply(c(1:nBoot),function(x)sample(tid,replace=TRUE))

	row.names(bootSamples ) <- NULL
	bootSamples <- data.frame(bootSamples )

	bESOutput <- NULL

	for(i in 1:ncol(bootSamples )){

		bData <- data[bootSamples[,i],]

		bmf <- model.frame(formula=formula, data=bData)
		bfixedDesignMatrix <- as.matrix(data.frame(model.matrix(attr(bmf, "terms"), data=bData)))
		bposttest <- model.response(bmf)
		btrt <- bData[,which(colnames(bData)==intervention)]
		boutput <-srt(posttest=bposttest,fixedDesignMatrix=bfixedDesignMatrix,intervention=intervention,trt=btrt)
	
		balpha_ITT_E <- boutput$ES[1,1]

		bcomp <- bData[,which(colnames(bData)==compliance)]
		bd <- sapply(theta,function(x)ifelse(bcomp > x,1,0))
		bibit <- bData[,tmp3]
      	btc.pp <- apply(bd,2,function(x) table(x,bibit))
      	bc.tpp1 <- btc.pp[2,]/(btc.pp[1,]+btc.pp[2,])
		bt.tpp1 <- btc.pp[4,]/(btc.pp[3,]+btc.pp[4,])
		bp1 <- bt.tpp1 - bc.tpp1 
		bcaceES  <- balpha_ITT_E/(bp1)
		bESOutput <- rbind(bESOutput ,bcaceES)  
	}

	bootES <- apply(bESOutput ,2,function(x)quantile(x,prob=c(0.025,0.975)))
	output <- data.frame(Compliance=paste("P>",as.character(theta)),ES=round(caceES,2),LB=round(bootES[1,],2),UB=round(bootES[2,],2))
	output2 <- list(CACE=output,Compliers=p1Table )
	return(output2 ) 

}


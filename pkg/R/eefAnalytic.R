# Project: eefAnalytics
# 
# Author: lucp8394
###############################################################################

#' Analytical Package for Education Trials
#'
#' \code{eefAnalytics} provides both frequentist and Bayesian multilevel models for analyzing education
#'trials. Multilevel model or random effect model is an "omnibus" model that can be used to
#'analyse multi-site trials, cluster randomised controlled trials and simple randomised trials
#'across schools. The model reduces to ordinary linear regression when intra-cluster correlation is zero and 
#'number of pupils is the same in each school . The frequentist method relies on lme4 package and 
#'the implementation of Hedges'effect size based on unequal cluster size for 'within' and 'total' variance (Hedges, 2007).
#'The Bayesian method relies on the packages geoR for scaled-inverse-Chi-square and mvtnorm for
#'multivariate Normal distribution. The estimation are based on Gibb's sampling from the
#'full conditional posterior distributions as discussed by Wang et al. (1993). To guarantee
#'convergence of the mcmc chains, a minimum of 10,000 iterations is recommended. All parameters are based on vague
#'priors and in most cases the results should be similar to the frequentist method.
#' 
#' @author Adetayo Kasim \email{a.s.kasim@@durham.ac.uk}
#' @author ZhiMin Xiao \email{zhimin.xiao@@durham.ac.uk}
#' @author Steve Higgins \email{s.e.higgins@@durham.ac.uk}
#' 
#' @references G. Verbeke and G. Molenberghs, Linear Mixed Models for Longitudinal Data, Springer, New York, NY, USA, 2000. 
#'
#' @references Hedges, L. V. (2007). Effect Sizes in Cluster-Randomized Designs. Journal of Educational and
#'Behavioral Statistics, 32 (4), 341-370.
#'
#' @references Torgerson, D. J., & Torgerson, C. J. (2008). Designing Randomised Trials in Health, Education and the Social Sciences: An Introduction. London: Palgrave Macmillan.
#'
#' @references Wang, C., Rutledge, J., & Gianola, D. (1993). Marginal inferences about variance components in
#'a mixed linear model using Gibbs sampling. Genetics Selection Evolution, 25 , 41-62.
#'
#' @references Weinberger, M., Oddone, E., Henderson, W., Smith, D., Huey, J., Giobbie-Hurder, A., & Feussner,
#'J. (2001). Multisite Randomized Controlled Trials in Health Services Research: Scientific
#'Challenges and Operational Issues. Medical Care, 39 (6), 627-634.
#'
#' @docType package
#' @name eefAnalytics
NULL


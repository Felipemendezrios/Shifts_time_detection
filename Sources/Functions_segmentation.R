#-----------------------------------------
# Module : Function of segmentation ------
#-----------------------------------------

#===============================================================
#' Initial guess to prior information of tau :
#'
#' @param XP Data of analysis period 
#' @param nS Number of segments
#' @return Initial guess to prior information of tau 
#===============================================================
prior_tau_ini	<- function(XP,nS) {
  XP_half=0
  tstart = 0
    for (i in 1:length(XP)) {
      XP_half[i] = (XP[i+1]+XP[i])/2
    }
    if ((length(XP)-nS) <= 2) {
      tstart=0
      for (i in 1: (nS-1)) {
        tstart[i] = XP_half[i]
      }
    } else {
      step = trunc((length(XP)-1)/(nS-1), digits = 0)
      init = trunc((length(XP)-1)/2)
      tstart=0
      aaa = 0
      for (i in 1: (nS-1)) {
        if (i %% 2 == 0){
          tstart[i] = XP_half[init + step*aaa]
        } else {
          tstart[i] = XP_half[init - step*aaa]
          aaa=aaa+1
        }
      }
      tstart =sort(tstart)
    }
  return(tstart)
}

#===============================================================
#' Computation of criterion DIC for optimal model selection :  
#'
#' @param mcmc.segm MCMC simulation of each segment 
#' @param nS Number of segments
#' @return Data frame with prior and posterior information to 
#'         estimate DIC and DIC value
#===============================================================

DIC_estimation <- function(mcmc.segm,nS){
  maxpost = 0; varLogpost=0; loglikelihood = 0;  maxLikelihood = 0; varLogLikelihood=0
  
  # LogPosterior simulations:
  logpost           = as.vector(unlist(rev(mcmc.segm)[1]))
  len.mcmc          = length(mcmc.segm[,1])
  
  # LogPriors of segments means parameters "mu", change point times "tau" and Remant uncertainty : Case FlatPrior+
  
  priors.mu         = matrix(NA, nrow = len.mcmc, ncol = nS)
  priors.gamma  = rep(0,len.mcmc)
  for (i in 1:nS) { 
    priors.mu[,i]  = 0
  }
  
  if (nS > 1) {
    priors.tau = matrix(0, nrow =len.mcmc, ncol = nS-1)
  } else {
    priors.tau = matrix(0, nrow = len.mcmc, ncol = 1)
  }
  
  # Compute LogPrior as sum of LogPriors of tau parameters (change point times) and mu parameters (segments means):
  logprior = 0
  
  for (ll in 1:len.mcmc){
    logprior[ll] = sum(priors.mu[ll,]) + sum(priors.tau[ll,])  + priors.gamma[ll]
  }
  
  # Compute LogLikelihood as difference between LogPosterior and LogPrior (Bayes' theorem):
  loglikelihood = logpost - logprior
  df.mcmc.LL = data.frame(loglikelihood = loglikelihood, logprio = logprior, logposterior  = logpost)
  
  # Compute the Maximum LogLikelihood and the maximum LogPosterior, and their variances:
  maxpost[nS]          = max(logpost)
  maxLikelihood[nS]    = max(loglikelihood)
  Likelihood.maxpost   = loglikelihood[which.max(logpost)]
  varLogpost[nS]       = var(logpost)
  varLogLikelihood[nS] = var(loglikelihood)
  MeanDev              = -2*mean(loglikelihood)
  
  #ref: Gelman 2004 "Bayesian data analysis"
  DIC =  MeanDev + 2*varLogLikelihood[nS]   
  
  DIC_info <- list(df.mcmc.LL,DIC)
  return(DIC_info)
}

#===============================================================
#' Computation of criterion DIC for optimal model selection :  
#'
#' @param mcmc.segm MCMC simulation of each segment 
#' @param nS Number of segments
#' @return Data frame with prior and posterior information to 
#'         estimate DIC and DIC value
#===============================================================

# DIC_estimation <- function(mcmc.segm,nS){
#   
# }



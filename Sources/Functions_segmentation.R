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
#' Classification of gaugings for each defined period :  
#'
#' @param ts.real Tau real values according to user's choice
#' @param nS Number of segments
#' @param X Initial Dataset, before segmentation
#' @param XP Data of analysis period 
#' @return Classification of gauging and identify the boundaries
#'         limits of shift-time detected
#===============================================================

clas_gauging <- function(ts.real,nS,X,XP){
  if (!is.null(ts.real[1])) {
    tss = 0
    for (j in 1:(nS-1)) {
      for (i in 1:(length(X)-1)) {
        if((  X[i+1] >= ts.real[j]) & ((X[i] <= ts.real[j] ))) {
          tss[j] = i+1
        }
      }
    }
    t.shift.plus <- XP[tss+1]
    t.shift.before <- XP[tss]
  } else {
    tss <- NULL
    t.shift.plus <- NULL
    t.shift.before <- NULL
  }
  
  return(list(tss[j],t.shift.plus,t.shift.before))
}

#===============================================================
#' Definition of the "stable" periods of data:  
#'
#' @param tss_tot_ns Classification of all gauging through segmentations
#' @param tss Classification of gauging 
#' @param X Initial Dataset, before segmentation
#' @return Vectors with initial and final data respectively
#===============================================================

stable_period <- function(tss_tot_ns,tss,X){
  tss_tot    =  c(tss_tot_ns, tss)
  tss_tot    =  tss_tot[c(TRUE, !tss_tot[-length(tss_tot)] == tss_tot[-1])]
  tss_tot    =  sort(tss_tot)
  tss_tot_ns =  tss_tot
  i_init     =  0;
  i_final    =  0;
  for (i in 1:(length(tss_tot)-1)) {
    i_init[i]  = tss_tot[i]
    i_final[i] = tss_tot[i+1]-1
  }
  if (tss_tot[length(tss_tot)] == length(X)) {
    i_init[length(tss_tot)]  = tss_tot[length(tss_tot)]
    i_final[length(tss_tot)] = tss_tot[length(tss_tot)]
  } else {
    i_init[length(tss_tot)]  = tss_tot[length(tss_tot)]
    i_final[length(tss_tot)] = length(X)
  }
  
  stable_period_ini_final <- list(i_init,i_final,tss_tot_ns)
  
  return(stable_period_ini_final)
}

#===============================================================
#' Update iteration and period index:  
#'
#' @param nS Number of segments
#' @param i_final Index of the final of stable period 
#' @param seg.period Index of period
#' @param seg.iter Index of iteration
#' @param Y Initial observation of Dataset, before segmentation
#' @return Vectors with initial and final data respectively
#===============================================================

update_seg <- function(nS,i_final,seg.period,seg.iter,Y){
  # update iteration and period index:
  if ((nS ==1) & (i_final[seg.period] != (length(Y)))){
    seg.period = seg.period +1
    seg.iter = seg.iter + 1
    end.end = FALSE
    # segmentation is finished for this current period!!! only one segment has been found in this period!
    # we pass to another period
    
    } else if ((nS != 1)) {
    seg.period = seg.period
    seg.iter = seg.iter + 1
    end.end = FALSE
    
    } else if ((nS ==1) & (i_final[seg.period] == (length(Y)))) { #all periods have been completely segmentated ) {
    end.end = TRUE
    } else {
    seg.period = seg.period +1
    seg.iter   = seg.iter + 1    
    end.end = FALSE
  }
  
  update_seg_p <- list(seg.period,seg.iter,end.end)
  return(update_seg_p)
}

#===============================================================
#' Assign the index of the segmentation iteration and where 
#' it is located in the tree structure :  
#' ex: tree structure - c(it1=1,it2=1,it3=0,it4=0,it5=1,it6=0,it7=0)
#' c(it1=1) so, segmentation detected -> 2 segments but it analyze just one until not segmentation 
#' c(it2=1) so, segmentation detected -> analyse first part of tree structure
#' c(it3=0) so, any segmentation detected. STOP and analyze the other segment
#' c(it4=0) so, any segmentation detected. STOP and go back in the tree structure
#' c(it5=1) so, segmentation detected -> this is the another part of first segmentation. Analyse first part again
#' c(it6=0) so, any segmentation detected. STOP and analyze the other segment
#' c(it7=0) so, any segmentation detected. STOP. Not more segmentation
#' Finally, the goal is to associate the segments to a iteration.
#' So, in this example, the recursive processes will proceed as follows:
#' "*" means the connextion in the tree structure
#' firs at all, we detect all zeros (stable periods),then we go through them 
#' (iit1=1,*it2=1*,*it3=0*,it4=0,it5=1,it6=0,it7=0) -> (1,*0*,*0*,0,1,0,0) -> iter = 3 link to iter_abov = 2 and it becomes *0*
#' (*1*,0,0,*0*,1,0,0) -> (*0*,0,0,*0*,1,0,0) -> iter = 4 link to iter_abov = 1
#' (0,0,0,0,*1*,*0*,0) -> (0,0,0,0,*1*,*0*,0) -> iter = 6 link to iter_abov = 5 
#' With not more number 1 is detected -> processus finished! - link to the same iteration
#' (0,0,0,0,*1*,0,*0*) -> (0,0,0,0,*0*,0,*0*) -> iter = 7 link to iter_abov = 5
#'  
#' @param final.period Tree structure with specific writing (automatically)
#' @return Index to associate a segment with an iteration 
#===============================================================

match_seg_iter <- function(final.period){
  final.period_test <- final.period
  id_iter_above=NULL
  for(i in 1:length(which(final.period==0))){
    id_stable_mu_period <- which(final.period==0)[i]
    j <- 1
    end.end.2 <- FALSE
    while(end.end.2==FALSE){
      if(final.period_test[id_stable_mu_period-j]==1){
        id_iter_above <- rbind(id_iter_above,data.frame(iter_above=id_stable_mu_period-j,
                                                        iter=id_stable_mu_period))
        if(length(which(final.period_test==1))>1){
          final.period_test[id_stable_mu_period-j] <- 0
        }
        end.end.2 = TRUE
      }else{
        j <- j+1
      }
    }
  }
  return(id_iter_above)
}


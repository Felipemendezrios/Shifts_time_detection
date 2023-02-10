#-----------------------------------------
# Module : Function of segmentation ------
#-----------------------------------------

#============================================
# Initial guess to prior information of tau
#============================================

prior_ini_guess_tau	<- function(XP,nS) {
  XP_half=0
    for (i in 1:length(XP)) {
      XP_half[j] = (XP[i+1]+XP[i])/2
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


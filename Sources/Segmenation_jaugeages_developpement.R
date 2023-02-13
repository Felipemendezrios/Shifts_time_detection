cat("\014") # Clear console
rm(list=ls())# Clean workspace

set.seed(2023)  #Set seed for reproducibility

#Install packages

list_packages <- c("RBaM", "ggplot2")
new_packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Libraries:
library(RBaM)
library(ggplot2)

# Directories : 
dir_proj <- c('C:/Users/famendezrios/Documents/Felipe_MENDEZ/GitHub/Shifts_time_detection')
dir_data <- c('C:/Users/famendezrios/Documents/Felipe_MENDEZ/GitHub/Shifts_time_detection/Datasets')
dir_BaM <- c('BaM_segmentation')
dir_results <- c('Segmentation_results')
dir_sources <- c('Sources')

# create directories 

dir.create(file.path(dir_proj, dir_BaM),showWarnings = FALSE)
dir.create(file.path(dir_proj, dir_results),showWarnings = FALSE)

# create directories for BaM results:
dir.segmentation      = file.path(dir_proj,dir_BaM)
dir.segment.gaug      = file.path(dir_proj,dir_results) # dir. with the results of gauging segmentation
dir.sources      = file.path(dir_proj,dir_sources)
  
setwd(dir = dir.sources)

# Functions:
source("Functions_segmentation.R")

setwd(dir_data)

load(file='datasets.RData')
# dataset_ben <- datasets[[3]]

# Input data

nX=1
nY=1

nSmax = 2  # Maximum number of segments in the series at each iteration 
# nS <- 2
nCycles <- 100
burn <- 0.5
nSlim <- 10
run_option_calibration <- T
remnant_unc <- c(0,100) # prior remnant uncertainty c(min,max)

# Config_Xtra
tmin_xtra= 0
nmin_xtra= 1
option_xtra= 1



mcmc_temp=mcmcOptions(nCycles=nCycles)
cook_temp=mcmcCooking(burn=burn,nSlim=nSlim)
# remnant_temp <- list(remnantErrorModel(funk = "Constant",
#                                        par = list(parameter(name="gamma1",
#                                                             init=1 ,
#                                                             prior.dist = "Uniform",
#                                                             prior.par = c(remnant_unc[1],remnant_unc[2])))))
remnant_temp <- list(remnantErrorModel(funk = "Constant",
                                       par = list(parameter(name="gamma1",
                                                            init=1 ,
                                                            prior.dist = "FlatPrior+"))))

###################
# start iterations:
###################

for(id_dataset in 1:length(datasets)){

  # Monitor computing time 
  T1<-Sys.time()
  
  # for(id_dataset in 2:2){
  end.end <- FALSE
  
  # plots:
  seg.plot = NULL
  plot_p <- NULL
  
  ################
  #Initialisation:
  ################
  # iteration indexes: "depth search tree approach"
  i_init = i_final = level = 0;  
  i = seg.period =  seg.iter = i_init[1] = 1; 
  iteration.list = list()
  end.end = FALSE;
  tss_tot_ns = c(1); 
  final.period=NULL
  # i_final[1] = length(Q_Gaug); tss_tot_ns = c(1); 
  # final.period = acf_resid = Pacf_resid = autocorr_lag_resid = NULL; 
  ## segments means:
  mean.of.segments = mu.results.df = NULL;
  ##shift times:
  times.of.shift.MAP = t.q10 = t.q90 = t.q2 = t.q97 = ts.all.real = ts.all.real.2 = NULL; 
  ts.all.MAP = ts.all.q2 = ts.all.q10 = ts.all.q90 =  ts.all.q97 = ts.morpho.real = NULL; 
  ts.morpho.MAP = ts.morpho.q2 = ts.morpho.q97 = tau.results.df = NULL; 
  pdf.ts = cbind();           
  shift.results.df =  data.frame(tMAP  = double(), 
                                 treal = double(), 
                                 t2    = double(),
                                 t10   = double(), 
                                 t90   = double(),
                                 t97   = double())
  #dataset
  dataset_P <- datasets[[id_dataset]]
  X <- dataset_P$t
  Y <- dataset_P$obs
  Yu <- dataset_P$u
  i_final[1] = length(Y)
  seg.iter <- 1
  
  #folder creation for the current dataset:
  dir.create(paste0(dir.segment.gaug,'/data_set',id_dataset),showWarnings = FALSE)
  dir_data_set <- paste0(dir.segment.gaug,'/data_set',id_dataset)
  
  while(end.end == FALSE) {
    
    # Monitor computing time 
    T3<-Sys.time()
    
    #folder creation for the current iteration:
    dir.create(paste0(dir_data_set,"/it",seg.iter),showWarnings = FALSE)
    
    #directories saving:
    
    dir.seg.iter    <- paste0(dir_data_set,"/it",seg.iter)
    
    #observations of the current period "P":
    XP             = c(X[i_init[seg.period] : i_final[seg.period]])
    YP             = c(Y[i_init[seg.period] : i_final[seg.period]])
    YuP            = c(Yu[i_init[seg.period] : i_final[seg.period]])
    
    data_P = data.frame("X"      = XP,
                        "Y"      = YP,
                        "Yu"     = YuP,
                        "Period" = 1)
    nobs.P    = length(XP)
    
    #**************************************************************************************************************     
    # Segmentation of residuals (iterative segmentation with increasing number of segments):
    #**************************************************************************************************************
    
    # initialisation:   DIC = Deviance Information Criterion
    
    # criteria for model selection:
    DIC=0;
    npar =0; 
   
     #Number of observations:
    N = length(YP) 
    
    # Perform the segmentation only if there are at least 3 data,
    # otherwise skip: 
    ##########
    if (N>2) {
      
      i = 1 # i = number of segments
      message("Applying Segmentation of residuals:")
      message("Trying increasing number of segments. Wait ...")
      
      while ((i <= nSmax) & ((N-i) >=1)){
        nS         = i
        
        dir.create(paste0(dir.seg.iter,"/nS",nS),showWarnings = FALSE)
        dir.nS    <- paste0(dir.seg.iter,"/nS",nS)
        #=======================#
        #segmentation
        #=======================#
        npar = nS + nS - 1
        N = length(XP)
        
        residual_seg <- vector(mode = 'list',length = npar)
        
        if(nS>1){
          tstart <- prior_tau_ini(XP=XP,nS=nS)
        }
        
        for(j in 1:npar){
          if (j<=nS){
            residual_seg [[j]] <- parameter(name=paste0('res',j),
                                            init=0,
                                            prior.dist = 'FlatPrior+' ,
                                            prior.par = NULL)
          }else{
            for (z in 1:length(tstart)) {
              residual_seg [[j]] <- parameter(name=paste0('tau',z),
                                              init=tstart[z],    # a revoir
                                              prior.dist = 'FlatPrior+' ,
                                              prior.par = NULL)
              j=j+1
            }
            break
          }
          
        }
        
        # dataset object
        
        data=dataset(X=data.frame(XP),
                     Y=data.frame(YP),
                     Yu=data.frame(YuP),
                     data.dir=dir.segmentation)
        
        # Config_Xtra 
        
        xtra=xtraModelInfo(fname="Config_Xtra.txt",
                           object=c(nS,
                                    tmin_xtra,
                                    nmin_xtra,
                                    option_xtra)
                           )
        
        mod=model(fname='Config_Model.txt',
                  ID='Segmentation',
                  nX=nX,
                  nY=nY,
                  par=residual_seg,
                  xtra=xtra)
        
        
        # Run BaM executable
        BaM(mod=mod,
            data=data,
            workspace=dir.segmentation,
            run=run_option_calibration,
            mcmc=mcmc_temp,
            cook = cook_temp,
            remnant = remnant_temp
        )
        
        # Save results of segmentation:
        list.of.files.segment = c(
          paste0(dir.segmentation,"/Config_model.txt"),
          paste0(dir.segmentation,"/CalibrationData.txt"),
          paste0(dir.segmentation,"/Results_Cooking.txt"),
          paste0(dir.segmentation,"/Results_Residuals.txt"),
          paste0(dir.segmentation,"/Results_Summary.txt")
        )
        for (ll in 1:length(list.of.files.segment)) {
          file.copy(list.of.files.segment[ll], dir.nS, overwrite = TRUE)
        }
        
        mcmc.segm    <- read.table(file=paste0(dir.segmentation,"/Results_Cooking.txt"),header=TRUE)
        resid.segm   <- read.table(file=paste0(dir.segmentation,"/Results_Residuals.txt"),header=TRUE)
        summary.segm <- read.table(file=paste0(dir.segmentation,"/Results_Summary.txt"),header=TRUE)

        #=======================#
        #end of segmentation
        #=======================#
        #*****************************************************
        # COMPUTATION OF CRITERIA FOR OPTIMAL MODEL SELECTION:
        #*****************************************************

        DIC_calcul <- DIC_estimation(mcmc.segm,nS)
        df.mcmc.LL <- DIC_calcul[[1]]
        DIC [i] <- DIC_calcul[[2]]
        
        write.table(df.mcmc.LL, file=paste0(dir.nS,"/likelihood_","it",seg.iter,"_","nS",nS,".txt"), sep="\t",row.names = F)
        
        i = i+1
      }
      
      DIC_export <- rbind(DIC)
      colnames(DIC_export) <- c("ns1","nS2") 
      write.table(DIC_export,file=file.path(dir.seg.iter,'DIC.txt'),sep="\t",col.names = TRUE, row.names = F)
     
      #*******************************************************************************************************
      # Analysis of segmentation results for this period:
      #*******************************************************************************************************
      # min values:
      
      DICmin = which.min(DIC);
      
      # Optimal number of segments according to criterion:
      nS = DICmin
      
      print(paste0("=========> Optimal number of segments (considering the minimum = DIC) = ", nS))
      dir.nS.ok <- paste0(dir.seg.iter,"/nS",nS)
      
      # Selecting only results of the optimal segmentation in dir.nS.ok:
      
      Residuals       <- read.table(file=file.path(dir.nS.ok,'Results_Residuals.txt'),sep="",header=TRUE)
      mu.s            <- as.numeric(Residuals[,5])
      Results.seg     <- read.table(file=file.path(dir.nS.ok,"Results_Summary.txt"),sep="",header=TRUE)
      mcmc.segment    <- read.table(file=file.path(dir.nS.ok,"Results_Cooking.txt"),sep="",header=TRUE)
      
      #initialisation of results:
      Q2.ts = NULL; Q97.ts = NULL; Q2.mu = NULL; Q97.mu = NULL;
      
      if ( nS > 1) { # if at least one change point:
        # change point times "tau":
        #*************************
        ts.res       <- as.numeric(c(Results.seg[16,(nS+1):(2*nS-1)]))  #the maxpost of all mcmc
        ts.mean      <- as.numeric(c(Results.seg[5,(nS+1):(2*nS-1)]))  #the maxpost of all mcmc
        ts.median    <- as.numeric(c(Results.seg[6,(nS+1):(2*nS-1)]))  #the maxpost of all mcmc
        ts.res.plus  <- as.numeric(c(Results.seg[16,(nS+1):(2*nS-1)], tail(XP,1)))
        ts.res.before <- as.numeric(c(XP[1],Results.seg[16,((nS+1):(2*nS-1))]))
        stdev.ts     <- as.numeric(c(Results.seg[11,(nS+1):(2*nS-1)]))
        for (i in 1:(nS-1)) {
          Q2.ts[i]  <- c(quantile(mcmc.segment[,(nS+i)], p = c(0.025)))
          Q97.ts[i] <- c(quantile(mcmc.segment[,(nS+i)], p = c(0.975)))
        }
        Q10.ts       <- as.numeric(c(Results.seg[7,(nS+1):(2*nS-1)]))
        Q90.ts       <- as.numeric(c(Results.seg[10,(nS+1):(2*nS-1)]))
        
        # segments mean "mu":
        #********************
        mu.res       <- as.numeric(Results.seg[16,1:nS])
        mu.mean      <- as.numeric(c(Results.seg[5,1:nS]))  #the maxpost of all mcmc
        mu.median    <- as.numeric(c(Results.seg[6,1:nS]))  #the maxpost of all mcmc
        for (j in 1:nS) {
          Q2.mu[j]  <- c(quantile(mcmc.segment[,j], p = c(0.025)))
          Q97.mu[j] <- c(quantile(mcmc.segment[,j], p = c(0.975)))
        }
        Q10.mu.res   <- as.numeric(Results.seg[7,1:(nS)])
        Q90.mu.res   <- as.numeric(Results.seg[10,1:(nS)])
        stdev.mu     <- as.numeric(Results.seg[11,1:(nS)])
        seg.iter.shift <- seg.iter
        #structural error parameter "gamma":
        #**********************************
        gamma.MAP    <- as.numeric(Results.seg[16,2*nS])
        gamma.stdev  <- as.numeric(Results.seg[11,2*nS])
        gamma.mean   <- as.numeric(Results.seg[5,2*nS])
        final.period = c(final.period,1)
        
        ######################################################   
      } else { # if no change points ==> only one segment !!
        ######################################################
        ts.res = ts.res.before = ts.res.plus = ts.mean = ts.median = NULL;
        stdev.ts = NULL; seg.iter.shift = NULL
        Q2.ts = Q97.ts = Q10.ts = Q90.ts = NULL;
        #-------------------------------------------------------
        mu.res       <- as.numeric(Results.seg[16,1])
        mu.mean      <- as.numeric(Results.seg[5,1])  #the maxpost of all mcmc         
        mu.median    <- as.numeric(Results.seg[6,1])  #the maxpost of all mcmc
        Q2.mu        <- quantile(mcmc.segment[,1], p = c(0.025))
        Q97.mu       <- quantile(mcmc.segment[,1], p = c(0.975))
        Q10.mu.res   <- as.numeric(Results.seg[7,1])
        Q90.mu.res   <- as.numeric(Results.seg[10,1])
        stdev.mu     <- as.numeric(Results.seg[11,1])
        #-------------------------------------------------------
        gamma.MAP    <- as.numeric(Results.seg[16,2])
        gamma.stdev  <- as.numeric(Results.seg[11,2])
        gamma.mean   <- as.numeric(Results.seg[5,2])
        final.period = c(final.period,0)
      }
      
      #saving shift times and segments statistics to two data.frames:
      tau.results.df[[seg.iter]] =  data.frame(# change point times "tau":
        tau.MAP    = ts.res, 
        tau.q2     = Q2.ts, 
        tau.q10    = Q10.ts, 
        tau.q90    = Q90.ts, 
        tau.q97    = Q97.ts, 
        tau.stdev  = stdev.ts, 
        tau.mean   = ts.mean, 
        tau.median = ts.median,
        tau.before = ts.res.before[1],
        tau.shift  = ts.res.before[2],
        tau.plus   = ts.res.plus[2],
        iteration  = seg.iter.shift)
      mu.results.df[[seg.iter]]  =  data.frame(# Segment mean "mu":
        mu.MAP     = mu.res, 
        mu.q2      = Q2.mu, 
        mu.q10     = Q10.mu.res,
        mu.q90     = Q90.mu.res, 
        mu.q97     = Q97.mu,
        mu.stdev   = stdev.mu, 
        mu.mean    = mu.mean,  
        mu.median  = mu.median,
        mu.iteration  = seg.iter)
      
      # Saving to vectors (updating the existing ones):     # THIS PART NEEDS IMPROVEMENT !!!
      times.of.shift.MAP <- c(times.of.shift.MAP, ts.res)
      mean.of.segments   <- c(mean.of.segments,mu.res)
      t.q2   <- c(t.q2, Q2.ts)
      t.q10  <- c(t.q10, Q10.ts)
      t.q90  <- c(t.q90, Q90.ts)
      t.q97  <- c(t.q97, Q97.ts)
      t.q2   <- t.q2[c(TRUE, !t.q2[-length(t.q2)] == t.q2[-1])]
      t.q10  <- t.q10[c(TRUE, !t.q10[-length(t.q10)] == t.q10[-1])]
      t.q90  <- t.q90[c(TRUE, !t.q90[-length(t.q90)] == t.q90[-1])]
      t.q97  <- t.q97[c(TRUE, !t.q97[-length(t.q97)] == t.q97[-1])]
      times.of.shift.MAP <- times.of.shift.MAP[c(TRUE, !times.of.shift.MAP[-length(times.of.shift.MAP)] == times.of.shift.MAP[-1])]
      mean.of.segments   <-  mean.of.segments[c(TRUE, !mean.of.segments[-length( mean.of.segments)] ==  mean.of.segments[-1])]
      
      #*********************************************************************
      # Shift times adjustment ==> Identification of the "true" shift times:
      #*********************************************************************
      ts.real  = NULL; 
      # Looking inside the u95% interval of the shift times.
      #
      # Options:  1) if there is a flood, assign the shift time to the max peak !
      #           2) if the shift is due to other causes (vegetation, works, ice ...)
      #              then assign the shift time to the maxpost
      #           3) if the shift time is known then fix it.
      
      i =0
      if (nS >1) {  
        
        while (i < length(tau.results.df[[seg.iter]]$tau.MAP)) {
          i = i +1
          message(paste0("shift time ", i, " = Tau MAP = ", tau.results.df[[seg.iter]]$tau.MAP[i]))
          ts.real[i] = tau.results.df[[seg.iter]]$tau.MAP[i] 
          
        } 
        
        # Saving shift times dataframes:
        df.shift.times = data.frame(Q2.ts, Q97.ts, ts.res, ts.real)
        df.shift.times.plus = data.frame(ts.res.before, ts.res.plus, Q2.mu, Q10.mu.res, Q90.mu.res, Q97.mu,  mu.res)
        
        # save results into txt files:
        write.table(mcmc.segment, file.path(dir.seg.iter,"mcmc_segmentation.txt"), sep ="\t", row.names=FALSE)
        write.table(tau.results.df[[seg.iter]],paste0(dir.seg.iter,"/df_tau_it", seg.iter,".txt"),sep ="\t", row.names=FALSE)
        write.table(mu.results.df[[seg.iter]], paste0(dir.seg.iter,"/df_mu_it", seg.iter,".txt"), sep ="\t", row.names=FALSE)
        write.table(df.shift.times, paste0(dir.seg.iter,"/df.shift.times_it", seg.iter,".txt"), sep ="\t", row.names=FALSE)
        write.table(df.shift.times.plus, paste0(dir.seg.iter,"/df.shift.times.plus_it", seg.iter,".txt"), sep ="\t", row.names=FALSE)
        
      }else{
        # Only one segment (no shift detected!)
        # Save results into txt files:
        write.table(mcmc.segment, paste0(dir.seg.iter,"/mcmc_segmentation.txt"), sep ="\t", row.names=FALSE)
        write.table(mu.results.df[[seg.iter]], paste0(dir.seg.iter,"/df_mu_it", seg.iter,".txt"), sep ="\t", row.names=FALSE)
        
      }
      
      #------------------------------------------------------------------------------------------
      # Arrange the detected shift times (this part needs to be cleaned!):
      if (!is.null(ts.real[1])) {
        ts.real                    = as.numeric(ts.real)
        ts.all.real                = c(ts.all.real, ts.real)
        ts.all.real                = sort(ts.all.real)
        ts.all.real                = ts.all.real[c(TRUE, !ts.all.real[-length(ts.all.real)] == ts.all.real[-1])]
        ts.all.real.2              = c(ts.all.real.2, ts.real)
        ts.all.real.2              = sort(ts.all.real.2)
        tau.results.df[[seg.iter]] = cbind(tau.results.df[[seg.iter]],  treal = ts.real)         
        shift.results.df           = rbind(shift.results.df , tau.results.df[[seg.iter]])
      } else {
        tau.results.df[[seg.iter]] = FALSE;
      }
      
      
      # #****************************************************
      # # Classification of gaugings for each defined period
      # #****************************************************
      
      classification_gaug <- clas_gauging(ts.real,nS,X,XP)
      
      tss <- classification_gaug[[1]]
      t.shift.plus <- classification_gaug[[2]] 
      t.shift.before <- classification_gaug[[3]]
     
      #***********************************************
      #definition of the "stable" periods of data:
      #***********************************************
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
      #***********************
      # end of segmentation:
      #**********************
      
      #**********************************************************************
      # Plot residuals time series with segmentation results :
      #**********************************************************************
      if (!is.null(ts.res[1])) {
        seg.plot[[seg.iter]] <- ggplot()+
          geom_point(data=data_P, aes(x=X,y=Y))+
          geom_errorbar(data=data_P,
                        aes(x =X,
                            ymin=Y-Yu,
                            ymax=Y+Yu),width=0.2)+
          geom_rect( mapping= aes(xmin = Q2.ts,
                                  xmax = Q97.ts,
                                  ymin = -Inf, 
                                  ymax = Inf),fill='#0b53c1', 
                     alpha=0.1 )+
          geom_vline(xintercept = ts.res,colour='#0b53c1',linetype='solid')+
          geom_vline(xintercept = ts.real,colour='#0b53c1',linetype='longdash')+
          geom_rect(mapping = aes(xmin = ts.res.before, 
                                  xmax = ts.res.plus, 
                                  ymin = Q2.mu ,
                                  ymax = Q97.mu),
                    fill = '#CB2027', alpha=0.3) + 
          geom_segment(mapping=aes(x    = ts.res.before, 
                                   y    = mu.res, 
                                   xend = ts.res.plus, 
                                   yend = mu.res),
                       col  = '#CB2027')+
          coord_cartesian(clip = 'off')+
          labs(x = 'Time [ ]', y='Measurement [ ]')+
          theme_bw()
      }else{
        seg.plot[[seg.iter]] <- ggplot()+
          geom_point(data=data_P, aes(x=X,y=Y))+
          geom_errorbar(data=data_P,
                        aes(x =X,
                            ymin=Y-Yu,
                            ymax=Y+Yu),width=0.2)+
          geom_rect(mapping = aes(xmin = data_P$X[1], 
                                  xmax = rev(data_P$X)[1], 
                                  ymin = Q2.mu ,
                                  ymax = Q97.mu),
                    fill = '#CB2027', alpha=0.3) + 
          geom_segment(mapping=aes(x = data_P$X[1], 
                                   xend = rev(data_P$X)[1], 
                                   y    = mu.res, 
                                   yend = mu.res),
                       col  = '#CB2027')+
          coord_cartesian(clip = 'off')+
          labs(x = 'Time [ ]', y='Measurement [ ]')+
          theme_bw()
        }
      
        ggsave(file=paste0(dir.seg.iter,"/segmentation_it",seg.iter,".png"))
    
        # update iteration and period index:
        if ((nS ==1) & (i_final[seg.period] != (length(Y)))){
          seg.period = seg.period +1
          seg.iter = seg.iter + 1
          # segmentation is finished for this current period!!! only one segment has been found in this period!
          # we pass to another period
        } else if ((nS != 1)) {
          seg.period = seg.period
          seg.iter = seg.iter + 1
          
        } else if ((nS ==1) & (i_final[seg.period] == (length(Y)))) { #all periods have been completely segmentated ) {
          end.end = TRUE
          #final.period = 1
          #segmentation is finished. stop !!!!!
          print("segmentation finished!")
          
        } else {
          seg.period = seg.period +1
          seg.iter   = seg.iter + 1         
          
        }
        # Monitor computing time 
        T4<-Sys.time()
        Tdiff=difftime(T4,T3)
        write.table(Tdiff,paste0(dir.seg.iter,"/computing_time.txt"),row.names = F, col.names = F)
    }else{
      print(paste0("No segmentation: there are only one or two points to segment for this period !!"))
      final.period = c(final.period,0)
      if (tail(XP,1) == (tail(X,1))) {
        end.end = TRUE
        mu.results.df[[seg.iter]]  =  data.frame(# Segment mean "mu":
          mu.MAP     = mean(data_P$Y), 
          mu.q2      = mean(data_P$Y)-mean(data_P$Yu), 
          mu.q10     = 0,
          mu.q90     = 0, 
          mu.q97     = mean(data_P$Y)+mean(data_P$Yu),
          mu.iteration = seg.iter)
        
        write.table(mu.results.df[[seg.iter]], paste0(dir.seg.iter,"/df_mu_it", seg.iter,".txt"), sep ="\t", row.names=FALSE)
        
        seg.plot[[seg.iter]] <- ggplot()+
          geom_point(data=data_P, aes(x=X,y=Y))+
          geom_errorbar(data=data_P,
                        aes(x =X,
                            ymin=Y-Yu,
                            ymax=Y+Yu),width=0.2)+
          geom_rect(mapping = aes(xmin = 0.9*data_P$X[1], 
                                  xmax = 1.1*rev(data_P$X)[1], 
                                  ymin = mu.results.df[[seg.iter]]$mu.q2 ,
                                  ymax = mu.results.df[[seg.iter]]$mu.q97),
                    fill = '#CB2027', alpha=0.1) + 
          geom_segment(mapping=aes(x = 0.9*data_P$X[1], 
                                   xend = 1.1*rev(data_P$X)[1], 
                                   y    = mu.results.df[[seg.iter]]$mu.MAP, 
                                   yend = mu.results.df[[seg.iter]]$mu.MAP),
                       col  = '#CB2027')+
          coord_cartesian(clip = 'off')+
          labs(x = 'Time [ ]', y='Measurement[ ]')+
          theme_bw()
      
      ggsave(file=paste0(dir.seg.iter,"/segmentation_it",seg.iter,".png"))
      
        #segmentation is finished. stop !!!!!
        print("segmentation finished!")
        
      } else {
        mu.results.df[[seg.iter]]  =  data.frame(# Segment mean "mu":
          mu.MAP     = mean(data_P$Y), 
          mu.q2      = mean(data_P$Y)-mean(data_P$Yu), 
          mu.q10     = 0,
          mu.q90     = 0, 
          mu.q97     = mean(data_P$Y)+mean(data_P$Yu),
          mu.iteration = seg.iter)
        
        write.table(mu.results.df[[seg.iter]], paste0(dir.seg.iter,"/df_mu_it", seg.iter,".txt"), sep ="\t", row.names=FALSE)
        
        seg.plot[[seg.iter]] <- ggplot()+
          geom_point(data=data_P, aes(x=X,y=Y))+
          geom_errorbar(data=data_P,
                        aes(x =X,
                            ymin=Y-Yu,
                            ymax=Y+Yu),width=0.2)+
          geom_rect(mapping = aes(xmin = 0.9*data_P$X[1], 
                                  xmax = 1.1*rev(data_P$X)[1], 
                                  ymin = mu.results.df[[seg.iter]]$mu.q2 ,
                                  ymax = mu.results.df[[seg.iter]]$mu.q97),
                    fill = '#CB2027', alpha=0.1) + 
          geom_segment(mapping=aes(x = 0.9*data_P$X[1], 
                                   xend = 1.1*rev(data_P$X)[1], 
                                   y    = mu.results.df[[seg.iter]]$mu.MAP, 
                                   yend = mu.results.df[[seg.iter]]$mu.MAP),
                       col  = '#CB2027')+
          coord_cartesian(clip = 'off')+
          labs(x = 'Time [ ]', y='Measurement[ ]')+
          theme_bw()
        
        ggsave(file=paste0(dir.seg.iter,"/segmentation_it",seg.iter,".png"))
        
        seg.period = seg.period +1
        seg.iter   = seg.iter + 1         
      }
      # Monitor computing time 
      T4<-Sys.time()
      Tdiff=difftime(T4,T3)
      write.table(Tdiff,paste0(dir.seg.iter,"/computing_time.txt"),row.names = F, col.names = F)
    }
  }
  
  if(seg.iter!=1){
    # times.of.shift.final <- c(times.of.shift.MAP, tail(Shift.Q$alpha_t,1))
    shift.times.gaugings <- data.frame(tMAP  = times.of.shift.MAP, 
                                       treal = ts.all.real, 
                                       t2    = t.q2,
                                       t10   = t.q10, 
                                       t90   = t.q90,
                                       t97   = t.q97)
    
    shift.results.df = shift.results.df[order(shift.results.df$treal),]
    shift.times.gaugings2 <- data.frame(tMAP  = shift.results.df$tau.MAP, 
                                        treal = shift.results.df$treal, 
                                        t2    = shift.results.df$tau.q2,
                                        t10   = shift.results.df$tau.q10, 
                                        t90   = shift.results.df$tau.q90,
                                        t97   = shift.results.df$tau.q97,
                                        t.before = shift.results.df$tau.before,
                                        t.shift = shift.results.df$tau.shift,
                                        t.plus   = shift.results.df$tau.plus,
                                        iter  = shift.results.df$iteration)
    
    write.table(shift.times.gaugings2, paste0(dir_data_set,"/shift_times.txt"), sep ="\t", row.names=FALSE)
    # tau.results.df
    stable_mu_period = NULL

    for(i in 1:length(which(final.period==0))){
      id_stable_mu_period <- which(final.period==0)[i]
      
      stable_mu_period <- rbind(stable_mu_period,data.frame(muMAP  = mu.results.df[[id_stable_mu_period]]$mu.MAP,
                                                            mu2    = mu.results.df[[id_stable_mu_period]]$mu.q2,
                                                            mu10   = mu.results.df[[id_stable_mu_period]]$mu.q10,
                                                            mu90   = mu.results.df[[id_stable_mu_period]]$mu.q90,
                                                            mu97   = mu.results.df[[id_stable_mu_period]]$mu.q97,  
                                                            iter   = mu.results.df[[id_stable_mu_period]]$mu.iteration))
    }
    
    # final.period <- c(1,1,0,1,1,0,1,0,0,0,0)
    # final.period <- c(1,1,1,0,0,1,1,0,0,0,1,0,0)
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
      # if (j >= 2){
      #   final.period_test[id_stable_mu_period-j] <- 0
      # }else if(final.period_test[id_stable_mu_period+1]==1){
      #   final.period_test[id_stable_mu_period-1] <- 0
      #   
      # }
    }
    
    stable_mu_period_t <- merge(stable_mu_period,id_iter_above,by='iter')
    data_result <- merge(stable_mu_period_t,shift.times.gaugings2, by.x = 'iter_above', by.y = 'iter', sort=F)
    
    write.table(data_result, paste0(dir_data_set,"/tau_mu_result.txt"), sep ="\t", row.names=FALSE)
    
    
    ### Plots
    
    
    plot_p <- ggplot(dataset_P, aes(x=t,y=obs))+
      geom_point()+
      geom_errorbar(aes(ymin=obs-u,ymax=obs+u),width=0.2)+
      geom_vline(xintercept = shift.times.gaugings2$treal,colour='#0b53c1',linetype='solid')+
      theme_bw()
    
    for(i in 1:nrow(shift.times.gaugings2)){
      plot_p <- plot_p+geom_rect(mapping = aes_string(xmin = shift.times.gaugings2$t2[i], 
                              xmax = shift.times.gaugings2$t97[i], 
                              ymin = -Inf, 
                              ymax = Inf),fill='#0b53c1',alpha=0.005 )
    }
    id_column_before <- which(colnames(data_result)=='t.before')
    id_column_plus <- which(colnames(data_result)=='t.plus')
    id_column_real <- which(colnames(data_result)=='t.shift')
    for (i in 1:nrow(data_result)){
      if (i==1){
        xlim_plot_min <- data_result[i,id_column_before]
        xlim_plot_max <- data_result[i,id_column_real]
      }else if(i==nrow(data_result)){
        xlim_plot_min <- data_result[i,id_column_real]
        xlim_plot_max <- data_result[i,id_column_plus]
      }else{
        xlim_plot_min <- data_result[i-1,id_column_real]
        xlim_plot_max <- data_result[i,id_column_real]
      }
      
      
      plot_p <- plot_p +
        geom_segment(mapping=aes_string(x    = xlim_plot_min,
                                        xend = xlim_plot_max, 
                                        y    = data_result$muMAP[i], 
                                        yend = data_result$muMAP[i]))+
        geom_rect(mapping = aes_string(xmin = xlim_plot_min, 
                                xmax = xlim_plot_max, 
                                ymin = data_result$mu2[i],
                                ymax = data_result$mu97[i]),
                  fill = '#CB2027', alpha=0.005)
    }
    
    ggsave(file=paste0(dir_data_set,"/segmentation.png"))
  
    }else{
      
      stable_mu_period  <- rbind(data.frame(muMAP  = mu.results.df[[1]]$mu.MAP,
                                            mu2    = mu.results.df[[1]]$mu.q2,
                                            mu10   = mu.results.df[[1]]$mu.q10,
                                            mu90   = mu.results.df[[1]]$mu.q90,
                                            mu97   = mu.results.df[[1]]$mu.q97,  
                                            iter   = mu.results.df[[1]]$mu.iteration,
                                            iter_above = 1))
      
      write.table(stable_mu_period, paste0(dir_data_set,"/mu_result.txt"), sep ="\t", row.names=FALSE)
      
      ### Plots
      
      plot_p <- ggplot(dataset_P, aes(x=t,y=obs))+
        geom_point()+
        geom_errorbar(aes(ymin=obs-u,ymax=obs+u),width=0.2)+
        geom_segment(mapping=aes_string(x    = sort(dataset_P$t)[1],
                                        xend = sort(dataset_P$t,decreasing =T)[1], 
                                        y    = stable_mu_period$muMAP, 
                                        yend = stable_mu_period$muMAP))+
        geom_rect(mapping = aes_string(xmin    = sort(dataset_P$t)[1],
                                       xmax = sort(dataset_P$t,decreasing =T)[1], 
                                       ymin = stable_mu_period$mu2,
                                       ymax = stable_mu_period$mu97),
                  fill = '#CB2027', alpha=0.005)+
        theme_bw()
      
      ggsave(file=paste0(dir_data_set,"/segmentation.png"))
      
    }
  
  T2<-Sys.time()
  Tdiff= difftime(T2,T1)
  write.table(Tdiff, file=paste0(dir_data_set,"/computing_time.txt"),row.names = F, col.names = F)
  
}





cat("\014") # Clear console
rm(list=ls())# Clean workspace
dev.off()

library(ggplot2);library(patchwork)

# Directories : 
dir_proj <- c('C:/Users/famendezrios/Documents/Felipe_MENDEZ/GitHub/Shifts_time_detection')
dir_results <- c('Datasets')

setwd(file.path(dir_proj,dir_results))

set.seed(76543210)
nChange=rep(c(0,1,5,10),3)
nObs=rep(c(30,60,100),each=4)
duration=10
m0=100
mChange=m0/10
sChange=mChange/2
mU=sChange/4
sU=0.4

datasets=truth=g=g0=vector('list',length(nChange))

for(k in 1:length(nChange)){
  tObs=sort(runif(nObs[k],0,duration))
  tChange=c(0,sort(runif(nChange[k],0,duration)))
  change=c(0,(-1)^rbinom(nChange[k],size=1,prob=0.5)*rnorm(nChange[k],mean=mChange,sd=sChange))
  mu=m0+cumsum(change)
  nPeriod=nChange[k]+1
  
  period=rep(0,nObs[k])
  for(i in 1:nPeriod){
    mask=(tObs>tChange[i])
    period[mask]=period[mask]+1
  }
  
  obs=rep(NA,nObs[k])
  uobs=rep(NA,nObs[k])
  for(i in 1:nObs[k]){
    uobs[i]=rlnorm(1,meanlog=log(mU),sdlog=sU)
    obs[i]=rnorm(1,mean=mu[period[i]],sd=uobs[i])
  }
  
  datasets[[k]]=data.frame(t=tObs,obs=obs,u=uobs)
  truth[[k]]=data.frame(t=tObs,obs=obs,u=uobs,period=as.factor(period))
  
  g[[k]]=ggplot(truth[[k]],aes(x=t))+geom_point(aes(y=obs,color=period))+
    geom_errorbar(aes(ymin=obs-1.96*u,ymax=obs+1.96*u,color=period))+
    theme_bw()+
    ggtitle(paste0('Dataset ',k))+
    theme(plot.title = element_text(hjust = 0.5))
  
  if(length(tChange)>0){
    g[[k]]=g[[k]]+
      geom_vline(xintercept = tChange[-1],colour='#0b53c1',linetype='solid')
  }
  
  g0[[k]]=ggplot(truth[[k]],aes(x=t))+geom_point(aes(y=obs))+
    geom_errorbar(aes(ymin=obs-1.96*u,ymax=obs+1.96*u))+
    theme_bw()
}

pdf(file='truth.pdf',width=16,height=9,useDingbats=F)
wrap_plots(g,ncol=4)
dev.off()

tiff(file="truth.tiff",
     res = 300,
     width = 400,
     height = 200,
     units = c('mm'),)
wrap_plots(g,ncol=4)
dev.off()

pdf(file='datasets.pdf',width=16,height=9,useDingbats=F)
wrap_plots(g0,ncol=4)
dev.off()

save(truth,g,file='truth.RData')
save(datasets,file='datasets.RData')
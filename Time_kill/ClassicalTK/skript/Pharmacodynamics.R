###################################################################################################################################
##################################fit pharmacodynamic function and impute missing data#############################################
###################################################################################################################################

#set working directory  
setwd("C:/Users/sfoerster/ABinteractions/Time_kill/GrowthRate")

# Negative log-likelihood
nll <- function(a, b, sigma, data = data2regress) {
  m <- a + b * data$x
  m <- ifelse(m == 0, 1e-5, m) # to avoid divisions by 0 (not really clean, but this should happen very rarely)
  #v <- (sigma / m)^2           # to avoid issues if m < 0 in dnorm, pnorm calls below  
  #v <- (sigma)^2              # constant error (additive on log scale)
  v <- (sigma/data$y_cens)^2             # constant error (additive on log scale)
  ll_i <- (1 - data$cens) * dnorm(data$y_cens, mean = m, sd = sqrt(v), log = TRUE) + 
  data$cens * pnorm(data$y_cens, mean = m, sd = sqrt(v), log.p = TRUE) 
  ll <- sum(ll_i)
  return(-ll)
}

#load libraries
library("drc")
library("XLConnect")
library("metafor")
#initialize some empty vectors and empty variables
parmlist=vector()
append2=vector()
append3=vector()
mydata=vector()
means=vector()
i <- 1
j <- 1
k <- 1

#define patterns to be recognized
sourcefiledir<-"input"
verzeichnis<-"input"
sourcefilepattern<-".txt"

#define some colors
mypalette<- c("red","lightblue","lightblue","lightblue","#C6DBEF", "#9ECAE1" ,"#6BAED6" ,"#4292C6" ,"#2171B5" ,"#08519C", "#08306B","black","black","black","black","black")

# list of all files and pattern to be recognized in these files
files <- list.files(path =sourcefiledir, pattern = sourcefilepattern, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
ablist<-c("WHOA_AZ","WHOA_CIP","WHOA_Tet","WHOA_CT","F89_AZ","F89_CIP","F89_Tet","F89_CT")


#do action for the dataset with the recognized pattern
for (ab in ablist){
  replicatespattern=ab
  replicates <- list.files(path =sourcefiledir, pattern = ab, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)

  #read data from file with the recognized pattern and append to a dataframe
  for (replicate in replicates){
    append1 <- read.table(paste(verzeichnis,"/",replicate,sep=""),header=TRUE)
    append2 <- cbind(append1,replicate,ab,i)
    mydata<-rbind(append2,mydata) 
    i=i+1
  }
}
problem<-subset(mydata,mydata$hours==0&mydata$count==0)
rows_to_garbage1 = which(mydata$sample==toString(problem$sample[1]), useNames = TRUE)
mydata <- mydata[-c(rows_to_garbage1),]
rows_to_garbage2 = which(mydata$sample==toString(problem$sample[2]), useNames = TRUE)
mydata <- mydata[-c(rows_to_garbage2),]
for(j in unique(mydata$replicate)){
  print("----------------------------filename------------------------------------------")
  data=mydata[mydata$replicate==j,]
  name=sub(".txt",".pdf",j)
  pdf(paste(name,sep=""))
  hours=data$hours
  time=data$hours
  CFU=data$count*data$dilution*100
  CFU[CFU==0&time==0]<-NA
  #cens=CFU
  #cens[cens==0]<-1
  #cens[cens>1]<-0
  CFU[CFU==0]<-100
  #options(scipen = 10)
  antibiotic=data$antibiotic.conc
  ab.conc<-factor(antibiotic)
  conc=round(data$antibiotic.conc,3)
  
  #conc=format(conc,digits=3,nsmall=5)

  conc2=mydata$antibiotic.conc
  conc2[conc2==min(conc)]=0
  conc[conc==min(conc)]=0
  conc=factor(conc)
  conc=levels(conc)
  #dat=data.frame(cbind(x=time, cens=cens, y_cens=CFU,ab=antibiotic))
  #dat[, c('x', 'cens', 'y_cens')]
  tmin<-0 
  tmax<-6 
  # define a vector for the decline rates:
  slopes<-vector("numeric",length=nlevels(ab.conc))
  #Calculate the linear regressions
 # for(idose in 1:nlevels(ab.conc)){
  #  dose<-as.numeric(levels(ab.conc)[idose])
   # data2regress<-data.frame(time=hours[hours>=tmin & hours<=tmax &antibiotic==dose],log.bacteria.rel=log(CFU[hours>=tmin&hours<=tmax & antibiotic==dose]/CFU[hours==tmin&antibiotic==dose]))
    #slopes[idose]<-coef(lm(log.bacteria.rel~time-1,data2regress))
  #}
  
  #start values a=log(CFU) at time 0, b=0
  slopes<-vector()

  
  for(idose in 1:nlevels(ab.conc)){
    dose<-as.numeric(levels(ab.conc)[idose])
    #data2regress=data.frame(cbind(x=time, cens=cens, y_cens=CFU,ab=dose)) 
    cens=CFU[hours>=tmin&hours<=tmax & antibiotic==dose]
    cens[cens==0]<-1
    cens[cens>1]<-0
    y_cens=log(CFU[hours>=tmin & hours<=tmax &antibiotic==dose])
    #data2regress<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=log(CFU[hours>=tmin&hours<=tmax & antibiotic==dose]/CFU[hours==tmin&antibiotic==dose]),cens=cens)
    data2points<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=CFU[hours>=tmin & hours<=tmax &antibiotic==dose],cens=cens)
    data2regress<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=y_cens,cens=cens)
    
    #names(data2regress)<-c("x")

    #slopes[idose]<-coef(lm(log.bacteria.rel~time-1,data2regress))
    fit <- mle2(nll, start=list(a=max(y_cens), b=-5, sigma = 0.1), method="Nelder-Mead", control=list(maxit=1e3))
    print("---------------------------------------------------------------------------------------------------")
    slope <- data.frame(coef(fit)[[2]],dose)
    predict<-coef(fit)[[1]]+data2regress$x*coef(fit)[[2]]  
    #points(data2points$x,data2points$y_cens,col=c(idose),pch=1,cex=1.5)
    #lines(data2points$x,data2points$y_cens,col=c(idose),pch=1,cex=1.5)
    #lines(data2regress$x,exp(predict),col=c(idose),lty=2) 
    slopes <- rbind(slope,slopes)
  }
  names(slopes)<-c("slope","dose")
  dose=as.numeric(levels(ab.conc))
  library("drc")
  model=drm(slopes$slope~slopes$dose,fct=LL.4())
  plot(model,log="x", xlim=c(min(conc2),max(conc2)),ylim=c(-2.5,1.5),xlab="Antibiotic concentration [mg/L]",ylab="Growth rate [per hour]",pch=21, bg=mypalette, cex=1.4,lwd=1.3,cex.lab=1.3)
  legend("bottomleft",conc, bty="n", text.col="gray47",pt.bg=mypalette,pch=21,cex=0.8,title="conc.[mg/L]")
  #legend("bottomleft",legend=i,cex=2,bty="n")
  library(plotrix)
  axis.break(axis=1,breakpos=min(data$antibiotic.conc)*2,style="slash") 
  #extract parameters from model summary
  kappa=coefficients(model)[1]
  kappaerr=coef(summary(model))[1,2]
  upper=coefficients(model)[3]
  uppererr=coef(summary(model))[3,2]
  lower=coefficients(model)[2]
  lowererr=coef(summary(model))[2,2]
  inflection=coefficients(model)[4]
  inflectionerr=coef(summary(model))[4,2]
  MIC=ED(model,type="absolute",interval="none",method="grid",respLev=0)
  zMIC=MIC[1]
  zMICerr=MIC[2]
  experiment=mydata$ab[mydata$replicate==j]
#   print(experiment)
  Liste=data.frame(j,toString(experiment[1]),kappa,kappaerr,upper,uppererr,lower,lowererr,inflection,inflectionerr,zMIC,zMICerr)
  parmlist=rbind(Liste,parmlist)
  dev.off()
  k=k+1
}
#names(parmlist)<-c("Identifier","Experiment","kappa","kappaerr","upper","uppererr","lower","lowererr","inflection","inflectionerr","zMIC","zMICerr")
#mean=data.frame(experiment[1],mean(parmlist$upper),mean(parmlist$upper),mean(parmlist$lower),mean(parmlist$inflection),mean(parmlist$zMIC),mean(parmlist$MICU),mean(parmlist$MICL))  
#means <-aggregate(parmlist, by=list(parmlist$experiment.1.),FUN=mean, na.rm=TRUE)
#print(aggdata)
#means=c()
# for (l in unique(parmlist$Experiment)){
#   parmlistsub <-subset(parmlist,parmlist$Experiment==l)
#   kappamodel <- rma(yi = parmlistsub$kappa, sei = parmlistsub$kappaerr, data = parmlistsub, method = 'FE')
#   kappasum <- data.frame(cbind(l,"kappa",summary(kappamodel)[1],summary(kappamodel)[2],summary(kappamodel)[3],summary(kappamodel)[4],summary(kappamodel)[5],summary(kappamodel)[6]))
#   uppermodel <- rma(yi = parmlistsub$upper, sei = parmlistsub$uppererr, data = parmlistsub, method = 'FE')
#   uppersum <- data.frame(cbind(l,"upper",summary(uppermodel)[1],summary(uppermodel)[2],summary(uppermodel)[3],summary(uppermodel)[4],summary(uppermodel)[5],summary(uppermodel)[6]))
#   lowermodel <- rma(yi = parmlistsub$lower, sei = parmlistsub$lowererr, data = parmlistsub, method = 'FE')
#   lowersum <- data.frame(cbind(l,"lower",summary(lowermodel)[1],summary(lowermodel)[2],summary(lowermodel)[3],summary(lowermodel)[4],summary(lowermodel)[5],summary(lowermodel)[6]))
#   inflectionmodel <- rma(yi = parmlistsub$inflection, sei = parmlistsub$inflectionerr, data = parmlistsub, method = 'FE')
#   inflectionsum <- data.frame(cbind(l,"inflection",summary(inflectionmodel)[1],summary(inflectionmodel)[2],summary(inflectionmodel)[3],summary(inflectionmodel)[4],summary(inflectionmodel)[5],summary(inflectionmodel)[6]))
#   zMICmodel <- rma(yi = parmlistsub$zMIC, sei = parmlistsub$zMICerr, data = parmlistsub, method = 'FE')
#   zMICsum <- data.frame(cbind(l,"zMIC",summary(zMICmodel)[1],summary(zMICmodel)[2],summary(zMICmodel)[3],summary(zMICmodel)[4],summary(zMICmodel)[5],summary(zMICmodel)[6]))
#   
#   mean=rbind(kappasum,uppersum,lowersum,inflectionsum,zMICsum)
#   means=rbind(mean,means)
# } 
#  names(means) <- c("experiment","parameter","estimate","se","zval","pval","ci.lb","ci.ub")
#  #write data to files, parameters from singleexperiments and mean of all parameters
#  data=data.frame(parmlist)
#  allplots <- loadWorkbook("summary_statistics5.xlsx", create = TRUE)  
#  createSheet(allplots, name = "parameterlist")
#  writeWorksheet(allplots, means, sheet = "parameterlist", startRow = 1, startCol = 1)
#  saveWorkbook(allplots)
# # # 
#   singleplots <- loadWorkbook("singleplot_statistics5.xlsx", create = TRUE)  
#   createSheet(singleplots, name = "parameterlist")
#   writeWorksheet(singleplots,parmlist, sheet = "parameterlist", startRow = 1, startCol = 1)
#   saveWorkbook(singleplots)
#  
# # rm(list=ls())
 graphics.off()

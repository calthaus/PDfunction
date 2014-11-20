###################################################################################################################################
##################################growth rates with censoring and pharmacodynamics#################################################
###################################################################################################################################
#set working directory
setwd("C:/Users/sfoerster/ABinteractions/Time_kill/ClassicalTK")
#load libraries
library("drc")
library("XLConnect")
library("metafor")
library("bbmle")
library("drc")
library("plotrix")
# define some functions:
# Negative log-likelihood
nll <- function(a, b, sigma, data = data2regress) {
  m <- a + b * data$x
  m <- ifelse(m == 0, 1e-5, m)
  v <- (sigma/data$y_cens)^2
  ll_i <- (1 - data$cens) * dnorm(data$y_cens, mean = m, sd = sqrt(v), log = TRUE) +
    data$cens * pnorm(data$y_cens, mean = m, sd = sqrt(v), log.p = TRUE)
  ll <- sum(ll_i)
  return(-ll)
}
#function to catch errors resulting from integer(0)
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

#function to remove data when the antibiotic killed cells immediately (count of zero at 0 hours after adding the compound)
remove <- function (x){ 
  for (i in 1:x){
    mydatasub1 = subset(mydata,mydata$sample==problem$sample[1]&mydata$replicate==problem$replicate[1])
    mydatasub2 = subset(mydata,mydata$sample==problem$sample[2]&mydata$replicate==problem$replicate[2])
    mydatasub3 = subset(mydata,mydata$sample==problem$sample[3]&mydata$replicate==problem$replicate[3])
    sub<-rbind(mydatasub1,mydatasub2,mydatasub3)
  }
  return(sub)
}
#initialize some empty vectors and empty variables
parmlist=vector()
append2=vector()
append3=vector()
mydata=vector()
means=vector()
slopes<-vector()
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
#loop for the dataset with the recognized pattern, appends them to a big dataset containing all experiments (mydata)
print("---------------------------------run started---------------------------------")
for (ab in ablist){
  replicatespattern=ab
  replicates <- list.files(path =sourcefiledir, pattern = ab, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
  for (replicate in replicates){
    append1 <- read.table(paste(verzeichnis,"/",replicate,sep=""),header=TRUE)
    append2 <- cbind(append1,replicate,ab,i)
    mydata<-rbind(append2,mydata)
    i=i+1
  }
}
mydata=cbind(mydata,index=rownames(mydata))
#remove values when already at timepoint 0 no colonies couldnt be counted anymore
problem<-subset(mydata,mydata$hours==0&mydata$count==0)
data_to_remove<-remove(length(problem$sample))
cleandata<-as.integer(rownames(data_to_remove))
if(is.integer0(cleandata)==F){
  mydata <- mydata[-c(cleandata),]
}
#make a plot and summary statistics for each experiment
for(j in unique(mydata$replicate)){
  data=mydata[mydata$replicate==j,]
  name=sub(".txt",".pdf",j)
  pdf(paste(name,sep=""))
  hours=data$hours
  time=data$hours
  #set CFU to 100 CFU (detection limit) when zero or 1 is counted
  CFU=data$count*data$dilution*100
  CFU[CFU==0]<-100
  antibiotic=data$antibiotic.conc
  ab.conc<-factor(antibiotic)
  slopes<-vector("numeric")
  #conc is needed to plot the points and for the regressions only for measured points conc2 to plot the lines until 6 hours
  #define concentration in different formats for plotting in legend and as levels for regression
  antibiotic=data$antibiotic.conc
  ab.conc<-factor(antibiotic)
  options(scipen=1000)
  conc=levels(ab.conc)
  conc[conc==min(conc)]= 0
  conc=as.vector(signif(as.numeric(conc),digits=2))
  conc2=mydata$antibiotic.conc
  conc2[conc2==min(conc)]=0
  #x values for linear regression
  tmin<-0
  tmax<-6
  #make regressions to caclulate growth rates including censoring, all data smaller or equal than 1 (Limit of detection) are censored
  for(idose in 1:nlevels(ab.conc)){
    dose<-as.numeric(levels(ab.conc)[idose])
    cens=CFU[hours>=tmin&hours<=tmax & antibiotic==dose]
    cens[cens<=1]<-1
    cens[cens>1]<-0
    y_cens=log(CFU[hours>=tmin & hours<=tmax &antibiotic==dose])
    data2points<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=CFU[hours>=tmin & hours<=tmax &antibiotic==dose],cens=cens)
    data2regress<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=y_cens,cens=cens)
    fit <- mle2(nll, start=list(a=log(1000000), b=-1, sigma = 0.1), method="Nelder-Mead", control=list(maxit=1e3))
    slope <- data.frame(coef(fit)[[2]],dose)
    predict<-coef(fit)[[1]]+data2regress$x*coef(fit)[[2]]
    slopes <- rbind(slope,slopes)
  }
  #fit pharmacodynamic function (Emax with four parameters) to growth rates using drc (based on optim)
  names(slopes)<-c("slope","dose")
  model=drm(slopes$slope~slopes$dose,fct=LL.4())
  #plot the pharmacodynamic model
  options(scipen=100)
  par(mar=c(4.7,5,2,1))
  par(oma=c(0,0,0,0) )
  plot(model,log="x",xlim=c(min(antibiotic),max(antibiotic)*10),ylim=c(-3,2.5),xlab="Antibiotic concentration [mg/L]", ylab=expression("Bacterial growth rate[h"^-1*"]"),pch=21, bg=mypalette, cex=1.4,lwd=1.3,cex.lab=1.5)
  legend("topright",legend=conc, bty="n",pt.bg=mypalette,pch=21,cex=0.85,title="conc.[mg/L]")
  axis.break(axis=1,breakpos=min(antibiotic)*1.5,style="slash")
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
  #make list with parameters for each plot
  Liste=data.frame(k,j,toString(experiment[1]),kappa,kappaerr,upper,uppererr,lower,lowererr,inflection,inflectionerr,zMIC,zMICerr)
  parmlist=rbind(Liste,parmlist)
  dev.off()
  k=k+1
}
names(parmlist)<-c("Number","Identifier","Experiment","kappa","kappaerr","upper","uppererr","lower","lowererr","inflection","inflectionerr","zMIC","zMICerr")
for (l in unique(parmlist$Experiment)){
  parmlistsub <-subset(parmlist,parmlist$Experiment==l)
  kappamodel <- rma(yi = parmlistsub$kappa, sei = parmlistsub$kappaerr, data = parmlistsub, method = 'FE')
  kappasum <- data.frame(cbind(l,"kappa",summary(kappamodel)[1],summary(kappamodel)[2],summary(kappamodel)[3],summary(kappamodel)[4],summary(kappamodel)[5],summary(kappamodel)[6]))
  uppermodel <- rma(yi = parmlistsub$upper, sei = parmlistsub$uppererr, data = parmlistsub, method = 'FE')
  uppersum <- data.frame(cbind(l,"upper",summary(uppermodel)[1],summary(uppermodel)[2],summary(uppermodel)[3],summary(uppermodel)[4],summary(uppermodel)[5],summary(uppermodel)[6]))
  lowermodel <- rma(yi = parmlistsub$lower, sei = parmlistsub$lowererr, data = parmlistsub, method = 'FE')
  lowersum <- data.frame(cbind(l,"lower",summary(lowermodel)[1],summary(lowermodel)[2],summary(lowermodel)[3],summary(lowermodel)[4],summary(lowermodel)[5],summary(lowermodel)[6]))
  inflectionmodel <- rma(yi = parmlistsub$inflection, sei = parmlistsub$inflectionerr, data = parmlistsub, method = 'FE')
  inflectionsum <- data.frame(cbind(l,"inflection",summary(inflectionmodel)[1],summary(inflectionmodel)[2],summary(inflectionmodel)[3],summary(inflectionmodel)[4],summary(inflectionmodel)[5],summary(inflectionmodel)[6]))
  zMICmodel <- rma(yi = parmlistsub$zMIC, sei = parmlistsub$zMICerr, data = parmlistsub, method = 'FE')
  zMICsum <- data.frame(cbind(l,"zMIC",summary(zMICmodel)[1],summary(zMICmodel)[2],summary(zMICmodel)[3],summary(zMICmodel)[4],summary(zMICmodel)[5],summary(zMICmodel)[6]))
  mean=rbind(kappasum,uppersum,lowersum,inflectionsum,zMICsum)
  means=rbind(mean,means)
}
names(means) <- c("experiment","parameter","estimate","se","zval","pval","ci.lb","ci.ub")
#write summary statistics to excel file
options(scipen=-1000)
data=data.frame(parmlist)
allplots <- loadWorkbook("summary_statistics.xlsx", create = TRUE)
createSheet(allplots, name = "parameterlist")
writeWorksheet(allplots, means, sheet = "parameterlist", startRow = 1, startCol = 1)
saveWorkbook(allplots)
#write parameters for each plot to excel file
singleplots <- loadWorkbook("singleplot_statistics.xlsx", create = TRUE)
createSheet(singleplots, name = "parameterlist")
writeWorksheet(singleplots,parmlist, sheet = "parameterlist", startRow = 1, startCol = 1)
saveWorkbook(singleplots)
print("---------------------------------finished---------------------------------")
#rm(list=ls())
graphics.off()
###################################################################################################################################
##################################growth rates with censoring from time-kill data and fitting pharmacodynamic Hill function########
###################################################################################################################################
#set working directory
setwd("C:/Users/sunny/Desktop/PDfunction/SupplementaryMaterial/TableS1/output")
#define patterns to be recognized
sourcefiledir<-"C:/Users/sunny/Desktop/PDfunction/SupplementaryMaterial/TableS1/data"
verzeichnis<-"C:/Users/sunny/Desktop/PDfunction/SupplementaryMaterial/TableS1/data"
#read in additional info about strains and antibiotics of the dataset to analyze
samplenames<-read.table("C:/Users/sunny/Desktop/PDfunction/SupplementaryMaterial/TableS1/data/samplenames_wt.txt",sep="\t",header=T)
ablist<-as.vector(samplenames$sampleid2)
#load libraries
library("drc")
library("XLConnect")
library("metafor")
library("plotrix")
library("magicaxis")
library("censReg")
#initialize some empty vectors and empty variables
parmlist=vector()
append2=vector()
append3=vector()
mydata=vector()
means=vector()
slopes<-vector()
i <- 1
k <- 1
m <- 1
l <- 1
sourcefilepattern<-".txt"
#define some colors
mypalette<- c("red","lightblue","lightblue","lightblue","#C6DBEF", "#9ECAE1" ,"#6BAED6" ,"#4292C6" ,"#2171B5" ,"#08519C", "#08306B","black","black","black","black","black")
# list of all files and pattern to be recognized in these files
files <- list.files(path =sourcefiledir, pattern = sourcefilepattern, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
#loop for the dataset with the recognized pattern, appends them to a big dataset containing all experiments (mydata)
print("---------------------------------run started---------------------------------")
for (ab in ablist){
  replicatespattern=ab
  replicates <- list.files(path =sourcefiledir, pattern = ab, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
  for (replicate in replicates){
    etest<-samplenames[samplenames$sampleid2==ab,]
    append1 <- read.table(paste(verzeichnis,"/",replicate,sep=""),header=TRUE)
    append2 <- cbind(append1,replicate,ab,etest,i,row.names = NULL)
    mydata<-rbind(append2,mydata)
    i=i+1
  }
}
mydata=cbind(mydata,index=rownames(mydata))
for(j in unique(mydata$replicate)){
  data=mydata[mydata$replicate==j,]
  tzero=subset(data,data$hours==0)
  data$count[data$hours==0]<-exp(mean(log(tzero$count)))
  data$dilution[data$hours==0]<-exp(mean(log(tzero$dilution)))
  sample=toString(data$sampleid[1])
  strainame=toString(data$strainname[1])
  letter=toString(data$letter[1])
  experiment=toString(data$replicate[1])
  antibiotics=toString(data$name[1])
  strains= toString(data$strain[1])
  etest=toString(data$etest[1])
  name=sub(".txt",".pdf",j)
  pdf(paste(k,name,sep="_"),width = 9, height = 8)
  hours=data$hours
  time=data$hours
  #set CFU to 100 CFU (detection limit) when zero or 1 is counted
  CFU=data$count*data$dilution*100
  CFU[CFU==0]<-100
  antibiotic=data$antibiotic.conc
  ab.conc<-factor(antibiotic)
  slopes<-vector("numeric")
  #define concentration in different formats for plotting in legend and as levels for regression
  data$antibiotic.conc[data$antibiotic.conc == min(data$antibiotic.conc)] 
  antibiotic=data$antibiotic.conc
  ab.conc<-factor(antibiotic)
  options(scipen=1)
  conc=levels(ab.conc)
  conc=as.vector(signif(as.numeric(conc),digits=2))
  conc[conc==min(conc)]= 0
  #x values for linear regression
  tmin<-0
  tmax<-6
  #make regressions to caclulate growth rates including censoring, all data smaller or equal than 1 (Limit of detection) are censored
  for(idose in 1:nlevels(ab.conc)){
    dose<-as.numeric(levels(ab.conc)[idose])
    cens=CFU[hours>=tmin&hours<=tmax & antibiotic==dose]
    cfumax=CFU[hours=6 & antibiotic==0]
    cens[cens<=100]<-1
    cens[cens>100]<-0
    y_cens=log(CFU[hours>=tmin & hours<=tmax &antibiotic==dose])
    data2points<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=CFU[hours>=tmin & hours<=tmax &antibiotic==dose],cens=cens)
    data2regress<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=y_cens,cens=cens)
    fit<-censReg((data2regress$y_cens) ~ data2regress$x , left = log(100), data = data2regress)
    slope <- data.frame(coef(fit)[[2]],dose)
    predict<-coef(fit)[[1]]+data2regress$x*coef(fit)[[2]]
    slopes <- rbind(slope,slopes)
  }
  # fit pharmacodynamic function (Emax with four parameters) to growth rates using drc (based on optim)
  names(slopes)<-c("slope","dose")
  slopes <- slopes[order(slopes$dose),] 
  model=drm(slopes$slope~slopes$dose,fct=LL.4())
  options(scipen=-1)
  par(mar=c(4.7,5,2,7))
  par(oma=c(0,0,0,0))
  plot(c(0.0000001,1000),c(-15,3),axes=F,bg=mypalette,bty="l",xlab=paste("ETX0914 concentration [mg/L]"), ylab=expression("Bacterial growth rate [h"^-1*"]"),log="x",cex.lab=1.8,col="white")
  plot(model,log="x",xlim=c(0.0000001,1000),axes=F,bty="l",xlab=paste(toString(antibiotics[1]),"concentration [mg/L]"), ylab=expression("Bacterial growth rate [h"^-1*"]"),pch=21,lty=1, bg=mypalette, cex=1.8,lwd=2,cex.lab=1.8,add=TRUE)
  points(slopes$dose,slopes$slope,pch=21, bg=mypalette,cex=1.8,lwd=2)
  box()
  options(scipen=1)
  magaxis(side=c(1),logpretty=TRUE,majorn=10, minorn=10,cex.axis=1.3)
  magaxis(side=c(2),cex.axis=1.3)
  legend("bottomleft",legend=conc, bty="n",pt.bg=mypalette,pch=21,cex=1.2,inset=0.02,title="conc.[mg/L]",xpd=T)
  legend("topleft",legend=k,bty="n",inset=c(-0.1,-0.05), cex=3.6) 
  legend("topright",legend=strainame,bty="n",inset=0, cex=2) 
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
  #make list with parameter statistics for each dataset
  Liste=data.frame(k,j,strains,antibiotics,sample,kappa,kappaerr,upper,uppererr,lower,lowererr,inflection,inflectionerr,zMIC,zMICerr,etest[1],max(conc))
  parmlist=rbind(Liste,parmlist)
  dev.off()
  k=k+1
}
#save the statistics for each file
singleplots <- loadWorkbook("singleplot_statistics_DOGK18.xlsx", create = TRUE)
createSheet(singleplots, name = "parameterlist")
writeWorksheet(singleplots,parmlist, sheet = "parameterlist", startRow = 1, startCol = 1)
saveWorkbook(singleplots)
print("---------------------------------finished---------------------------------")
#rm(list=ls())
graphics.off()
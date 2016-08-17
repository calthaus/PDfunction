###################################################################################################################################
##################################growth rates with censoring and pharmacodynamics#################################################
###################################################################################################################################
#set working directory

setwd("C:/Users/sunny/Desktop/PDfunction/Figure4/Figure4B/output/")
samplenames<-read.table("C:/Users/sunny/Desktop/PDfunction/Figure4/Figure4B/data/samplenames_off-target.txt",sep="\t",header=T)
#define patterns to be recognized
sourcefiledir<-"C:/Users/sunny/Desktop/PDfunction/Figure4/Figure4B/data"
verzeichnis<-"C:/Users/sunny/Desktop/PDfunction/Figure4/Figure4B/data"

sourcefilepattern<-".txt"
ablist<-as.vector(samplenames$sampleid)
foldername<-"off-target"
letter<-""
#load libraries
library("drc")
library("XLConnect")
library("metafor")
library("drc")
library("plotrix")
library("magicaxis")
library("censReg")
#initialize some empty vectors and empty variables
model=NULL
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
#define patterns to be recognized
# sourcefiledir<-paste("C:/Documents and Settings/mollab/Mina dokument/Dropbox/Time_kill/ClassicalTK/input/",foldername,sep="")
# verzeichnis<-paste("C:/Documents and Settings/mollab/Mina dokument/Dropbox/Time_kill/ClassicalTK/input/",foldername,sep="")
# sourcefilepattern<-".txt"
#define some colors
mypalette<- c("red","lightblue","lightblue","lightblue","#C6DBEF", "#9ECAE1" ,"#6BAED6" ,"#4292C6" ,"#2171B5" ,"#08519C", "#08306B","black","black","black","black","black")
# list of all files and pattern to be recognized in these files
files <- list.files(path =sourcefiledir, pattern = sourcefilepattern, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
pdf(paste(foldername,"B.pdf",sep=""),width = 9, height = 8)
options(scipen=100)
par(mar=c(4.7,5,2,7))
par(oma=c(0,0,0,0))
plot(model,log="x",type="n",xlim=c(0.0000001,0.1),ylim=c(-15,3),axes=F,bty="l",xlab=paste("antimicrobial concentration [mg/l]"), ylab=expression("Bacterial growth rate [h"^-1*"]"),pch=21, bg=mypalette, cex=1.8,lwd=5,cex.lab=1.8,col="black")
legend("bottomleft",legend=conc, bty="n",pt.bg=mypalette,pch=21,cex=1.2,inset=0.02,title="conc.[mg/L]",xpd=T)
strainname<-c("azithromycin","cefixime","ceftriaxone","chloramphenicol","ciprofloxacin","gentamycin","penicillin","spectinomycin","tetracycline")
box()
options(scipen=1)
magaxis(side=c(1,2),logpretty=TRUE,cex.axis=1.3)
#loop for the dataset with the recognized pattern, appends them to a big dataset containing all experiments (mydata)
print("---------------------------------run started---------------------------------")
for (ab in ablist){
  replicatespattern=ab
  replicates <- list.files(path =sourcefiledir, pattern = ab, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
  for (replicate in replicates){
    etest<-samplenames[samplenames$sampleid==(ab),]
    append1 <- read.table(paste(verzeichnis,"/",replicate,sep=""),header=TRUE)
    append2 <- cbind(append1,replicate,etest,i,row.names = NULL)
    mydata<-rbind(append2,mydata)
    i=i+1
  }
}
mydata=cbind(mydata,index=rownames(mydata))
print(mydata)
#remove values when already at timepoint 0 no colonies couldnt be counted anymore
# problem<-subset(mydata,mydata$hours==1&mydata$count==0)
# data_to_remove<-remove(length(problem$sample))
# cleandata<-as.integer(rownames(data_to_remove))
# if(is.integer0(cleandata)==F){
#   mydata <- mydata[-c(cleandata),]
# }
#make a plot and summary statistics for each experiment
for(j in unique(mydata$replicate)){
  data=mydata[mydata$replicate==j,]
  tzero=subset(data,data$hours==0)
  data$count[data$hours==0]<-exp(mean(log(tzero$count)))
  data$dilution[data$hours==0]<-exp(mean(log(tzero$dilution)))
  sample=toString(data$sampleid[1])
  experiment=toString(data$replicate[1])
  antibiotics=toString(data$name[1])
  strains= toString(data$strain[1])
  etest=toString(data$etest[1])
  color=toString(data$straincol[1])
  linet=data$linetype[1]
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
  options(scipen=1000)
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
  #fit pharmacodynamic function (Emax with four parameters) to growth rates using drc (based on optim)
  names(slopes)<-c("slope","dose")
  model=drm(slopes$slope~slopes$dose,fct=LL.4())

  #start parameters
  weights=slopes$slope[11]/2+slopes$slope[10]/4+slopes$slope[9]/8+slopes$slope[8]/16
  K=weights/2
  larger=subset(slopes,slopes$slope>K)
  smaller=subset(slopes,slopes$slope<K)
  cn=max(larger$dose)
  ci=min(smaller$dose)
  ck=(ci+cn)/2
  biflist=NULL
  for (i in 1:length(smaller$dose)){
    ybif=smaller$slope[i]-(smaller$slope[i+1])
    #xbif=(smaller$dose[i]-smaller$dose[i+1])/(log(smaller$dose[i])-log(smaller$dose[i+1]))
    xbif=(smaller$dose[i]+smaller$dose[i+1])/2
    
    bif=cbind(xbif,ybif)
    biflist=rbind(bif,biflist)
    print(ybif)
    }
  distance=data.frame(biflist)
  wendepunkt=subset(distance,distance$ybif==max(distance$ybif,na.rm=TRUE))
  firstdf=subset(slopes,slopes$dose<=wendepunkt$xbif)
  model=drm(firstdf$slope~firstdf$dose,fct=LL.4())
  secondf=subset(slopes,slopes$dose>=0.0002)
  df=slopes$dose
  #start parameters
  einphasic<-function(b,c,d,e,x) {
    #c+((d-c)/1+exp(b*(log(x)-log(e))))
    1+((c-1)/(1+(e/x)^b))
  } 
  uniphasic<-function(psimax1,psimin1,kappa1,MIC1,conc) {
    (psimax1-(psimax1-psimin1)*(conc/MIC1)^kappa1/((conc/MIC1)^kappa1-psimin1/psimax1))
  }
  biphasic<-function(psimax1,psimin1,kappa1,MIC1,psimax2,psimin2,kappa2,MIC2,conc) {
    (psimax1-(psimax1-psimin1)*(conc/MIC1)^kappa1/((conc/MIC1)^kappa1-psimin1/psimax1))*(psimax2-(psimax2-psimin2)*(conc/MIC2)^kappa2/((conc/MIC2)^kappa2-psimin2/psimax2))
  }
  x=rev(firstdf$dose)
  curve(einphasic( 1.95234614,-0.67667978, 0.90417939, 0.00008599,x),2e-7,80000,add=TRUE,col="red",lty=1,lwd=2)
  x=rev(secondf$dose)
  curve(uniphasic(-0.67667978,-2.33377, 2.50278, 0.00593,rev(x)),2e-7,80000,add=TRUE,col="green",lty=1,lwd=2)
  x=df
  curve(biphasic(0.90417939,-0.67667978, 1.95234614, 0.00008599,-0.67667978,-2.33377, 2.50278, 0.00593,x),2e-7,1,add=TRUE,col="black",lty=1,lwd=2)
  plot(model,log="x",xlim=c(0.0000001,100000),ylim=c(-15,3),col=color,lty=linet,axes=F,bty="l",main="",xlab=paste(antibiotics,"concentration [mg/L]"), ylab=expression("Bacterial growth rate [h"^-1*"]"),pch=21, bg=mypalette, cex=1.8,lwd=2,cex.lab=1.8,add=T)
  legend("topleft",legend="(B)",bty="n",inset=c(-0.1,-0.05), cex=3.6) 
  legend("topright",legend="strain: WT 103",bty="n",inset=0, cex=2) 
  
  
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
  Liste=data.frame(k,j,strains,antibiotics,toString(experiment[1]),kappa,kappaerr,upper,uppererr,lower,lowererr,inflection,inflectionerr,zMIC,zMICerr,etest[1])
  parmlist=rbind(Liste,parmlist)
  k=k+1
}
singleplots <- loadWorkbook("singleplot_statistics_overlay.xlsx", create = TRUE)
createSheet(singleplots, name = "parameterlist")
writeWorksheet(singleplots,parmlist, sheet = "parameterlist", startRow = 1, startCol = 1)
saveWorkbook(singleplots)
print("---------------------------------finished---------------------------------")
#rm(list=ls())
graphics.off()
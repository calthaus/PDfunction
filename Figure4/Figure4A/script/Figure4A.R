###################################################################################################################################
##################################fit linear regressions and plot lines############################################################
###################################################################################################################################
# set working directory
setwd("C:/Users/sunny/Desktop/PDfunction/Figure4/Figure4A/output")
samplenames<-read.table("C:/Users/sunny/Desktop/PDfunction/Figure4/Figure4A/data/samplenames_off-target.txt",sep="\t",header=T)
#define patterns to be recognized
sourcefiledir<-"C:/Users/sunny/Desktop/PDfunction/Figure4/Figure4A/data"
verzeichnis<-"C:/Users/sunny/Desktop/PDfunction/Figure4/Figure4A/data"
sourcefilepattern<-".txt"
ablist<-as.vector(samplenames$sampleid)
options(scipen=1)
# load libraries
library("Hmisc")
library("drc")
library("XLConnect")
library("drc")
library("magicaxis")
library("censReg")
#initialize some empty vectors and empty variables
savefont <- par(font=2) 
sub=vector() 
parmlist=vector()
append2=vector()
append3=vector()
mydata=vector()
means=vector()
i <- 1
j <- 1
k <- 1
b=-4
#define some colors
mypalette<- c("red","lightblue","lightblue","lightblue","#C6DBEF", "#9ECAE1" ,"#6BAED6" ,"#4292C6" ,"#2171B5" ,"#08519C", "#08306B","black","black","black","black","black")
# list of all files and pattern to be recognized in these files
files <- list.files(path =sourcefiledir, pattern = sourcefilepattern, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
#loop for the dataset with the recognized pattern, appends them to a big dataset containing all experiments (mydata)
print("---------------------------------run started---------------------------------")
for (ab in ablist){
  replicatespattern=ab
  replicatelist<- list.files(path =sourcefiledir, pattern = ab, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
  replicates<-Filter(length,replicatelist)
  #read data from file with the recognized pattern and append to a dataframe
  for (replicate in replicates){
    etest<-samplenames[samplenames$sampleid==ab,]
    etest$abname <- toString(etest$abname)
    append1 <- read.table(paste(verzeichnis,"/",replicate,sep=""),header=TRUE)
    append2 <- cbind(append1,replicate,etest,i,row.names = NULL)
    mydata<-rbind(append2,mydata)
    i=i+1
  }
}
mydata=cbind(mydata,index=rownames(mydata))

#plot time against CFU and make the summary statistics for all experiments done for a strain treated with an antibiotic (j)
for(j in unique(mydata$replicate)){
  data=mydata[mydata$replicate==j,]
  tzero=subset(data,data$hours==0)
  data$count[data$hours==0]<-exp(mean(log(tzero$count)))
  data$dilution[data$hours==0]<-exp(mean(log(tzero$dilution)))
  sample=toString(data$sample[1])
  strainame=toString(data$strainname[1])
  letter=toString(data$letter[1])
  experiment=toString(data$replicate[1])
  antibiotics=data$abname[1]
  strains= toString(data$strain[1])
  etest=toString(data$etest[1])
  pdf("Figure4A.pdf",width = 9, height = 8 )
  hours=data$hours 
  time=data$hours
  hours2=data$hours
  #Calculate the CFU from counts, set hundred as lower limit of detection for all zero counts
  CFU=data$count*data$dilution*100
  CFU2=data$count*data$dilution*100
  CFU[CFU==0]<-100
  #define concentration in different formats for plotting in legend and as levels for regression
  antibiotic=data$antibiotic.conc
  ab.conc<-factor(antibiotic)
  #options(scipen=1000)
  conc=levels(ab.conc)
  conc=as.vector(signif(as.numeric(conc),digits=2))
  conc[conc==min(conc)]= 0
  conc2=mydata$antibiotic.conc
  conc2[conc2==min(conc)]=0
  #time range for regression
  tmin<-0
  tmax<-6
  slopes<-vector("numeric")
  options(scipen=1)
  #plot empty frame and legend so lines and points for each concentration can be added later on
  par(mar=c(4.7,5,2,7))
  par(oma=c(0,0,0,0))
  plot(c(-4,6),c(10,2e9),type="n",xlab="Time [h]",axes=F,log="y",bty="l",main="", ylab="Bacteria [CFU/ml]",cex.lab=1.8,cex=1.8,lty=2,lwd=2)
  box()
  magaxis(side=c(1,2),logpretty=TRUE,cex.axis=1.3)
  par(savefont)
  legend(("bottomleft" ),legend=conc,lty=1,inset=0.02,col=mypalette,lwd=2,pch=1,cex=1.2,bty="n",title="conc.[mg/L]",xpd=T)
  legend(-5.5,9000000000,legend=letter,bty="n",inset=0, cex=3.6) 
  for(idose in 1:nlevels(ab.conc)){
    dose<-as.numeric(levels(ab.conc)[idose])
    cens=CFU[hours>=tmin&hours<=tmax & antibiotic==dose]
    #data for plotting the lines contains all data from -2 to 6 hours, data below detection are plottet as CFU/ml=100
    y_lines=CFU[hours>=-4&hours<=6 & antibiotic==dose]
    x_lines=hours[hours>=-4&hours<=6 & antibiotic==dose]
    #data for plotting the points contains only REAL measured data, data below detection limit are assigned NA 
    CFU2[CFU2<=1]<-NA
    hours2[CFU2<=1]<-NA
    x_points=hours2[hours2>=-4&hours2<=6 & antibiotic==dose]
    y_points=CFU2[hours2>=-4&hours2<=6 & antibiotic==dose]
    doses=antibiotic[hours>=-4&hours<=6 & antibiotic==dose]
    #calculate censoring and log likelihood
    cens[cens<=100]<-1
    cens[cens>100]<-0
    data2points<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=CFU[hours>=tmin & hours<=tmax &antibiotic==dose],cens=cens)
    data2regress<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=log(CFU[hours>=tmin & hours<=tmax &antibiotic==dose]),cens=cens)
    fit<- censReg((data2regress$y_cens) ~ data2regress$x , left = log(100), data = data2regress)
    predict<-coef(fit)[[1]]+data2regress$x*coef(fit)[[2]]
    curve(exp(coef(fit)[1])*exp(coef(fit)[2]*x),0,6,add=TRUE,col=mypalette[idose],lty=2)
    lines(x_lines,y_lines,col=mypalette[idose],lty=1,cex=1,lwd=2)
    points(x_points,y_points,col=mypalette[idose],lty=1, pch=1, cex=1.5)
    slope <- data.frame(coef(fit)[[2]],dose)
    slopes <- rbind(slope,slopes)
  }
  k=k+1
  dev.off()
}
#clean up workspace
print("---------------------------------run finished--------------------------------")
rm(list=ls())
graphics.off()

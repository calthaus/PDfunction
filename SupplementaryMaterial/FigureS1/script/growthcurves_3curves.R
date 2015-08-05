########################################################################################################
########################## Fitting the Gompertz Model to growth data####################################
########################################################################################################
#load libraries
library("cellGrowth")
library("XLConnect")
library("drc")
detach("package:drc", unload=TRUE)
library("censReg")
library("magicaxis")
#read in the data
setwd("C:/Users/sunny/Desktop/PDfunction/SupplementaryMaterial/FigureS1/output")
curves=read.table("C:/Users/sunny/Desktop/PDfunction/SupplementaryMaterial/FigureS1/data/Growth_GW_short.txt",header=T)
#calculate colony forming units from counts and make new data frame
CFU <- curves$Count*curves$Dilution*100
datatemp <- data.frame(curves$Strain,curves$Count,curves$Dilution,curves$Time,CFU,curves$Experiment)
data <- data.frame(datatemp)
names(data) <- c("strain","count","dilution","time","CFU","exp")
strains <- sort(levels(data$strain))
experiment<-levels(data$experiment)
regdata <- subset(data,data$time <= 20 & data$time >= 2)
#Symbols and colors
colList<-c("grey39","grey","grey10")
pchList <-c(21,24,23,NA)
ltyList <-c(2,2,2,1)
#initialize vectors
parameters <- c()
#fit model and calculate growth rate

  
for (i in strains){
  pdf(paste(i,".pdf",sep=""))
  y=data$CFU[data$strain==i]
  x=data$time[data$strain==i]
  y=log(y)
  #fit linear regression from t4-t12
  x1=regdata$time[regdata$strain==i]
  y1=log(regdata$CFU[regdata$strain==i])
  #regression <- lm(y1~x1)
  
  
  regression<- censReg(y1 ~ x1 , left = log(100), data = regdata)
  
  #x_reg <- data.frame(x1=seq(from=-10,to=50,length.out=100))
  #y_reg <- predict(regression,x_reg,se.fit = TRUE)
  
  
  
  #guess sensible start parameters for Gompertz model
  guess <- guessCellGrowthParams(x,y,relative.height.at.lag=0.2)
  #pars<-c(mu=guess$mu,l=guess$l,z0=guess$z0,zmax=guess$zmax)
  fit <- nls(y~gompertz(x,mu,l,z0,zmax),start=guess)
  #fit <- optim(pars,gompertz(x,mu=guess$mu,l=guess$l,z0=guess$z0,zmax=guess$zmax),method="Nelder-Mead")
  xnew.pooled=list(x=seq(min(x),max(x),length.out=1000))
  ypredict.pooled=exp(predict(fit,xnew.pooled))
  options("scipen"=1)
  #plot legend and frame to add lines and points after

  par(mar=c(4.7,5,2,1))
  par(oma=c(0,0,0,0) )
  plot(data$time[data$strain==i],data$CFU[data$strain==i],axes=FALSE, ylim=c(100,100000000000),xlim=c(0,36), xlab = "", ylab = "",type="n",log="y",cex=3,cex.lab=1.5)
  mtext(side = 1, text = "time[h]", cex=1.8,line = 3)
  mtext(side = 2, text = expression("CFU [ml"^-1*"]"),cex=1.8, line = 3)
  magaxis(side=c(1,2),logpretty=TRUE,cex.axis=1.5)
#   axis(1, c(0,4,8,12,16,20,24,28,32,36,40),cex.axis=1.5)
#   options("scipen"=-100, "digits"=4)
#   axis(2, c(1000,100000,10000000,1000000000),cex.axis=1.5)
  box(lty = 1)
  #plot points in different symbols for each experiment
  for(j in unique(data$exp)){
    newdata <- subset(data, exp == j& strain==i,select=c(CFU,time,strain,exp))
    y3=log(newdata$CFU)
    x3=newdata$time

#     #guess sensible start parameters for Gompertz model
    guess.exp <- guessCellGrowthParams(as.vector(x3),as.vector(y3),relative.height.at.lag=0.1)
    fit.exp <- nls(y3~gompertz(x3,mu,l,z0,zmax),start=guess.exp)
    xnew.exp=list(x3=seq(min(x3),max(x3),length.out=1000))
    ypredict.exp=exp(predict(fit.exp,xnew.exp))
    lines(xnew.exp$x3,ypredict.exp,lty=2)
    
    points(newdata$time,newdata$CFU,pch=pchList[j],bg=colList[j],cex=1.7)
  }
  # add regression line and gompertz fit
  #lines(x_reg$x1,exp(y_reg$fit),lty=2,lwd=2)
  lines(xnew.pooled$x,ypredict.pooled)
  # calculate parameters
  growth=coef(summary(regression))[2]
  mu=round(coef(fit)[[1]],2)
  zmax=exp(round(coef(fit)[[4]],0))
  zmax=formatC(zmax,digits=2,format="E")
  parameter=data.frame(i,mu,growth,zmax)
  parameters=rbind(parameter,parameters)
  # plot legend
  legend("bottomright",c("experiment 1","experiment 2","experiment 3","pooled data"),pch=pchList,lty=ltyList,cex=1.3, pt.bg=colList, col="black",bty="n")
  legend("topleft",i,cex=3,bty="n")
  
  dev.off()
}
names(parameters)<-c("strain","max. growth rate","growth rate","asymptote")
allplots <- loadWorkbook("growth_parameter.xlsx", create = TRUE)
createSheet(allplots, name = "parameterlist")
writeWorksheet(allplots, parameters, sheet = "parameterlist", startRow = 1, startCol = 1)
saveWorkbook(allplots)
#clean workspace
#rm(list=ls())
graphics.off()
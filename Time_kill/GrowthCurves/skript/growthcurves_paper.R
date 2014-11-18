########################################################################################################
########################## Fitting the Gompertz Model to growth data####################################
########################################################################################################

#options


#load packages
library("cellGrowth", lib.loc="C:/Program Files/R/R-3.1.1/library")
library("XLConnect")

#read in the data
setwd("C:/users/sfoerster/ABinteractions/Time_kill/GrowthCurves")
curves=read.table("input/Growth_GW_short.txt",header=T)

#calculate colony forming units from counts and make new data frame
CFU <- curves$Count*curves$Dilution*100
datatemp <- data.frame(curves$Strain,curves$Count,curves$Dilution,curves$Time,CFU,curves$Experiment)
data <- data.frame(datatemp)
names(data) <- c("strain","count","dilution","time","CFU","exp")
strains <- sort(levels(data$strain))
regdata <- subset(data,data$time <= 20 & data$time >= 4)

#Symbols and colors 
colList<-c("grey39","grey","grey10")
pchList <-c(21,24,23,NA)
ltyList <-c(1,1,1,2)

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
  regression <- lm(y1~x1)
  x_reg <- data.frame(x1=seq(from=-10,to=50,length.out=100))
  y_reg <- predict(regression,x_reg,se.fit = TRUE)
  
  #guess sensible start parameters for Gompertz model
  guess <- guessCellGrowthParams(x,y,relative.height.at.lag=0.2)
  pars<-c(mu=guess$mu,l=guess$l,z0=guess$z0,zmax=guess$zmax)
  fit <- nls(y~gompertz(x,mu,l,z0,zmax),start=guess)
  print(fit)
  xnew=list(x=seq(min(x),max(x),length.out=1000))
  ypredict=exp(predict(fit,xnew))
  options("scipen"=100, "digits"=4)
  plot(data$time[data$strain==i],data$CFU[data$strain==i],axes=FALSE, xlab = "", ylab = "",type="n",log="y",cex=1.5)
  mtext(side = 1, text = "time[h]", cex=1.3,line = 2.5)
  mtext(side = 2, text = expression("CFU [ml"^-1*"]"),cex=1.3, line = 2.5)
  axis(1, c(0,4,8,12,16,20,24,28,32,36,40))
  options("scipen"=-100, "digits"=4)
  axis(2, c(1000,100000,10000000,1000000000))
  box(lty = 1)

 #plot points in different symbols for each experiment
  for(j in unique(data$exp)){
    newdata <- subset(data, exp == j& strain==i,select=c(CFU,time,strain,exp)) 
    points(newdata$time,newdata$CFU,pch=pchList[j],bg=colList[j],cex=1.7) 
  }
  
  # add regression line and gompertz fit
  lines(x_reg$x1,exp(y_reg$fit),lty=2,lwd=2)
  lines(xnew$x,ypredict)
  
  # calculate parameters
  growth=coef(summary(regression))[2]
  mu=round(coef(fit)[[1]],2)
  zmax=exp(round(coef(fit)[[4]],0))
  zmax=formatC(zmax,digits=2,format="E")
  parameter=data.frame(i,mu,growth,zmax)
  parameters=rbind(parameter,parameters)
  # plot legend
  legend("topleft",c("experiment 1","experiment 2","experiment 3","growth rate"),pch=pchList,lty=ltyList, pt.bg=colList, col="black",bty="n")
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


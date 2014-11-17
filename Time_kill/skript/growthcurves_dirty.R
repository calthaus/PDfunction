### Fitting Gompertz growth curves###
#read in data
setwd("C:/users/sfoerster/ABinteractions/Time_kill")
curves=read.table("input/Growth_GW_short.txt",header=T)
#calculate colony forming units from counts
CFU <- curves$Count*curves$Dilution*100
datatemp <- data.frame(curves$Strain,curves$Count,curves$Dilution,curves$Time,CFU,curves$Experiment)
data <- data.frame(datatemp)

names(data) <- c("strain","count","dilution","time","CFU","exp")
strains <- sort(levels(data$strain))
#data <- subset(data, exp == i,select=c(CFU, time,strain)) 

for (i in strains){
  pdf(paste(i,"short.pdf",sep=""))
  y=data$CFU[data$strain==i]
  print("-------------------------")
  print(y)
  x=data$time[data$strain==i]
  print(x)
  y=log(y)

  guess <- guessCellGrowthParams(x,y,relative.height.at.lag=0.2)
  pars<-c(mu=guess$mu,l=guess$l,z0=guess$z0,zmax=guess$zmax)
  #gomfit<-gompertz(x=x,mu=guess$mu,l=guess$l,z0=guess$z0,zmax=guess$zmax)
  #fit2 <- optim(pars,ssr,gr=NULL,dat,hessian=TRUE)
  fit <- nls(y~gompertz(x,mu,l,z0,zmax),start=guess)
  print(fit)
  xnew=list(x=seq(min(x),max(x),length.out=1000))
  ypredict=exp(predict(fit,xnew))
  plot(data$time[data$strain==i],data$CFU[data$strain==i],main=paste(i,"Gompertz growth model"),xlab="time[h]",ylab="bacteria[CFU/ml]",type="n",log="y")
  #points(x,y,col=c("black","red","green","blue"))
  for(j in unique(data$exp)){
    newdata <- subset(data, exp == j& strain==i,select=c(CFU,time,strain,exp)) 
    points(newdata$time,newdata$CFU,col=c(j+1),cex=1.5) 
  }
  lines(xnew$x,ypredict)
  mu=round(coef(fit)[[1]],2)
  zmax=exp(round(coef(fit)[[4]],0))
  zmax=formatC(zmax,digits=2,format="E")
  legend("topleft",c(paste("growth rate [per hour]:",mu),paste("asymptote [CFU/ml]:",zmax)))
  #Plotting the fitted function to the data

  dev.off()
}




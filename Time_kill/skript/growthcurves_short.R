### Fitting logistic growth curve to data from time-kill assay


j=1
### Example: Sample No. 3 (Pre FB minus)
setwd("C:/users/sfoerster/ABinteractions/Time_kill")
curves=read.table("input/Growth_GW_short.txt",header=T)
vector=c()
counts=c()
CFUs=c()
growthcurves=c()
CFU <- curves$Count*curves$Dilution*100
datatemp <- data.frame(curves$Strain,curves$Count,curves$Dilution,curves$Time,CFU,j)
#counts<-rbind(count,counts)
data <- data.frame(datatemp)

names(data) <- c("strain","count","dilution","time","CFU","exp")
print(data)
#meandata=aggregate(CFU~strain+time, data, mean)

for (j in unique(data$strain)){
  time=data$time[data$strain==j]
  foo=c(0 , 2 , 4 , 6 , 8 ,10 ,12 ,20, 22 ,24, 26 ,28 ,30 ,32 ,34, 36)
  print("foo-------------------------------------------------")
  print(time)
  CFU=data$CFU[data$strain==j]
  CFU[CFU==0]<-1
  pdf(paste(j,"short.pdf",sep=""))
  x=0:36
  y=data$CFU[data$strain==j]
  print(y)
  y=log(y)
  #plot(x,y)
  guess <- guessCellGrowthParams(x,y,relative.height.at.lag=0.5)
  print(guess)
  fit <- nls(y~gompertz(x,mu,l,z0,zmax),start=guess)
  print(fit) 
  x <- seq(0,36,length=100)
  ypredict=predict(fit,newdata=list(x))

  #plot(data$time[data$strain==j],log(data$CFU[data$strain==j]),col=c(1,2,3),main=j,xlab="time[h]",ylab="bacteria[CFU/ml]")
  plot(x,ypredict,lwd=2,col="red",type="l")  
  points(x,ypredict)
  points(data$time[data$strain==j],log(data$CFU[data$strain==j]),col=c(1,2,3),main=j,xlab="time[h]",ylab="bacteria[CFU/ml]")
  
  #m1 <- drm(log(data$CFU[data$strain==j])~data$time[data$strain==j], logDose=NULL, fct=gompGrowth.1())
  #m1 <- smooth.spline(data$time[data$strain==j],data$CFU[data$strain==j])
  #foo<-predict(m1,data$time[data$strain==j])
  #plot(m1,log="y")
  # Define your growth model (logistic or any other)

  # Plotting the fitted function to the data
  #x <- seq(min(data$time),max(data$time),length=1000)
  #lines(x,model(fit$par,x),col="blue")
  dev.off()
}




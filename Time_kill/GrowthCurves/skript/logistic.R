### Fitting logistic growth curve to data from time-kill assay

samplename<-function(j){
  switch(j,
         "1" = {
           name<- "AZ"
         },                 
         "2" = {
           name<- "F89"
         },
         "3" = {
           name<- "QC18"
         },
         "4" = {
           name<- "Pre_FB_minus"
         },
         "5" = {
            name<- "Pre_FB_plus"
         },
         "6" = {
           name<- "Pre_GCBL"
         },
  )
  return(Name=name)
}

### Example: Sample No. 3 (Pre FB minus)
setwd("C:/users/sfoerster/ABinteractions/Time_kill")
curves=read.table("input/Growth_GW.txt",header=T)
vector=c()
for (i in unique(curves$sample_no)){
	count1 <- curves$count1[curves$sample_no==i]
	count2 <- curves$count2[curves$sample_no==i]
	count3 <- curves$count3[curves$sample_no==i]
	dilution <- curves$dilution[curves$sample_no==i]
  expresp <- ((count1+count2+count3)/3)*dilution
	expdose <- curves$time[curves$sample_no==i]
	expdata <- cbind(as.data.frame(expdose), as.data.frame(expresp), as.data.frame(i))
	colnames(expdata) = c("time",paste("CFU.",i,sep=""),"i")
	vector=c(expdata,vector)
}
vector=as.data.frame(vector)
loop=c(1,2,3)
for(j in loop){
	run=paste("CFU",j,sep=".")
	print(run)
	titelname=samplename(j)
	png(paste("test-growthcurves-logistic-fit-FB",run,".png",sep=""))
	# Define the data and put them into a data frame

	time <- vector$time
	count <- vector[j+(j*2-1)]

	data <- as.data.frame(cbind(time,count))
	colnames(data)=c("time","count")
	print(data)
	# Plot the data
	plot(data$time,data$count,log="y",ylim=c(min(count),max(count)*10),main=titelname,xlab="Time (hours)",ylab="CFU (per ml)",frame=FALSE)
	
	# Define your growth model (logistic or any other)
	model <- function(pars,t) {
		with(as.list(pars), N*K/(N+(K-N)*exp(-r*t)))
	}

	# The function to be minimized over (sum of squared residuals of the log-transformed values)
	ssr <- function(pars,data) {
		predicted <- model(pars,data$time)
		sum((log(predicted)-log(data$count))^2)
	}

	# Initial values of your parameters
	pars <- c(r = 0.5, K = max(data$count), N = min(data$count))

	# Fitting the model to the data
	fit <- optim(pars,ssr,gr=NULL,data)

  
	#fit<-nls(log(count)~log(N*K/(N+(K-N)*exp(-r*time))),data,start=list(r=0.5,N=min(data$count),K=max(data$count)))
	print(fit)
	# Plotting the fitted function to the data
	x <- seq(min(data$time),max(data$time),length=100)
	lines(x,model(fit$par,x),col="blue")

	# Add the linear regression from 2 to 8 hours
	dat <- data[2:5,]
	lines(dat$time,exp(predict(lm(log(as.matrix(count)) ~ time, dat))),col="red")
	modelreg=lm(log(count)~time,dat)
	legend("topright",legend=c(paste("slope_red=",round(modelreg$coefficients[2],3)),paste("slope_blue=",round(fit$par[1],3))))
}


graphics.off()

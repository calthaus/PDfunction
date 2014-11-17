###################################################################################################################################
##################################fit pharmacodynamic function and impute missing data#############################################
###################################################################################################################################

#set working directory  
setwd("C:/Users/sfoerster/ABinteractions/Time_kill/GrowthRate")
options(scipen=1)
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
library("bbmle")
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
#titellist <- c("WHOA Azithromycin","WHOA Ciprofloxacin","WHOA Tetracyclin","WHOA Ceftriaxone","F89 Azithromycin","F89 Ciprofloxacin","F89 Tetracyclin","F89 Ceftriaxone")

#do action for the dataset with the recognized pattern
for (ab in ablist){
  replicatespattern=ab
  replicates <- list.files(path =sourcefiledir, pattern = ab, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
  print(replicates)
  #read data from file with the recognized pattern and append to a dataframe
  for (replicate in replicates){
    append1 <- read.table(paste(verzeichnis,"/",replicate,sep=""),header=TRUE)
    append2 <- cbind(append1,replicate,ab,i)
    mydata<-rbind(append2,mydata) 
    i=i+1
  }
}

#remove values when already at timepoint 0 no colonies couldnt be counted anymore
# mydata1 <- subset(mydata,mydata$hours==0)
# mydata1 <- subset(mydata1,mydata1$count>0)
# mydata2 <- subset(mydata,hours>0)
# mydata <- rbind(mydata1,mydata2)
#mean
# mydata1<-subset(mydata,mydata$hours==0)
# mean0 <- mean(mydata1$count)
# mydata1$count<-mean0
# mydata1 <- subset(mydata1,mydata1$count>0)
# mydata2 <- subset(mydata,hours>0)
# mydata <- rbind(mydata1,mydata2)
problem<-subset(mydata,mydata$hours==0&mydata$count==0)
rows_to_garbage1 = which(mydata$sample==toString(problem$sample[1]), useNames = TRUE)
mydata <- mydata[-c(rows_to_garbage1),]
rows_to_garbage2 = which(mydata$sample==toString(problem$sample[2]), useNames = TRUE)
mydata <- mydata[-c(rows_to_garbage2),]

print(rows_to_garbage1)
print(rows_to_garbage2)

for(j in unique(mydata$replicate)){
  print("----------------------------filename------------------------------------------")
  data=mydata[mydata$replicate==j,]
  name=sub(".txt","regression.pdf",j)
  pdf(paste(name,sep=""))
  hours=data$hours
  time=data$hours
  hours2=data$hours
  CFU=data$count*data$dilution*100
  CFU2=data$count*data$dilution*100
  #CFU[CFU==0&time==0]<- NA
  #cens=CFU
  #cens[cens==0]<-1
  #cens[cens>1]<-0
  CFU[CFU==0]<-100
  antibiotic=round(data$MIC,3)
  ab.conc<-factor(antibiotic)
  conc=levels(ab.conc)
  conc[conc==min(conc)]= 0
  
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
  #par(mar=c(5,5,1,1))
  #plot(c(-2,6),c(10,1e9),type="n",xlab="Time [h]",main=j,log="y",ylab="Bacteria [CFU/ml]",cex.lab=1.3,cex=1.4,lty=2,lwd=8)
  library("drc")
  plot(c(-2,6),c(10,1e10),type="n",xlab="Time [h]",log="y",ylab="Bacteria [CFU/ml]",cex.lab=1.3,cex=1.4,lty=2,lwd=8)
  legend(x=3,y=100,"Limit of detection: 100 CFU/ml", bty="n",text.col="gray47",cex=0.8)
  legend("topleft",conc,col=mypalette,bty="n",text.col="gray47",pch=1,cex=0.8,title="conc.[mg/L]")
  axis.break(axis=1,breakpos=-1,style="slash") 
  abline(h =100, untf = FALSE,col="gray47",lty=2)
  for(idose in 1:nlevels(ab.conc)){
    dose<-as.numeric(levels(ab.conc)[idose])
    #data2regress=data.frame(cbind(x=time, cens=cens, y_cens=CFU,ab=dose)) 
    cens=CFU[hours>=tmin&hours<=tmax & antibiotic==dose]
    y_lines=CFU[hours>=-2&hours<=6 & antibiotic==dose]
    x_lines=hours[hours>=-2&hours<=6 & antibiotic==dose]    
    CFU2[CFU2==0]<-NA
    hours2[CFU2==0]<-NA
    
    x_points=hours2[hours2>=-2&hours2<=6 & antibiotic==dose]
    y_points=CFU2[hours2>=-2&hours2<=6 & antibiotic==dose]  
    doses=antibiotic[hours>=-2&hours<=6 & antibiotic==dose]
    
    cens[cens==0]<-1
    cens[cens>1]<-0
    y_cens=log(CFU[hours>=tmin & hours<=tmax &antibiotic==dose])
    #data2regress<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=log(CFU[hours>=tmin&hours<=tmax & antibiotic==dose]/CFU[hours==tmin&antibiotic==dose]),cens=cens)
    data2points<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=CFU[hours>=tmin & hours<=tmax &antibiotic==dose],cens=cens)
    data2regress<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=log(CFU[hours>=tmin & hours<=tmax &antibiotic==dose]),cens=cens)
    #data2lines<-data.frame(x_lines=hours[hours>=-2 & hours<=6 &antibiotic==dose],y_lines=CFU[hours>=-2 & hours<=6 &antibiotic==dose],cens=cens)
    
    #names(data2regress)<-c("x")
    print(data2regress)
    #slopes[idose]<-coef(lm(log.bacteria.rel~time-1,data2regress))
    #fit <- mle2(nll, start=list(a =log(1000000), b=-5, sigma = 1), method="Nelder-Mead", control=list(maxit=1e3))
    fit <- mle2(nll, start=list(a =log(1000000), b=-1,sigma=0.1), method="Nelder-Mead", control=list(maxit=1e3))
    
    print("---------------------------------------------------------------------------------------------------")
    print(fit)
    b <- coef(fit)[[2]]
    fita <- c(coef(fit)[[1]],coef(fit)[[2]])
    predict<-coef(fit)[[1]]+data2regress$x*coef(fit)[[2]]  
    #points(data2points$x,data2points$y_cens,col=mypalette[idose],pch=1,cex=1.5)
    #lines(data2lines$x,data2lines$y,col=mypalette[idose],pch=1,cex=1.5,lty=2)
    lines(x_lines,y_lines,col=mypalette[idose],lty=1) 
    points(x_points,y_points,col=mypalette[idose],lty=1, pch=1, cex=1.5)
    slopes <- c(fita,slopes)
  }
  dev.off()
  k=k+1
}


#rm(list=ls())
graphics.off()
    
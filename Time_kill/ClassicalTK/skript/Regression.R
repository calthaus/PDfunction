###################################################################################################################################
##################################fit linear regressions and plot lines############################################################
###################################################################################################################################

#set working directory  
setwd("C:/Users/sfoerster/ABinteractions/Time_kill/ClassicalTK")
options(scipen=1)

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

#load libraries
library("drc")
library("XLConnect")
library("bbmle")
library("drc")

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
problem<-subset(mydata,mydata$hours==0&mydata$count==0)
rows_to_garbage1 = which(mydata$sample==toString(problem$sample[1]), useNames = TRUE)
mydata <- mydata[-c(rows_to_garbage1),]
rows_to_garbage2 = which(mydata$sample==toString(problem$sample[2]), useNames = TRUE)
mydata <- mydata[-c(rows_to_garbage2),]


for(j in unique(mydata$replicate)){
  data=mydata[mydata$replicate==j,]
  name=sub(".txt","regression.pdf",j)
  pdf(paste(name,sep=""))
  hours=data$hours
  time=data$hours
  hours2=data$hours
  CFU=data$count*data$dilution*100
  CFU2=data$count*data$dilution*100
  CFU[CFU==0]<-100
  antibiotic=data$antibiotic.conc
  ab.conc<-factor(antibiotic)
  options(scipen=-1000)
  conc=levels(ab.conc)
  conc[conc==min(conc)]= 0
  conc=as.vector(round((as.numeric(conc)),3))
  tmin<-0 
  tmax<-6 
  # define a vector for the decline rates:
  slopes<-vector("numeric",length=nlevels(ab.conc))


  plot(c(-2,6),c(10,1e10),type="n",xlab="Time [h]",log="y",ylab="Bacteria [CFU/ml]",cex.lab=1.3,cex=1.4,lty=2,lwd=8)
  legend(x=3,y=100,"Limit of detection: 100 CFU/ml", bty="n",cex=0.8)
  options(scipen=1)

  legend("topleft",legend=conc,col=mypalette,bty="n",pch=1,cex=0.8,title="conc.[mg/L]")
  axis.break(axis=1,breakpos=-1,style="slash") 
  abline(h =100, untf = FALSE,lty=2)
  
  for(idose in 1:nlevels(ab.conc)){
    dose<-as.numeric(levels(ab.conc)[idose])

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
    data2points<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=CFU[hours>=tmin & hours<=tmax &antibiotic==dose],cens=cens)
    data2regress<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=log(CFU[hours>=tmin & hours<=tmax &antibiotic==dose]),cens=cens)

    fit <- mle2(nll, start=list(a =log(1000000), b=-1,sigma=0.1), method="Nelder-Mead", control=list(maxit=1e3))
    
    print(fit)
    b <- coef(fit)[[2]]
    fita <- c(coef(fit)[[1]],coef(fit)[[2]])
    predict<-coef(fit)[[1]]+data2regress$x*coef(fit)[[2]]  
    lines(x_lines,y_lines,col=mypalette[idose],lty=1) 
    points(x_points,y_points,col=mypalette[idose],lty=1, pch=1, cex=1.5)
    slopes <- c(fita,slopes)
  }
  dev.off()
  k=k+1
}


#rm(list=ls())
graphics.off()
    
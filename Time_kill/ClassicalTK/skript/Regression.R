###################################################################################################################################
##################################fit linear regressions and plot lines############################################################
###################################################################################################################################
# set working directory
setwd("C:/Users/sfoerster/ABinteractions/Time_kill/ClassicalTK")
options(scipen=1)
# load libraries
library("drc")
library("XLConnect")
library("bbmle")
library("drc")
# define some functions:
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
#function to catch errors resulting from integer(0)
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}
#function to calculate decimalplaces for the antibiotic concentrations
decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
#function to replace filename with antibioticname
pattern <- c(".+AZ.*",".+CIP.*",".+Tet.*",".+CT.*")
replacement <- c("Azithromycin","Ciprofloxacin","Tetracyclin","Ceftriaxone")
replace <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}
#function to remove data when the antibiotic killed cells immediately (count of zero at 0 hours after adding the compound)
remove <- function (x){ 
  for (i in 1:x){
    mydatasub1 = subset(mydata,mydata$sample==problem$sample[1]&mydata$replicate==problem$replicate[1])
    mydatasub2 = subset(mydata,mydata$sample==problem$sample[2]&mydata$replicate==problem$replicate[2])
    mydatasub3 = subset(mydata,mydata$sample==problem$sample[3]&mydata$replicate==problem$replicate[3])
    sub<-rbind(mydatasub1,mydatasub2,mydatasub3)
  }
  return(sub)
}
#function to replace filename with antibioticname
# replace <- function (j){
#   titel=c(gsub(".+AZ.*","Azithromycin",j),gsub(".+CT.*","Ceftriaxon",j))
#   
# return(titel)
# }
#initialize some empty vectors and empty variables
sub=vector() 
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
#loop for the dataset with the recognized pattern, appends them to a big dataset containing all experiments (mydata)
print("---------------------------------run started---------------------------------")
for (ab in ablist){
  replicatespattern=ab
  replicates <- list.files(path =sourcefiledir, pattern = ab, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
  #read data from file with the recognized pattern and append to a dataframe
  for (replicate in replicates){
    append1 <- read.table(paste(verzeichnis,"/",replicate,sep=""),header=TRUE)
    append2 <- cbind(append1,replicate,ab,i)
    mydata<-rbind(append2,mydata)
    i=i+1
  }
}
mydata=cbind(mydata,index=rownames(mydata))
#remove values when already at timepoint 0 no colonies couldnt be counted anymore
problem<-subset(mydata,mydata$hours==0&mydata$count==0)
data_to_remove<-remove(length(problem$sample))
cleandata<-as.integer(rownames(data_to_remove))
 if(is.integer0(cleandata)==F){
   mydata <- mydata[-c(cleandata),]
 }
#plot time against CFU and make the summary statistics for all experiments done for a strain treated with an antibiotic (j)
for(j in unique(mydata$replicate)){
  data=mydata[mydata$replicate==j,]
  name=sub(".txt","regression.pdf",j)
  pdf(paste(name,sep=""))
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
  options(scipen=-1000)
  conc=levels(ab.conc)
  conc[conc==min(conc)]= 0
  conc=as.vector(round((as.numeric(conc)),decimalplaces(as.numeric(min(conc)))))
  #time range for regression
  tmin<-0
  tmax<-6
  slopes<-vector("numeric")
  options(scipen=1)
  #plot empty frame and legend so lines and points for each concentration can be added later on
  plot(c(-2,6),c(10,1e11),type="n",xlab="Time [h]",log="y",ylab="Bacteria [CFU/ml]",cex.lab=1.45,cex=1.4,lty=2,lwd=8)
  legend(x=3,y=100,"Limit of detection: 100 CFU/ml", bty="n",cex=0.8)
  options(scipen=1000,digits=5)
  legend("topleft",legend=conc,col=mypalette,bty="n",pch=1,cex=0.7,title="conc.[mg/L]")
  legend("topright",legend=replace(pattern,replacement,j),bty="n",cex=2)
  axis.break(axis=1,breakpos=-1,style="slash")
  abline(h =100, untf = FALSE,lty=2)
  #calculate the censoring (all values smaller or equal than 1 are censored)
  for(idose in 1:nlevels(ab.conc)){
    dose<-as.numeric(levels(ab.conc)[idose])
    cens=CFU[hours>=tmin&hours<=tmax & antibiotic==dose]
    #data for plotting the lines contains all data from -2 to 6 hours, data below detection are plottet as CFU/ml=100
    y_lines=CFU[hours>=-2&hours<=6 & antibiotic==dose]
    x_lines=hours[hours>=-2&hours<=6 & antibiotic==dose]
    #data for plotting the points contains only REAL measured data, data below detection limit are assigned NA 
    CFU2[CFU2<=1]<-NA
    hours2[CFU2<=1]<-NA
    x_points=hours2[hours2>=-2&hours2<=6 & antibiotic==dose]
    y_points=CFU2[hours2>=-2&hours2<=6 & antibiotic==dose]
    doses=antibiotic[hours>=-2&hours<=6 & antibiotic==dose]
    #calculate censoring and log likelihood
    cens[cens<=1]<-1
    cens[cens>1]<-0
    y_cens=log(CFU[hours>=tmin & hours<=tmax &antibiotic==dose])
    data2points<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=CFU[hours>=tmin & hours<=tmax &antibiotic==dose],cens=cens)
    data2regress<-data.frame(x=hours[hours>=tmin & hours<=tmax &antibiotic==dose],y_cens=log(CFU[hours>=tmin & hours<=tmax &antibiotic==dose]),cens=cens)
    fit <- mle2(nll, start=list(a =log(1000000), b=-1,sigma=0.1), method="Nelder-Mead", control=list(maxit=1e3))
    b <- coef(fit)[[2]]
    fita <- c(coef(fit)[[1]],coef(fit)[[2]])
    predict<-coef(fit)[[1]]+data2regress$x*coef(fit)[[2]]
    lines(x_lines,y_lines,col=mypalette[idose],lty=1)
    points(x_points,y_points,col=mypalette[idose],lty=1, pch=1, cex=1.5)
    slopes <- c(fita,slopes)
  }
  k=k+1
  dev.off()
}
#clean up workspace
print("---------------------------------run finished--------------------------------")
#rm(list=ls())
graphics.off()

#input
setwd("C:/Users/sfoerster/ABinteractions/Fluorescence/Auto-Time-kill")
sourcefiledir<-"input"
sourcefilepattern<-".txt"
files<-list.files(path =sourcefiledir, pattern = sourcefilepattern, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
print(files)
verzeichnis<-"input"
for (file in files) {
  outfilename <- paste(sub(".txt",".pdf", file), sep="/")
  df<-read.table(paste(verzeichnis,"/",file,sep=""),header=TRUE)
  outputfile <- paste("pharmacodynamics", "-",outfilename, sep="")
  pdf(outputfile)


xmax<-3
#plot(c(10,xmax),c(200,4000),type="n",log="y",xlab="time[min]",ylab="Fluorescence[AU]",lty=2,lwd=8,main="WHOA Azithromycin")
# legend("bottomright", 1.9, c("2xMIC", "1xMIC", "0.5xMIC","0.25xMIC","0.125xMIC","0xMIC ctr"), col = c(6,5,4,3,2,1,0),
#        text.col = "black", lty = c(1), pch = c(6,5,4,3,2,1,0),
#        merge = TRUE, bg = "white")
ab.conc=factor(df$conc)
antibiotic.conc=df$conc
hours=df$time
CFU.per.ml=df$Fluorescence



tmin<-20 # lower bound of the time interval over which to calculate the decline
tmax<-120 # upper bound of the time interval over which to calculate the decline

# define a vector for the decline rates:
slopes<-vector("numeric",length=nlevels(ab.conc))

print(hours)
# NET GROWTH RATES ARE CALCULATED BY LINEARLY REGRESSING 
# BACTERIAL CONCENTRATIONS OVER THE CHOSEN TIME INTERVAL
# THE DECLINE RATE IS THE COEFFICIENT OF REGRESSION

for(idose in 1:nlevels(ab.conc)){
  dose<-as.numeric(levels(ab.conc)[idose])
  data2regress<-data.frame(time=hours[hours>=tmin & hours<=tmax &antibiotic.conc==dose],log.bacteria.rel=log(CFU.per.ml[hours>=tmin&hours<=tmax & antibiotic.conc==dose]/CFU.per.ml[hours==tmin&antibiotic.conc==dose]))
  slopes[idose]<-coef(lm(log.bacteria.rel~time-1,data2regress))

}

dose=as.numeric(levels(ab.conc))
library("drc")
slope=round(slopes*60,2)
print(slope)
print(dose)
model=drm(slope~dose,fct=LL.4())
MIC=ED(model,99)
kappa=coefficients(model)[1]
  print(file)
print("kappa")
print(kappa)

#par(mar= c(5, 6, 6, 3) + 0.1)
#par(mgp= c(4, 1, 0))
par(mar=c(5,5,1,1))
mypalette<-c("red","lightblue","lightblue","lightblue","#C6DBEF", "#9ECAE1" ,"#6BAED6" ,"#4292C6" ,"#2171B5" ,"#08519C", "#08306B","black","black","black","black","black")
plot(model,ylim=c(-0.4,0.6),xlab="Antibiotic concentration[fold change MIC]",ylab="Change in Fluorescence [per hour]",lwd=2,pch=21,bg=mypalette,cex.lab=1.3,cex=1.4,lty=1)
    # legend("topright", 1.9 ,c("Control","0.2","0.4","0.8","1.56","3.125","6.25","12.5","25","50","100"), pt.bg = mypalette,
     #        text.col = "black", pch = 21,cex=0.75,title="conc.[xMIC]",
      #       merge = TRUE)
  print(file)
#print(summary(model))
MIC=ED(model,type="absolute",interval="fls",method="grid",respLev=0)
zMIC=MIC[1]

  #legend("topright",legend=c(paste("kappa=",round(kappa,4))),bty = "n")
library(plotrix)
axis.break(axis=1,breakpos=0.03,style="slash")   

#  print(data2fit)
# FITTING PROCEDURE

#fittedmodel<-nls(growth~psi(psimax,psimin,kappa,MIC,dose),
#                data=data2fit,
#               start=startpar)

#print(summary(fittedmodel)) # print a summary of the fit


# PLOT THE FITTED FUNCTION ONTO THE PREVIOUS PLOT
# (AS BEFORE, WE PLOT A CURVE SHIFTED SLIGHTLY,
# BY A TENTH OF THE LOWEST ANTIBIOTIC CONCENTRATION)

# curve(psi(coef(fittedmodel)[[1]],
#          coef(fittedmodel)[[2]],
#         coef(fittedmodel)[[3]],
#        coef(fittedmodel)[[4]],
#       x-as.numeric(levels(ab.conc))[2]/10),
#      add=TRUE,col=2)
dev.off()
}

graphics.off()

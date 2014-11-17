#input
library(drc)
library(plotrix)
setwd("C:/Users/sfoerster/ABinteractions/Fluorescence/Auto-Time-kill")
sourcefiledir<-"input"
sourcefilepattern<-".txt"
files<-list.files(path =sourcefiledir, pattern = sourcefilepattern, all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
print(files)
verzeichnis<-"input"
for (file in files) {
  outfilename <- paste(sub(".txt",".pdf", file), sep="/")
  df<-read.table(paste(verzeichnis,"/",file,sep=""),header=TRUE)
  outputfile <- paste("linear", "-",outfilename, sep="")
  pdf(outputfile)


xmax<-120
  # LET'S LOOK AT THE DATA:
  mypalette<-c("red","lightblue","lightblue","lightblue","#C6DBEF", "#9ECAE1" ,"#6BAED6" ,"#4292C6" ,"#2171B5" ,"#08519C", "#08306B","black","black","black","black","black")
  par(mar=c(5,5,1,1))
  #par(mar=c(6,6,6,6))
 # plot(c(0,xmax),c(300,4000),xlab="Time [h]",ylab="Bacteria [CFU/ml]",lwd=8)
  plot(c(20,xmax),c(500,4000),log="y",type="n",xlab="Time [min]",ylab="Fluorescence [AU]",cex.lab=1.3,cex=1.4,lty=2,lwd=8)
     legend("topleft", 1.9 ,c("Control","0.2","0.4","0.8","1.56","3.125","6.25","12.5","25","100"), col = mypalette,
        text.col = "black", lty = c(1), pch = 1,cex=0.75,title="conc.[xMIC]",
       merge = TRUE, bg = "white")
  
ab.conc=factor(df$conc)
antibiotic.conc=df$conc
hours=df$time
CFU.per.ml=df$Fluorescence
  mypalette<-c("red","lightblue","lightblue","lightblue","#C6DBEF", "#9ECAE1" ,"#6BAED6" ,"#4292C6" ,"#2171B5" ,"#08519C", "#08306B","black")
  
  for(idose in 1:nlevels(ab.conc)){
    dose<-as.numeric(levels(ab.conc)[idose])    
    points(hours[antibiotic.conc==dose],CFU.per.ml[antibiotic.conc==dose],
           col=mypalette[idose],pch=1,cex=2)
    lines(hours[antibiotic.conc==dose],CFU.per.ml[antibiotic.conc==dose],
          col=mypalette[idose],lwd=2,cex=1.5)
  } 

}



dev.off()


graphics.off()

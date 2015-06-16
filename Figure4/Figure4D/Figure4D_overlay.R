###################################################################################################################################
##################################growth rates with censoring and pharmacodynamics#################################################
###################################################################################################################################
#set working directory
setwd("C:/Users/Sunny/Dropbox/Time_kill/ClassicalTK")
library("drc")
library("plotrix")
library("magicaxis")
pdf("figure4D.pdf",width = 9, height = 8)
savefont <- par(font=2) 
means<-c()
m<-1
## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

color=c("black","maroon2","blue","orange","red","darkgrey","purple3","cyan3","green")
colorAlpha=add.alpha(color,alpha=0.6)
data<-readWorksheetFromFile("C:/Users/Sunny/Dropbox/Time_kill/ClassicalTK/summary_statistics_paper1.xlsx", sheet = 2, header = TRUE, startCol = 0, startRow = 0)
par(mar=c(4.7,5,2,7))
par(oma=c(0,0,0,0))
plot(c(0.0000001,100),c(-10,3),log="x",axes=F,bty="l",xlab=paste("antimicrobial concentration [mg/L])"), ylab=expression("Bacterial growth rate [h"^-1*"]"),cex.lab=1.8,col="white")
box()
magaxis(side=c(1,2),logpretty=TRUE,cex.axis=1.3)

psi<-function(psimax,psimin,kappa,MIC,conc) {
  psimax-(psimax-psimin)*(conc/MIC)^kappa/((conc/MIC)^kappa-psimin/psimax)
}
for(j in unique(data$Antibiotic)){
curve(psi(data$psi_max[m],data$psi_min[m],data$Kappa[m],data$zMIC[m],x),1e-10,1000,add=TRUE,col=color[m],lwd=3)
m=m+1
print(m)
}
strainname<-c("azithromycin","cefixime","ceftriaxone","chloramphenicol","ciprofloxacin","gentamycin","penicillin","spectinomycin","tetracycline")
legend("topleft",legend=letter,bty="n",inset=c(-0.1,-0.05), cex=3.6) 
legend("bottomleft",legend=strainname, bty="n",col=c("black","maroon2","blue","orange","red","darkgrey","purple3","cyan3","green"),lty=c(1,1,1,1,1,1),lwd=3, cex=1.5)
legend("topright",legend="strain: WT 103",bty="n",inset=0, cex=2) 
dev.off()
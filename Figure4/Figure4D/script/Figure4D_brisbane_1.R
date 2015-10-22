###################################################################################################################################
##################################growth rates with censoring and pharmacodynamics#################################################
###################################################################################################################################
#set working directory
setwd("C:/Users/sunny/Desktop/PDfunction/Figure4/Figure4D/output")
library("drc")
library("plotrix")
library("magicaxis")
library("XLConnect")
pdf("figure4D_1.pdf",width = 9, height = 8)
means<-c()
m<-1
letter=""
## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

color=c("black","maroon2","blue")
colorAlpha=add.alpha(color,alpha=0.6)
data<-readWorksheetFromFile("C:/Users/sunny/Desktop/PDfunction/Figure4/Figure4D/data/summary_statistics_9AB_1.xlsx", sheet = 2, header = TRUE, startCol = 0, startRow = 0)
par(mar=c(4.7,5,2,7))
par(oma=c(0,0,0,0))
plot(c(0.0000001,10000),c(-10,3),log="x",axes=F,bty="l",xlab=paste("Antimicrobial concentration [mg/L]"), ylab=expression("Bacterial growth rate [h"^-1*"]"),cex.lab=1.8,col="white")
box()
magaxis(side=c(1,2),logpretty=TRUE,cex.axis=1.3)

psi<-function(psimax,psimin,kappa,MIC,conc) {
  psimax-(psimax-psimin)*(conc/MIC)^kappa/((conc/MIC)^kappa-psimin/psimax)
}
for(j in unique(data$Antibiotic)){
curve(psi(data$psi_max[m],data$psi_min[m],data$Kappa[m],data$zMIC[m],x),1e-7,10000,add=TRUE,col=colorAlpha[m],lwd=6)
m=m+1
print(m)
}
strainname<-c("ciprofloxacin","gentamicin","spectinomycin")
legend("topleft",legend=letter,bty="n",inset=c(-0.1,-0.05), cex=3.6) 
legend("bottomleft",legend=strainname, bty="n",col=c("#00000099" ,"#EE30A799", "#0000FF99"),lty=c(1,1,1),lwd=6, cex=1.5)
dev.off()
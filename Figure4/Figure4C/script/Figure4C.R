###################################################################################################################################
##################################growth rates with censoring and pharmacodynamics#################################################
###################################################################################################################################
#set working directory
setwd("C:/Users/sunny/Desktop/PDfunction/Figure4/Figure4C/output")
library("drc")
library("plotrix")
library("magicaxis")
library("XLConnect")
pdf("figure4C.pdf",width = 9, height = 8)
#savefont <- par(font=2) 
means<-c()
m<-1
letter<-"(C)"
## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

color=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","red")
colorAlpha=add.alpha(color,alpha=0.6)
data<-readWorksheetFromFile("C:/Users/sunny/Desktop/PDfunction/Figure4/Figure4C/data/singleplot_statistics_panel.xlsx", sheet = 1, header = TRUE, startCol = 0, startRow = 0)
par(mar=c(4.7,5,2,7))
par(oma=c(0,0,0,0))
plot(c(0.0000001,10000),c(-10,3),log="x",axes=F,bty="l",xlab=paste("Antimicrobial concentration [mg/L]"), ylab=expression("Bacterial growth rate [h"^-1*"]"),cex.lab=1.8,col="white")
box()
magaxis(side=c(1,2),logpretty=TRUE,cex.axis=1.3)

psi<-function(psimax,psimin,kappa,MIC,conc) {
  psimax-(psimax-psimin)*(conc/MIC)^kappa/((conc/MIC)^kappa-psimin/psimax)
}
for(j in unique(data$strains)){
curve(psi(data$upper[m],data$lower[m],data$kappa[m],data$zMIC[m],x),1e-7,10000,add=TRUE,col=color[m],lwd=3)
m=m+1
print(m)
}
# strainname<-c("azithromycin","cefixime","ceftriaxone","chloramphenicol","ciprofloxacin","gentamycin","penicillin","spectinomycin","tetracycline")
strainname<-c("WHO G (LLR)","WHO K (HLR)","WHO L (HLR)","WHO M (R)","WHO N (R)")
legend("topleft",legend=letter,bty="n",inset=c(-0.1,-0.05), cex=3.6) 
legend("bottomleft",legend=strainname, bty="n",col=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E"),lty=c(1,1,1,1,1,1),lwd=3, cex=1.5)
dev.off()
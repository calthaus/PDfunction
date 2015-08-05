#load libraries
library("XLConnect")
library("drc")
library("XLConnect")
library("metafor")
library("drc")
library("plotrix")
library("magicaxis")
library("censReg")
means<-c()
m<-0
#read in data
setwd("C:/Users/sunny/Desktop/PDfunction/Table1/output")
data<-readWorksheetFromFile("C:/Users/sunny/Desktop/PDfunction/Table1/output/singleplot_statistics_DOGK18.xlsx", sheet = 1, 
                      header = TRUE, startCol = 0, 
                      startRow = 0)
parmlist=data.frame(data,row.names=NULL)
names(parmlist)<-c("Number","Identifier","Strain","Antibiotic", "Experiment","kappa","kappaerr","upper","uppererr","lower","lowererr","inflection","inflectionerr","zMIC","zMICerr","ETestMIC")
#data=data.frame(parmlist)
for (l in unique(parmlist$Experiment)){
    parmlistsub <-subset(parmlist,parmlist$Experiment==l)
    Strain=toString(parmlistsub$Strain[1])
    Antibiotic=toString(parmlistsub$Antibiotic[1])
    kappamodel <- rma(yi = parmlistsub$kappa, sei = parmlistsub$kappaerr, data = parmlistsub,weighted=FALSE, method = 'FE')
    kappasum <- data.frame(cbind(Strain,Antibiotic,signif(t(coef(summary(kappamodel))[1:6]),3)))
    names(kappasum)<-c("Strain","Antibiotic","Kappa")
    uppermodel <- rma(yi = parmlistsub$upper, sei = parmlistsub$uppererr, data = parmlistsub,weighted=FALSE, method = 'FE')
    uppersum <- data.frame(cbind(signif(t(coef(summary(uppermodel))[1:6]),3)))
    names(uppersum)<-c("psi_max")
    lowermodel <- rma(yi = parmlistsub$lower, sei = parmlistsub$lowererr, data = parmlistsub,weighted=FALSE, method = 'FE')
    lowersum <-  data.frame(cbind(signif(t(coef(summary(lowermodel))[1:6]),3)))
    names(lowersum)<-c("psi_min")
    inflectionmodel <- rma(yi = parmlistsub$inflection, sei = parmlistsub$inflectionerr,weighted=FALSE, data = parmlistsub, method = 'FE')
    inflectionsum <-  data.frame(cbind(signif(t(coef(summary(inflectionmodel))[1:6]),3)))
    names(inflectionsum)<-c("EC50")
    zMICmodel <- rma(yi = parmlistsub$zMIC, sei = parmlistsub$zMICerr, weighted=FALSE, data = parmlistsub, method = 'FE')
    zMICsum <-  data.frame(cbind(signif(t(coef(summary(zMICmodel))[1:6]),3)))
    names(zMICsum)<-c("zMIC")
    mean=cbind(kappasum,uppersum,lowersum,inflectionsum,zMICsum)
    mean=mean[-c(2, 3, 4), ]
    means=rbind(mean,means)
    m=m+1
}
estimate=means[grepl("estimate*", rownames(means)),]
lowerCI=means[grepl("ci.lb*", rownames(means)),]
names(lowerCI) <-c("Strain","Antibiotic","Kappa_lower","psi_max_lower","psi_min_lower","EC50_lower","zMIC_lower")
upperCI=means[grepl("ci.ub*", rownames(means)),]
names(upperCI) <-c("Strain","Antibiotic","Kappa_upper","psi_max_upper","psi_min_upper","EC50_upper","zMIC_upper")
newdf=cbind(estimate,lowerCI,upperCI)
newmeans=newdf[-c(8,9,15,16)]
#write summary statistics to excel file
 allplots <- loadWorkbook("summary_statistics_Table1.xlsx", create = TRUE)
 createSheet(allplots, name = "parameterlist")
 createSheet(allplots, name = "estimates")
 writeWorksheet(allplots, newmeans, sheet = "parameterlist", startRow = 1, startCol = 1)
 writeWorksheet(allplots, estimate, sheet = "estimates", startRow = 1, startCol = 1)
 saveWorkbook(allplots)




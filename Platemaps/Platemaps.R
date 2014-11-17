#dev.new(width=6,height=4)
setwd("C:/Users/sfoerster/ABinteractions/Platemaps")
platelay=read.table("input.txt",header=T)
pdf("platemap_color2.pdf",width=10,height=8)
rown <- unique(platelay$rown)
coln <- unique(platelay$coln)

plot(NA,ylim=c(0.3,length(rown)+0.5),xlim=c(0.3,length(coln)+0.5),ann=FALSE,axes=FALSE)
#box()
box()
axis(2,at=seq_along(rown),labels=rev(rown),las=2)
axis(1,at=seq_along(coln),labels=coln)
mypalette<-c("red","lightblue","lightblue","lightblue","#C6DBEF", "#9ECAE1" ,"#6BAED6" ,"#4292C6" ,"#2171B5" ,"#08519C", "#08306B","black")
#colgrp <- findInterval(platelay$colorvar,seq(min(platelay$colorvar),max(platelay$colorvar),length.out=10))
#colfunc <- colorRampPalette(c("green", "blue"))
#collist <- colfunc(length(unique(colgrp))) 

symbols(platelay$coln,
        factor(platelay$rown, rev(levels(platelay$rown))),
        circles=rep(0.4,nrow(platelay)),
        add=TRUE,
        inches=FALSE,
        bg=mypalette)
dev.off()
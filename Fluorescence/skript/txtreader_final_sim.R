#!/usr/bin/Rscript
library("gplots")
library("bmd")
library("gplots")
library("splines")
library("plyr")
library("emdbook")
setwd("C:/users/sfoerster/ABinteractions/Fluorescence")
Sys.setlocale(locale="C") # Fixes "input string X is invalid in this locale"

outdir <- "output"
sourcefilepattern<-".txt$"
testfolders <- c("input_data")
testfolder <- c("input_data")

test <- function(){
      source("skript/4in1skript_ankervariation_backup.R")
      for (testfolder in testfolders) {
          files <- list.files(path=testfolder, pattern = sourcefilepattern, all.files = FALSE, recursive = FALSE,
                            ignore.case = FALSE, include.dirs = FALSE)
          print(files)
          for (thisfile in files){
          #thisfile = "Azithromycin_49.txt"
          #thisfile = "BIBF_k.txt"
              nr <- as.integer(sub(unlist(strsplit(thisfile, "_"))[5], pattern=".txt", replacement=""))
              myData <- read.table(paste(testfolder,"/",thisfile,sep=""),header=TRUE)
              myData <- myData[, c("dose","response","experiment")]
              doseunit <- "Concentration"
              od <- paste(outdir, sep="/")
              dir.create(od, recursive=TRUE, showWarnings=FALSE)
              plotname <- paste(od, sub(".txt",".png", thisfile), sep="/") # for the plots
              outname <- sub(".txt","",thisfile) # for the tables
              outname <- sub("_", " ",outname)
              outname <- paste(outname, sep="")
              #outname<-("blub.txt")
              doseunit<-paste("Concentration [AU]",sep="")
              simulation = FALSE
              nsteps = 400
              figures = TRUE
              run = 4
              # Simulation
              if (simulation) {
                  for (plotnr in 1:nsteps) {
                      try(processData(myData,outname,doseunit,plotname,run=4,plotnr,nsteps,figures=figures))
                  }
              } else {
                  try(processData(myData,outname,doseunit,plotname,run,plotnr=100,figures=figures))
              }
      }
  }
}


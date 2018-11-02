#######################
#gating.R       
#
#started: 10/23/2016
#
#author1: D Gresham
#modified: Nathan Brandt
#Purpose: To produce gates for the GAP1::GFP Screen
######################

###########################################################################################################################
#This script is intended to read in .fcs files and perform manual gating for i) single cells, ii) debris and iii) fluorescence 
#
#Gating is performed with untransformed data
#
#Individual gates are saved in a file gates.Rdata for use with the Gresham Lab Flow Cytometry Analysis.Rmd pipeline
###########################################################################################################################

##To be run the first time if packages are not installed.
#source("http://bioconductor.org/biocLite.R")
#biocLite("flowViz")
#biocLite("flowCore")

#Load libraries
library(flowCore)
library(flowViz)

#Read in the data
flowData <- read.flowSet(path = "/Users/Brandt/Google Drive/Transfer_Folder/GAP1_GFP_Screen/GAP1_GFP_Gating_Controls", pattern=".fcs", alter.names=TRUE)

#sampel to use in gating
s <- 7

#Check how many cells were counted in each fcs file
fsApply(flowData, each_col, length)

##############################
#1. Generate gate for singlet cells
#this gate is defined on the basis of the relationship between forward scatter height and area
plot(flowData[[s]], c('FSC.H','FSC.A'), xlim=c(0,3e6), ylim=c(0,3e6),smooth=F)
Agate <- locator(10)
gm.1 <- matrix(,10,2)
colnames(gm.1) <- c('FSC.H','FSC.A')
gm.1[,1] <- Agate$x
gm.1[,2] <- Agate$y
pg.singlets <- polygonGate(filterId="singlets",.gate=gm.1)

#test that the singlet gate looks reasonable for the sample
xyplot(FSC.A~FSC.H,data=flowData[[s]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)

#test that the gate looks reasonable over all the samples
xyplot(FSC.A~FSC.H, data=flowData, xlim=c(0,3e6), ylim=c(0,3e6), 
       smooth=F, filter=pg.singlets, outline=T, displayFilter=TRUE,
       stat=T, pos=0.5, abs=T)

##############################
#2. Generate Gate for debris based on forward scatter and side scatter
plot(flowData[[s]], c('FSC.A','SSC.A'), xlim=c(0,3e6), ylim=c(0,5e5), smooth=F)
Bgate <- locator(10)
gm.2 <- matrix(,10,2)
colnames(gm.2) <- c('FSC.A','SSC.A')
gm.2[,1] <- Bgate$x
gm.2[,2] <- Bgate$y
pg.nondebris <- polygonGate(filterId="nonDebris",.gate=gm.2)

#test that the debris gate looks reasonable for the sample
xyplot(SSC.A ~ FSC.A, data=flowData[[s]], displayFilter=TRUE, xlim=c(0,3e6), ylim=c(0,5e5), filter=pg.nondebris, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#test that the decris gate looks reasonable over all the samples
xyplot(SSC.A ~ FSC.A, data=flowData, displayFilter=TRUE, xlim=c(0,3e6), ylim=c(0,5e5), filter=pg.nondebris, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

##############################
#3. FL1
#  a) Generate Gate for non-fluorescencing cells
plot(flowData[[s]], c('FSC.A','FL1.A'), xlim=c(0,3e6), ylim=c(0,2e4), smooth=F)
Cgate <- locator(10)
gm.3 <- matrix(,10,2)
colnames(gm.3) <- c('FSC.A','FL1.A')
gm.3[,1] <- Cgate$x
gm.3[,2] <- Cgate$y
pg.FL1neg <- polygonGate(filterId="FL1neg",.gate=gm.3)

#test that the non-fluorescing gate looks reasonable for the sample
xyplot(FL1.A~FSC.A,data=flowData[[s]], xlim=c(0,5e6), ylim=c(0,5e4), filter=pg.FL1neg, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#test that the non-fluorescing gate looks reasonable over all the samples
xyplot(FL1.A~FSC.A,data=flowData, xlim=c(0,5e6), ylim=c(0,5e4), filter=pg.FL1neg, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#  b) Generate Gate for fluorescencing cells
plot(flowData[[s]], c('FSC.A','FL1.A'), xlim=c(0,5e6), ylim=c(0,5e5), smooth=F)
Dgate <- locator(10)
gm.4 <- matrix(,10,2)
colnames(gm.4) <- c('FSC.A','FL1.A')
gm.4[,1] <- Dgate$x
gm.4[,2] <- Dgate$y
pg.FL1pos <- polygonGate(filterId="FL1pos",.gate=gm.4)

#test that the fluorescing gate looks reasonable for the sample
xyplot(FL1.A~FSC.A,data=flowData[[s]], xlim=c(0,5e6), ylim=c(0,5e5), filter=pg.FL1pos, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#test that the fluorescing gate looks reasonable over all the samples
xyplot(FL1.A~FSC.A,data=flowData, xlim=c(0,5e6), ylim=c(0,5e5), filter=pg.FL1pos, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

################################

##############################
#4. FL2
#  a) Generate Gate for non-fluorescencing cells
plot(flowData[[s]], c('FSC.A','FL2.A'), xlim=c(0,3e6), ylim=c(0,2e4), smooth=F)
Egate <- locator(10)
gm.5 <- matrix(,10,2)
colnames(gm.5) <- c('FSC.A','FL2.A')
gm.5[,1] <- Egate$x
gm.5[,2] <- Egate$y
pg.FL2neg <- polygonGate(filterId="FL2neg",.gate=gm.5)

#test that the non-fluorescing gate looks reasonable for the sample
xyplot(FL2.A~FSC.A,data=flowData[[s]], xlim=c(0,5e6), ylim=c(0,5e4), filter=pg.FL2neg, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#test that the non-fluorescing gate looks reasonable over all the samples
xyplot(FL2.A~FSC.A,data=flowData, xlim=c(0,5e6), ylim=c(0,5e4), filter=pg.FL2neg, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#  b) Generate Gate for fluorescencing cells
plot(flowData[[s]], c('FSC.A','FL2.A'), xlim=c(0,5e6), ylim=c(0,5e5), smooth=F)
Fgate <- locator(10)
gm.6 <- matrix(,10,2)
colnames(gm.6) <- c('FSC.A','FL2.A')
gm.6[,1] <- Fgate$x
gm.6[,2] <- Fgate$y
pg.FL2pos <- polygonGate(filterId="FL2pos",.gate=gm.6)

#test that the fluorescing gate looks reasonable for the sample
xyplot(FL2.A~FSC.A,data=flowData[[s]], xlim=c(0,5e6), ylim=c(0,5e5), filter=pg.FL2pos, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#test that the fluorescing gate looks reasonable over all the samples
xyplot(FL2.A~FSC.A,data=flowData, xlim=c(0,5e6), ylim=c(0,5e5), filter=pg.FL2pos, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

################################

#6. Save the gate information to an R data file

rm(list=c("flowData"))
save.image(file="gates_Sample3.Rdata")



#######################
#gating.R       
#
#started: 01/07/2016
#
#author1: S Lauer, G Avecilla, D Gresham
######################

######################
#This script is specific for analyzing the data obtained in S. Lauer's experiment: SLL08_EE06_Mini02.
#For any other purpose, the script must be modified accordingly. 

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
flowData <- read.flowSet(path = ".", pattern=".fcs", alter.names=TRUE)

##Confirm the number of .fcs files in your folder. The script below is only accurate if there are 32 .fcs files
str(flowData) #it should say there are 32 observations. if not, rework the script

##############################
#1. Generate gate for singlet cells.. glucose needs a different gate
#this gate is defined on the basis of the relationship between forward scatter height and area
plot(flowData[[32]], c('FSC.H','FSC.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,3e6),smooth=F)
Agate <- locator(10)
gm.1 <- matrix(,10,2)
colnames(gm.1) <- c('FSC.H','FSC.A')
gm.1[,1] <- Agate$x
gm.1[,2] <- Agate$y
pg.singlets <- polygonGate(filterId="singlets",.gate=gm.1)

polygon(Agate$x, Agate$y, border='red')

#test that the singlet gate looks reasonable for some samples (these are all 1 copy controls)
xyplot(FSC.A~FSC.H,data=flowData[[32]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)
xyplot(FSC.A~FSC.H,data=flowData[[29]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)

#a two copy control
xyplot(FSC.A~FSC.H,data=flowData[[24]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)

####GLUCOSE###
plot(flowData[[3]], c('FSC.H','FSC.A'), xlim=c(0,4e6), xaxs = "i", yaxs = "i", ylim=c(0,4e6),smooth=F)
Aagate <- locator(10)
gm.1a <- matrix(,10,2)
colnames(gm.1a) <- c('FSC.H','FSC.A')
gm.1a[,1] <- Aagate$x
gm.1a[,2] <- Aagate$y
pg.singletsGLU <- polygonGate(filterId="singletsGLU",.gate=gm.1a)

polygon(Aagate$x, Aagate$y, border='red')
xyplot(FSC.A~FSC.H,data=flowData[[3]],xlim=c(0,4e6), ylim=c(0,4e6), smooth=F, filter=pg.singletsGLU, outline=T)
xyplot(FSC.A~FSC.H,data=flowData[[2]],xlim=c(0,4e6), ylim=c(0,4e6), smooth=F, filter=pg.singletsGLU, outline=T)
xyplot(FSC.A~FSC.H,data=flowData[[11]],xlim=c(0,4e6), ylim=c(0,4e6), smooth=F, filter=pg.singletsGLU, outline=T)

##############################
#2. Generate Gate for debris based on forward scatter and side scatter. 
#This needs to be done separately for each media condition.

######GLUCOSE#####
plot(flowData[[3]], c('FSC.A','SSC.A'), xlim=c(0,4e6), xaxs = "i", yaxs = "i", ylim=c(0,1e6), smooth=F)
Bgate <- locator(10)
gm.2 <- matrix(,10,2)
colnames(gm.2) <- c('FSC.A','SSC.A')
gm.2[,1] <- Bgate$x
gm.2[,2] <- Bgate$y
pg.nondebris.glu <- polygonGate(filterId="nonDebrisGlu",.gate=gm.2)

#test that the debris gate looks reasonable for some samples (these are all 1 copy controls)
xyplot(SSC.A ~ FSC.A, data=flowData[[3]], displayFilter=TRUE, xlim=c(0,4e6), ylim=c(0,3e6), filter=pg.nondebris.glu, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)
xyplot(SSC.A ~ FSC.A, data=flowData[[2]], displayFilter=TRUE, xlim=c(0,4e6), ylim=c(0,3e6), filter=pg.nondebris.glu, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)
xyplot(SSC.A ~ FSC.A, data=flowData[[11]], displayFilter=TRUE, xlim=c(0,4e6), ylim=c(0,3e6), filter=pg.nondebris.glu, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#######UREA######
plot(flowData[[29]], c('FSC.A','SSC.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,6e5), smooth=F)
Cgate <- locator(10)
gm.3 <- matrix(,10,2)
colnames(gm.3) <- c('FSC.A','SSC.A')
gm.3[,1] <- Cgate$x
gm.3[,2] <- Cgate$y
pg.nondebris.ur <- polygonGate(filterId="nonDebrisUr",.gate=gm.3)

#test that the debris gate looks reasonable for some samples (these are all 1 copy controls)
xyplot(SSC.A ~ FSC.A, data=flowData[[29]], displayFilter=TRUE, xlim=c(0,3e6), ylim=c(0,1e6), filter=pg.nondebris.ur, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)
xyplot(SSC.A ~ FSC.A, data=flowData[[6]], displayFilter=TRUE, xlim=c(0,3e6), ylim=c(0,1e6), filter=pg.nondebris.ur, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

######GLUTAMINE#####
plot(flowData[[32]], c('FSC.A','SSC.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,6e5), smooth=F)
Dgate <- locator(10)
gm.4 <- matrix(,10,2)
colnames(gm.4) <- c('FSC.A','SSC.A')
gm.4[,1] <- Dgate$x
gm.4[,2] <- Dgate$y
pg.nondebris.gln <- polygonGate(filterId="nonDebrisGln",.gate=gm.4)

#test that the debris gate looks reasonable for some samples (these are all 1 copy controls)
xyplot(SSC.A ~ FSC.A, data=flowData[[24]], displayFilter=TRUE, xlim=c(0,2e6), ylim=c(0,1e6), filter=pg.nondebris.gln, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)
xyplot(SSC.A ~ FSC.A, data=flowData[[32]], displayFilter=TRUE, xlim=c(0,2e6), ylim=c(0,1e6), filter=pg.nondebris.gln, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)

#############################
#3. FLUORESCENCE

##Generate gates for 0, 1, 2, and 3+ copies for GLUCOSE-limitation##

##Plot the control sample that has non-fluorescing cells (0 copy)
plot(flowData[[3]], c('FSC.A','FL1.A'), xlim=c(-1e4,3e6), xaxs = "i", yaxs = "i", ylim=c(-1e4,1e5), smooth=F)
htgate.0 <- locator(10, type='l', col='red')

##Plot the control sample that has 1 copy
plot(flowData[[3]], c('FSC.A','FL1.A'), xlim=c(0,5e6), xaxs = "i", yaxs = "i", ylim=c(0,2e6), smooth=F)

##Overlay your previous gate
polygon(htgate.0$x, htgate.0$y, border='red')

##Draw a new gate for the one copy
htgate.1 <- locator(10, type='l', col='blue')

##Overlay and check the new gate
plot(flowData[[3]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)
polygon(htgate.0$x, htgate.0$y, border='red')
polygon(htgate.1$x, htgate.1$y, border='blue')

##Plot the control sample that has 2 copies
plot(flowData[[11]], c('FSC.A','FL1.A'), xlim=c(0,5e6), xaxs = "i", yaxs = "i", ylim=c(0,2e6), smooth=F)
polygon(htgate.0$x, htgate.0$y, border='red')
polygon(htgate.1$x, htgate.1$y, border='blue')
htgate.2 <- locator(10, type='l', col='green')
htgate.3 <- locator(10, type="l", col='purple')

##Check how the gates look on the sample
plot(flowData[[11]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)
polygon(htgate.0$x, htgate.0$y, border='red')
polygon(htgate.1$x, htgate.1$y, border='blue')
polygon(htgate.2$x, htgate.2$y, border='green')
polygon(htgate.3$x, htgate.3$y, border='purple')

htgm.0 <- matrix(,10,2)
htgm.1 <- matrix(,10,2)
htgm.2 <- matrix(,10,2)
htgm.3 <- matrix(,10,2)

colnames(htgm.0) <- c('FSC.A','FL1.A')
colnames(htgm.1) <- c('FSC.A','FL1.A')
colnames(htgm.2) <- c('FSC.A','FL1.A')
colnames(htgm.3) <- c('FSC.A','FL1.A')

htgm.0[,1] <- htgate.0$x
htgm.1[,1] <- htgate.1$x
htgm.2[,1] <- htgate.2$x
htgm.3[,1] <- htgate.3$x

htgm.0[,2] <- htgate.0$y
htgm.1[,2] <- htgate.1$y
htgm.2[,2] <- htgate.2$y
htgm.3[,2] <- htgate.3$y

pg.zcht<- polygonGate(filterId="ZeroCopyGlucose",.gate=htgm.0)
pg.ocht<- polygonGate(filterId="OneCopyGlucose",.gate=htgm.1)
pg.twcht<- polygonGate(filterId="TwoCopyGlucose",.gate=htgm.2)
pg.thrcht<- polygonGate(filterId="ThreeCopyGlucose",.gate=htgm.3)


#################################################################
##Generate gates for 0, 1, 2, and 3+ copies for UREA-limitation##

##Plot the control sample that has non-fluorescing cells (0 copy)
plot(flowData[[29]], c('FSC.A','FL1.A'), xlim=c(-1e4,3e6), xaxs = "i", yaxs = "i", ylim=c(-1e4,1e5), smooth=F)
urgate.0 <- locator(10, type='l', col='red')

##Plot the control sample that has 1 copy
plot(flowData[[29]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)

##Overlay your previous gate
polygon(urgate.0$x, urgate.0$y, border='red')

##Draw a new gate for the one copy
urgate.1 <- locator(10, type='l', col='blue')

##Overlay and check the new gate
plot(flowData[[29]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)
polygon(urgate.0$x, urgate.0$y, border='red')
polygon(urgate.1$x, urgate.1$y, border='blue')

##Plot the control sample that has 2 copies
plot(flowData[[6]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,1e6), smooth=F)
polygon(urgate.0$x, urgate.0$y, border='red')
polygon(urgate.1$x, urgate.1$y, border='blue')
urgate.2 <- locator(10, type='l', col='green')
urgate.3 <- locator(10, type="l", col='purple')

##Check how the gates look on the sample
plot(flowData[[6]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)
polygon(urgate.0$x, urgate.0$y, border='red')
polygon(urgate.1$x, urgate.1$y, border='blue')
polygon(urgate.2$x, urgate.2$y, border='green')
polygon(urgate.3$x, urgate.3$y, border='purple')

urgm.0 <- matrix(,10,2)
urgm.1 <- matrix(,10,2)
urgm.2 <- matrix(,10,2)
urgm.3 <- matrix(,10,2)

colnames(urgm.0) <- c('FSC.A','FL1.A')
colnames(urgm.1) <- c('FSC.A','FL1.A')
colnames(urgm.2) <- c('FSC.A','FL1.A')
colnames(urgm.3) <- c('FSC.A','FL1.A')

urgm.0[,1] <- urgate.0$x
urgm.1[,1] <- urgate.1$x
urgm.2[,1] <- urgate.2$x
urgm.3[,1] <- urgate.3$x

urgm.0[,2] <- urgate.0$y
urgm.1[,2] <- urgate.1$y
urgm.2[,2] <- urgate.2$y
urgm.3[,2] <- urgate.3$y

ur.zcht<- polygonGate(filterId="ZeroCopyUrea",.gate=urgm.0)
ur.ocht<- polygonGate(filterId="OneCopyUrea",.gate=urgm.1)
ur.twcht<- polygonGate(filterId="TwoCopyUrea",.gate=urgm.2)
ur.thrcht<- polygonGate(filterId="ThreeCopyUrea",.gate=urgm.3)

#################################################################
##Generate gates for 0, 1, 2, and 3+ copies for GLUTAMINE-limitation##

##Plot the control sample that has non-fluorescing cells (0 copy)
plot(flowData[[32]], c('FSC.A','FL1.A'), xlim=c(-1e4,3e6), xaxs = "i", yaxs = "i", ylim=c(-1e4,1e5), smooth=F)
glngate.0 <- locator(10, type='l', col='red')

##Plot the control sample that has 1 copy
plot(flowData[[24]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)

##Overlay your previous gate
polygon(glngate.0$x, glngate.0$y, border='red')

##Draw a new gate for the one copy
glngate.1 <- locator(10, type='l', col='blue')

##Overlay and check the new gate
plot(flowData[[24]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)
polygon(glngate.0$x, glngate.0$y, border='red')
polygon(glngate.1$x, glngate.1$y, border='blue')

##Plot the control sample that has 2 copies
plot(flowData[[32]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,1e6), smooth=F)
polygon(glngate.0$x, glngate.0$y, border='red')
polygon(glngate.1$x, glngate.1$y, border='blue')
glngate.2 <- locator(10, type='l', col='green')
glngate.3 <- locator(10, type="l", col='purple')

##Check how the gates look on another sample
plot(flowData[[24]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)
polygon(glngate.0$x, glngate.0$y, border='red')
polygon(glngate.1$x, glngate.1$y, border='blue')
polygon(glngate.2$x, glngate.2$y, border='green')
polygon(glngate.3$x, glngate.3$y, border='purple')

glngm.0 <- matrix(,10,2)
glngm.1 <- matrix(,10,2)
glngm.2 <- matrix(,10,2)
glngm.3 <- matrix(,10,2)

colnames(glngm.0) <- c('FSC.A','FL1.A')
colnames(glngm.1) <- c('FSC.A','FL1.A')
colnames(glngm.2) <- c('FSC.A','FL1.A')
colnames(glngm.3) <- c('FSC.A','FL1.A')

glngm.0[,1] <- glngate.0$x
glngm.1[,1] <- glngate.1$x
glngm.2[,1] <- glngate.2$x
glngm.3[,1] <- glngate.3$x

glngm.0[,2] <- glngate.0$y
glngm.1[,2] <- glngate.1$y
glngm.2[,2] <- glngate.2$y
glngm.3[,2] <- glngate.3$y

gln.zcht<- polygonGate(filterId="ZeroCopyGln",.gate=glngm.0)
gln.ocht<- polygonGate(filterId="OneCopyGln",.gate=glngm.1)
gln.twcht<- polygonGate(filterId="TwoCopyGln",.gate=glngm.2)
gln.thrcht<- polygonGate(filterId="ThreeCopyGln",.gate=glngm.3)

#Save the gate information to an R data file
rm(list=c("flowData")) 
save.image(file="gates_05252017_MIDTPT.Rdata")

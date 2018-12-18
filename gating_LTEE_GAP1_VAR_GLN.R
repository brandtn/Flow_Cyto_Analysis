#######################
#gating.R       
#
#started: 01/07/2016
#modified: 12/18/2018
#
#author1: G Avecilla, S Lauer, D Gresham
#author2: N Brandt
######################

######################
#This script is specific for analyzing the data obtained in LTEE_GAP1_Variants in Gln,
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

#Set working directory to the folder in which you have stored your .fcs files
#Read in all the fcs files in the directory.


#working directory
dir = '.'

#file location
path.data = "/Users/Brandt/Google Drive/MiniStatRun_10_2018/"


list.folders <- c("LTEE_mCitrine_GAP1_Variants_T00", "LTEE_mCitrine_GAP1_Variants_T06", "LTEE_mCitrine_GAP1_Variants_T07", "LTEE_mCitrine_GAP1_Variants_T08.1", "LTEE_mCitrine_GAP1_Variants_T08.2", "LTEE_mCitrine_GAP1_Variants_T08.3", "LTEE_mCitrine_GAP1_Variants_T11.1", "LTEE_mCitrine_GAP1_Variants_T11.2", "LTEE_mCitrine_GAP1_Variants_T13.1", "LTEE_mCitrine_GAP1_Variants_T13.2", "LTEE_mCitrine_GAP1_Variants_T14", "LTEE_mCitrine_GAP1_Variants_T15", "LTEE_mCitrine_GAP1_Variants_T18")

#fcs run sample name
name <- "LTEE_mCitrine_GAP1_Variants_TON"

#name <- list.folders[1]

flowData <- read.flowSet(path = paste(path.data, name,"/", sep=""), pattern=".fcs", alter.names = TRUE)

sample.sheet <- read.csv(paste(path.data,"samplesheet_",name,".csv", sep=""))

#Adds a sample sheet data to the pData of the flowset

sampleNames(flowData) <- paste(gsub(" ","_",sample.sheet$Strain),"_",sub(" ","_",sample.sheet$Well), sep="")

#Need to determine samplesheet inputs
pData(flowData)$Well <- sample.sheet$Well
pData(flowData)$Strain <- sample.sheet$Strain
pData(flowData)$Genotype <- sample.sheet$Genotype
pData(flowData)$Ploidy <- sample.sheet$Ploidy
pData(flowData)$Media <- sample.sheet$Media
pData(flowData)$Experiment <- sample.sheet$Experiment

gateData <- flowData

##Confirm the number of .fcs files in your folder. The script below is only accurate if there are 32 .fcs files
str(gateData) #it should say there are 32 observations. if not, rework the script

zerocopy <- 1
onecopy <- 3
two copy <-4

##############################
#1. Generate gate for singlet cells####
#this gate is defined on the basis of the relationship between forward scatter height and area
plot(gateData[[zerocopy]], c('FSC.H','FSC.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,3e6),smooth=F)
Agate <- locator(10, type='l', col='red')
gm.1 <- matrix(,length(Agate$x),2)
colnames(gm.1) <- c('FSC.H','FSC.A')
gm.1[,1] <- Agate$x
gm.1[,2] <- Agate$y
pg.singlets <- polygonGate(filterId="singlets",.gate=gm.1)

#test that the singlet gate looks reasonable for some samples (these are all 0 copy controls)
xyplot(FSC.A~FSC.H,data=gateData[[16]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)
xyplot(FSC.A~FSC.H,data=gateData[[15]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)
xyplot(FSC.A~FSC.H,data=gateData[[14]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)


#a zero copy control
xyplot(FSC.A~FSC.H,data=gateData[[zerocopy]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)
#a one copy control
xyplot(FSC.A~FSC.H,data=gateData[[onecopy]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)
#a two copy control
xyplot(FSC.A~FSC.H,data=gateData[[twocopy]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.singlets, outline=T)


##############################
#2. Generate Gate for debris based on forward scatter and side scatter. ####
#This needs to be done separately for each media condition.

######GLUTAMINE#####
plot(gateData[[zerocopy]], c('FSC.A','SSC.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,6e5), smooth=F)
Bgate <- locator(10, type='l', col='red')
gm.2 <- matrix(,length(Bgate$x),2)
colnames(gm.2) <- c('FSC.A','SSC.A')
gm.2[,1] <- Bgate$x
gm.2[,2] <- Bgate$y
pg.nondebris <- polygonGate(filterId="nonDebris",.gate=gm.2)

#test that the debris gate looks reasonable for some samples (these are 1 and 2 copy controls)
xyplot(SSC.A ~ FSC.A, data=gateData[[13]], displayFilter=TRUE, xlim=c(0,2e6), ylim=c(0,1e6), filter=pg.nondebris, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)
xyplot(SSC.A ~ FSC.A, data=gateData[[12]], displayFilter=TRUE, xlim=c(0,2e6), ylim=c(0,1e6), filter=pg.nondebris, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)
xyplot(SSC.A ~ FSC.A, data=gateData[[11]], displayFilter=TRUE, xlim=c(0,2e6), ylim=c(0,1e6), filter=pg.nondebris, smooth=F, xbin=1024, stat=T, pos=0.5, abs=T)


#a zero copy control
xyplot(FSC.A~FSC.H,data=gateData[[zerocopy]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.nondebris, outline=T)
#a one copy control
xyplot(FSC.A~FSC.H,data=gateData[[onecopy]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.nondebris, outline=T)
#a two copy control
xyplot(FSC.A~FSC.H,data=gateData[[twocopy]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=pg.nondebris, outline=T)


#############################
#3. FLUORESCENCE####

####Generate gates for 0, 1, 2, and 3+ copies####

##Plot the control sample that has non-fluorescing cells (0 copy)
plot(gateData[[zerocopy]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(-1e4,1e5), smooth=F)
Cgate <- locator(10, type='l', col='red')
gm.3 <- matrix(,length(Cgate$x),2)
colnames(gm.3) <- c('FSC.A','SSC.A')
gm.3[,1] <- Cgate$x
gm.3[,2] <- Cgate$y
fl1gate.0 <- polygonGate(filterId="zeroFL1",.gate=gm.3)

#a zero copy control
xyplot(FSC.A~FSC.H,data=gateData[[zerocopy]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=fl1gate.0, outline=T)
#a one copy control
xyplot(FSC.A~FSC.H,data=gateData[[onecopy]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=fl1gate.0, outline=T)
#a two copy control
xyplot(FSC.A~FSC.H,data=gateData[[twocopy]],xlim=c(0,3e6), ylim=c(0,3e6), smooth=F, filter=fl1gate.0, outline=T)


##Draw a new gate for the one copy
plot(gateData[[onecopy]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(-1e4,1e5), smooth=F)
Dgate <- locator(10, type='l', col='blue')
gm.4 <- matrix(,length(Dgate$x),2)
colnames(gm.4) <- c('FSC.A','FL1.A')
gm.4[,4] <- Dgate$x
gm.4[,4] <- Dgate$y
fl1gate.1 <- polygonGate(filterId="oneCopyFL1",.gate=gm.4)


##Overlay and check the new gate
plot(gateData[[onecopy]], c('FSC.A','FL1.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,3e5), smooth=F)
polygon(fl1gate.0$x, fl1gate.0$y, border='red')
polygon(fl1gate.1$x, fl1gate.1$y, border='blue')

##Plot the control sample that has 2 copies
plot(gateData[[twocopy]], c('FSC.A','FL1.A'), xlim=c(0,2e6), xaxs = "i", yaxs = "i", ylim=c(0,4e5), smooth=F)
polygon(fl1gate.0$x, fl1gate.0$y, border='red')
polygon(fl1gate.1$x, fl1gate.1$y, border='blue')
fl1gate.2 <- locator(10, type='l', col='green')
fl1gate.3 <- locator(10, type="l", col='purple')

##Check how the gates look on the sample
plot(gateData[[onecopy]], c('FSC.A','FL1.A'), xlim=c(0,3e6), xaxs = "i", yaxs = "i", ylim=c(0,5e5), smooth=F)




polygon(glngate.0$x, glngate.0$y, border='red')
polygon(glngate.1$x, glngate.1$y, border='blue')
polygon(glngate.2$x, glngate.2$y, border='green')
polygon(glngate.3$x, glngate.3$y, border='purple')

glngm.0 <- matrix(,6,2)
glngm.1 <- matrix(,5,2)
glngm.2 <- matrix(,7,2)
glngm.3 <- matrix(,6,2)

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

gln.zero<- polygonGate(filterId="ZeroCopyGlutamine",.gate=glngm.0)
gln.one<- polygonGate(filterId="OneCopyGlutamine",.gate=glngm.1)
gln.two<- polygonGate(filterId="TwoCopyGlutamine",.gate=glngm.2)
gln.three<- polygonGate(filterId="ThreeCopyGlutamine",.gate=glngm.3)



#Save the gate information to an R data file
rm(list=c("gateData")) 
save(pg.singlets, pg.nondebris, gln.zero, gln.one, gln.two, gln.three, file=paste(name,"_gates"_,Sys.Date(),".Rdata",sep=""))

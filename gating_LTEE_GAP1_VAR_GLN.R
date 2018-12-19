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
library(ggcyto)


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
pData(flowData)$name <- sampleNames(flowData)
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
twocopy <-4

##############################
#1. Generate gate for singlet cells####
#this gate is defined on the basis of the relationship between forward scatter height and area

ggcyto(gateData[zerocopy], aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) 
Agate <- locator(10, type='l', col='red')
gm.1 <- matrix(,length(Agate$x),2)
colnames(gm.1) <- c('FSC.H','FSC.A')
gm.1[,1] <- Agate$x
gm.1[,2] <- Agate$y
pg.singlets <- polygonGate(filterId="singlets",.gate=gm.1)


#Look at the gating on the controls
ggcyto(gateData[c(zerocopy,onecopy,twocopy)], aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.singlets)


#test that the singlet gate looks reasonable for All samples
ggcyto(gateData, aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.singlets)




##############################
#2. Generate Gate for debris based on forward scatter and side scatter. ####
#This needs to be done separately for each media condition.

plot(flowData[[3]], c('FSC.A','SSC.A'), xlim=c(0,4e6), xaxs = "i", yaxs = "i", ylim=c(0,1e6), smooth=F)

ggcyto(gateData[zerocopy], aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + geom_gate(pg.nondebris)
Bgate <- locator(10, type='l', col='red')
gm.2 <- matrix(,length(Bgate$x),2)
colnames(gm.2) <- c('FSC.A','SSC.A')
gm.2[,1] <- Bgate$x
gm.2[,2] <- Bgate$y
pg.nondebris <- polygonGate(filterId="nonDebris",.gate=gm.2)

#Look at the gating on the controls
ggcyto(gateData[c(zerocopy,onecopy,twocopy)], aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.nondebris)


#test that the singlet gate looks reasonable for All samples
ggcyto(gateData, aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.nondebris)


#############################
#3. FLUORESCENCE####

####Generate gates for 0, 1, 2, and 3+ copies####

##Plot the control sample that has non-fluorescing cells (0 copy)
ggcyto(gateData[zerocopy], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512)

Cgate <- locator(10, type='l', col='red')
gm.3 <- matrix(,length(Cgate$x),2)
colnames(gm.3) <- c('FSC.A','FL1.A')
gm.3[,1] <- Cgate$x
gm.3[,2] <- Cgate$y
fl1gate.0 <- polygonGate(filterId="zeroFL1",.gate=gm.3)

#Look at the gating on the controls
ggcyto(gateData[c(zerocopy,onecopy,twocopy)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(fl1gate.0)


##Draw a new gate for the one copy
ggcyto(gateData[onecopy], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512)
Dgate <- locator(10, type='l', col='blue')
gm.4 <- matrix(,length(Dgate$x),2)
colnames(gm.4) <- c('FSC.A','FL1.A')
gm.4[,4] <- Dgate$x
gm.4[,4] <- Dgate$y
fl1gate.1 <- polygonGate(filterId="oneCopyFL1",.gate=gm.4)


##Overlay and check the new gate
ggcyto(gateData[onecopy], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + geom_gate(fl1gate.0) + geom_gate(fl1gate.1)
<<<<<<< HEAD
=======

>>>>>>> cdacfaf1f9316168a96e6af34b24d7a7d9100fb7

##Plot the control sample that has 2 copies and draw a new gate for two copy
ggcyto(gateData[twocopy], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + geom_gate(fl1gate.0) + geom_gate(fl1gate.1)
Egate <- locator(10, type='l', col='green')
gm.5 <- matrix(,length(Egate$x),2)
colnames(gm.5) <- c('FSC.A','FL1.A')
gm.5[,4] <- Egate$x
gm.5[,4] <- Egate$y
fl1gate.2 <- polygonGate(filterId="twoCopyFL1",.gate=gm.5)

<<<<<<< HEAD
##Plot the control sample that has 2 copies and draw a new gate for two copy
ggcyto(gateData[twocopy], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + geom_gate(fl1gate.0) + geom_gate(fl1gate.1)
Egate <- locator(10, type='l', col='green')
gm.5 <- matrix(,length(Egate$x),2)
colnames(gm.5) <- c('FSC.A','FL1.A')
gm.5[,4] <- Egate$x
gm.5[,4] <- Egate$y
fl1gate.2 <- polygonGate(filterId="twoCopyFL1",.gate=gm.5)

##Plot the control sample that has 2 copies and draw a new gate for more then 2 copies
ggcyto(gateData[twocopy], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + geom_gate(fl1gate.0) + geom_gate(fl1gate.1) + geom_gate(fl1gate.2)
Fgate <- locator(10, type='l', col='purple')
gm.6 <- matrix(,length(Fgate$x),2)
colnames(gm.5) <- c('FSC.A','FL1.A')
gm.6[,4] <- Fgate$x
gm.6[,4] <- Fgate$y
fl1gate.3 <- polygonGate(filterId="2plusCopyFL1",.gate=gm.6)




#Look at the gating on the controls
ggcyto(gateData[c(zerocopy,onecopy,twocopy)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(fl1gate.0) + geom_gate(fl1gate.1) + geom_gate(fl1gate.2) + geom_gate(fl1gate.3)

=======
##Plot the control sample that has 2 copies and draw a new gate for more then 2 copies
ggcyto(gateData[twocopy], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + geom_gate(fl1gate.0) + geom_gate(fl1gate.1) + geom_gate(fl1gate.2)
Fgate <- locator(10, type='l', col='purple')
gm.6 <- matrix(,length(Fgate$x),2)
colnames(gm.5) <- c('FSC.A','FL1.A')
gm.6[,4] <- Fgate$x
gm.6[,4] <- Fgate$y
fl1gate.3 <- polygonGate(filterId="2plusCopyFL1",.gate=gm.6)




#Look at the gating on the controls
ggcyto(gateData[c(zerocopy,onecopy,twocopy)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(fl1gate.0) + geom_gate(fl1gate.1) + geom_gate(fl1gate.2) + geom_gate(fl1gate.3)

>>>>>>> cdacfaf1f9316168a96e6af34b24d7a7d9100fb7
##Check how the gates look on the sample
ggcyto(gateData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(fl1gate.0) + geom_gate(fl1gate.1) + geom_gate(fl1gate.2) + geom_gate(fl1gate.3)



#Save the gate information to an R data file
rm(list=c("gateData")) 
save(pg.singlets, pg.nondebris, fl1gate.0, fl1gate.1, fl1gate.2, fl1gate.3, file=paste(name,"_gates"_,Sys.Date(),".Rdata",sep=""))

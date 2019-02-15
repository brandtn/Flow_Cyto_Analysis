#######################
#gating_LTEE_GAP1_Variants.R       
#
#started: 01/07/2016
#modified: 1/25/2019
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
# CytoRSuite development version on GitHub
#devtools::install_github("DillonHammill/CytoRSuite", build_vignettes = TRUE)

# CytoRSuiteData development version on GitHub
#devtools::install_github("DillonHammill/CytoRSuiteData")


#Load libraries
library(flowCore)
library(flowViz)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(ggforce)
#library(CytoRSuite)
devtools::install_github("DillonHammill/CytoRSuite", build_vignettes = TRUE)
library(CytoRSuiteData)




#Read in the data
#Set working directory to the folder in which you have stored your .fcs files
#Read in all the fcs files in the directory.

#working directory
dir = '.'

#file location
path.data = "/Users/Brandt/Google Drive/MiniStatRun_10_2018/"
#path.data = "/Users/nathanbrandt/Google Drive/MiniStatRun_10_2018/"


#set name of run to create gates for
list.folders <- c("LTEE_mCitrine_GAP1_Variants_T00", 
                  "LTEE_mCitrine_GAP1_Variants_T06", 
                  "LTEE_mCitrine_GAP1_Variants_T07", 
                  "LTEE_mCitrine_GAP1_Variants_T08.3", 
                  #"LTEE_mCitrine_GAP1_Variants_T11.1", 
                  "LTEE_mCitrine_GAP1_Variants_T11.2", 
                  #"LTEE_mCitrine_GAP1_Variants_T13.1", 
                  "LTEE_mCitrine_GAP1_Variants_T13.2", 
                  "LTEE_mCitrine_GAP1_Variants_T14", 
                  "LTEE_mCitrine_GAP1_Variants_T15", 
                  "LTEE_mCitrine_GAP1_Variants_T18", 
                  "LTEE_mCitrine_GAP1_Variants_T22", 
                  "LTEE_mCitrine_GAP1_Variants_T25", 
                  "LTEE_mCitrine_GAP1_Variants_T27", 
                  "LTEE_mCitrine_GAP1_Variants_T29", 
                  "LTEE_mCitrine_GAP1_Variants_T34")

name <- list.folders[1]

#load sample sheet
sample.sheet <- read.csv(paste(path.data,"samplesheet_",name,".csv", sep=""))

#read in fcs files in order presented in sample sheet (based on well identifier)
files <- paste(path.data,name,"/",sort(factor(list.files(paste(path.data,name,"/", sep=""),full.names=FALSE), levels = paste(sample.sheet$Well,".fcs",sep="" ), ordered=TRUE)),sep="")
flowData <- read.ncdfFlowSet(files=files, pattern=".fcs", alter.names = TRUE)

sample.ind <- which(paste(sample.sheet$Well,".fcs", sep="") %in% sampleNames(flowData))
sample.sheet <- sample.sheet[sample.ind,]
sample.sheet <- sample.sheet[order(sample.sheet$Well),]

#rename sample name of flow set to make it easier to identify
sampleNames(flowData) <- paste(gsub(" ","_",sample.sheet$Strain),"_",sub(" ","_",sample.sheet$Well), sep="")

#set copy number controls
zerocopy <- 1
onecopy <- 3
twocopy <- 4


# Add flowSet to GatingSet
gs <- GatingSet(flowData)

# Gate Cells
gate_draw(gs, 
          parent = "root",
          alias = "Cells",
          channels = c("FSC-A","SSC-A"),
          type = "polygon",
          gatingTemplate = "Compensation-gatingTemplate.csv")

# Gate Single Cells
gate_draw(gs, 
          parent = "Cells",
          alias = "Single Cells",
          channels = c("FSC-A","FSC-H"),
          type = "polygon",
          gatingTemplate = "Compensation-gatingTemplate.csv")





##############################
#1. Generate gate for singlet cells####
#this gate is defined on the basis of the relationship between forward scatter height and area
#******Please note you may need to adjust the x and y plot limits to properly visualize your data

plot(flowData[[zerocopy]], c('FSC.H','FSC.A'), xlim=c(0,3e6), ylim=c(0,3e6),smooth=T)


singlet.gate <- locator(100, type='l', col='red')
gm.1 <- matrix(,length(singlet.gate$x),2)
colnames(gm.1) <- c('FSC.H','FSC.A')
gm.1[,1] <- singlet.gate$x
gm.1[,2] <- singlet.gate$y
pg.singlets <- polygonGate(filterId="singlets",.gate=gm.1)

ggcyto(flowData[zerocopy], aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + geom_gate(pg.singlets)

#Look at the gating on the controls
ggcyto(flowData[c(zerocopy,onecopy,twocopy)], aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.singlets) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)

#test that the singlet gate looks reasonable for All samples
for(i in 1:round(length(flowData)/4,0)){
  plot <- ggcyto(flowData, aes(x = `FSC.H`, y =  `FSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(pg.singlets) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = i)
  print(plot)
}

#Filter out the doublets
flowData.singlets <- Subset(flowData,pg.singlets)

##############################
#2. Generate Gate for debris based on forward scatter and side scatter. ####
#This needs to be done separately for each media condition.
#******Please note you may need to adjust the x and y plot limits to properly visualize your data


plot(flowData.singlets[[zerocopy]], c('FSC.A','SSC.A'), xlim=c(0,3e6), ylim=c(0,1e6),smooth=T)

debris.gate <- locator(100, type='l', col='red')
gm.2 <- matrix(,length(debris.gate$x),2)
colnames(gm.2) <- c('FSC.A','SSC.A')
gm.2[,1] <- debris.gate$x
gm.2[,2] <- debris.gate$y
pg.nondebris <- polygonGate(filterId="nonDebris",.gate=gm.2)

#Look at the gating on the controls
ggcyto(flowData.singlets[c(zerocopy,onecopy,twocopy)], aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e6) + geom_gate(pg.nondebris) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)


#test that the singlet gate looks reasonable for All samples
for(i in 1:round(length(flowData)/4,0)){
  plot <- ggcyto(flowData.singlets, aes(x = `FSC.A`, y =  `SSC.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,1e6) + geom_gate(pg.nondebris) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = i)
  print(plot)
}

#Filter out the debris
flowData.nondebris <- Subset(flowData.singlets,pg.nondebris)

#############################
#3. FLUORESCENCE####

####Generate gates for 0, 1, 2, and 3+ copies####
#******Please note you may need to adjust the x and y plot limits to properly visualize your data

##Plot the control sample that has non-fluorescing cells (0 copy)
plot(flowData.nondebris[[zerocopy]], c('FSC.A','FL1.A'), xlim=c(0,3e6), ylim=c(0,5e4),smooth=T)


#ggcyto(flowData[zerocopy], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512)
zero.gate <- locator(100, type='l', col='red')
gm.3 <- matrix(,length(zero.gate$x),2)
colnames(gm.3) <- c('FSC.A','FL1.A')
gm.3[,1] <- zero.gate$x
gm.3[,2] <- zero.gate$y
fl1gate.0 <- polygonGate(filterId="zeroFL1",.gate=gm.3)

#Look at the gating on the controls
ggcyto(flowData.nondebris[c(zerocopy,onecopy,twocopy)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(fl1gate.0) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)

for(i in 1:round(length(flowData)/4,0)){
  plot <- ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,5e4) + geom_gate(fl1gate.0) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = i)
  print(plot)
}



##Draw a new gate for the one copy include the gate for zero copies

plot(flowData.nondebris[[onecopy]], c('FSC.A','FL1.A'), xlim=c(0,3e6), ylim=c(0,5e5),smooth=T)
polygon(zero.gate)

one.gate <- locator(100, type='l', col='blue')
gm.4 <- matrix(,length(one.gate$x),2)
colnames(gm.4) <- c('FSC.A','FL1.A')
gm.4[,1] <- one.gate$x
gm.4[,2] <- one.gate$y
fl1gate.1 <- polygonGate(filterId="oneCopyFL1",.gate=gm.4)

##Overlay and check the new gate
ggcyto(flowData.nondebris[onecopy], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) +  geom_gate(fl1gate.0) + geom_gate(fl1gate.1)


ggcyto(flowData.nondebris[c(onecopy,twocopy)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) +  geom_gate(fl1gate.0) + geom_gate(fl1gate.1)

##Plot the control sample that has 2 copies along with the one and zero copy gates and draw a new gate for two copy
plot(flowData.nondebris[[twocopy]], c('FSC.A','FL1.A'), xlim=c(0,2e6), ylim=c(0,5e5),smooth=T)
polygon(zero.gate)
polygon(one.gate)

two.gate <- locator(100, type='l', col='green')
gm.5 <- matrix(,length(two.gate$x),2)
colnames(gm.5) <- c('FSC.A','FL1.A')
gm.5[,1] <- two.gate$x
gm.5[,2] <- two.gate$y
fl1gate.2 <- polygonGate(filterId="twoCopyFL1",.gate=gm.5)

##Overlay and check the new gate
ggcyto(flowData.nondebris[twocopy], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + geom_gate(fl1gate.0) + geom_gate(fl1gate.1) + geom_gate(fl1gate.2)

##Plot the control sample that has 2 copies along with the two, one, and zero copy gates and draw a new gate for more then 2 copies
plot(flowData.nondebris[[twocopy]], c('FSC.A','FL1.A'), xlim=c(0,3e6), ylim=c(0,1e6), smooth=T)
polygon(zero.gate)
polygon(one.gate)
polygon(two.gate)

three.gate <- locator(10, type='l', col='purple')
gm.6 <- matrix(,length(three.gate$x),2)
colnames(gm.6) <- c('FSC.A','FL1.A')
gm.6[,1] <- three.gate$x
gm.6[,2] <- three.gate$y
fl1gate.3 <- polygonGate(filterId="2plusCopyFL1",.gate=gm.6)

#Look at the gating on the controls
ggcyto(flowData.nondebris[c(zerocopy,onecopy,twocopy)], aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(fl1gate.0) + geom_gate(fl1gate.1) + geom_gate(fl1gate.2) + geom_gate(fl1gate.3) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)


##Check how the gates look on all the samples
for(i in 1:round(length(flowData)/4,0)){
  plot <- ggcyto(flowData, aes(x = `FSC.A`, y =  `FL1.A`)) + geom_hex(bins = 512) + xlim(0,3e6) + ylim(0,3e6) + geom_gate(fl1gate.0) + geom_gate(fl1gate.1) + geom_gate(fl1gate.2) + geom_gate(fl1gate.3) + facet_wrap_paginate(~name, ncol = 2, nrow = 2, page = 1)
  print(plot)
}

#Save the gate information to an R data file
rm(list=c("flowData")) 
save(pg.singlets, pg.nondebris, fl1gate.0, fl1gate.1, fl1gate.2, fl1gate.3, file=paste(name,"_gates_",Sys.Date(),".Rdata",sep=""))

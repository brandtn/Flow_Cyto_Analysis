###############################
###CODE TO GENERATE FIGURE 2###       ### One slight fix you could do, if editing this figure- is make sure the colors are the same between Figure 2B and 2C (GLN)
###############################

###To get libraries from BioConductor.  This only needs to be run the first time
source("http://bioconductor.org/biocLite.R")
biocLite("flowViz")
biocLite("flowCore")

#Load libraries
library(flowCore)
library(flowViz)
library(ggplot2)
library(tidyverse)
library(reshape2) ##**Needed for melt() function

biocLite("devtools")
library(devtools)
install_github("clauswilke/ggjoy")
library(ggjoy)
biocLite("ggthemes")
library(ggthemes)
biocLite("hrbrthemes")
library(hrbrthemes)

#######################
#######FIGURE 2A#######   GGJOY PLOT, EXAMPLE OF ALL THREE MEDIA CONDITIONS, FLUOR DISTRIBUTIONS OVER TIME
#######################

#If desired, load previously generated gates to filter cell debris
#Load from any folder located in Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/
load("gates_05252017_MIDTPT.Rdata")

##
#To extract the data as a timecourse, I extract all the data for a given population sample 
#One population sample at a time

#Define the directory, or directories, containing your .fcs files using absolute path names 
dir1 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 01/"
dir2 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 02/"
dir3 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 03/"
dir4 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 04/"
dir5 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 05/"
dir6 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 06/"
dir7 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 07/"
dir8 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 08/"
dir9 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 09/"
dir10 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 10/"
dir11 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 11/"
dir12 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 12/"
dir13 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 13/"
dir14 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 14/"
dir15 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 15/"
dir16 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 16/"
dir17 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 17/"
dir18 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 18/"
dir19 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 19/"
dir20 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 20/"
dir21 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 21/"
dir22 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 22/"
dir23 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 23/"
dir24 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 24/"
dir25 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 25/"
dir26 <- "/Users/Steff/Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 26/"


#Read in all the fcs files in the directory, with alter.names changing "-" to "."

##Get the controls from your desired timepoint first.. recommend dir3:
flowData.1a <- read.flowSet(path = dir17, pattern="C10.fcs", alter.names=TRUE)
flowData.1b <- read.flowSet(path = dir17, pattern="C10.fcs", alter.names=TRUE)

##Then one population from all timepoints
#Change the pattern to reflect which population sample you are extracting
#Use Edit --> Find --> Replace all to easily change the filename
flowData.1 <- read.flowSet(path = dir1, pattern="C10.fcs", alter.names=TRUE)
flowData.2 <- read.flowSet(path = dir2, pattern="C10.fcs", alter.names=TRUE)
flowData.3 <- read.flowSet(path = dir3, pattern="C10.fcs", alter.names=TRUE)
flowData.4 <- read.flowSet(path = dir4, pattern="C10.fcs", alter.names=TRUE)
flowData.5 <- read.flowSet(path = dir5, pattern="C10.fcs", alter.names=TRUE)
flowData.6 <- read.flowSet(path = dir6, pattern="C10.fcs", alter.names=TRUE)
flowData.7 <- read.flowSet(path = dir7, pattern="C10.fcs", alter.names=TRUE)
flowData.8 <- read.flowSet(path = dir8, pattern="C10.fcs", alter.names=TRUE)
flowData.9 <- read.flowSet(path = dir9, pattern="C10.fcs", alter.names=TRUE)
flowData.10 <- read.flowSet(path = dir10, pattern="C10.fcs", alter.names=TRUE)
flowData.11 <- read.flowSet(path = dir11, pattern="C10.fcs", alter.names=TRUE)
flowData.12 <- read.flowSet(path = dir12, pattern="C10.fcs", alter.names=TRUE)
flowData.13 <- read.flowSet(path = dir13, pattern="C10.fcs", alter.names=TRUE)
flowData.14 <- read.flowSet(path = dir14, pattern="C10.fcs", alter.names=TRUE)
flowData.15 <- read.flowSet(path = dir15, pattern="C10.fcs", alter.names=TRUE)
flowData.16 <- read.flowSet(path = dir16, pattern="C10.fcs", alter.names=TRUE)
flowData.17 <- read.flowSet(path = dir17, pattern="C10.fcs", alter.names=TRUE)
flowData.18 <- read.flowSet(path = dir18, pattern="C10.fcs", alter.names=TRUE)
flowData.19 <- read.flowSet(path = dir19, pattern="C10.fcs", alter.names=TRUE)
flowData.20 <- read.flowSet(path = dir20, pattern="C10.fcs", alter.names=TRUE)
flowData.21 <- read.flowSet(path = dir21, pattern="C10.fcs", alter.names=TRUE)
flowData.22 <- read.flowSet(path = dir22, pattern="C10.fcs", alter.names=TRUE)
flowData.23 <- read.flowSet(path = dir23, pattern="C10.fcs", alter.names=TRUE)
flowData.24 <- read.flowSet(path = dir24, pattern="C10.fcs", alter.names=TRUE)
flowData.25 <- read.flowSet(path = dir25, pattern="C10.fcs", alter.names=TRUE)
flowData.26 <- read.flowSet(path = dir26, pattern="C10.fcs", alter.names=TRUE)


#Read in a sample sheet that has relevant information
#I simply used a sheet that had my generations listed.., I keep it in Directory 1
sample.sheet <- read.delim(paste(dir1, "SampleSheetA.txt", sep="/"))

#Change names of samples to those specified in the sample sheet

sampleNames(flowData.1) <- paste(sample.sheet[1,1], sep=" ")
sampleNames(flowData.2) <- paste(sample.sheet[2,1], sep=" ")
sampleNames(flowData.3) <- paste(sample.sheet[3,1], sep=" ")
sampleNames(flowData.4) <- paste(sample.sheet[4,1], sep=" ")
sampleNames(flowData.5) <- paste(sample.sheet[5,1], sep=" ")
sampleNames(flowData.6) <- paste(sample.sheet[6,1], sep=" ")
sampleNames(flowData.7) <- paste(sample.sheet[7,1], sep=" ")
sampleNames(flowData.8) <- paste(sample.sheet[8,1], sep=" ")
sampleNames(flowData.9) <- paste(sample.sheet[9,1], sep=" ")
sampleNames(flowData.10) <- paste(sample.sheet[10,1], sep=" ")
sampleNames(flowData.11) <- paste(sample.sheet[11,1], sep=" ")
sampleNames(flowData.12) <- paste(sample.sheet[12,1], sep=" ")
sampleNames(flowData.13) <- paste(sample.sheet[13,1], sep=" ")
sampleNames(flowData.14) <- paste(sample.sheet[14,1], sep=" ")
sampleNames(flowData.15) <- paste(sample.sheet[15,1], sep=" ")
sampleNames(flowData.16) <- paste(sample.sheet[16,1], sep=" ")
sampleNames(flowData.17) <- paste(sample.sheet[17,1], sep=" ")
sampleNames(flowData.18) <- paste(sample.sheet[18,1], sep=" ")
sampleNames(flowData.19) <- paste(sample.sheet[19,1], sep=" ")
sampleNames(flowData.20) <- paste(sample.sheet[20,1], sep=" ")
sampleNames(flowData.21) <- paste(sample.sheet[21,1], sep=" ")
sampleNames(flowData.22) <- paste(sample.sheet[22,1], sep=" ")
sampleNames(flowData.23) <- paste(sample.sheet[23,1], sep=" ")
sampleNames(flowData.24) <- paste(sample.sheet[24,1], sep=" ")
sampleNames(flowData.25) <- paste(sample.sheet[25,1], sep=" ")
sampleNames(flowData.26) <- paste(sample.sheet[26,1], sep=" ")

#Put all the data from one population's timecourse in a giant flow data set
flowData <- rbind2(flowData.1, flowData.2)
flowData <- rbind2(flowData, flowData.3)
flowData <- rbind2(flowData, flowData.4)
flowData <- rbind2(flowData, flowData.5)
flowData <- rbind2(flowData, flowData.6)
flowData <- rbind2(flowData, flowData.7)
flowData <- rbind2(flowData, flowData.8)
flowData <- rbind2(flowData, flowData.9)
flowData <- rbind2(flowData, flowData.10)
flowData <- rbind2(flowData, flowData.11)
flowData <- rbind2(flowData, flowData.12)
flowData <- rbind2(flowData, flowData.13)
flowData <- rbind2(flowData, flowData.14)
flowData <- rbind2(flowData, flowData.15)
flowData <- rbind2(flowData, flowData.16)
flowData <- rbind2(flowData, flowData.17)
flowData <- rbind2(flowData, flowData.18)
flowData <- rbind2(flowData, flowData.19)
flowData <- rbind2(flowData, flowData.20)
flowData <- rbind2(flowData, flowData.21)
flowData <- rbind2(flowData, flowData.22)
flowData <- rbind2(flowData, flowData.23)
flowData <- rbind2(flowData, flowData.24)
flowData <- rbind2(flowData, flowData.25)
flowData <- rbind2(flowData, flowData.26)

##Subset the data by applying sequential gates if desired##
#apply doublet gate
flowData.singlets <- Subset(flowData, pg.singlets) 
#apply debris gate
filteredGlnData <- Subset(flowData.singlets, pg.nondebris.gln) 


#This is a cool technique Darach taught me. Use this loop to generate a series of lists.
alldata <- list()
for (i in 1:length(filteredGlnData)){
  
  temp <- exprs(filteredGlnData[[i]]) #exprs() extracts a matrix of the values from the flowframe
  alldata[[(sampleNames(filteredGlnData[i]))]] <- data.frame(gens=sampleNames(filteredGlnData[i]),FSC=temp[,1],FL1=temp[,3],FL1_FSC=(temp[,3]/temp[,1]), logFL1=log(temp[,3]), logFL1_FSC=log(temp[,3]/temp[,1]))
  
}
#Put all the lists from the loop into one big data frame
big_data_GLN = do.call(rbind,alldata)



#######REPEAT FOR UREA LIMITATOIN######
##Subset the data by applying sequential gates if desired##
#apply doublet gate
flowData.singlets <- Subset(flowData, pg.singlets) 
#apply debris gate
filteredUrData <- Subset(flowData.singlets, pg.nondebris.ur) 


#This is a cool technique Darach taught me. Use this loop to generate a series of lists.
alldata <- list()
for (i in 1:length(filteredUrData)){
  
  temp <- exprs(filteredUrData[[i]]) #exprs() extracts a matrix of the values from the flowframe
  alldata[[(sampleNames(filteredUrData[i]))]] <- data.frame(gens=sampleNames(filteredUrData[i]),FSC=temp[,1],FL1=temp[,3],FL1_FSC=(temp[,3]/temp[,1]), logFL1=log(temp[,3]), logFL1_FSC=log(temp[,3]/temp[,1]))
  
}
#Put all the lists from the loop into one big data frame
big_data_UR = do.call(rbind,alldata)







#######REPEAT FOR GLUCOSE LIMITATOIN######
##Subset the data by applying sequential gates if desired##
#apply doublet gate
flowData.singlets <- Subset(flowData, pg.singletsGLU) 
#apply debris gate
filteredGLUData <- Subset(flowData.singlets, pg.nondebris.glu) 


#This is a cool technique Darach taught me. Use this loop to generate a series of lists.
alldata <- list()
for (i in 1:length(filteredGLUData)){
  
  temp <- exprs(filteredGLUData[[i]]) #exprs() extracts a matrix of the values from the flowframe
  alldata[[(sampleNames(filteredGLUData[i]))]] <- data.frame(gens=sampleNames(filteredGLUData[i]),FSC=temp[,1],FL1=temp[,3],FL1_FSC=(temp[,3]/temp[,1]), logFL1=log(temp[,3]), logFL1_FSC=log(temp[,3]/temp[,1]))
  
}
#Put all the lists from the loop into one big data frame
big_data_GLU = do.call(rbind,alldata)




#Make sure data is in the correct format.. ggjoy needs factors
big_data_GLN$gens<-factor(big_data_GLN$gens,levels=unique(big_data_GLN$gens))
big_data_UR$gens<-factor(big_data_UR$gens,levels=unique(big_data_UR$gens))
big_data_GLU$gens<-factor(big_data_GLU$gens,levels=unique(big_data_GLU$gens))

ggplot(big_data_GLN, aes(x = `FL1_FSC`, y = `gens`, fill = ..x.., height=..density..)) +
  geom_joy_gradient(scale = 3, rel_min_height = 0.03) +
  scale_x_continuous("Normalized Fluorescence", limits=c(0.05,1), expand = c(0.01, 0)) +
  scale_y_discrete("Generations", expand = c(0.01, 0)) +
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") +
  theme_joy(font_size = 11, grid = TRUE) + theme(legend.position = 'none', legend.title = element_blank(), axis.title.y = element_blank())

ggplot(big_data_UR, aes(x = `FL1_FSC`, y = `gens`, fill = ..x.., height=..density..)) +
  geom_joy_gradient(scale = 3, rel_min_height = 0.03) +
  scale_x_continuous("Normalized Fluorescence (a.u.)", limits=c(0.0,.4), expand = c(0.01, 0)) +
  scale_y_discrete("Generations", expand = c(0.01, 0)) +
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") +
  theme_joy(font_size = 11, grid = TRUE) + theme(legend.position = 'none', legend.title = element_blank(), axis.title.y = element_blank())

ggplot(big_data_GLU, aes(x = `FL1_FSC`, y = `gens`, fill = ..x.., height=..density..)) +
  geom_joy_gradient(scale = 3, rel_min_height = 0.03) +
  scale_x_continuous("Normalized Fluorescence (a.u.)", limits=c(0.05,1), expand = c(0.01, 0)) +
  scale_y_discrete("Generations", expand = c(0.01, 0)) +
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") +
  theme_joy(font_size = 11, grid = TRUE) + theme(legend.position = 'none', legend.title = element_blank(), axis.title.y = element_blank())



#Make sure data is in the correct format.. ggjoy needs factors
big_data_GLN_CONTROLS$gens<-factor(big_data_GLN_CONTROLS$gens,levels=unique(big_data_GLN_CONTROLS$gens))

##Make the namez list a factor and re-order by number
alldata <- data.frame(SingleCellDistributions)
alldata$gens<-factor(alldata$gens,levels=unique(alldata$gens))
gln_ctrl01 <- alldata[which(alldata$namez == 'gln_ctrl01'),]
gln_ctrl02 <- alldata[which(alldata$namez == 'gln_ctrl02'),]
gln_gap08 <- alldata[which(alldata$namez == 'gln_gap08'),]


ggplot(gln_08, aes(x = `FL1_FSC`, y = `gens`, fill = ..x.., height=..density..)) +
  geom_joy_gradient(scale = 2.25, rel_min_height = 0.03) +
  scale_x_continuous("Normalized Fluorescence", limits=c(0.05,1.25), expand = c(0.01, 0), breaks = c(0.1, 0.5, 1)) +
  scale_y_discrete("Generations", expand = c(0.01, 0)) +
  facet_grid(. ~ namez) + 
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") +
  theme_joy(font_size = 12, grid = TRUE) + theme(legend.position = 'none', legend.title = element_blank())

gln_ctr01 <- ggplot(gln_ctrl01, aes(x = FL1_FSC, y = gens, fill = ..x.., height=..density..)) +
  geom_joy_gradient(scale = 2.25, rel_min_height = 0.01) +
  scale_x_continuous(" ", limits=c(0.05,1.0), expand = c(0.01, 0), breaks = c(0.1, 0.55, 1.0)) +
  scale_y_discrete("Generations", expand = c(0.01, 0)) +
  facet_grid(. ~ namez) + 
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") +
  theme_joy(font_size = 12, grid = TRUE) + theme(legend.position = 'none', legend.title = element_blank())
gln_ctr01

gln_ctr02 <- ggplot(gln_ctrl02, aes(x = FL1_FSC, y = gens, fill = ..x.., height=..density..)) +
  geom_joy_gradient(scale = 2.25, rel_min_height = 0.01) +
  scale_x_continuous("Normalized Fluorescence", limits=c(0.05,1.0), expand = c(0.01, 0), breaks = c(0.1, 0.55, 1.0)) +
  scale_y_discrete(" ", expand = c(0.01, 0)) +
  facet_grid(. ~ namez) + 
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") +
  theme_joy(font_size = 12, grid = TRUE) + theme(legend.position = 'none', legend.title = element_blank())
gln_ctr02

gln_08 <- ggplot(gln_gap08, aes(x = FL1_FSC, y = gens, fill = ..x.., height=..density..)) +
  geom_joy_gradient(scale = 2.25, rel_min_height = 0.01) +
  scale_x_continuous(" ", limits=c(0.05,1.0), expand = c(0.01, 0), breaks = c(0.1, 0.55, 1.0)) +
  scale_y_discrete(" ", expand = c(0.01, 0)) +
  facet_grid(. ~ namez) + 
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") +
  theme_joy(font_size = 12, grid = TRUE) + theme(legend.position = 'none', legend.title = element_blank())
gln_08

plot_grid(gln_ctr01, gln_ctr02, gln_08, labels = c("A)"), nrow=1, align = "h")

##Write all the data to a datatable -- this table is currently living in Sample 25 directory
write.table(big_data_GLN_CONTROLS, file="SingleCellDistributionsCtrls.txt", row.names=FALSE, quote=F, sep="\t")




#######################
#######FIGURE 2B#######   LINE GRAPH, ALL CONDITIONS, MEDIAN FLUOR
#######################

#Data located here: SLL08/All Flow Cyto/FlowAnalysisSumm_MidTpt Gates_Fix
#You can download to your computer as an XLS, then import it to R
alldata <- data.frame(FlowAnalysisSumm_FINAL_2_)

alldata1 <- data.frame(FlowAnalysisSumm_FINAL_Sheet1_6_)
alldata1 <- alldata1[which(alldata1$Cells.counted > 70000),]

library(plotly)
gln_ctrl01 <- alldata1[which(alldata1$Unique_ID == 'gln_ctrl01'),]
gln_ctrl01IQR <- IQR(gln_ctrl01$normalizedGFP_median)
gln_ctrl01med <- median(gln_ctrl01$normalizedGFP_median)

gln_ctrl02 <- alldata1[which(alldata1$Unique_ID == 'gln_ctrl02'),]
gln_ctrl02IQR <- IQR(gln_ctrl02$normalizedGFP_median)
gln_ctrl02med <- median(gln_ctrl02$normalizedGFP_median)


glnplot <- ggplot(alldata1, aes(x=Generation,y=(normalizedGFP_median), colour = Unique_ID)) +
  geom_line(aes(linetype = Type), show.legend = FALSE, size=1.2) +
  scale_y_continuous(' ', expand = c(0, 0), limits=c(0,0.7)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0), limits=c(20,260), " ") +
  theme_classic() +
  scale_fill_manual(values=c("grey", "grey", "green")) + 
  scale_color_manual(name="Nutrient",values=c('#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b', '#08306b','#2171b5','#6baed6', 'grey', 'grey')) +
  scale_linetype_manual(values=c('twodash', 'twodash', "solid")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme(legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
glnplot <- glnplot %+% subset(alldata1, Nutrient %in% c("Glutamine"))
glnplot <- glnplot + aes(group=rev(Unique_ID))

urplot <- ggplot(alldata1, aes(x=Generation,y=(normalizedGFP_median), colour = Unique_ID)) +
  #geom_point(size=1, show.legend = FALSE) +
  geom_line(aes(linetype = Type), show.legend = FALSE, size=1.2) +
  scale_y_continuous('Normalized median fluorescence (a.u.)', expand = c(0, 0), limits=c(0,0.7)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0), limits=c(20,260), " ") +
  theme_classic() +
  scale_color_manual(name="Nutrient",values=c('#662506','#cc4c02','#fed976','#feb24c','#e31a1c','#fc4e2a', '#fd8d3c','#bd0026','#800026', 'grey', 'grey')) +
  scale_linetype_manual(values=c('twodash', 'twodash', "solid")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme(legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
urplot <- urplot %+% subset(alldata1, Nutrient %in% c("Urea"))
urplot <- urplot + aes(group=rev(Unique_ID))

glUplot <- ggplot(alldata1, aes(x=Generation,y=(normalizedGFP_median), colour = Unique_ID)) +
  #geom_point(size=1, show.legend = FALSE) +
  geom_line(aes(linetype = Type), show.legend = FALSE, size=1.2) +
  scale_y_continuous(' ', expand = c(0, 0), limits=c(0,0.7)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0), limits=c(20,260), "Generations") +
  theme_classic() +
  scale_color_manual(name="Nutrient",values=c('#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a', 'gray', 'gray')) +
  scale_linetype_manual(values=c('twodash', 'twodash', "solid")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme(legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
glUplot <- glUplot %+% subset(alldata1, Nutrient %in% c("Glucose"))
glUplot <- glUplot + aes(group=rev(Unique_ID))

plot_grid(glUplot, urplot, glnplot, labels = c("B)"), ncol=1, align = "v")


#######################
#######FIGURE 2C#######   BARPLOT, EXAMPLE OF ONE POPULATION, PROPORTION 0, 1, 2, 3+ COPIES 
#######################

##For this plot, we actually need to change the format of the data. You can do that using melt() .. it's an awesome command!!
new <- melt(alldata1, id.vars = c("Generation", "Unique_ID", "Population", "Nutrient"), measure.vars = c("PopProp_0copy","PopProp_1copy","PopProp_CNV"))
new <- data.frame(Generation=as.character(new[,1]), Unique_ID=as.character(new[,2]), Population=as.character(new[,3]), Nutrient=as.character(new[,4]), variable=as.character(new[,5]), value=as.numeric(new[,6])) 

new$Generation <- factor(new$Generation, levels = c(25, 33, 41, 54, 62, 70, 79, 87, 95, 103, 116, 124, 132, 145, 153, 161, 174, 182, 190, 211, 219, 232, 244, 257, 267))
new$variable <- factor(new$variable, levels = c("PopProp_CNV", "PopProp_1copy", "PopProp_0copy"))


theplot <- ggplot(new, aes(Generation, value, fill=variable)) + 
  theme_classic() +
  theme(legend.position = "top", legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  geom_bar(colour = "black", position = "fill", stat = "identity", width = 0.9) +
  scale_y_continuous("Proportion", expand = c(0, 0), limits=c(0,1)) + 
  scale_fill_manual(values = c('#006d2c', '#a1d99b','#e5f5e0'))

theplot %+% subset(new, Population %in% c("C10"))



#######################
#######FIGURE 2D#######   LINE GRAPH, ALL GLN POPS, PROPORTION OF CNV CELLS 
#######################

myplot <- ggplot(alldata1, aes(x=Generation,y=(PopProp_CNV), colour=Unique_ID)) +
  geom_point(size=1.1, show.legend = FALSE) +
  geom_line(size=1, show.legend = FALSE) +
  scale_y_continuous(expand = c(0, 0), 'Proportion of cells with CNV', limits=c(0,1)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), "Generations") +
  theme_classic() +
  scale_color_manual(values = c('#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b', '#08306b','#2171b5','#6baed6')) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme(legend.position = 'bottom',legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))

myplot %+% subset(alldata1, Nutrient %in% c("Glutamine") & Type %in% c("Experimental"))



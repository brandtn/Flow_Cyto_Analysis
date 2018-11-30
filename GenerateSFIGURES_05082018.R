######################## FINAL SCRIPTS USED TO GENERATE SUPPLEMENTARY PLOTS ##########################
###########
##Figure S1
###########

#Data to load
data1 <- data.frame(SL_clone_fitness_csv_SL_clone_fitness_csv_2_)

##Plot
base <- ggplot(data1, aes(x=clone, y=fitness_coef+1)) +
  geom_bar(stat = 'identity', width = 0.9) + 
  facet_grid(~Replicate) + 
  geom_errorbar(aes(ymin=1+ci2.5, ymax=1+ci97.5), width=.3) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 11, angle=75, vjust=0.6), legend.position = 'none', legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=11), axis.title.x = element_text(face="bold", size=12)) + 
  scale_fill_manual(values = c('#74c476','#00441b')) +
  scale_x_discrete('Strain') +
  geom_hline(yintercept=1.0, linetype="dotted", color="black", size=.8) +
  scale_y_continuous('Relative Fitness', expand = c(0, 0), limits = c(0,1.3))
base %+% subset(data1, Type %in% c("Control") & Replicate %in% c("Replicate 2", "Replicate 3"))

ggsave(filename="Figure S1 Test.pdf", scale = 1.5, width=3, height=2, units="in", dpi=300)
###########
##Figure S2
###########

data2 <- read_delim("Google Drive/Gresham Lab_Steff/SLL08_EE06_Mini02/ALL FLOW CYTO/Sample 25/SingleCellDistributions.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
data2 <- data.frame(data2)
data2$gens<-factor(data2$gens,levels=unique(data2$gens))

ggplot(data2, aes(x = FL1_FSC, y = gens, fill = ..x.., height=..density..)) +
  geom_joy_gradient(scale = 1.25, rel_min_height = 0.02) +
  scale_x_continuous("Normalized Fluorescence (a.u.)", limits=c(0.05,1.4), expand = c(0.01, 0), breaks = c(0.1, 0.7, 1.2)) +
  scale_y_discrete("Generations", expand = c(0.01, 0)) +
  facet_grid(. ~ namez) + 
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") +
  theme_bw() + 
  theme(text=element_text(size=11),legend.position = 'none', legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=10), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=10))

ggsave(filename="Figure 2 Test.pdf", scale = 1.25, width=7.5, height=6, units="in", dpi=300)

###########
##Figure S3
###########

data3 <- data.frame(FlowAnalysisSumm_FINAL_Sheet1_6_)

glnplot <- ggplot(data3, aes(x=Generation,y=(FL1_median), colour = Unique_ID)) +
  geom_line(aes(linetype = Type), show.legend = FALSE, size=1.2) +
  scale_y_continuous(' ', expand = c(0, 0), limits=c(0,3e5)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0), limits=c(20,260), " ") +
  theme_classic() +
  annotate("text", x = 50, y = 3e4, label = 'Glutamine', parse = TRUE) + 
  scale_fill_manual(values=c("grey", "grey", "green")) + 
  scale_color_manual(name="Nutrient",values=c('#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b', '#08306b','#2171b5','#6baed6', 'grey', 'grey')) +
  scale_linetype_manual(values=c('twodash', 'twodash', "solid")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme(text=element_text(size=12), legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
glnplot <- glnplot %+% subset(data3, Nutrient %in% c("Glutamine"))
glnplot <- glnplot + aes(group=rev(Unique_ID))
glnplot

urplot <- ggplot(data3, aes(x=Generation,y=(FL1_median), colour = Unique_ID)) +
  #geom_point(size=1, show.legend = FALSE) +
  geom_line(aes(linetype = Type), show.legend = FALSE, size=1.2) +
  scale_y_continuous('Median fluorescence (a.u.)', expand = c(0, 0), limits=c(0,3e5)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0), limits=c(20,260), " ") +
  theme_classic() +
  annotate("text", x = 50, y = 3e4, label = 'Urea', parse = TRUE) +
  scale_color_manual(name="Nutrient",values=c('#662506','#cc4c02','#fed976','#feb24c','#e31a1c','#fc4e2a', '#fd8d3c','#bd0026','#800026', 'grey', 'grey')) +
  scale_linetype_manual(values=c('twodash', 'twodash', "solid")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme(text=element_text(size=12), legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
urplot <- urplot %+% subset(data3, Nutrient %in% c("Urea"))
urplot <- urplot + aes(group=rev(Unique_ID))
urplot

glUplot <- ggplot(data3, aes(x=Generation,y=(FL1_median), colour = Unique_ID)) +
  #geom_point(size=1, show.legend = FALSE) +
  geom_line(aes(linetype = Type), show.legend = FALSE, size=1.2) +
  scale_y_continuous(' ', expand = c(0, 0), limits=c(0,3e5)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0), limits=c(20,260), "Generations") +
  theme_classic() +
  scale_color_manual(name="Nutrient",values=c('#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a', 'gray', 'gray')) +
  scale_linetype_manual(values=c('twodash', 'twodash', "solid")) +
  annotate("text", x = 50, y = 3e4, label = 'Glucose', parse = TRUE) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme(text=element_text(size=12), legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
glUplot <- glUplot %+% subset(data3, Nutrient %in% c("Glucose"))
glUplot <- glUplot + aes(group=rev(Unique_ID))
glUplot

fluor <- plot_grid(glnplot,urplot, glUplot, labels = c("A)"), ncol=1, align = "v")
fluor

glnplot <- ggplot(data3, aes(x=Generation,y=(FSC_median), colour = Unique_ID)) +
  geom_line(aes(linetype = Type), show.legend = FALSE, size=1.2) +
  scale_y_continuous(' ', expand = c(0, 0), limits=c(2e5,14e5)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0), limits=c(20,260), " ") +
  theme_classic() +
  scale_fill_manual(values=c("grey", "grey", "green")) + 
  scale_color_manual(name="Nutrient",values=c('#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b', '#08306b','#2171b5','#6baed6', 'grey', 'grey')) +
  scale_linetype_manual(values=c('twodash', 'twodash', "solid")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme(text=element_text(size=12), legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
glnplot <- glnplot %+% subset(data3, Nutrient %in% c("Glutamine"))
glnplot <- glnplot + aes(group=rev(Unique_ID))
glnplot

urplot <- ggplot(data3, aes(x=Generation,y=(FSC_median), colour = Unique_ID)) +
  #geom_point(size=1, show.legend = FALSE) +
  geom_line(aes(linetype = Type), show.legend = FALSE, size=1.2) +
  scale_y_continuous('Forward scatter (a.u.)', expand = c(0, 0), limits=c(2e5,14e5)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0), limits=c(20,260), " ") +
  theme_classic() +
  scale_color_manual(name="Nutrient",values=c('#662506','#cc4c02','#fed976','#feb24c','#e31a1c','#fc4e2a', '#fd8d3c','#bd0026','#800026', 'grey', 'grey')) +
  scale_linetype_manual(values=c('twodash', 'twodash', "solid")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme(text=element_text(size=12), legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
urplot <- urplot %+% subset(data3, Nutrient %in% c("Urea"))
urplot <- urplot + aes(group=rev(Unique_ID))
urplot

glUplot <- ggplot(data3, aes(x=Generation,y=(FSC_median), colour = Unique_ID)) +
  #geom_point(size=1, show.legend = FALSE) +
  geom_line(aes(linetype = Type), show.legend = FALSE, size=1.2) +
  scale_y_continuous(' ', expand = c(0, 0), limits=c(2e5,14e5)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), expand = c(0, 0), limits=c(20,260), "Generations") +
  theme_classic() +
  scale_color_manual(name="Nutrient",values=c('#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a', 'gray', 'gray')) +
  scale_linetype_manual(values=c('twodash', 'twodash', "solid")) +
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme(text=element_text(size=12), legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
glUplot <- glUplot %+% subset(data3, Nutrient %in% c("Glucose"))
glUplot <- glUplot + aes(group=rev(Unique_ID))
glUplot

size <- plot_grid(glnplot, urplot, glUplot, labels = c("B)"), ncol=1, align = "v")
size

plot_grid(fluor, size, ncol = 2, align = 'h')
ggsave(filename="Figure 3 Test.tiff", scale = 1.25, width=6.5, height=5, units="in", dpi=300)
###########
##Figure S4
###########

new <- melt(data3, id.vars = c("Generation", "Unique_ID", "Population", "Nutrient", "Type"), measure.vars = c("PopProp_0copy","PopProp_1copy","PopProp_2copy", "PopProp_3plus"))
new <- data.frame(Generation=as.character(new[,1]), Unique_ID=as.character(new[,2]), Population=as.character(new[,3]), Nutrient=as.character(new[,4]), Type=as.factor(new[,5]), variable=as.character(new[,6]), value=as.numeric(new[,7])) 

new$Generation <- factor(new$Generation, levels = c(25, 33, 41, 54, 62, 70, 79, 87, 95, 103, 116, 124, 132, 145, 153, 161, 174, 182, 190, 211, 219, 232, 244, 257, 267))
new$variable <- factor(new$variable, levels = c("PopProp_3plus", "PopProp_2copy", "PopProp_1copy", "PopProp_0copy"))


theplot <- ggplot(new, aes(Generation, value, fill=variable)) + 
  theme_classic() +
  theme(legend.position = "top", legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  geom_bar(colour = "black", position = "fill", stat = "identity", width = 0.9) +
  scale_y_continuous("Proportion", expand = c(0, 0), breaks = c(0.5, 1.0)) + 
  facet_wrap('Unique_ID', ncol = 1, strip.position = c("right")) + 
  scale_fill_manual(values = c('#006d2c', 'green', '#a1d99b', '#e5f5e0')) + 
  theme(text=element_text(size=12), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=11), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=11))

glnonly <- theplot %+% subset(new, Nutrient %in% c("Glutamine") & Type %in% c("Experimental"))
glnonly

###########
##Figure S5
###########

data5 <- read_csv("Downloads/MASTER statistics spreadsheet - PopCorrPlot (1).csv")
data5 <- data.frame(MASTER_statistics_spreadsheet_PopCorrPlot_1_)
data5$Fluorescence=as.numeric(as.character(data5$Fluorescence))

equation = function(x) {
  lm_coef <- list(a = round(coef(x)[1], digits = 2),
                  b = round(coef(x)[2], digits = 2),
                  r2 = round(summary(x)$r.squared, digits = 2));
  lm_eq <- substitute(~italic(R)^2~"="~r2,lm_coef)
  as.character(as.expression(lm_eq));                 
}

fit <- lm(Copy..at.GAP1 ~ Fluorescence, data = subset(data5, Media %in% c("Glucose", "Glutamine", "Urea")))
summary(fit)

myplot <- ggplot(data5, aes(x=Fluorescence,y=(Copy..at.GAP1))) +
  geom_point(size=1.5, aes(colour=Media, shape = as.factor(Generation))) +
  scale_y_continuous(expand = c(0, 0), 'Relative GAP1 read depth', limits=c(-.2,4.5)) +
  scale_x_continuous(expand = c(0, 0), "Median normalized fluorescence", limits=c(-.01,.75)) +
  theme_classic() +
  #annotate("rect", xmin = 0.4, xmax = .8, ymin = .4, ymax = .8, fill="white", colour="red") +
  annotate("text", x = .4, y = .6, label = equation(fit), parse = TRUE) + 
  scale_color_manual(values = c('#ae017e','#2171b5','#fd8d3c')) +
  #guides(colour = guide_legend(override.aes = list(size=2))) + 
  geom_smooth(method=lm, se=FALSE, color = 'black') + 
  theme(text=element_text(size=12), legend.position = c(.1,.8), legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))

myplot %+% subset(data5, Media %in% c("Glucose", "Glutamine", "Urea"))
ggsave(filename="Figure 5 Test.tiff", scale = 4, width=1.5, height=1.5, units="in", dpi=300)


'#dd3497','#ae017e' <- glucose
'#feb24c', '#fd8d3c' <- urea
'#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b','#2171b5','#6baed6' <- gln



###########
##Figure S6
###########
data1 <- data.frame(Fitness_Final_csv_SL_clone_fitness_pval_NEW_csv_1_)
data1$Generation<-factor(data1$Generation,levels=unique(data1$Generation))

##Other fitness plots
base <- ggplot(data1, aes(x=as.factor(Clone), y=fitness_coef+1)) +
  geom_bar(aes(fill = as.factor(Generation)), stat = 'identity', width = 0.9) + 
  facet_wrap(~Population, nrow=1, strip.position = c("bottom"), scales = "free_x") + 
  geom_errorbar(aes(ymin=1+ci2.5, ymax=1+ci97.5), width=.4) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 11, angle=75, vjust=0.6), legend.position = 'none', legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=11), axis.title.x = element_text(face="bold", size=12)) + 
  scale_fill_manual(values = c('grey','#74c476','#00441b')) +
  scale_x_discrete('Clone') +
  geom_hline(yintercept=1.0, linetype="dotted", color="black", size=.8) +
  scale_y_continuous('Relative Fitness', expand = c(0, 0), limits = c(0,1.5))
base
base %+% subset(data1, Type %in% c("Experimental") & Generation %in% c("150", "250"))
base %+% subset(data1, Keep.Remove %in% c("Keep"))
base %+% subset(data1, Keep.Remove %in% c("Keep"))

ggsave(filename="Figure S6 Test.tiff", scale = 1.5, width=5, height=4, units="in", dpi=300)


###########
##Figure S7
###########
clonedata <- read_csv("Downloads/MASTER Ministat Summary Table - Copy of SUMMARY-CLONES STATS (2).csv")
clonedata <- data.frame(MASTER_Ministat_Summary_Table_Copy_of_SUMMARY_CLONES_STATS_2_)
glnplot <- ggplot(clonedata, aes(x=as.character(Generation),y=Copy..at.GAP1)) +
  #geom_dotplot(aes(fill = Population), binaxis = 'y', stackdir='center', dotsize = .5) + 
  theme_classic() + 
  scale_y_continuous('Estimated GAP1 Copy Number', expand = c(0, 0), limits=c(1,7)) +
  geom_jitter(aes(color = Population), size = 2, shape=16, show.legend = FALSE, position=position_jitter(0.15)) + 
  scale_x_discrete("Generation") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") +
  scale_color_manual(name="Nutrient",values=c('#a1d99b','#74c476','#41ab5d', 'green', '#238b45','#006d2c','#00441b', '#08306b','#2171b5', 'blue')) + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="black") +
  theme(legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))

glnplot <- glnplot %+% subset(clonedata, Condition %in% c("Glutamine") & Type %in% c("Clone"))
glnplot

ureaplot <- ggplot(clonedata, aes(x=as.character(Generation),y=as.numeric(Copy..at.DUR3))) +
  #geom_dotplot(aes(fill = Population), width = .5, binaxis = 'y', stackdir='center', dotsize = .5) + 
  theme_classic() + 
  scale_y_continuous('Estimated DUR3 Copy Number', expand = c(0, 0), limits=c(1,7)) +
  geom_jitter(aes(color = Population), size = 2, shape=16, show.legend = FALSE, position=position_jitter(0.25)) + 
  scale_x_discrete("Generation") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="black") +
  theme(legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))

ureaplot <- ureaplot %+% subset(clonedata, Condition %in% c("Urea") & Type %in% c("Clone"))
ureaplot

ureapopz <- ggplot(clonedata, aes(x=as.character(Generation),y=as.numeric(Copy..at.DUR3))) +
  #geom_dotplot(aes(fill = Population), width = .5, binaxis = 'y', stackdir='center', dotsize = .5) + 
  theme_classic() + 
  scale_y_continuous('Estimated DUR3 Copy Number', expand = c(0, 0), limits=c(1,7)) +
  geom_jitter(aes(color = Population), size = 2, shape=16, show.legend = FALSE, position=position_jitter(0.25)) + 
  scale_x_discrete("Generation") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="black") +
  #scale_color_manual(name="Nutrient",values=c('#662506','#cc4c02','#fed976','#feb24c','#e31a1c','#fc4e2a', '#fd8d3c','#bd0026','#800026', 'grey', 'grey')) +
  theme(legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
ureapopz <- ureapopz %+% subset(clonedata, Condition %in% c("Urea") & Type %in% c("Pop") & Population %in% c('ure_01', 'ure_02', 'ure_05', 'ure_07', 'ure_ctrl01', 'ure_ctrl02'))
ureapopz

glnsize <- ggplot(clonedata, aes(x=as.character(Generation),y=as.numeric(Size..KB.)/1000)) +
  #geom_dotplot(aes(fill = Population), width = .5, binaxis = 'y', stackdir='center', dotsize = .5) + 
  theme_classic() + 
  scale_y_continuous('Estimated CNV allele size (kb)', expand = c(0, 0), limits=c(-10,300)) +
  geom_jitter(aes(color = Population), size = 2, shape=16, show.legend = TRUE, position=position_jitter(0.15)) + 
  scale_x_discrete("Generation") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="black") +
  scale_color_manual(name="Nutrient",values=c('#a1d99b','#74c476','#41ab5d', 'green', '#238b45','#006d2c','#00441b', '#08306b','#2171b5','#6baed6', 'blue')) +
  theme(legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
glnsize <- glnsize %+% subset(clonedata, Condition %in% c("Glutamine", "Control") & Type %in% c("Clone"))
glnsize

ursize <- ggplot(clonedata, aes(x=as.character(Generation),y=as.numeric(Size..KB.)/1000)) +
  #geom_dotplot(aes(fill = Population), width = .5, binaxis = 'y', stackdir='center', dotsize = .5) + 
  theme_classic() + 
  scale_y_continuous('Estimated CNV allele size (kb)', expand = c(0, 0), limits=c(-10,300)) +
  geom_jitter(aes(color = Population), size = 2, shape=16, show.legend = TRUE, position=position_jitter(0.25)) + 
  scale_x_discrete("Generation") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="black") + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="black") +
  #scale_color_manual(name="Nutrient",values=c('#662506','#cc4c02','#fed976','#feb24c','#e31a1c','#fc4e2a', '#fd8d3c','#bd0026','#800026', 'grey', 'grey')) +
  theme(legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))
ursize <- ursize %+% subset(clonedata, Condition %in% c("Urea") & Type %in% c("Clone"))
ursize

CN <- plot_grid(glnplot, ureaplot, labels = c("B)", "C)"), align = "h", nrow = 1)
sizes <- plot_grid(glnsize, ursize, labels = c("E)", "F)"), align = "h", nrow = 1)

plot_grid(CN, sizes, align = "h", nrow = 2)
ggsave(filename="Figure S7B.pdf", scale = 1.5, width=7.5, height=5, units="in", dpi=300)

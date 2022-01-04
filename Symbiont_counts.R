setwd("~/Documents/Bioinformatics_scripts/R_scripts/VTK/VTK_2022/")

#library(car)
library(drc)
library(tidyr)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
P4=c("#2E33D1", "#FFEE32","#D37D47", "#F43535") #27, 29, 32, 34
P3=c("#6a994e", "#4ea8de","#bc4749")
  
input=read.table("./developing/symbiont_counts.txt", header = T, sep = "\t") 

#long version of the table (easier to handle for ggplot)
input_long=melt(input, id.vars=c("Treatment","Temperature", "Replicate"))
input_long$Temperature=as.factor(input_long$Temperature) #make sure temperature is numeric
input_long$Treatment=factor(input_long$Treatment, levels = c("Native", "Antibiotic", "Pathogenic")) #make sure temperature is numeric


#check for normality 
hist(summary1$mean) #Gaussian shape?
shapiro.test(summary1$mean) #normality?


#ANOVAs overall comparisons
anova_results=aov(value ~ Temperature*Treatment, data = input_long)
summary(anova_results)


#t-tests multiple comparisons between temperatures
native=subset(input_long, Treatment== "Native" & !is.na(input_long$value))
pairwise.t.test(native$value, native$Temperature, p.adj = "fdr")

antibiotic=subset(input_long, Treatment== "Antibiotic" & !is.na(input_long$value))
pairwise.t.test(antibiotic$value, antibiotic$Temperature, p.adj = "fdr")

pato=subset(input_long, Treatment== "Pathogenic"& !is.na(input_long$value))
pairwise.t.test(pato$value, pato$Temperature, p.adj = "fdr")

summary1=input_long %>% drop_na() %>% group_by(Treatment, Temperature, Replicate) %>% dplyr::summarise(mean=mean(value), sd=sd(value))
positions=summary1 %>% group_by(Treatment, Temperature) %>% dplyr::summarise(position=max(mean)+25000)

ggplot(summary1, aes(x=Temperature, y=mean, fill=Temperature)) +
  geom_boxplot()+ labs(x="", y= "Symbiodiniaceae density \n(cells/ml)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme_bw() + facet_grid(~Treatment) + scale_fill_manual(values=P4) +
  # annotate(geom="text", x=1, y=positions$position[1], label= "A") + 
  # annotate(geom="text", x=2, y=positions$position[2], label= "B") + 
  # annotate(geom="text", x=3, y=positions$position[3], label= "C") +
  # annotate(geom="text", x=4, y=positions$position[4], label= "D") +
  # annotate(geom="text", x=5, y=positions$position[5], label= "A") + 
  # annotate(geom="text", x=6, y=positions$position[6], label= "A") + 
  # annotate(geom="text", x=7, y=positions$position[7], label= "B") +
  # annotate(geom="text", x=8, y=positions$position[8], label= "C") +
  # annotate(geom="text", x=9, y=positions$position[9], label= "A") + 
  # annotate(geom="text", x=10, y=positions$position[10], label= "B") + 
  # annotate(geom="text", x=11, y=positions$position[11], label= "C") +
  # annotate(geom="text", x=12, y=positions$position[12], label= "D") +
  # theme(legend.position = "none")  

#t-tests multiple comparisons between treatments

summary2=input_long %>% drop_na() %>% group_by(Treatment, Temperature, Replicate) %>% dplyr::summarise(mean=mean(value), sd=sd(value))
ggplot(summary1, aes(x=Treatment, y=mean, fill=Treatment)) +
  geom_boxplot()+ labs(x="", y= "Symbiodiniaceae density \n(cells/ml)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  theme_bw() + facet_grid(~Temperature) + scale_fill_manual(values=P3) 
  


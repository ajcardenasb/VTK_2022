setwd("~/Documents/Bioinformatics_scripts/R_scripts/VTK/VTK_2022/")

  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(reshape2)

#read data in -> table should contain 6 columns (ID,	Strain,	Temperature,	Treatment,	Replicate, symbiont_count_adj) and be a 
# simplified text file with no formatting! 
input<-read.table("Input_data/VTK2022_Symbiont_counts.txt", header = T, sep = "\t") 

input$Temperature<-as.factor(input$Temperature) #convert temperature to factors
input$Treatment<-factor(input$Treatment, levels = c("native", "antibiotic", "pathogenic")) #re-order level of factors

#check for normality 
hist(sqrt(input$symbiont_count_adj)) #Gaussian shape?
shapiro.test(sqrt(input$symbiont_count_adj)) #normality?

hist(log10(input$symbiont_count_adj)) #Gaussian shape?
shapiro.test(log10(input$symbiont_count_adj)) #normality?

#transform variables befor doing the anova
input$sqrt_counts<-sqrt(input$symbiont_count_adj)

#ANOVAs overall comparisons
anova_results=aov(sqrt_counts ~ Strain+Temperature+Treatment+Strain:Treatment:Temperature, data = input)
summary(anova_results)

#subset your dataframe by strains
H2_counts=subset(input, Strain == "H2")
FOO3_counts=subset(input, Strain == "FOO3")

anova_results_H2=aov(sqrt_counts ~ Treatment:Temperature, data = H2_counts)
summary(anova_results_H2)

anova_results_FOO3=aov(sqrt_counts ~ Treatment:Temperature, data = FOO3_counts)
summary(anova_results_FOO3)

#t-tests multiple comparisons between temperatures

#h2
h2_native=subset(input, Strain == "H2" & Treatment== "native" & !is.na(input$sqrt_counts))
pairwise.t.test(h2_native$sqrt_counts,h2_native$Temperature, p.adj = "fdr")

h2_antibiotic=subset(input, Strain == "H2" & Treatment== "antibiotic" & !is.na(input$sqrt_counts))
pairwise.t.test(h2_antibiotic$sqrt_counts, h2_antibiotic$Temperature, p.adj = "fdr")

h2_pato=subset(input, Strain == "H2" & Treatment== "pathogenic" & !is.na(input$sqrt_counts))
pairwise.t.test(h2_pato$sqrt_counts, h2_pato$Temperature, p.adj = "fdr")

#foo3
foo_native=subset(input, Strain == "FOO3" & Treatment== "native" & !is.na(input$sqrt_counts))
pairwise.t.test(foo_native$sqrt_counts,foo_native$Temperature, p.adj = "fdr")

foo_antibiotic=subset(input, Strain == "FOO3" & Treatment== "antibiotic" & !is.na(input$sqrt_counts))
pairwise.t.test(foo_antibiotic$sqrt_counts, foo_antibiotic$Temperature, p.adj = "fdr")

foo_pato=subset(input, Strain == "FOO3" & Treatment== "pathogenic" & !is.na(input$sqrt_counts))
pairwise.t.test(foo_pato$sqrt_counts, foo_pato$Temperature, p.adj = "fdr")

#plotting
#defining color palette values
P4<-c("#2E33D1", "#FFEE32","#D37D47", "#F43535") #27, 29, 32, 34
P3<-c("#C3BF6D", "#4ea8de","#bc4749")


#plots by ttemperature
ggplot(H2_counts, aes(x=Temperature, y=sqrt_counts, fill=Temperature)) +
  geom_boxplot()+ labs(x="", y= "Square root of Symbiodiniaceae density \n(cells/mm host oral disc)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme_bw() + facet_grid(~Treatment) + scale_fill_manual(values=P4) 

ggplot(FOO3_counts, aes(x=Temperature, y=sqrt_counts, fill=Temperature)) +
  geom_boxplot()+ labs(x="", y= "Square root of Symbiodiniaceae density \n(cells/mm host oral disc)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme_bw() + facet_grid(~Treatment) + scale_fill_manual(values=P4) 


#t-tests multiple comparisons between treeatments
h2_30=subset(input, Strain == "H2" & Temperature== "30" & !is.na(input$sqrt_counts))
pairwise.t.test(h2_30$sqrt_counts, h2_30$Treatment, p.adj = "fdr")

h2_33=subset(input, Strain == "H2" & Temperature== "33" & !is.na(input$sqrt_counts))
pairwise.t.test(h2_33$sqrt_counts, h2_33$Treatment, p.adj = "fdr")

h2_36=subset(input, Strain == "H2" & Temperature== "36" & !is.na(input$sqrt_counts))
pairwise.t.test(h2_36$sqrt_counts, h2_36$Treatment, p.adj = "fdr")

h2_39=subset(input, Strain == "H2" & Temperature== "39" & !is.na(input$sqrt_counts))
pairwise.t.test(h2_39$sqrt_counts, h2_39$Treatment, p.adj = "fdr")

#foo3
foo_30=subset(input, Strain == "FOO3" & Temperature== "30" & !is.na(input$sqrt_counts))
pairwise.t.test(foo_30$sqrt_counts, foo_30$Treatment, p.adj = "fdr")

foo_33=subset(input, Strain == "FOO3" & Temperature== "33" & !is.na(input$sqrt_counts))
pairwise.t.test(foo_33$sqrt_counts, foo_33$Treatment, p.adj = "fdr")

foo_36=subset(input, Strain == "FOO3" & Temperature== "36" & !is.na(input$sqrt_counts))
pairwise.t.test(foo_36$sqrt_counts, foo_36$Treatment, p.adj = "fdr")

foo_39=subset(input, Strain == "FOO3" & Temperature== "39" & !is.na(input$sqrt_counts))
pairwise.t.test(foo_39$sqrt_counts, foo_39$Treatment, p.adj = "fdr")

#plots by treatment
ggplot(H2_counts, aes(x=Treatment, y=sqrt_counts, fill=Treatment)) +
  geom_boxplot()+ labs(x="", y= "Square root of Symbiodiniaceae density \n(cells/mm host oral disc)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme_bw() + facet_grid(~Temperature) + scale_fill_manual(values=P3) 


ggplot(FOO3_counts, aes(x=Treatment, y=sqrt_counts, fill=Treatment)) +
  geom_boxplot()+ labs(x="", y= "Square root of Symbiodiniaceae density \n(cells/mm host oral disc)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme_bw() + facet_grid(~Temperature) + scale_fill_manual(values=P3) 

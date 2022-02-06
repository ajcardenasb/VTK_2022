library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)

#read data in -> table should contain 6 columns (ID,	Strain,	Temperature,	Treatment,	Replicate, symbiont_count_adj) and be a 
# simplified text file with no formatting! 
input<-read.table(???) 

input$Temperature<-as.factor(input$Temperature) #convert temperature to factors
input$Treatment<-factor(input$Treatment, levels = c("native", "antibiotic", "pathogenic")) #re-order level of factors

#check for normality 
hist(???) #Gaussian shape?
shapiro.test(???) #normality?

#transform variables before doing the anova
input$???_counts<-sqrt(input$symbiont_count_adj)

#ANOVAs overall comparisons
anova_results=aov(sqrt_counts ~ Strain+Temperature+Treatment+Strain:Treatment:Temperature, data = input)
summary(anova_results)

#subset your dataframe by strains
H2_counts=subset(input, ???) #only H2
FOO3_counts=subset(input, ???) #only FOO3

anova_results_H2=aov(??? ~ Treatment:Temperature, data = ???)
summary(anova_results_H2)

anova_results_FOO3=aov(??? ~ Treatment:Temperature, data = ???)
summary(anova_results_FOO3)

#t-tests multiple comparisons between temperatures

#example
h2_native=subset(input, Strain == "H2" & Treatment== "native" & !is.na(input$sqrt_counts))
pairwise.t.test(h2_native$sqrt_counts,h2_native$Temperature, p.adj = "fdr")


#plotting
#defining color palette values
P4<-c("#2E33D1", "#FFEE32","#D37D47", "#F43535") #27, 29, 32, 34
P3<-c("#6a994e", "#4ea8de","#bc4749")


#plots by temperature
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


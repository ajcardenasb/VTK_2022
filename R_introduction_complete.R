
#1. set your working directory using the "setwd" function
setwd("~/Documents/Bioinformatics_scripts/R_scripts/VTK/VTK_2022/")
#2. load ggplot2, dplyr and reshape2 libraries using the "library" function
library(ggplot2)
library(dplyr)
library(reshape2)

#3. import and inspect data using the read.table function

input=read.table("Input_data/VTK2022_colony_taxonomy.txt", header = T, sep = "\t")
input$closest_related_species=gsub("A.", "Alteromonas", input$closest_related_species)
input$Isolate_ID=gsub("F003", "FOO3", input$Isolate_ID)

#4. filter dataframe to keep only Vibrios 
library(data.table)
subset(input,   Isolate_ID  %like% "*H2*" )
input %>% filter(Genus == "Vibrio")

input[input$Genus == "Vibrio",]
input[grepl("Vibrio",input$Genus ),]

#5. create a new column with the strain 
input$Strain=gsub(".*_H2_.*", "H2" , input$Isolate_ID)
input$Strain=gsub(".*_FOO3_.*", "FOO3" , input$Strain)

input$Strain2=ifelse(grepl("H2", input$Isolate_ID),"H2","FOO3")

#6 filter dataframe to keep only strains isolated from H2, without using subset
input[grepl("H2",input$Isolate_ID ),]

#7 calculate the frecuencies per Genus per strain
input %>% group_by(Strain, Genus) %>% tally()

#8. calculate the total number of different species per strain
#hint use the unique command to remove duplicates (will remove rows that with ALL identical values)
unique(input[,c(7:8)]) %>% group_by(Strain) %>% tally()

#9. create a new column with (Phylum - Genus)
input$tax_new= paste(input$Phylum, input$Genus, sep = " - ")

#10 claculate the overlap of genera between the two strains using the function "intersect"
intersect(input[input$Strain == "H2",]$Genus, input[input$Strain == "FOO3",]$Genus)

# 11 - create stack bar plots of your isolates uwing the following formula:
#ggplots needs to have a colum for values in x and values in y. 
#Additionally you can assign colors to another categorical/discrete variable
frec_genus=input %>% group_by(Strain, Genus) %>% tally()
ggplot(frec_genus, aes(fill=Genus, y=n, x=Strain)) +
  geom_bar(position="stack", stat="identity")+
  labs(x="", y = "Number of isolates") + theme_bw()

frec_fam=input %>% group_by(Strain, Family) %>% tally()
ggplot(frec_fam, aes(fill=Family, y=n, x=Strain)) +
  geom_bar(position="stack", stat="identity")  + theme_classic() +
  labs(x="", y = "Number of isolates")


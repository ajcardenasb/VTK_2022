
#1. set your working directory using the "setwd" function

#2. load dplyr and reshape2 libraries using the "library" function

#3. import and inspect data using the read.table function

#4. filter dataframe to keep only strains classified as Vibrios 

#5. create a new column with the strain 

#6 filter dataframe to keep only strains isolated from H2, without using subset

#7 calculate the frecuencies per Genus per strain

#8. calculate the total number of different species per strain
#hint use the unique command to remove duplicates (will remove rows that with ALL identical values)

#9. create a new column with (Phylum - Genus)

#10 claculate the overlap of genera between the two strains using the function "intersect"

# 11 - create stack bar plots of your isolates uwing the following formula:
#ggplots needs to have a colum for values in x and values in y. 
#Additionally you can assign colors to another categorical/discrete variable
ggplot(???, aes(fill=???, y=???, x=???) + 
  geom_bar(position="stack", stat="identity")

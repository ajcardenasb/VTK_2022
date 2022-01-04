setwd("~/Documents/Bioinformatics_scripts/R_scripts/VTK/VTK_2022/")

#library(car)
library(drc)
library(tidyr)
library(RColorBrewer)
library(dplyr)
library(ggplot2)

input=read.table("./developing/Input_PAM.txt", header = T, sep = "\t") 

#break sample names into sample variables
input=tidyr::separate(input,Sample,into =c("Strain", "Treatment", "Temperature", "Replicate" ),sep = "_",remove = FALSE,extra = "merge")
input$Replicate=as.factor(paste("R", sprintf("%02d",as.numeric(input$Replicate)), sep = ""))
input$Temperature=as.numeric(input$Temperature)

#Making groups
input$Group1=as.factor(paste(input$Strain,input$Treatment, input$Replicate, sep = "_")) # replicate per treatment per strain
levels(input$Group1)
input$Group2=as.factor(paste(input$Strain,input$Treatment, sep = "_")) # treatment per strain
levels(input$Group2)
#input$Group2=paste(input$Site,input$Species, input$Colony, sep = "_") # species per site
#input$Group3=paste(input$Species, sep = "_") # across species
#input$Group4=paste(input$Site, sep = "_") # across sites

#if PAM values are > 1, then convert to 0-1 values
#input$PAM=ifelse(input$PAM >= 1, input$PAM/1000,input$PAM )
temp_x<- seq(30, 40.5, length = 100)

################################# drm  Demo ####################################

drm(PAM ~ Temperature, data=input[input$Group1=="H2_Antibiotic_R01",],
    fct = LL.3(names = c('hill', 'max', 'ed50')), 
    upperl = c(120, 0.72, 40), lowerl = c(10, 0.55, 30))

################################# drm in loops ####################################

#####################################################
## fitting the curves for each replicate (group 1) ## 
#####################################################

#identify and remove replicates with missing data
missing=input$Group1[is.na(input$PAM)] 
input2=subset(input, !Group1 %in% missing)

mod1=lapply(unique(input2$Group1), 
            function(x) drm(PAM ~ Temperature, data=input2[input2$Group1==x,],
                            fct = LL.3(names = c('hill', 'max', 'ed50')), 
                            upperl = c(120, 0.72, 40), lowerl = c(10, 0.55, 30)))
pred1=lapply(mod1, function(x) predict(x, data.frame(Temperature = temp_x)))
pred_df=as.data.frame(do.call(rbind, pred1))
colnames(pred_df)= round(temp_x, digits = 2)
pred_df$Sample=unique(input2$Group1)
pred_df_long=reshape2::melt(pred_df)
colnames(pred_df_long)[2:3]=c("Temperature", "Fv/Fm")

#extract ED50 
ed50_list=lapply(c(1:length(mod1)), function(x) mod1[[x]][["coefficients"]][["ed50:(Intercept)"]])
ed50_df=as.data.frame(do.call(rbind, ed50_list))
ed50_df$Sample=unique(input2$Group1)
colnames(ed50_df)[1]="ED50"
#write.table(ed50_df[,c(2,1)], "outputs/CBASS_ED50s", row.names = F, sep = "\t")

#prepare for plotting
pred_df_long$ED50=ed50_df$ED50[match(pred_df_long$Sample, ed50_df$Sample)]
pred_df_long$Temperature=as.numeric(as.character(pred_df_long$Temperature))
pred_df_long=tidyr::separate(pred_df_long,Sample,into =c("Strain", "Treatment", "Replicate" ),sep = "_",remove = FALSE,extra = "merge")
pred_df_long$Replicate_ED50=paste(pred_df_long$Replicate, " - ",round(pred_df_long$ED50, 2),"ºC", sep=" ")

#define colors
n <- length(unique(input$Group1))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
Colours = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(pred_df_long, aes(x=Temperature, y=`Fv/Fm`, group=Sample)) +
  geom_line(aes(color=Replicate)) +
  geom_segment(data=pred_df_long, aes(x = ED50, y = 0, xend = ED50, yend = 0.65,color=Replicate), linetype=3) +
  geom_text(data=pred_df_long, mapping=aes(x=ED50, y=0.68, label=Replicate_ED50), size=2.5, angle=45, hjust=0, check_overlap = T) +
  guides(fill=guide_legend())+
  theme(panel.grid = element_blank()) +
  scale_color_manual(values=Colours) +  theme_classic() +
  scale_x_continuous(breaks=c(30, 33, 36, 39), limits = c(30,40), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,0.9), expand = c(0, 0)) + labs(color='') +
  facet_grid(Strain~Treatment) 


################################################################
## fitting the curves for each treament and strain (group 2) ## 
################################################################

#identify and remove replicates with missing data

mod2=lapply(unique(input2$Group2), 
            function(x) drm(PAM ~ Temperature, data=input2[input2$Group2==x,],
                            fct = LL.3(names = c('hill', 'max', 'ed50')), 
                            upperl = c(120, 0.72, 40), lowerl = c(10, 0.55, 30)))
pred2=lapply(mod2, function(x) predict(x, data.frame(Temperature = temp_x)))
pred2_df=as.data.frame(do.call(rbind, pred2))
colnames(pred2_df)= round(temp_x, digits = 2)
pred2_df$Sample=unique(input2$Group2)
pred2_df_long=reshape2::melt(pred2_df)
colnames(pred2_df_long)[2:3]=c("Temperature", "Fv/Fm")

#extract ED50 
ed50_list2=lapply(c(1:length(mod2)), function(x) mod2[[x]][["coefficients"]][["ed50:(Intercept)"]])
ed50_df2=as.data.frame(do.call(rbind, ed50_list2))
ed50_df2$Sample=unique(input2$Group2)
colnames(ed50_df2)[1]="ED50"
#write.table(ed50_df2[,c(2,1)], "outputs/CBASS_ED50s", row.names = F, sep = "\t")

#prepare for plotting
pred2_df_long$ED50=ed50_df2$ED50[match(pred2_df_long$Sample, ed50_df2$Sample)]
pred2_df_long$Temperature=as.numeric(as.character(pred2_df_long$Temperature))
pred2_df_long=tidyr::separate(pred2_df_long,Sample,into =c("Strain", "Treatment" ),sep = "_",remove = FALSE,extra = "merge")
pred2_df_long$treatment_ED50=paste(pred2_df_long$Treatment, " - ",round(pred2_df_long$ED50, 2),"ºC", sep=" ")

#define colors
n <- length(unique(input$Group1))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
Colours = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(pred2_df_long, aes(x=Temperature, y=`Fv/Fm`, group=Treatment)) +
  geom_line(aes(color=Treatment)) +
  geom_segment(data=pred2_df_long, aes(x = ED50, y = 0, xend = ED50, yend = 0.65,color=Treatment), linetype=3) +
  geom_text(data=pred2_df_long, mapping=aes(x=ED50, y=0.68, label=treatment_ED50), size=2.5, angle=45, hjust=0, check_overlap = T) +
  guides(fill=guide_legend())+
  theme(panel.grid = element_blank()) +
  scale_color_manual(values=Colours) +  theme_classic() +
  scale_x_continuous(breaks=c(30, 33, 36, 39), limits = c(30,40), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,0.9), expand = c(0, 0)) + labs(color='') +
  facet_grid(Strain~.) 


summary(mod2)
compParm(mod2, 'ed50')
compParm(DRCpam, 'ed50', "-")

#### statistics

pop=drm(PAM ~ Temperature, data = input2, curveid=Group2,
    fct = LL.3(names = c('hill', 'max', 'ed50')), 
    upperl = c(120, 0.72, 40), lowerl = c(10, 0.55, 30))

summary(pop)
compParm(pop, 'ed50', "-")
plot(pop)
points(input2$Temperature, input2$PAM)
ED(pop, c(50))

col=drm(PAM ~ Temperature, data = input2, curveid=Group1,
        fct = LL.3(names = c('hill', 'max', 'ed50')), 
        upperl = c(120, 0.72, 40), lowerl = c(10, 0.55, 30))

summary(col)
compParm(col, 'ed50')
plot(col)
points(input2$Temperature, input2$PAM)
ED(col, c(50))


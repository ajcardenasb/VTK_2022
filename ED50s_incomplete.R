library(drc)
library(tidyr)
library(dplyr)
library(ggplot2)

input<-???

???=as.numeric(???) # we need temperature as numeric

#############################################
## fitting the curves for each individual ## 
#############################################
#remove rows with missing data
input<-input[complete.cases(input), ]

#create a group to define colony ID
input$Colony=paste(???) # replicate per treatment per strain

#fit the PAM data from each colony to the response curve

# drm  Demo to one colony
drm(PAM ~ Temperature, data=input[input$Colony=="H2_antibiotic_1",],
    fct = LL.3(names = c('Slope', 'Max', 'ED50')), 
    upperl = c(120, 0.72, 40), lowerl = c(10, 0.55, 30))

### now to all colonies in a loop
mod1<-lapply(unique(input$Colony), 
             function(x) drm(PAM ~ Temperature, data=input[input$Colony==x,],
                             fct = LL.3(names = c('Slope', 'Max', 'ED50')),
                             upperl = c(120, 0.72, 40), lowerl = c(10, 0.55, 30)))

#############################
## extract and plot ED50s ## 
#############################

#extract ED50 
ed50_list<-lapply(c(1:length(mod1)), function(x) mod1[[x]][["coefficients"]][["ED50:(Intercept)"]])
ed50_df<-as.data.frame(do.call(rbind, ed50_list))
ed50_df$Sample=unique(input$Colony)
ed50_df=tidyr::separate(ed50_df,Sample,into =c("Strain", "Treatment" , "Replicate"),sep = "_",remove = FALSE,extra = "merge")
ed50_df$Treatment=factor(ed50_df$Treatment, levels=c("native", "antibiotic", "pathogenic"))
colnames(ed50_df)[1]="ED50"
#write.table(ed50_df, "VTK2022_ED50.txt", row.names = F, sep = "\t")

#plot 1 ED50 dsitributions x treatment 
treatment_colours<-c("#C3BF6D", "#4ea8de","#bc4749")
ggplot(???, aes(x=???, y=???, fill =???)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, lwd=0.8) + theme_bw() +
  scale_fill_manual(values=treatment_colours) + facet_grid(~Strain)

#### statistic comparison on ED50s
#check for normality 
???

#transform variables befor doing the anova
???

#Overall anova.
???

#t-tests multiple comparisons between treatments
???
  
#######################
## plot ED50 curves ## 
#######################

# temperature ranges 
temp_x<- seq(30, 40.5, length = 100)

# prediction of the fitted values corresponding to the range of temperatures above
pred1<-lapply(mod1, function(x) predict(x, data.frame(Temperature = temp_x)))
pred_df<-as.data.frame(do.call(rbind, pred1))
colnames(pred_df)= round(temp_x, digits = 2)
pred_df$Sample=unique(input$Colony)
pred_df_long=reshape2::melt(pred_df, id.vars=c("Sample"))
colnames(pred_df_long)[2:3]<-c("Temperature", "Fv/Fm")
pred_df_long=tidyr::separate(pred_df_long,Sample,into =c("Strain", "Treatment" , "Replicate"),sep = "_",remove = FALSE,extra = "merge")
pred_df_long$Treatment=factor(pred_df_long$Treatment, levels=c("native", "antibiotic", "pathogenic"))
pred_df_long$Temperature=as.numeric(as.character(pred_df_long$Temperature))
pred_df_long$group=paste(pred_df_long$Strain, pred_df_long$Treatment)
#calculate mean ED50 per treatment
ED50_means<-ed50_df %>% 
  group_by(Strain, Treatment) %>%
  summarise(mean=mean(ED50), sd=sd(ED50)) %>%
  unite(Group, c(Treatment), sep = "-", remove = FALSE)
ED50_means$group=paste(ED50_means$Strain, ED50_means$Treatment)

pred_df_long$meanED50=round(ED50_means$mean[match(pred_df_long$group,ED50_means$group)], 2)

curve_H2_input<-pred_df_long[pred_df_long$Strain=="H2",]

ggplot(curve_H2_input, aes(x=Temperature, y=`Fv/Fm`, group=Sample)) +
  geom_line(aes(color=Treatment)) +
  geom_segment(data=curve_H2_input, aes(x = meanED50, y = 0, xend = meanED50, yend = 0.65,color=Treatment), linetype=3) +
  geom_text(data=curve_H2_input, mapping=aes(x=meanED50, y=0.68, label=meanED50, color=Treatment), size=3, angle=90, hjust=0, check_overlap = T) +
  guides(fill=guide_legend())+
  theme(panel.grid = element_blank()) +
  theme_classic() +
  scale_color_manual(values=treatment_colours)+
  scale_x_continuous(breaks=c(30, 33, 36, 39), limits = c(30,41), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,0.75), expand = c(0, 0)) + labs(color='') 

curve_FOO3_input<-pred_df_long[pred_df_long$Strain=="FOO3",]

ggplot(curve_FOO3_input, aes(x=Temperature, y=`Fv/Fm`, group=Sample)) +
  geom_line(aes(color=Treatment)) +
  geom_segment(data=curve_FOO3_input, aes(x = meanED50, y = 0, xend = meanED50, yend = 0.65,color=Treatment), linetype=3) +
  geom_text(data=curve_FOO3_input, mapping=aes(x=meanED50, y=0.68, label=meanED50, color=Treatment), size=3, angle=90, hjust=0, check_overlap = T) +
  guides(fill=guide_legend())+
  theme(panel.grid = element_blank()) +
  theme_classic() +
  scale_color_manual(values=treatment_colours)+
  scale_x_continuous(breaks=c(30, 33, 36, 39), limits = c(30,41), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,0.75), expand = c(0, 0)) + labs(color='')



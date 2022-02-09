library(reshape2)
library(ggplot2)
library(dplyr)
library(emmeans)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/VTK/VTK_2022/")

### read the table
ct_vals <- read.table("./Input_data/VTK2022_qPCR.txt", header = T)

#rename 16S to bacteria and b-actin to Aiptasia
ct_vals$Gene <-gsub("16S", "Bacteria",ct_vals$Gene)
ct_vals$Gene <-gsub("b-actin", "Aptasia",ct_vals$Gene)

############################################
### Average the Ct values for replicates ###
############################################
# re-arrange the table by averaging ct values by gene (16S and b-actin)
ct_vals_wide <-reshape2::dcast(formula = Sample_name~Gene, fun.aggregate = mean,data = ct_vals)

ct_vals_wide2 <- ct_vals %>% group_by(Sample_name, Gene) %>% summarise(Mean=mean(Ct)) %>% 
  spread(key=Gene, value=Mean)

# remove the rows with uncomplete values (NAs)
ct_vals_wide <-ct_vals_wide[complete.cases(ct_vals_wide), ]

#############################
### calculate ∆ct values ###
############################
#break sample names into sample variables
ct_vals_wide_samp=tidyr::separate(ct_vals_wide,Sample_name,into =c("Strain",  "Temperature", "Treatment","Replicate" ),
                                  sep = "_",remove = FALSE,extra = "merge")

#rename treatments
ct_vals_wide_samp$Treatment <-gsub("N", "Native",ct_vals_wide_samp$Treatment)
ct_vals_wide_samp$Treatment <-gsub("A", "Antibiotic",ct_vals_wide_samp$Treatment)
ct_vals_wide_samp$Treatment <-gsub("P", "Pathogen",ct_vals_wide_samp$Treatment)

#order treatmens in the way we want them on the plot
ct_vals_wide_samp$Treatment <- factor(ct_vals_wide_samp$Treatment, levels = c("Native", "Antibiotic", "Pathogen"))


#calculate ∆ct values
ct_vals_wide_samp$dCT <- ct_vals_wide_samp$Bacteria - ct_vals_wide_samp$Aptasia

#split tables per strain
H2=ct_vals_wide_samp[ct_vals_wide_samp$Strain == "H2",]
FOO3=ct_vals_wide_samp[ct_vals_wide_samp$Strain == "FOO3",]




#############################
### calculate ∆∆ct values ###
############################

############################# H2 #############################
#calculate the average dct values for our reference == Native at 30ºC
mean_native_30C_H2 <- mean(H2[H2$Treatment == "Native" &  H2$Temperature=="30",]$dCT)

#calculate ∆∆ct values
H2$ddCT <- H2$dCT - mean_native_30C_H2

#change signs for a more intuitive representation on the plots
H2$ddCT<- -(H2$ddCT)

#calculate mean and standard error per group
summary_h2 <- H2 %>% group_by (Temperature, Treatment) %>% 
  summarise(Abu_ratio_mean=mean(ddCT), Abu_ratio_se=sd(ddCT)/sqrt(length(ddCT)))


pdf("Outputs/qPCR_bacAbundance_H2.pdf" , height = 4, width = 7, pointsize = 12 )
summary_h2 %>% ggplot(aes(x = Temperature, y = Abu_ratio_mean)) +
  geom_col(width=0.5) +
  geom_errorbar(aes(ymin=Abu_ratio_mean-Abu_ratio_se, 
                    ymax=Abu_ratio_mean+Abu_ratio_se),width=.2, 
                position=position_dodge(.9), 
                color="red") +
  facet_grid(~Treatment) + 
  geom_hline(yintercept=mean(summary_h2[summary_h2$Treatment == "Native" & summary_h2$Temperature=="30",]$Abu_ratio_mean),
             linetype="dashed", color = "red") +
  theme_bw() + labs(x= "Temperature (ºC)", y= "-ddct values")
  dev.off()

############################# FOO3 #############################
#calculate the average dct values for our reference == Native at 30ºC
mean_native_30C_FOO3 <- mean(FOO3[FOO3$Treatment == "Native" &  FOO3$Temperature=="30",]$dCT)

#calculate ∆∆ct values
FOO3$ddCT <- FOO3$dCT - mean_native_30C_FOO3

#change signs for a more intuitive representation on the plots
FOO3$ddCT<- -(FOO3$ddCT)

#calculate mean and standard error per group
summary_foo3 <- FOO3 %>% group_by (Temperature, Treatment) %>% 
  summarise(Abu_ratio_mean=mean(ddCT), Abu_ratio_se=sd(ddCT)/sqrt(length(ddCT)))

pdf("Outputs/qPCR_bacAbundance_FOO3.pdf" , height = 4, width = 7, pointsize = 12 )
summary_foo3 %>% ggplot(aes(x = Temperature, y = Abu_ratio_mean)) +
  geom_col(width=0.5)+
  geom_errorbar(aes(ymin=Abu_ratio_mean-Abu_ratio_se, 
                    ymax=Abu_ratio_mean+Abu_ratio_se),width=.2, 
                position=position_dodge(.9), 
                color="red") +
  facet_grid(~Treatment) + 
  geom_hline(yintercept=mean(summary_foo3[summary_foo3$Treatment == "Native" & summary_foo3$Temperature=="30",]$Abu_ratio_mean),
             linetype="dashed", color = "red") +
  theme_bw() + labs(x= "Temperature (ºC)", y= "-ddct values")
dev.off()


###################
### statistics ###
###################
plot(density(H2$ddCT))
plot(density(FOO3$ddCT))
shapiro.test(H2$ddCT)
shapiro.test(FOO3$ddCT)

####### H2

lm_h2=lm(ddCT ~ Treatment*Temperature , data = H2)
summary(lm_h2)
anova(lm_h2)

#residuals
plot(lm_h2$fitted.values,rstudent(lm_h2),
     main="Residuals",
     xlab="", ylab="Residuals")
abline(h=0,lty=2)

#Histogram of Residuals
hist(lm_h2$residuals, main="Histogram of Residuals",ylab="Residuals")
shapiro.test(lm_h2$residuals)

#Q-Q Plot
qqnorm(lm_h2$residuals)
qqline(lm_h2$residuals)

#test for pairwise differences
pair_h2 = emmeans(lm_h2, pairwise ~ Treatment*Temperature, 
                     weights = "proportional", adjust="none", pbkrtest.limit = 10000)
rbind(pair_h2$contrasts, adjust="fdr")


##### FOO3 #####
lm_foo3=lm(ddCT ~ Treatment*Temperature , data = FOO3)
summary(lm_foo3)
anova(lm_foo3)
#residuals
plot(lm_foo3$fitted.values,rstudent(lm_foo3),
     main="Residuals",
     xlab="", ylab="Residuals")
abline(h=0,lty=2)

#Histogram of Residuals
hist(lm_foo3$residuals, main="Histogram of Residuals",
     ylab="Residuals")
shapiro.test(lm_foo3$residuals)

#Q-Q Plot
qqnorm(lm_foo3$residuals)
qqline(lm_foo3$residuals)

#test for pairwise differences
pair_foo3 = emmeans(lm_foo3, pairwise ~ Treatment*Temperature, 
                  weights = "proportional", adjust="none", pbkrtest.limit = 10000)
rbind(pair_foo3$contrasts, adjust="fdr")


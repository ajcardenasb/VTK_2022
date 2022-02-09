
### read the table

#rename 16S to "Bacteria" and b-actin to "Aiptasia"

############################################
### Average the Ct values for replicates ###
############################################
# re-arrange the table by averaging ct values by gene (16S and b-actin)

# remove the rows with uncomplete values (NAs)

#############################
### calculate ∆ct values ###
############################
#break sample names into sample variables

#rename treatments (A, N, P)

#order treatmens in the way we want them on the plot (native first)

#calculate ∆ct values

#split tables per strain

#############################
### calculate ∆∆ct values ###
############################

############################# H2 #############################

#calculate the average dct values for our reference == Native at 30ºC

#calculate ∆∆ct values

#change signs for a more intuitive representation on the plots (-x)

#calculate ddCT mean and standard error per group

summary_h2 <- ??? %>% group_by (???, ???) %>% 
  summarise(Abu_ratio_mean=mean(???), Abu_ratio_se=sd(???)/sqrt(length(???)))


#################################################################
### plot bacterial ratios across treatments and temperatures ###
#################################################################

summary_h2 %>% ggplot(aes(x = Temperature, y = Abu_ratio_mean)) +
  geom_col(width=0.5)+
  geom_errorbar(aes(ymin=Abu_ratio_mean-Abu_ratio_se, 
                    ymax=Abu_ratio_mean+Abu_ratio_se),width=.2, 
                position=position_dodge(.9), 
                color="red") +
  facet_grid(~Treatment) + 
  geom_hline(yintercept=mean(summary_h2[summary_h2$Treatment == "Native" & summary_h2$Temperature=="30",]$Abu_ratio_mean),
             linetype="dashed", color = "red") +
  theme_bw() + labs(x= "Temperature (ºC)", y= "-∆∆ct values")


############################# FOO3 #############################

#############################
### calculate ∆∆ct values ###
############################
#calculate the average dct values for our reference == Native at 30ºC

#calculate ∆∆ct values

#change signs for a more intuitive representation on the plots

#calculate mean and standard error per group
summary_foo3 <- ??? %>% group_by (???, ???) %>% 
  summarise(Abu_ratio_mean=mean(???), Abu_ratio_se=sd(???)/sqrt(length(???)))


#################################################################
### plot bacterial ratios across treatments and temperatures ###
#################################################################

summary_foo3 %>% ggplot(aes(x = Temperature, y = Abu_ratio_mean)) +
  geom_col(width=0.5)+
  geom_errorbar(aes(ymin=Abu_ratio_mean-Abu_ratio_se, 
                    ymax=Abu_ratio_mean+Abu_ratio_se),width=.2, 
                position=position_dodge(.9), 
                color="red") +
  facet_grid(~Treatment) + 
  geom_hline(yintercept=mean(summary_foo3[summary_foo3$Treatment == "Native" & summary_foo3$Temperature=="30",]$Abu_ratio_mean),
             linetype="dashed", color = "red") +
  theme_bw() + labs(x= "Temperature (ºC)", y= "-∆∆ct values")


###################
### statistics ###
###################
library(emmeans)
lm_h2=lm(ddCT ~ Treatment*Temperature , data = H2)
summary(lm_h2)

#residuals
plot(lm_h2$fitted.values,rstudent(lm_h2),
     main="Residuals",
     xlab="", ylab="Residuals")
abline(h=0,lty=2)

#Histogram of Residuals
hist(lm_h2$residuals, main="Histogram of Residuals",
     ylab="Residuals")
shapiro.test(lm_h2$resid)

#Q-Q Plot
qqnorm(lm_h2$resid)
qqline(lm_h2$resid)

#test for differences
anova(lm_h2)

pair_h2 = emmeans(lm_h2, pairwise ~ Treatment*Temperature, 
                  weights = "proportional", adjust="none", pbkrtest.limit = 10000)
rbind(pair_h2$contrasts, adjust="fdr")


##### FOO3 #####
lm_foo3=lm(ddCT ~ Treatment*Temperature , data = FOO3)
summary(lm_foo3)

#residuals
plot(lm_foo3$fitted.values,rstudent(lm_foo3),
     main="Residuals",
     xlab="", ylab="Residuals")
abline(h=0,lty=2)

#Histogram of Residuals
hist(lm_foo3$residuals, main="Histogram of Residuals",
     ylab="Residuals")
shapiro.test(lm_foo3$resid)

#Q-Q Plot
qqnorm(lm_foo3$resid)
qqline(lm_foo3$resid)

#test for differences
anova(lm_foo3)
pair_foo3 = emmeans(lm_foo3, pairwise ~ Treatment*Temperature, 
                    weights = "proportional", adjust="none", pbkrtest.limit = 10000)
rbind(pair_foo3$contrasts, adjust="fdr")

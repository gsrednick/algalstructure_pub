#### Spring 2019 - Algal structure work - extended from CSUN M.S. work ####
# Clear Workspace 
#rm(list=ls())

library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)

# Algal volume and surface area were calculated as the values predicted from the power equations when applied to mean height and multiplied by density for each species on each transect, and then summed across all species on the transect.

#### (1) Part 1 - Building models ####
#### (2) Data upload ####
setwd("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/")
alg <- read.csv("algalstructure_R.csv")


#### (3) Data Wrangling ####
# Height
MAPY_H <- subset(alg, subset = Species == "MAPY", select = Height)
SAHO_H <- subset(alg, subset = Species == "SAHO", select = Height)
SAPA_H <- subset(alg, subset = Species == "SAPA", select = Height)
SAMU_H <- subset(alg, subset = Species == "SAMU", select = Height)
STOS_H <- subset(alg, subset = Species == "STOS", select = Height)

# Volume
MAPY_V <- subset(alg, subset = Species == "MAPY", select = Vol)
SAHO_V <- subset(alg, subset = Species == "SAHO", select = Vol)
SAPA_V <- subset(alg, subset = Species == "SAPA", select = Vol)
SAMU_V <- subset(alg, subset = Species == "SAMU", select = Vol)
STOS_V <- subset(alg, subset = Species == "STOS", select = Vol)

# Surface area
MAPY_SA <- subset(alg, subset = Species == "MAPY", select = SA)
SAHO_SA <- subset(alg, subset = Species == "SAHO", select = SA)
SAPA_SA <- subset(alg, subset = Species == "SAPA", select = SA)
SAMU_SA <- subset(alg, subset = Species == "SAMU", select = SA)
STOS_SA <- subset(alg, subset = Species == "STOS", select = SA)


#### (4) Step 1: Run this for data wrangling - All together ####
MAPY<-subset(alg, subset = Species == "MAPY")
SAHO<-subset(alg, subset = Species == "SAHO")
SAPA<-subset(alg, subset = Species == "SAPA")
SAMU<-subset(alg, subset = Species == "SAMU")
STOS<-subset(alg, subset = Species == "STOS")



#### Satisfy assumptions ####
# Linear distribution
qqnorm(MAPY_H$Height)  # NO
qqnorm(SAHO_H$Height)  # NO
qqnorm(SAPA_H$Height)  # NO
qqnorm(SAMU_H$Height)  # NO
qqnorm(STOS_H$Height)  # Yes, ish - maybe linear fit

qqnorm(MAPY_V$Vol)  # NO
qqnorm(SAHO_V$Vol)  # NO
qqnorm(SAPA_V$Vol)  # NO
qqnorm(SAMU_V$Vol)  # NO
qqnorm(STOS_V$Vol)  # NO

qqnorm(MAPY_SA$SA)  # NO
qqnorm(SAHO_SA$SA)  # NO
qqnorm(SAPA_SA$SA)  # NO
qqnorm(SAMU_SA$SA)  # NO
qqnorm(STOS_SA$SA)  # NO

# Homogenity of variances
hov(Height~Species, data=alg)   # Satisfied
hov(Vol~Species, data=alg)      # Satisfied
hov(SA~Species, data=alg)       # Satisfied


#### ANOVA - Do metrics differ among species ####
height_anova <- aov(Height~Species,data=alg)   # Yes
summary(height_anova)      

vol_anova <- aov(Vol~Species,data=alg)   # Yes
summary(vol_anova)      

SA_anova <- aov(SA~Species,data=alg)   # Yes
summary(vol_anova)      


#### Initial plot data ####

# Height v. Vol

# M. pyrifera
ggplot(subset(alg,Species %in% "MAPY")) + 
  geom_point(aes(Height,Vol,group=Species)) +
  geom_smooth(aes(Height,Vol,group=Species),method = "lm")

# S. horneri
ggplot(subset(alg,Species %in% "SAHO")) + 
  geom_point(aes(Height,Vol,group=Species)) +
  geom_smooth(aes(Height,Vol,group=Species),method = "lm")

# S. palmeri
ggplot(subset(alg,Species %in% "SAPA")) + 
  geom_point(aes(Height,Vol,group=Species)) +
  geom_smooth(aes(Height,Vol,group=Species),method = "lm")

# S. muticum
ggplot(subset(alg,Species %in% "SAMU")) + 
  geom_point(aes(Height,Vol,group=Species)) +
  geom_smooth(aes(Height,Vol,group=Species),method = "lm")

# S. osmundacea ---> Must remove high volume (300 mL sample)
ggplot(subset(alg,Species %in% "STOS")) + 
  geom_point(aes(Height,Vol,group=Species)) +
  geom_smooth(aes(Height,Vol,group=Species),method = "lm")

#
# Height v. Sfc

# M. pyrifera
ggplot(subset(alg,Species %in% "MAPY")) + 
  geom_point(aes(Height,SA,group=Species)) +
  geom_smooth(aes(Height,SA,group=Species),method = "lm")

# S. horneri
ggplot(subset(alg,Species %in% "SAHO")) + 
  geom_point(aes(Height,SA,group=Species)) +
  geom_smooth(aes(Height,SA,group=Species),method = "lm")

# S. palmeri
ggplot(subset(alg,Species %in% "SAPA")) + 
  geom_point(aes(Height,SA,group=Species)) +
  geom_smooth(aes(Height,SA,group=Species),method = "lm")

# S. muticum
ggplot(subset(alg,Species %in% "SAMU")) + 
  geom_point(aes(Height,SA,group=Species)) +
  geom_smooth(aes(Height,SA,group=Species),method = "lm")

# S. osmundacea
ggplot(subset(alg,Species %in% "STOS")) + 
  geom_point(aes(Height,SA,group=Species)) +
  geom_smooth(aes(Height,SA,group=Species),method = "lm")


#### (5) Part 2 - Actual Models ####
# coefficents used for start
# MAPY
lm(log(Vol)~log(Height),data=subset(alg, Species=="MAPY")) 
lm(log(SA)~log(Height),data=subset(alg, Species=="MAPY")) 

# SAHO
lm(log(Vol)~log(Height),data=SAHO) 
lm(log(SA)~log(Height),data=SAHO)

# SAPA
lm(log(Vol)~log(Height),data=SAPA) 
lm(log(SA)~log(Height),data=SAPA)

# STOS
lm(log(Vol)~log(Height),data=STOS) 
lm(log(SA)~log(Height),data=STOS)

## power model -- coefficients from section above  
# MAPY
MAPY_vol<- nls(Vol ~ a * (Height)^b, start=list(a=3.007,b=0.631), data= MAPY)
MAPY_sfc<- nls(SA ~ a * Height^b, start=list(a=5.8353,b=0.3443), data= MAPY)

# SAHO
SAHO_vol<- nls(Vol ~ a * Height^b, start=list(a=-3.682,b=1.851), data= SAHO)
SAHO_sfc<- nls(SA ~ a * Height^b, start=list(a=2.711,b=1.013), data= SAHO)

# SAPA
SAPA_vol<- nls(Vol ~ a * Height^b, start=list(a=0.7099,b=1.3142), data= SAPA)
SAPA_sfc<- nls(SA ~ a * Height^b, start=list(a=2.213,b=1.462), data= SAPA)

# STOS
STOS_vol<- nls(Vol ~ a * Height^b, start=list(a=2.7161,b=0.3443), data= STOS)
STOS_sfc<- nls(SA ~ a * Height^b, start=list(a=6.0245,b=0.1717), data= STOS)

## Summaries
# MAPY
summary(MAPY_vol)
summary(MAPY_sfc)

# SAHO
summary(SAHO_vol)
summary(SAHO_sfc)

# SAPA
summary(SAPA_vol)
summary(SAPA_sfc)

# STOS
summary(STOS_vol)
summary(STOS_sfc)

## Goodness of fit
# MAPY
MAPY_VOL_GOF<-cor(MAPY$Vol,predict(MAPY_vol))
MAPY_SFC_GOF<-cor(MAPY$SA,predict(MAPY_sfc))

# SAHO
SAHO_VOL_GOF<-cor(SAHO$Vol,predict(SAHO_vol))
SAHO_SFC_GOF<-cor(SAHO$SA,predict(SAHO_sfc))

# SAPA
SAPA_VOL_GOF<-cor(SAPA$Vol,predict(SAPA_vol))
SAPA_SFC_GOF<-cor(SAPA$SA,predict(SAPA_sfc))

# STOS
STOS_VOL_GOF<-cor(STOS$Vol,predict(STOS_vol))
STOS_SFC_GOF<-cor(STOS$SA,predict(STOS_sfc))

models_vol_r2<-data.frame(rbind(MAPY_VOL_GOF,SAHO_VOL_GOF,SAPA_VOL_GOF,STOS_VOL_GOF))
models_sfc_r2<-data.frame(rbind(MAPY_SFC_GOF,SAHO_SFC_GOF,SAPA_SFC_GOF,STOS_SFC_GOF))

## plots 
# MAPY
MAPY_V_plot<-ggplot(MAPY, aes(x = Height,y = Vol)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_continuous(limits = c(0,2000)) +
  geom_smooth(method = 'nls', se = FALSE,
              method.args = list(formula = y ~ a * (x^b), 
                                 start = list(a=3.01,b=0.63)))

MAPY_SA_plot<-ggplot(MAPY, aes(x = Height,y = SA)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_continuous(limits = c(0,8000)) +
  geom_smooth(method = 'nls', se = FALSE,
              method.args = list(formula = y ~ a * (x^b), 
                                 start = list(a=5.8353,b=0.3443)))

# new per-reviewer
MAPY_SAV_plot<-ggplot(MAPY, aes(x = SA,y = Vol)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_continuous(limits = c(0,3000)) +
  scale_x_continuous(limits = c(0,6000)) +
  ggtitle('MAPY')



# SAHO
SAHO_V_plot<-ggplot(SAHO, aes(x = Height,y = Vol)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_continuous(limits = c(0,2000)) +
  geom_smooth(method = 'nls', se = FALSE,
              method.args = list(formula = y ~ a * (x^b), 
                                 start = list(a=-3.682,b=1.851)))

SAHO_SA_plot<-ggplot(SAHO, aes(x = Height,y = SA)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_continuous(limits = c(0,8000)) +
  geom_smooth(method = 'nls', se = FALSE,
              method.args = list(formula = y ~ a * (x^b), 
                                 start = list(a=2.711,b=1.013)))

# new per-reviewer
SAHO_SAV_plot<-ggplot(SAHO, aes(x = SA,y = Vol)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title.y=element_blank()) +
  scale_y_continuous(limits = c(0,1500)) +
  ggtitle('SAHO')



# SAPA
SAPA_V_plot<-ggplot(SAPA, aes(x = Height,y = Vol)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_continuous(limits = c(0,2000)) +
  geom_smooth(method = 'nls', se = FALSE,
              method.args = list(formula = y ~ a * (x^b), 
                                 start = list(a=0.7099,b=1.3142)))

SAPA_SA_plot<-ggplot(SAPA, aes(x = Height,y = SA)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_continuous(limits = c(0,8000)) +
  geom_smooth(method = 'nls', se = FALSE,
              method.args = list(formula = y ~ a * (x^b), 
                                 start = list(a=2.213,b=1.462)))

# new per-reviewer
SAPA_SAV_plot<-ggplot(SAPA, aes(x = SA,y = Vol)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title.y=element_blank()) +
  scale_y_continuous(limits = c(0,600)) +
  ggtitle('SAPA')



# STOS
STOS_V_plot<-ggplot(STOS, aes(x = Height,y = Vol)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_continuous(limits = c(0,2000)) +
  geom_smooth(method = 'nls', se = FALSE,
              method.args = list(formula = y ~ a * (x^b), 
                                 start = list(a=2.7161,b=0.3443)))

STOS_SA_plot<-ggplot(STOS, aes(x = Height,y = SA)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_continuous(limits = c(0,8000)) +
  geom_smooth(method = 'nls', se = FALSE,
              method.args = list(formula = y ~ a * (x^b), 
                                 start = list(a=6.0245,b=0.1717)))

# new per-reviewer
STOS_SAV_plot<-ggplot(STOS, aes(x = SA,y = Vol)) + 
  geom_point() + 
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.y=element_blank()) +
  scale_y_continuous(limits = c(0,400)) +
  ggtitle('STOS')






structure_models<-plot_grid(MAPY_V_plot,SAHO_V_plot,SAPA_V_plot,STOS_V_plot,MAPY_SA_plot,SAHO_SA_plot,SAPA_SA_plot,STOS_SA_plot, labels = c("MAPY","SAHO","SAPA","STNE","MAPY","SAHO","SAPA","STNE"), ncol = 4, nrow = 2, align = "v")

structure_models_SAV <- MAPY_SAV_plot + SAHO_SAV_plot +SAPA_SAV_plot+STOS_SAV_plot+  plot_layout(ncol = 4)


ggsave("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Figures/structure_models_SAV.eps", plot = structure_models_SAV, width = 10, height = 2.7)


#### (6) Predict Volume -- HAVE TO MODIFY COLUMN NAMES TO MATCH MODELS --- looks to be something weird with Mapy ####
field <- read.csv("algal_height_field.csv")

# MAPY
MAPY_Height <- data.frame(Height = field$MAPY_H_plant) ## make a new vector
MAPY_VOL_PREDICT <- data.frame(predict(MAPY_vol, newdata = MAPY_Height))
colnames(MAPY_VOL_PREDICT)[1] <- "MAPY_VOL"


# SAHO
SAHO_Height <- data.frame(Height = field$SAHO_H) ## make a new vector
SAHO_VOL_PREDICT <- data.frame(predict(SAHO_vol, newdata = SAHO_Height))
colnames(SAHO_VOL_PREDICT)[1] <- "SAHO_VOL"

# SAPA
SAPA_Height <- data.frame(Height = field$SAPA_H) ## make a new vector
SAPA_VOL_PREDICT <- data.frame(predict(SAPA_vol, newdata = SAPA_Height))
colnames(SAPA_VOL_PREDICT)[1] <- "SAPA_VOL"

# STOS
STOS_Height <- data.frame(Height = field$STOS_H) ## make a new vector
STOS_VOL_PREDICT <- data.frame(predict(STOS_vol, newdata = STOS_Height))
colnames(STOS_VOL_PREDICT)[1] <- "STOS_VOL"


#### (7) Predict Surface area -- HAVE TO MODIFY COLUMN NAMES TO MATCH MODELS --- looks to be something weird with Mapy ####
# MAPY
MAPY_Height <- data.frame(Height = field$MAPY_H_plant) ## make a new vector
MAPY_sfc_PREDICT <- data.frame(predict(MAPY_sfc, newdata = MAPY_Height))
colnames(MAPY_sfc_PREDICT)[1] <- "MAPY_SA"

# SAHO
SAHO_Height <- data.frame(Height = field$SAHO_H) ## make a new vector
SAHO_sfc_PREDICT <- data.frame(predict(SAHO_sfc, newdata = SAHO_Height))
colnames(SAHO_sfc_PREDICT)[1] <- "SAHO_SA"

# SAPA
SAPA_Height <- data.frame(Height = field$SAPA_H) ## make a new vector
SAPA_sfc_PREDICT <- data.frame(predict(SAPA_sfc, newdata = SAPA_Height))
colnames(SAPA_sfc_PREDICT)[1] <- "SAPA_SA"

# STOS
STOS_Height <- data.frame(Height = field$STOS_H) ## make a new vector
STOS_sfc_PREDICT <- data.frame(predict(STOS_sfc, newdata = STOS_Height))
colnames(STOS_sfc_PREDICT)[1] <- "STOS_SA"


## Write to new csv
`algal_volume`<-cbind(field$Season,field$Site,field$Tran,MAPY_VOL_PREDICT,SAHO_VOL_PREDICT,SAPA_VOL_PREDICT,STOS_VOL_PREDICT)
algal_sfc<-cbind(field$Season,field$Site,field$Tran,MAPY_sfc_PREDICT,SAHO_sfc_PREDICT,SAPA_sfc_PREDICT,STOS_sfc_PREDICT)


#### (8) Step 4 -- USE ME INSTEAD OF STEP 3 -- run this to build a new dataframe to run analyses with  ####
#setwd("~/Documents/Research/M.S. Thesis/Catalina Data/Data/Main Setup & Analysis/R/Ecology/algal structure - new models")
field <- read.csv("algal_height_field.csv")

colnames(MAPY_Height)[1] <- "MAPY_H"
colnames(SAHO_Height)[1] <- "SAHO_H"
colnames(SAPA_Height)[1] <- "SAPA_H"
colnames(STOS_Height)[1] <- "STOS_H"

algal_all_structure<-cbind(field$Season,field$Site,field$Tran,MAPY_Height,SAHO_Height,SAPA_Height,STOS_Height,MAPY_VOL_PREDICT,SAHO_VOL_PREDICT,SAPA_VOL_PREDICT,STOS_VOL_PREDICT,MAPY_sfc_PREDICT,SAHO_sfc_PREDICT,SAPA_sfc_PREDICT,STOS_sfc_PREDICT)

# Sum Volume across species
vol_colnms=c("MAPY_VOL", "SAHO_VOL", "SAPA_VOL", "STOS_VOL")
algal_all_structure$Total_volume<-rowSums(algal_all_structure[,vol_colnms],na.rm = T)

# Sum surface area across species
sa_colnms=c("MAPY_SA", "SAHO_SA", "SAPA_SA", "STOS_SA")
algal_all_structure$Total_SA<-rowSums(algal_all_structure[,sa_colnms],na.rm = T)

# mean height across species
H_colnms=c("MAPY_H", "SAHO_H", "SAPA_H", "STOS_H")
algal_all_structure$Total_H<-rowMeans(algal_all_structure[,H_colnms],na.rm = T)

setwd("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/")

write.csv(algal_all_structure,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/algal_all_structure.csv")

#### (9) Wide to longform -- shooting for converting individual thallus height to SA/VOL, instead of transect averages ####

setwd("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/Data Tidying")
SAHO_WIDE <- read.csv("SAHO_WIDE.csv",check.names=FALSE)
SAPA_WIDE <- read.csv("SAPA_WIDE.csv",check.names=FALSE)
STOS_WIDE <- read.csv("STOS_WIDE.csv",check.names=FALSE)
MAPY_WIDE <- read.csv("MAPY_WIDE.csv",check.names=FALSE)


SAHO_LONG <- SAHO_WIDE %>% gather(Height, Count, c(`3`:`210`),convert = TRUE)
MAPY_LONG <- MAPY_WIDE %>% gather(Height, Count, c(`50`:`900`),convert = TRUE)
SAPA_LONG <- SAPA_WIDE %>% gather(Height, Count, c(`7`:`80`),convert = TRUE)
STOS_LONG <- STOS_WIDE %>% gather(Height, Count, c(`7`:`60`),convert = TRUE)

#### (10) New predictions ####
# (1) Have to predict for individual heights
# (2) Multiply conversions by abundance ("SA_2")
# (3) Spread back to wide format
# (4) Sum up predictions, divide by abundance ("SA_f")
# (5) MUST SUM ALL VOLUME AND SURFACE AREA ACROSS SPECIES 

## Predict Surface area and volume
# MAPY
MAPY_LONG$SA <- NA
MAPY_LONG$VOL <- NA
MAPY_LONG$SA<-predict(MAPY_sfc, newdata = MAPY_LONG) # Step 1 - sfc
MAPY_LONG$VOL<-predict(MAPY_vol, newdata = MAPY_LONG) # Step 1 - vol

MAPY_LONG$SA_2 <- MAPY_LONG$SA * MAPY_LONG$Count # Step 2 - sfc
MAPY_LONG$VOL_2 <- MAPY_LONG$VOL * MAPY_LONG$Count # Step 2 - vol

MAPY_LONG_F<-MAPY_LONG %>%  # Step 4 - both sfc and vol
  dplyr::group_by(SEASON,SITE,TRAN) %>% 
  dplyr::summarize(SA_3 = sum(SA_2), VOL_3 = sum(VOL_2),MAPY_DENS = sum(Count)) %>%
  dplyr::mutate(SA_f = SA_3/MAPY_DENS, VOL_f = VOL_3/MAPY_DENS)

# SAHO
SAHO_LONG$SA <- NA
SAHO_LONG$VOL <- NA
SAHO_LONG$SA<-predict(SAHO_sfc, newdata = SAHO_LONG) # Step 1 - sfc
SAHO_LONG$VOL<-predict(SAHO_vol, newdata = SAHO_LONG) # Step 1 - vol

SAHO_LONG$SA_2 <- SAHO_LONG$SA * SAHO_LONG$Count # Step 2 - sfc
SAHO_LONG$VOL_2 <- SAHO_LONG$VOL * SAHO_LONG$Count # Step 2 - vol

SAHO_LONG_F<-SAHO_LONG %>%  # Step 4 - both sfc and vol
  dplyr::group_by(SEASON,SITE,TRAN) %>% 
  dplyr::summarize(SA_3 = sum(SA_2), VOL_3 = sum(VOL_2),SAHO_DENS = sum(Count)) %>%
  dplyr::mutate(SA_f = SA_3/SAHO_DENS, VOL_f = VOL_3/SAHO_DENS)

# SAPA
SAPA_LONG$SA <- NA
SAPA_LONG$VOL <- NA
SAPA_LONG$SA<-predict(SAPA_sfc, newdata = SAPA_LONG) # Step 1 - sfc
SAPA_LONG$VOL<-predict(SAPA_vol, newdata = SAPA_LONG) # Step 1 - vol

SAPA_LONG$SA_2 <- SAPA_LONG$SA * SAPA_LONG$Count # Step 2 - sfc
SAPA_LONG$VOL_2 <- SAPA_LONG$VOL * SAPA_LONG$Count # Step 2 - vol

SAPA_LONG_F<-SAPA_LONG %>%  # Step 4 - both sfc and vol
  dplyr::group_by(SEASON,SITE,TRAN) %>% 
  dplyr::summarize(SA_3 = sum(SA_2), VOL_3 = sum(VOL_2),SAPA_DENS = sum(Count)) %>%
  dplyr::mutate(SA_f = SA_3/SAPA_DENS, VOL_f = VOL_3/SAPA_DENS)

# STOS
STOS_LONG$SA <- NA
STOS_LONG$VOL <- NA
STOS_LONG$SA<-predict(STOS_sfc, newdata = STOS_LONG) # Step 1 - sfc
STOS_LONG$VOL<-predict(STOS_vol, newdata = STOS_LONG) # Step 1 - vol

STOS_LONG$SA_2 <- STOS_LONG$SA * STOS_LONG$Count # Step 2 - sfc
STOS_LONG$VOL_2 <- STOS_LONG$VOL * STOS_LONG$Count # Step 2 - vol

STOS_LONG_F<-STOS_LONG %>%  # Step 4 - both sfc and vol
  dplyr::group_by(SEASON,SITE,TRAN) %>% 
  dplyr::summarize(SA_3 = sum(SA_2), VOL_3 = sum(VOL_2),STOS_DENS = sum(Count))%>%
  dplyr::mutate(SA_f = SA_3/STOS_DENS, VOL_f = VOL_3/STOS_DENS)

#### (11) bring all together #### 

all<-data.frame(MAPY_LONG_F$SEASON,MAPY_LONG_F$SITE,MAPY_LONG_F$TRAN,MAPY_LONG_F$SA_3,MAPY_LONG_F$VOL_3,MAPY_LONG_F$MAPY_DENS,SAHO_LONG_F$SA_3,SAHO_LONG_F$VOL_3,SAHO_LONG_F$SAHO_DENS,SAPA_LONG_F$SA_3,SAPA_LONG_F$VOL_3,SAPA_LONG_F$SAPA_DENS,STOS_LONG_F$SA_3,STOS_LONG_F$VOL_3,STOS_LONG_F$STOS_DENS)
colnames(all)<-c("SEASON","SITE","TRAN","MAPY_SFA","MAPY_VOL","MAPY_DENS","SAHO_SFA","SAHO_VOL","SAHO_DENS","SAPA_SFA","SAPA_VOL","SAPA_DENS","STOS_SFA","STOS_VOL","STOS_DENS")

# sum across species
SFA_colnms=c("MAPY_SFA", "SAHO_SFA", "SAPA_SFA", "STOS_SFA")
DENS_colnms=c("MAPY_DENS", "SAHO_DENS", "SAPA_DENS", "STOS_DENS")
VOL_colnms=c("MAPY_VOL", "SAHO_VOL", "SAPA_VOL", "STOS_VOL")
H_colnms=c("MAPY_H", "SAHO_H", "SAPA_H", "STOS_H")

all$ALL_DENS<-(rowSums(all[,DENS_colnms],na.rm = T))
all$ALL_SFC<-(rowSums(all[,SFA_colnms],na.rm = T)/(all$ALL_DENS))
all$ALL_VOL<-(rowSums(all[,VOL_colnms],na.rm = T)/(rowSums(all[,DENS_colnms],na.rm = T)))
all$ALL_H<-(rowSums(algal_all_structure[,H_colnms],na.rm = T)/(rowSums(all[,DENS_colnms],na.rm = T)))

all[is.na(all)] <- 0

# put all together
write.csv(all,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/algal_structure_data.csv")


# END #


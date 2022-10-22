#### Analysis code 
# Macroalgal physical structure predicts variation in some attributes of temperate fish assemblages better than macroalgal species composition
# in Marine Biology

# Author script for analyses -- univariate and multivariate

# packages 

library(dotwhisker)
library(tidyverse)
library(broom.mixed)
library(piecewiseSEM)
library(arm)
library(vegan)
library(patchwork)
library(MuMIn)
library(ggExtra)
library(egg)
library(lmerTest)
library(ggrepel)


# data import 
# full MLR datasets with fish
MLR_data<-read.csv("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/MLR_data_updated.csv")

allspp_data<-read.csv("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/fish_main_tran.csv")

# Log transform fish densities
allspp_data_log<-allspp_data %>% 
  mutate_at(c(4:39),funs(log(.+1)))

# algae data
algae<-read.csv("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/algae_density.csv")

MLR_data_ext<-merge(MLR_data,allspp_data_log)

tran_id<-read.csv("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/tran_ID.csv")

MLR_data_mod<-merge(tran_id,MLR_data)
MLR_data_ext_mod<-merge(tran_id,MLR_data_ext)




# transformations 

algae_sqrt<-sqrt(algae+0.5) # for S. horneri, S. palmeri, S. osmundacea
algae_log<-log(algae+1) # for M.pyrifera 

# bring new biomass data; merge with MLR and run
biomass<-read.csv("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/biomass.csv")


# diversity data
divers<-read.csv("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/diversity.csv")


MLR_data_mod<-merge(MLR_data_mod,divers)
MLR_data_mod<-merge(MLR_data_mod,biomass)


### ready to test !!! 






#### Multiple linear regression models ####

## (1) Fish abundance -- no bbgoby ####

FISH2_SP_m1<-lmer(FISH2 ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS, 
                  data = MLR_data_mod,REML=FALSE)
FISH2_SP_m0<-lm(FISH2 ~ Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS,
                data = MLR_data_mod)
FISH2_SP_rand<-anova(FISH2_SP_m1,FISH2_SP_m0)

FISH2_SP_r2<-rsquared(FISH2_SP_m1)

anova(FISH2_SP_m1)

# Height
FISH2_H_m1<-lmer(FISH2 ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + Height, 
                 data = MLR_data_mod,REML=FALSE)
FISH2_H_m0<-lm(FISH2 ~ Site + Habitat_PC1 + Habitat_PC2 + Height,
               data = MLR_data_mod)
FISH2_H_rand<-anova(FISH2_H_m1,FISH2_H_m0)

FISH2_H_r2<-rsquared(FISH2_H_m1)

anova(FISH2_H_m1)

#Volume
FISH2_V_m1<-lmer(FISH2 ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + Volume, 
                 data = MLR_data_mod,REML=FALSE)
FISH2_V_m0<-lm(FISH2 ~ Site + Habitat_PC1 + Habitat_PC2 + Volume,
               data = MLR_data_mod)
FISH2_V_rand<-anova(FISH2_V_m1,FISH2_V_m0)

FISH2_V_r2<-rsquared(FISH2_V_m1)

anova(FISH2_V_m1)


# Surface area
FISH2_SA_m1<-lmer(FISH2 ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + SA, 
                  data = MLR_data_mod,REML=FALSE)

FISH2_SA_m0<-lm(FISH2 ~ Site + Habitat_PC1 + Habitat_PC2 + SA,
                data = MLR_data_mod)
FISH2_SA_rand<-anova(FISH2_SA_m1,FISH2_SA_m0)

FISH2_SA_r2<-rsquared(FISH2_SA_m1)

anova(FISH2_SA_m1)

# Habitat
FISH2_HAB_m1<-lmer(FISH2 ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2, 
                   data = MLR_data_mod,REML=FALSE)
FISH2_HAB_m0<-lm(FISH2 ~ Site + Habitat_PC1 + Habitat_PC2,
                 data = MLR_data_mod)
FISH2_HAB_rand<-anova(FISH2_HAB_m1,FISH2_HAB_m0)

FISH2_HAB_r2<-rsquared(FISH2_HAB_m1)

anova(FISH2_HAB_m1)


# results in matrix
ms_RM_FISH2   <- data.frame(model.sel(FISH2_SP_m1,
                                   FISH2_H_m1,
                                   FISH2_V_m1,
                                   FISH2_SA_m1,
                                   FISH2_HAB_m1))


# clean
FISH2_model_comp<-ms_RM_FISH2[-c(1:12,14)]


# join rsquared
model_ls_FISH2   <- list(FISH2_SP_m1,
                                      FISH2_H_m1,
                                      FISH2_V_m1,
                                      FISH2_SA_m1,
                                      FISH2_HAB_m1)

FISH2_model_r2<-rsquared(model_ls_FISH2)
FISH2_model_r2$model <-c("SP","H","V","SA","HAB")

# join rsquared and model

# standardize coefficients
FISH2_SP_m1_Z<-standardize(FISH2_SP_m1,standardize.y = FALSE, unchanged = NULL)
FISH2_H_m1_Z<-standardize(FISH2_H_m1,standardize.y = FALSE, unchanged = NULL)
FISH2_V_m1_Z<-standardize(FISH2_V_m1,standardize.y = FALSE, unchanged = NULL)
FISH2_SA_m1_Z<-standardize(FISH2_SA_m1,standardize.y = TRUE, unchanged = NULL)
FISH2_HAB_m1_Z<-standardize(FISH2_HAB_m1,standardize.y = TRUE, unchanged = NULL)



# stats in table
FISH2_SP_tidy <-tidy(FISH2_SP_m1_Z)
FISH2_H_tidy <-tidy(FISH2_H_m1_Z)
FISH2_V_tidy <-tidy(FISH2_V_m1_Z)
FISH2_SA_tidy <-tidy(FISH2_SA_m1_Z)
FISH2_HAB_tidy <-tidy(FISH2_HAB_m1_Z)

# all together
FISH2_ME<-rbind(FISH2_H_tidy, 
                FISH2_V_tidy, 
                FISH2_SA_tidy,
                FISH2_SP_tidy,
                FISH2_HAB_tidy)


# remove coefficients for  and intercept
FISH2_GO<-FISH2_ME %>%
  filter(!grepl('Intercept|Habitat|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))


FISH2_print<-FISH2_ME %>%
  filter(!grepl('Intercept|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

#write.csv(FISH2_print,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure/tables/fish2.csv",row.names = FALSE)

# make siginificance column 
FISH2_GO_SIG<-FISH2_GO %>% mutate(SIG = ifelse(p.value < 0.05, ifelse(estimate > 0, "positive","negative"), "no"))



FISH2_GO_SIG <- transform(FISH2_GO_SIG,
                          term=plyr::revalue(term,c("algae_sqrt$HORN_DENS"="S. horneri",
                                                    "algae_sqrt$SAPA_DENS"="S. palmeri",
                                                    "algae_sqrt$STOS_DENS"="S. osmundacea",
                                                    "algae_log$MAPY_DENS"="M. pyrifera",
                                                    "z.Height" = "Height",
                                                    "z.Volume" = "Volume",
                                                    "z.SA" = "Surface area")))


# figure 
FISH2<-dwplot(FISH2_GO_SIG,dot_args = list(size = 5,aes(fill = SIG),colour="black",pch=21), whisker_args = list(size = 1, color = "black")) + 
  scale_fill_manual(values=c(negative="red", 
                             no="white", 
                             positive="blue")) +
  scale_y_discrete(labels=c(expression(paste(italic("S. neglecta"))),
                            expression(paste(italic("M. pyrifera"))),
                            expression(paste(italic("S. palmeri"))),
                            expression(paste(italic("S. horneri"))),
                            "Surface-area", "Volume","Height")) +
  geom_vline(xintercept = 0, 
             colour = "grey60", 
             linetype = 2) +
  theme(text = element_text(size=15,color="black"),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))+
  scale_x_continuous(limits = c(-1,1))

#ggsave("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Figures/FISH2_dot.pdf", plot = FISH2, width = 7, height = 6.5)






## (2) Fish abundance -- no blacksmith or bbgoby ####
FISH3_SP_m1<-lmer(FISH3 ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS, 
                        data = MLR_data_mod,REML=FALSE)
FISH3_SP_m0<-lm(FISH3 ~ Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS,
                           data = MLR_data_mod)
FISH3_SP_rand<-anova(FISH3_SP_m1,FISH3_SP_m0)
anova(FISH3_SP_m1)

# Height
FISH3_H_m1<-lmer(FISH3 ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + Height, 
                      data = MLR_data_mod,REML=FALSE)
FISH3_H_m0<-lm(FISH3 ~ Site + Habitat_PC1 + Habitat_PC2 + Height,
             data = MLR_data_mod)
FISH3_H_rand<-anova(FISH3_H_m1,FISH3_H_m0)
anova(FISH3_H_m1)


#Volume
FISH3_V_m1<-lmer(FISH3 ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + Volume, 
                 data = MLR_data_mod,REML=FALSE)
FISH3_V_m0<-lm(FISH3 ~ Site + Habitat_PC1 + Habitat_PC2 + Volume,
               data = MLR_data_mod)
FISH3_V_rand<-anova(FISH3_V_m1,FISH3_V_m0)
anova(FISH3_V_m1)


# Surface area
FISH3_SA_m1<-lmer(FISH3 ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + SA, 
                 data = MLR_data_mod,REML=FALSE)
FISH3_SA_m0<-lm(FISH3 ~ Site + Habitat_PC1 + Habitat_PC2 + SA,
               data = MLR_data_mod)
FISH3_SA_rand<-anova(FISH3_SA_m1,FISH3_SA_m0)
anova(FISH3_SA_m1)


# Habitat
FISH3_HAB_m1<-lmer(FISH3 ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2, 
                  data = MLR_data_mod,REML=FALSE)
FISH3_HAB_m0<-lm(FISH3 ~ Site + Habitat_PC1 + Habitat_PC2,
                data = MLR_data_mod)
FISH3_HAB_rand<-anova(FISH3_HAB_m1,FISH3_HAB_m0)
anova(FISH3_HAB_m1)


# join rsquared
model_ls_FISH3   <- list(FISH3_SP_m1,
                         FISH3_H_m1,
                         FISH3_V_m1,
                         FISH3_SA_m1,
                         FISH3_HAB_m1)

FISH3_model_r2<-rsquared(model_ls_FISH3)
FISH3_model_r2$model <-c("SP","H","V","SA","HAB")

# join rsquared and model


# results in matrix
ms_RM_FISH3 <- data.frame(model.sel(FISH3_SP_m1,
                                   FISH3_H_m1,
                                   FISH3_V_m1,
                                   FISH3_SA_m1,
                                   FISH3_HAB_m1))

# clean
FISH3_model_comp<-ms_RM_FISH3[-c(1:11,13)]



# standardize coefficients
FISH3_SP_m1_Z<-standardize(FISH3_SP_m1,standardize.y = FALSE, unchanged = NULL)
FISH3_H_m1_Z<-standardize(FISH3_H_m1,standardize.y = FALSE, unchanged = NULL)
FISH3_V_m1_Z<-standardize(FISH3_V_m1,standardize.y = FALSE, unchanged = NULL)
FISH3_SA_m1_Z<-standardize(FISH3_SA_m1,standardize.y = TRUE, unchanged = NULL)
FISH3_HAB_m1_Z<-standardize(FISH3_HAB_m1,standardize.y = TRUE, unchanged = NULL)



# stats in table
FISH3_SP_tidy <-tidy(FISH3_SP_m1_Z)
FISH3_H_tidy <-tidy(FISH3_H_m1_Z)
FISH3_V_tidy <-tidy(FISH3_V_m1_Z)
FISH3_SA_tidy <-tidy(FISH3_SA_m1_Z)
FISH3_HAB_tidy <-tidy(FISH3_HAB_m1_Z)

# all together
FISH3_ME<-rbind(FISH3_H_tidy, 
                  FISH3_V_tidy, 
                  FISH3_SA_tidy,
                FISH3_SP_tidy,
                FISH3_HAB_tidy)


# remove coefficients for  and intercept
FISH3_GO<-FISH3_ME %>%
  filter(!grepl('Intercept|Habitat|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

FISH3_print<-FISH3_ME %>%
  filter(!grepl('Intercept|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

write.csv(FISH3_print,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Tables/fish3.csv",row.names = FALSE)

# make siginificance column 
FISH3_GO_SIG<-FISH3_GO %>% mutate(SIG = ifelse(p.value < 0.05, ifelse(estimate > 0, "positive","negative"), "no"))



FISH3_GO_SIG <- transform(FISH3_GO_SIG,
                          term=plyr::revalue(term,c("algae_sqrt$HORN_DENS"="S. horneri",
                                                    "algae_sqrt$SAPA_DENS"="S. palmeri",
                                                    "algae_sqrt$STOS_DENS"="S. neglecta",
                                                    "algae_log$MAPY_DENS"="M. pyrifera",
                                                    "z.Height" = "Height",
                                                    "z.Volume" = "Volume",
                                                    "z.SA" = "Surface area")))


# figure 
FISH3<-dwplot(FISH3_GO_SIG,dot_args = list(size = 5,aes(fill = SIG),colour="black",pch=21), whisker_args = list(size = 1, color = "black")) + 
  scale_fill_manual(values=c(negative="red", 
                             no="white", 
                             positive="blue")) +
  scale_y_discrete(labels=c(expression(paste(italic("S. neglecta"))),
                            expression(paste(italic("M. pyrifera"))),
                            expression(paste(italic("S. palmeri"))),
                            expression(paste(italic("S. horneri"))),
                            "Surface-area", "Volume","Height")) +
  geom_vline(xintercept = 0, 
             colour = "grey60", 
             linetype = 2) +
  theme(text = element_text(size=15,color="black"),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))+
  scale_x_continuous(limits = c(-1,1))


ggsave("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Figures/FISH3_dot.pdf", plot = FISH3, width = 7, height = 6.5)

















## (3) Species richness ####
SPP_SP_m1<-lmer(H ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS, 
                  data = MLR_data_mod,REML=FALSE)
SPP_SP_m0<-lm(H ~ Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS,
                data = MLR_data_mod)
SPP_SP_rand<-anova(SPP_SP_m1,SPP_SP_m0)
anova(SPP_SP_m1)

# Height
SPP_H_m1<-lmer(H ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + Height, 
                 data = MLR_data_mod,REML=FALSE)
SPP_H_m0<-lm(H ~ Site + Habitat_PC1 + Habitat_PC2 + Height,
               data = MLR_data_mod)
SPP_H_rand<-anova(SPP_H_m1,SPP_H_m0)

anova(SPP_H_m1)

#Volume
SPP_V_m1<-lmer(H ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + Volume, 
                 data = MLR_data_mod,REML=FALSE)
SPP_V_m0<-lm(H ~ Site + Habitat_PC1 + Habitat_PC2 + Volume,
               data = MLR_data_mod)
SPP_V_rand<-anova(SPP_V_m1,SPP_V_m0)
anova(SPP_V_m1)


# Surface area
SPP_SA_m1<-lmer(H ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + SA, 
                  data = MLR_data_mod,REML=FALSE)
SPP_SA_m0<-lm(H ~ Site + Habitat_PC1 + Habitat_PC2 + SA,
                data = MLR_data_mod)
SPP_SA_rand<-anova(SPP_SA_m1,SPP_SA_m0)
anova(SPP_SA_m1)


# Habitat
SPP_HAB_m1<-lmer(H ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2, 
                   data = MLR_data_mod,REML=FALSE)
SPP_HAB_m0<-lm(H ~ Site + Habitat_PC1 + Habitat_PC2,
                 data = MLR_data_mod)
SPP_HAB_rand<-anova(SPP_HAB_m1,SPP_HAB_m0)
anova(SPP_HAB_m1)


# join rsquared
model_ls_SPP   <- list(SPP_SP_m1,
                         SPP_H_m1,
                         SPP_V_m1,
                         SPP_SA_m1,
                         SPP_HAB_m1)

SPP_model_r2<-rsquared(model_ls_SPP)
SPP_model_r2$model <-c("SP","H","V","SA","HAB")



# results in matrix
ms_RM_SPP   <- data.frame(model.sel(SPP_SP_m1,
                                      SPP_H_m1,
                                      SPP_V_m1,
                                      SPP_SA_m1,
                                      SPP_HAB_m1))

# clean
SPP_model_comp<-ms_RM_SPP[-c(1:11,13)]



# standardize coefficients
SPP_SP_m1_Z<-standardize(SPP_SP_m1,standardize.y = FALSE, unchanged = NULL)
SPP_H_m1_Z<-standardize(SPP_H_m1,standardize.y = FALSE, unchanged = NULL)
SPP_V_m1_Z<-standardize(SPP_V_m1,standardize.y = FALSE, unchanged = NULL)
SPP_SA_m1_Z<-standardize(SPP_SA_m1,standardize.y = TRUE, unchanged = NULL)
SPP_HAB_m1_Z<-standardize(SPP_HAB_m1,standardize.y = TRUE, unchanged = NULL)



# stats in table
SPP_SP_tidy <-tidy(SPP_SP_m1_Z)
SPP_H_tidy <-tidy(SPP_H_m1_Z)
SPP_V_tidy <-tidy(SPP_V_m1_Z)
SPP_SA_tidy <-tidy(SPP_SA_m1_Z)
SPP_HAB_tidy <-tidy(SPP_HAB_m1_Z)

# all together
SPP_ME<-rbind(SPP_H_tidy, 
                SPP_V_tidy, 
                SPP_SA_tidy,
              SPP_SP_tidy,
              SPP_HAB_tidy)


# remove coefficients for  and intercept
SPP_GO<-SPP_ME %>%
  filter(!grepl('Intercept|Habitat|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

SPP_print<-SPP_ME %>%
  filter(!grepl('Intercept|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

write.csv(SPP_print,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Tables/spp.csv",row.names = FALSE)

# make siginificance column 
SPP_GO_SIG<-SPP_GO %>% mutate(SIG = ifelse(p.value < 0.05, ifelse(estimate > 0, "positive","negative"), "no"))




SPP_GO_SIG <- transform(SPP_GO_SIG,
                          term=plyr::revalue(term,c("algae_sqrt$HORN_DENS"="S. horneri",
                                                    "algae_sqrt$SAPA_DENS"="S. palmeri",
                                                    "algae_sqrt$STOS_DENS"="S. neglecta",
                                                    "algae_log$MAPY_DENS"="M. pyrifera",
                                                    "z.Height" = "Height",
                                                    "z.Volume" = "Volume",
                                                    "z.SA" = "Surface area")))


# figure 
SPP<-dwplot(SPP_GO_SIG,dot_args = list(size = 5,aes(fill = SIG),colour="black",pch=21), whisker_args = list(size = 1, color = "black")) + 
  scale_fill_manual(values=c(negative="red", 
                             no="white", 
                             positive="blue")) +
  scale_y_discrete(labels=c(expression(paste(italic("S. neglecta"))),
                            expression(paste(italic("M. pyrifera"))),
                            expression(paste(italic("S. palmeri"))),
                            expression(paste(italic("S. horneri"))),
                            "Surface-area", "Volume","Height")) +
  geom_vline(xintercept = 0, 
             colour = "grey60", 
             linetype = 2) +
  theme(text = element_text(size=15,color="black"),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))+
  scale_x_continuous(limits = c(-0.5,0.5))

ggsave("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Figures/SPP_dot.pdf", plot = SPP, width = 7, height = 6.5)













## (4) Water column fishes ####
WC_SP_m1<-lmer(WC ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS, 
                  data = MLR_data_mod,REML=FALSE)
WC_SP_m0<-lm(WC ~ Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS,
                data = MLR_data)
WC_SP_rand<-anova(WC_SP_m1,WC_SP_m0)
anova(WC_SP_m1)

# Height
WC_H_m1<-lmer(WC ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + Height, 
                 data = MLR_data_mod,REML=FALSE)
WC_H_m0<-glm(WC ~ Site + Habitat_PC1 + Habitat_PC2 + Height,
               data = MLR_data_mod)
WC_H_rand<-anova(WC_H_m1,WC_H_m0)
anova(WC_H_m1)


#Volume
WC_V_m1<-lmer(WC ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + Volume, 
                 data = MLR_data_mod,REML=FALSE)
WC_V_m0<-lm(WC ~ Site + Habitat_PC1 + Habitat_PC2 + Volume,
               data = MLR_data_mod)
WC_V_rand<-anova(WC_V_m1,WC_V_m0)
anova(WC_V_m1)


# Surface area
WC_SA_m1<-lmer(WC ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + SA, 
                  data = MLR_data_mod,REML=FALSE)
WC_SA_m0<-lm(WC ~ Site + Habitat_PC1 + Habitat_PC2 + SA,
                data = MLR_data_mod)
WC_SA_rand<-anova(WC_SA_m1,WC_SA_m0)
anova(WC_SA_m1)


# Habitat
WC_HAB_m1<-lmer(WC ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2, 
                   data = MLR_data_mod,REML=FALSE)
WC_HAB_m0<-lm(WC ~ Site + Habitat_PC1 + Habitat_PC2,
                 data = MLR_data_mod)
WC_HAB_rand<-anova(WC_HAB_m1,WC_HAB_m0)
anova(WC_HAB_m1)

# join rsquared
model_ls_WC   <- list(WC_SP_m1,
                         WC_H_m1,
                         WC_V_m1,
                         WC_SA_m1,
                         WC_HAB_m1)

WC_model_r2<-rsquared(model_ls_WC)
WC_model_r2$model <-c("SP","H","V","SA","HAB")


# results in matrix
ms_RM_WC   <- data.frame(model.sel(WC_SP_m1,
                                      WC_H_m1,
                                      WC_V_m1,
                                      WC_SA_m1,
                                      WC_HAB_m1))

# clean
WC_model_comp<-ms_RM_WC[-c(1:11,13)]



# standardize coefficients
WC_SP_m1_Z<-standardize(WC_SP_m1,standardize.y = FALSE, unchanged = NULL)
WC_H_m1_Z<-standardize(WC_H_m1,standardize.y = FALSE, unchanged = NULL)
WC_V_m1_Z<-standardize(WC_V_m1,standardize.y = FALSE, unchanged = NULL)
WC_SA_m1_Z<-standardize(WC_SA_m1,standardize.y = TRUE, unchanged = NULL)
WC_HAB_m1_Z<-standardize(WC_HAB_m1,standardize.y = TRUE, unchanged = NULL)



# stats in table
WC_SP_tidy <-tidy(WC_SP_m1_Z)
WC_H_tidy <-tidy(WC_H_m1_Z)
WC_V_tidy <-tidy(WC_V_m1_Z)
WC_SA_tidy <-tidy(WC_SA_m1_Z)
WC_HAB_tidy <-tidy(WC_HAB_m1_Z)

# all together
WC_ME<-rbind(WC_H_tidy, 
                WC_V_tidy, 
                WC_SA_tidy,
                WC_SP_tidy,
             WC_HAB_tidy)


# remove coefficients for  and intercept
WC_GO<-WC_ME %>%
  filter(!grepl('Intercept|Habitat|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

WC_print<-WC_ME %>%
  filter(!grepl('Intercept|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

write.csv(WC_print,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Tables/WC.csv",row.names = FALSE)


# make siginificance column 
WC_GO_SIG<-WC_GO %>% mutate(SIG = ifelse(p.value < 0.05, ifelse(estimate > 0, "positive","negative"), "no"))



WC_GO_SIG <- transform(WC_GO_SIG,
                          term=plyr::revalue(term,c("algae_sqrt$HORN_DENS"="S. horneri",
                                                    "algae_sqrt$SAPA_DENS"="S. palmeri",
                                                    "algae_sqrt$STOS_DENS"="S. neglecta",
                                                    "algae_log$MAPY_DENS"="M. pyrifera",
                                                    "z.Height" = "Height",
                                                    "z.Volume" = "Volume",
                                                    "z.SA" = "Surface area")))


# figure 
WC<-dwplot(WC_GO_SIG,dot_args = list(size = 5,aes(fill = SIG),colour="black",pch=21), whisker_args = list(size = 1, color = "black")) + 
  scale_fill_manual(values=c(negative="red", 
                             no="white", 
                             positive="blue")) +
  scale_y_discrete(labels=c(expression(paste(italic("S. neglecta"))),
                            expression(paste(italic("M. pyrifera"))),
                            expression(paste(italic("S. palmeri"))),
                            expression(paste(italic("S. horneri"))),
                            "Surface-area", "Volume","Height")) +
  geom_vline(xintercept = 0, 
             colour = "grey60", 
             linetype = 2) +
  theme(text = element_text(size=15,color="black"),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_x_continuous(limits = c(-1,2.5))

ggsave("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Figures/WC_dot.pdf", plot = WC, width = 7, height = 6.5)







## (5) total fish size -- all species ####
SIZE1_SP_m1<-lmer(log_total_biomass ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS, 
                  data = MLR_data_mod,REML=FALSE)
SIZE1_SP_m0<-lm(log_total_biomass ~ Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS,
                data = MLR_data_mod)
SIZE1_SP_rand<-anova(SIZE1_SP_m1,SIZE1_SP_m0)
anova(SIZE1_SP_m1)

# Height
SIZE1_H_m1<-lmer(log_total_biomass ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + Height, 
                 data = MLR_data_mod,REML=FALSE)
SIZE1_H_m0<-lm(log_total_biomass ~ Site + Habitat_PC1 + Habitat_PC2 + Height,
               data = MLR_data_mod)
SIZE1_H_rand<-anova(SIZE1_H_m1,SIZE1_H_m0)
anova(SIZE1_H_m1)


#Volume
SIZE1_V_m1<-lmer(log_total_biomass ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + Volume, 
                 data = MLR_data_mod,REML=FALSE)
SIZE1_V_m0<-lm(log_total_biomass ~ Site + Habitat_PC1 + Habitat_PC2 + Volume,
               data = MLR_data_mod)
SIZE1_V_rand<-anova(SIZE1_V_m1,SIZE1_V_m0)
anova(SIZE1_V_m1)


# Surface area
SIZE1_SA_m1<-lmer(log_total_biomass ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + SA, 
                  data = MLR_data_mod,REML=FALSE)
SIZE1_SA_m0<-lm(log_total_biomass ~ Site + Habitat_PC1 + Habitat_PC2 + SA,
                data = MLR_data_mod)
SIZE1_SA_rand<-anova(SIZE1_SA_m1,SIZE1_SA_m0)
anova(SIZE1_SA_m1)


# Habitat
SIZE1_HAB_m1<-lmer(log_total_biomass ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2, 
                   data = MLR_data_mod,REML=FALSE)
SIZE1_HAB_m0<-lm(log_total_biomass ~ Site + Habitat_PC1 + Habitat_PC2,
                 data = MLR_data_mod)
SIZE1_HAB_rand<-anova(SIZE1_HAB_m1,SIZE1_HAB_m0)
anova(SIZE2_HAB_m1)


# join rsquared
model_ls_SIZE1   <- list(SIZE1_SP_m1,
                         SIZE1_H_m1,
                         SIZE1_V_m1,
                         SIZE1_SA_m1,
                         SIZE1_HAB_m1)

SIZE1_model_r2<-rsquared(model_ls_SIZE1)
SIZE1_model_r2$model <-c("SP","H","V","SA","HAB")



# results in matrix
ms_RM_SIZE1   <- data.frame(model.sel(SIZE1_SP_m1,
                                      SIZE1_H_m1,
                                      SIZE1_V_m1,
                                      SIZE1_SA_m1,
                                      SIZE1_HAB_m1))

# clean
SIZE1_model_comp<-ms_RM_SIZE1[-c(1:11,13)]



# standardize coefficients
SIZE1_SP_m1_Z<-standardize(SIZE1_SP_m1,standardize.y = FALSE, unchanged = NULL)
SIZE1_H_m1_Z<-standardize(SIZE1_H_m1,standardize.y = FALSE, unchanged = NULL)
SIZE1_V_m1_Z<-standardize(SIZE1_V_m1,standardize.y = FALSE, unchanged = NULL)
SIZE1_SA_m1_Z<-standardize(SIZE1_SA_m1,standardize.y = TRUE, unchanged = NULL)
SIZE1_HAB_m1_Z<-standardize(SIZE1_HAB_m1,standardize.y = TRUE, unchanged = NULL)



# stats in table
SIZE1_SP_tidy <-tidy(SIZE1_SP_m1_Z)
SIZE1_H_tidy <-tidy(SIZE1_H_m1_Z)
SIZE1_V_tidy <-tidy(SIZE1_V_m1_Z)
SIZE1_SA_tidy <-tidy(SIZE1_SA_m1_Z)
SIZE1_HAB_tidy <-tidy(SIZE1_HAB_m1_Z)

# all together
SIZE1_ME<-rbind(SIZE1_H_tidy, 
                SIZE1_V_tidy, 
                SIZE1_SA_tidy,
                SIZE1_SP_tidy,
                SIZE1_HAB_tidy)


# remove coefficients for  and intercept
SIZE1_GO<-SIZE1_ME %>%
  filter(!grepl('Intercept|Habitat|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

SIZE1_print<-SIZE1_ME %>%
  filter(!grepl('Intercept|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

write.csv(SIZE1_print,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Tables/size1.csv",row.names = FALSE)


# make siginificance column 
SIZE1_GO_SIG<-SIZE1_GO %>% mutate(SIG = ifelse(p.value < 0.05, ifelse(estimate > 0, "positive","negative"), "no"))




SIZE1_GO_SIG <- transform(SIZE1_GO_SIG,
                          term=plyr::revalue(term,c("algae_sqrt$HORN_DENS"="S. horneri",
                                                    "algae_sqrt$SAPA_DENS"="S. palmeri",
                                                    "algae_sqrt$STOS_DENS"="S. neglecta",
                                                    "algae_log$MAPY_DENS"="M. pyrifera",
                                                    "z.Height" = "Height",
                                                    "z.Volume" = "Volume",
                                                    "z.SA" = "Surface area")))


# figure 
SIZE1<-dwplot(SIZE1_GO_SIG,dot_args = list(size = 5,aes(fill = SIG),colour="black",pch=21), whisker_args = list(size = 1, color = "black")) + 
  scale_fill_manual(values=c(negative="red", 
                             no="white", 
                             positive="blue")) +
  scale_y_discrete(labels=c(expression(paste(italic("S. neglecta"))),
                            expression(paste(italic("M. pyrifera"))),
                            expression(paste(italic("S. palmeri"))),
                            expression(paste(italic("S. horneri"))),
                            "Surface-area", "Volume","Height")) +
  geom_vline(xintercept = 0, 
             colour = "grey60", 
             linetype = 2) +
  theme(text = element_text(size=15,color="black"),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))#+
  scale_x_continuous(limits = c(-1,1))

ggsave("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Figures/SIZE1_dot.pdf", plot = SIZE1, width = 7, height = 6.5)










## (6) fish size -- no blacksmith or bbgoby ####
SIZE2_SP_m1<-lmer(log_biomass_reduced ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS, 
                  data = MLR_data_mod,REML=FALSE)
SIZE2_SP_m0<-lm(log_biomass_reduced ~ Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS,
                data = MLR_data_mod)
SIZE2_SP_rand<-anova(SIZE2_SP_m1,SIZE2_SP_m0)
anova(SIZE2_SP_m1)

# Height
SIZE2_H_m1<-lmer(log_biomass_reduced ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + Height, 
                 data = MLR_data_mod,REML=FALSE)
SIZE2_H_m0<-lm(log_biomass_reduced ~ Site + Habitat_PC1 + Habitat_PC2 + Height,
               data = MLR_data_mod)
SIZE2_H_rand<-anova(SIZE2_H_m1,SIZE2_H_m0)
anova(SIZE2_H_m1)


#Volume
SIZE2_V_m1<-lmer(log_biomass_reduced ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + Volume, 
                 data = MLR_data_mod,REML=FALSE)
SIZE2_V_m0<-lm(log_biomass_reduced ~ Site + Habitat_PC1 + Habitat_PC2 + Volume,
               data = MLR_data_mod)
SIZE2_V_rand<-anova(SIZE2_V_m1,SIZE2_V_m0)
anova(SIZE2_V_m1)


# Surface area
SIZE2_SA_m1<-lmer(log_biomass_reduced ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2 + SA, 
                  data = MLR_data_mod,REML=FALSE)
SIZE2_SA_m0<-lm(log_biomass_reduced ~ Site + Habitat_PC1 + Habitat_PC2 + SA,
                data = MLR_data_mod)
SIZE2_SA_rand<-anova(SIZE2_SA_m1,SIZE2_SA_m0)
anova(SIZE2_SA_m1)


# Habitat
SIZE2_HAB_m1<-lmer(log_biomass_reduced ~ Site + (1|Tran_ID) + Habitat_PC1 + Habitat_PC2, 
                   data = MLR_data_mod,REML=FALSE)
SIZE2_HAB_m0<-lm(log_biomass_reduced ~ Site + Habitat_PC1 + Habitat_PC2,
                 data = MLR_data_mod)
SIZE2_HAB_rand<-anova(SIZE2_HAB_m1,SIZE2_HAB_m0)
anova(SIZE2_HAB_m1)


# join rsquared
model_ls_SIZE2   <- list(SIZE2_SP_m1,
                         SIZE2_H_m1,
                         SIZE2_V_m1,
                         SIZE2_SA_m1,
                         SIZE2_HAB_m1)

SIZE2_model_r2<-rsquared(model_ls_SIZE2)
SIZE2_model_r2$model <-c("SP","H","V","SA","HAB")

# results in matrix
ms_RM_SIZE2   <- data.frame(model.sel(SIZE2_SP_m1,
                                      SIZE2_H_m1,
                                      SIZE2_V_m1,
                                      SIZE2_SA_m1,
                                      SIZE2_HAB_m1))

# clean
SIZE2_model_comp<-ms_RM_SIZE2[-c(1:12,14)]



# standardize coefficients
SIZE2_SP_m1_Z<-standardize(SIZE2_SP_m1,standardize.y = FALSE, unchanged = NULL)
SIZE2_H_m1_Z<-standardize(SIZE2_H_m1,standardize.y = FALSE, unchanged = NULL)
SIZE2_V_m1_Z<-standardize(SIZE2_V_m1,standardize.y = FALSE, unchanged = NULL)
SIZE2_SA_m1_Z<-standardize(SIZE2_SA_m1,standardize.y = TRUE, unchanged = NULL)
SIZE2_HAB_m1_Z<-standardize(SIZE2_HAB_m1,standardize.y = TRUE, unchanged = NULL)



# stats in table
SIZE2_SP_tidy <-tidy(SIZE2_SP_m1_Z)
SIZE2_H_tidy <-tidy(SIZE2_H_m1_Z)
SIZE2_V_tidy <-tidy(SIZE2_V_m1_Z)
SIZE2_SA_tidy <-tidy(SIZE2_SA_m1_Z)
SIZE2_HAB_tidy <-tidy(SIZE2_HAB_m1_Z)

# all together
SIZE2_ME<-rbind(SIZE2_H_tidy, 
                SIZE2_V_tidy, 
                SIZE2_SA_tidy,
                SIZE2_SP_tidy,
                SIZE2_HAB_tidy)


# remove coefficients for  and intercept
SIZE2_GO<-SIZE2_ME %>%
  filter(!grepl('Intercept|Habitat|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

SIZE2_print<-SIZE2_ME %>%
  filter(!grepl('Intercept|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

write.csv(SIZE2_print,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Tables/size2.csv",row.names = FALSE)


# make siginificance column 
SIZE2_GO_SIG<-SIZE2_GO %>% mutate(SIG = ifelse(p.value < 0.05, ifelse(estimate > 0, "positive","negative"), "no"))


SIZE2_GO_SIG <- transform(SIZE2_GO_SIG,
                          term=plyr::revalue(term,c("algae_sqrt$HORN_DENS"="S. horneri",
                                        "algae_sqrt$SAPA_DENS"="S. palmeri",
                                        "algae_sqrt$STOS_DENS"="S. neglecta",
                                        "algae_log$MAPY_DENS"="M. pyrifera",
                                        "z.Height" = "Height",
                                        "z.Volume" = "Volume",
                                        "z.SA" = "Surface area")))

SIZE2_GO_SIG$term <- factor(SIZE2_GO_SIG$term , levels = c("Height","Volume","Surface area","S. horneri","S. palmeri","M. pyrifera","S. neglecta"),
                            labels = c("Height","Volume","Surface area","S. horneri","S. palmeri","M. pyrifera","S. neglecta"))

# figure 
SIZE2<-dwplot(SIZE2_GO_SIG,
              dot_args = list(size = 5,aes(fill = SIG),colour="black",pch=21), 
              whisker_args = list(size = 1, color = "black")) + 
  scale_fill_manual(values=c(negative="red", 
                             no="white", 
                             positive="blue")) +
  scale_y_discrete(labels=c(expression(paste(italic("S. neglecta"))),
                            expression(paste(italic("M. pyrifera"))),
                            expression(paste(italic("S. palmeri"))),
                            expression(paste(italic("S. horneri"))),
                            "Surface-area", "Volume","Height")) +
  geom_vline(xintercept = 0, 
             colour = "grey60", 
             linetype = 2) +
  theme(text = element_text(size=15,color="black"),
        axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2))+
  coord_cartesian(x = c(-1,1))

ggsave("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Figures/SIZE2_dot.pdf", plot = SIZE2, width = 7, height = 6.5)
















#### MLR Figure ####

library(egg)

SPP_group <- SPP + theme(axis.title.y=element_blank(),
                          axis.text.y=element_blank())

SIZE1_group <- SIZE1 + theme(axis.title.y=element_blank(),
                         axis.text.y=element_blank())

SIZE2_group <- SIZE2 + theme(axis.title.y=element_blank(),
                             axis.text.y=element_blank())


dotplot<-egg::ggarrange(tag_facet(FISH3,tag_pool = "(A) Fish density", open = "", close = "",
                                  hjust = 0.01),
                        tag_facet(SPP_group,tag_pool = "(B) Diversity (H')", open = "", close = "",
                                  hjust = 0.01),
                        tag_facet(WC,tag_pool = "(C) Vertical distribution ", open = "", close = "",
                                  hjust = 0.01),
                        tag_facet(SIZE1_group,tag_pool = "(D) Total biomass", open = "", close = "",
                                  hjust = 0.01))

library(ggpubr)

dotplot_acutal<-annotate_figure(dotplot, bottom = text_grob("Standarized coefficient", hjust = 0.1, vjust = -1))



ggsave("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Figures/dotplot_actual.pdf", plot = dotplot_acutal, width = 10, height = 11)










#### Coefficients table ####
FISH2_GO[c(3,8,4)]
FISH3_GO[c(3,8,4)]
SPP_GO[c(3,8,4)]
WC_GO[c(3,8,4)]
SIZE1_GO[c(3,8,4)]
SIZE2_GO[c(3,8,4)]

LYDA_GO[c(3,8,4)]


#### Model Comparison ####

all_models<-rbind(FISH2_model_comp,
                  FISH3_model_comp,
                  SPP_model_comp,
                  WC_model_comp,
                  SIZE1_model_comp,
                  SIZE2_model_comp)










#### Individual Spp #### 


# L. dalli
# P. clathratus
# H. rubicundus
# H. semicinctus
# O. californica



# C. punctipinnis ####

# Species
CHPU_SP<-lmer(BLACKSMITH ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS, 
                  data = MLR_data_ext_mod,REML=FALSE)
anova(CHPU_SP)

# height
CHPU_H<-lmer(BLACKSMITH ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + Height, 
              data = MLR_data_ext_mod,REML=FALSE)
anova(CHPU_H)

# volume
CHPU_V<-lmer(BLACKSMITH ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + Volume, 
             data = MLR_data_ext_mod,REML=FALSE)
anova(CHPU_V)

# sfc
CHPU_SA<-lmer(BLACKSMITH ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + SA, 
             data = MLR_data_ext_mod,REML=FALSE)
anova(CHPU_SA)

# site and habitat
CHPU_HAB<-lmer(BLACKSMITH ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2, 
             data = MLR_data_ext_mod,REML=FALSE)

anova(CHPU_HAB)

# join rsquared
model_ls_CHPU   <- list(CHPU_SP,
                        CHPU_H,
                        CHPU_V,
                        CHPU_SA,
                        CHPU_HAB)

CHPU_model_r2<-rsquared(model_ls_CHPU)
CHPU_model_r2$model <-c("SP","H","V","SA","HAB")

# results in matrix
ms_RM_CHPU   <- data.frame(model.sel(CHPU_SP,
                                      CHPU_H,
                                      CHPU_V,
                                      CHPU_SA,
                                      CHPU_HAB))

# clean
CHPU_model_comp<-ms_RM_CHPU[-c(1:11,13)]


# standardize coefficients
CHPU_SP_Z<-standardize(CHPU_SP,standardize.y = FALSE, unchanged = NULL)
CHPU_H_Z<-standardize(CHPU_H,standardize.y = FALSE, unchanged = NULL)
CHPU_V_Z<-standardize(CHPU_V,standardize.y = FALSE, unchanged = NULL)
CHPU_SA_Z<-standardize(CHPU_SA,standardize.y = TRUE, unchanged = NULL)
CHPU_HAB_Z<-standardize(CHPU_HAB,standardize.y = TRUE, unchanged = NULL)



# stats in table
CHPU_SP_tidy <-tidy(CHPU_SP_Z)
CHPU_H_tidy <-tidy(CHPU_H_Z)
CHPU_V_tidy <-tidy(CHPU_V_Z)
CHPU_SA_tidy <-tidy(CHPU_SA_Z)
CHPU_HAB_tidy <-tidy(CHPU_HAB_Z)

# all together
CHPU_ME<-rbind(CHPU_H_tidy, 
                CHPU_V_tidy, 
                CHPU_SA_tidy,
                CHPU_SP_tidy,
               CHPU_HAB_tidy)


# remove coefficients for  and intercept
CHPU_GO<-CHPU_ME %>%
  filter(!grepl('Intercept|Habitat|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

CHPU_print<-CHPU_ME %>%
  filter(!grepl('Intercept|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

write.csv(CHPU_print,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Tables/CHPU.csv",row.names = FALSE)



# L. dalli ####
# Species
LYDA_SP<-lmer(BLUE.BAND.GOBY ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS, 
              data = MLR_data_ext_mod,REML=FALSE)
anova(LYDA_SP)

# height
LYDA_H<-lmer(BLUE.BAND.GOBY ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + Height, 
             data = MLR_data_ext_mod,REML=FALSE)
anova(LYDA_H)

# volume
LYDA_V<-lmer(BLUE.BAND.GOBY ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + Volume, 
             data = MLR_data_ext_mod,REML=FALSE)
anova(LYDA_V)

# sfc
LYDA_SA<-lmer(BLUE.BAND.GOBY ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + SA, 
              data = MLR_data_ext_mod,REML=FALSE)
anova(LYDA_SA)

# site and habitat
LYDA_HAB<-lmer(BLUE.BAND.GOBY ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2, 
               data = MLR_data_ext_mod,REML=FALSE)
anova(LYDA_HAB)



# join rsquared
model_ls_LYDA   <- list(LYDA_SP,
                        LYDA_H,
                        LYDA_V,
                        LYDA_SA,
                        LYDA_HAB)

LYDA_model_r2<-rsquared(model_ls_LYDA)
LYDA_model_r2$model <-c("SP","H","V","SA","HAB")

# results in matrix
ms_RM_LYDA   <- data.frame(model.sel(LYDA_SP,
                                     LYDA_H,
                                     LYDA_V,
                                     LYDA_SA,
                                     LYDA_HAB))

# clean
LYDA_model_comp<-ms_RM_LYDA[-c(1:11,13)]


# standardize coefficients
LYDA_SP_Z<-standardize(LYDA_SP,standardize.y = FALSE, unchanged = NULL)
LYDA_H_Z<-standardize(LYDA_H,standardize.y = FALSE, unchanged = NULL)
LYDA_V_Z<-standardize(LYDA_V,standardize.y = FALSE, unchanged = NULL)
LYDA_SA_Z<-standardize(LYDA_SA,standardize.y = TRUE, unchanged = NULL)
LYDA_HAB_Z<-standardize(LYDA_HAB,standardize.y = TRUE, unchanged = NULL)



# stats in table
LYDA_SP_tidy <-tidy(LYDA_SP_Z)
LYDA_H_tidy <-tidy(LYDA_H_Z)
LYDA_V_tidy <-tidy(LYDA_V_Z)
LYDA_SA_tidy <-tidy(LYDA_SA_Z)
LYDA_HAB_tidy <-tidy(LYDA_HAB_Z)

# all together
LYDA_ME<-rbind(LYDA_H_tidy, 
               LYDA_V_tidy, 
               LYDA_SA_tidy,
               LYDA_SP_tidy,
               LYDA_HAB_tidy)


# remove coefficients for  and intercept
LYDA_GO<-LYDA_ME %>%
  filter(!grepl('Intercept|Habitat|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

LYDA_print<-LYDA_ME %>%
  filter(!grepl('Intercept|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

write.csv(LYDA_print,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Tables/LYDA.csv",row.names = FALSE)






# P. clathratus ####
# Species
PACL_SP<-lmer(KELP.BASS.COUNT ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS, 
              data = MLR_data_ext_mod,REML=FALSE)
anova(PACL_SP)

# height
PACL_H<-lmer(KELP.BASS.COUNT ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + Height, 
             data = MLR_data_ext_mod,REML=FALSE)
anova(PACL_H)

# volume
PACL_V<-lmer(KELP.BASS.COUNT ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + Volume, 
             data = MLR_data_ext_mod,REML=FALSE)
anova(PACL_V)

# sfc
PACL_SA<-lmer(KELP.BASS.COUNT ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + SA, 
              data = MLR_data_ext_mod,REML=FALSE)
anova(PACL_SA)

# site and habitat
PACL_HAB<-lmer(KELP.BASS.COUNT ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2, 
               data = MLR_data_ext_mod,REML=FALSE)
anova(PACL_HAB)



# join rsquared
model_ls_PACL   <- list(PACL_SP,
                        PACL_H,
                        PACL_V,
                        PACL_SA,
                        PACL_HAB)

PACL_model_r2<-rsquared(model_ls_PACL)
PACL_model_r2$model <-c("SP","H","V","SA","HAB")

# results in matrix
ms_RM_PACL   <- data.frame(model.sel(PACL_SP,
                                     PACL_H,
                                     PACL_V,
                                     PACL_SA,
                                     PACL_HAB))

# clean
PACL_model_comp<-ms_RM_PACL[-c(1:11,13)]


# standardize coefficients
PACL_SP_Z<-standardize(PACL_SP,standardize.y = FALSE, unchanged = NULL)
PACL_H_Z<-standardize(PACL_H,standardize.y = FALSE, unchanged = NULL)
PACL_V_Z<-standardize(PACL_V,standardize.y = FALSE, unchanged = NULL)
PACL_SA_Z<-standardize(PACL_SA,standardize.y = TRUE, unchanged = NULL)
PACL_HAB_Z<-standardize(PACL_HAB,standardize.y = TRUE, unchanged = NULL)



# stats in table
PACL_SP_tidy <-tidy(PACL_SP_Z)
PACL_H_tidy <-tidy(PACL_H_Z)
PACL_V_tidy <-tidy(PACL_V_Z)
PACL_SA_tidy <-tidy(PACL_SA_Z)
PACL_HAB_tidy <-tidy(PACL_HAB_Z)

# all together
PACL_ME<-rbind(PACL_H_tidy, 
               PACL_V_tidy, 
               PACL_SA_tidy,
               PACL_SP_tidy,
               PACL_HAB_tidy)


# remove coefficients for  and intercept
PACL_GO<-PACL_ME %>%
  filter(!grepl('Intercept|Habitat|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

PACL_print<-PACL_ME %>%
  filter(!grepl('Intercept|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

write.csv(PACL_print,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Tables/PACL.csv",row.names = FALSE)






# H. rubicundus ####
# Species
HYRU_SP<-lmer(GARIBALDI ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS, 
              data = MLR_data_ext_mod,REML=FALSE)
anova(HYRU_SP)

# height
HYRU_H<-lmer(GARIBALDI ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + Height, 
             data = MLR_data_ext_mod,REML=FALSE)
anova(HYRU_H)

# volume
HYRU_V<-lmer(GARIBALDI ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + Volume, 
             data = MLR_data_ext_mod,REML=FALSE)
anova(HYRU_V)

# sfc
HYRU_SA<-lmer(GARIBALDI ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + SA, 
              data = MLR_data_ext_mod,REML=FALSE)
anova(HYRU_SA)

# site and habitat
HYRU_HAB<-lmer(GARIBALDI ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2, 
               data = MLR_data_ext_mod,REML=FALSE)
anova(HYRU_HAB)


# join rsquared
model_ls_HYRU   <- list(HYRU_SP,
                        HYRU_H,
                        HYRU_V,
                        HYRU_SA,
                        HYRU_HAB)

HYRU_model_r2<-rsquared(model_ls_HYRU)
HYRU_model_r2$model <-c("SP","H","V","SA","HAB")

# results in matrix
ms_RM_HYRU   <- data.frame(model.sel(HYRU_SP,
                                     HYRU_H,
                                     HYRU_V,
                                     HYRU_SA,
                                     HYRU_HAB))

# clean
HYRU_model_comp<-ms_RM_HYRU[-c(1:11,13)]


# standardize coefficients
HYRU_SP_Z<-standardize(HYRU_SP,standardize.y = FALSE, unchanged = NULL)
HYRU_H_Z<-standardize(HYRU_H,standardize.y = FALSE, unchanged = NULL)
HYRU_V_Z<-standardize(HYRU_V,standardize.y = FALSE, unchanged = NULL)
HYRU_SA_Z<-standardize(HYRU_SA,standardize.y = TRUE, unchanged = NULL)
HYRU_HAB_Z<-standardize(HYRU_HAB,standardize.y = TRUE, unchanged = NULL)



# stats in table
HYRU_SP_tidy <-tidy(HYRU_SP_Z)
HYRU_H_tidy <-tidy(HYRU_H_Z)
HYRU_V_tidy <-tidy(HYRU_V_Z)
HYRU_SA_tidy <-tidy(HYRU_SA_Z)
HYRU_HAB_tidy <-tidy(HYRU_HAB_Z)


# all together
HYRU_ME<-rbind(HYRU_H_tidy, 
               HYRU_V_tidy, 
               HYRU_SA_tidy,
               HYRU_SP_tidy,
               HYRU_HAB_tidy)


# remove coefficients for  and intercept
HYRU_GO<-HYRU_ME %>%
  filter(!grepl('Intercept|Habitat|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))


HYRU_print<-HYRU_ME %>%
  filter(!grepl('Intercept|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

write.csv(HYRU_print,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Tables/HYRU.csv",row.names = FALSE)




# H. semicinctus ####
# Species
HASE_SP<-lmer(ROCK.WRASSE ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS, 
              data = MLR_data_ext_mod,REML=FALSE)
anova(HASE_SP)
# height
HASE_H<-lmer(ROCK.WRASSE ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + Height, 
             data = MLR_data_ext_mod,REML=FALSE)
anova(HASE_H)

# volume
HASE_V<-lmer(ROCK.WRASSE ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + Volume, 
             data = MLR_data_ext_mod,REML=FALSE)
anova(HASE_V)

# sfc
HASE_SA<-lmer(ROCK.WRASSE ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + SA, 
              data = MLR_data_ext_mod,REML=FALSE)
anova(HASE_SA)

# site and habitat
HASE_HAB<-lmer(ROCK.WRASSE ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2, 
               data = MLR_data_ext_mod,REML=FALSE)
anova(HASE_HAB)


# join rsquared
model_ls_HASE   <- list(HASE_SP,
                        HASE_H,
                        HASE_V,
                        HASE_SA,
                        HASE_HAB)

HASE_model_r2<-rsquared(model_ls_HASE)
HASE_model_r2$model <-c("SP","H","V","SA","HAB")

# results in matrix
ms_RM_HASE   <- data.frame(model.sel(HASE_SP,
                                     HASE_H,
                                     HASE_V,
                                     HASE_SA,
                                     HASE_HAB))

# clean
HASE_model_comp<-ms_RM_HASE[-c(1:11,13)]


# standardize coefficients
HASE_SP_Z<-standardize(HASE_SP,standardize.y = FALSE, unchanged = NULL)
HASE_H_Z<-standardize(HASE_H,standardize.y = FALSE, unchanged = NULL)
HASE_V_Z<-standardize(HASE_V,standardize.y = FALSE, unchanged = NULL)
HASE_SA_Z<-standardize(HASE_SA,standardize.y = TRUE, unchanged = NULL)
HASE_HAB_Z<-standardize(HASE_HAB,standardize.y = TRUE, unchanged = NULL)



# stats in table
HASE_SP_tidy <-tidy(HASE_SP_Z)
HASE_H_tidy <-tidy(HASE_H_Z)
HASE_V_tidy <-tidy(HASE_V_Z)
HASE_SA_tidy <-tidy(HASE_SA_Z)
HASE_HAB_tidy <-tidy(HASE_HAB_Z)


# all together
HASE_ME<-rbind(HASE_H_tidy, 
               HASE_V_tidy, 
               HASE_SA_tidy,
               HASE_SP_tidy,
               HASE_HAB_tidy)


# remove coefficients for  and intercept
HASE_GO<-HASE_ME %>%
  filter(!grepl('Intercept|Habitat|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))


HASE_print<-HASE_ME %>%
  filter(!grepl('Intercept|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

write.csv(HASE_print,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Tables/HASE.csv",row.names = FALSE)


# O. californica ####
# Species
OXCA_SP<-lmer(SENORITA.COUNT ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS, 
              data = MLR_data_ext_mod,REML=FALSE)
anova(OXCA_SP)

# height
OXCA_H<-lmer(SENORITA.COUNT ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + Height, 
             data = MLR_data_ext_mod,REML=FALSE)
anova(OXCA_H)

# volume
OXCA_V<-lmer(SENORITA.COUNT ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + Volume, 
             data = MLR_data_ext_mod,REML=FALSE)
anova(OXCA_V)

# sfc
OXCA_SA<-lmer(SENORITA.COUNT ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2 + SA, 
              data = MLR_data_ext_mod,REML=FALSE)
anova(OXCA_SA)

# site and habitat
OXCA_HAB<-lmer(SENORITA.COUNT ~ (1|Tran_ID) + Site + Habitat_PC1 + Habitat_PC2, 
               data = MLR_data_ext_mod,REML=FALSE)
anova(OXCA_HAB)

# join rsquared
model_ls_OXCA   <- list(OXCA_SP,
                        OXCA_H,
                        OXCA_V,
                        OXCA_SA,
                        OXCA_HAB)

OXCA_model_r2<-rsquared(model_ls_OXCA)
OXCA_model_r2$model <-c("SP","H","V","SA","HAB")

# results in matrix
ms_RM_OXCA   <- data.frame(model.sel(OXCA_SP,
                                     OXCA_H,
                                     OXCA_V,
                                     OXCA_SA,
                                     OXCA_HAB))

# clean
OXCA_model_comp<-ms_RM_OXCA[-c(1:11,13)]


# standardize coefficients
OXCA_SP_Z<-standardize(OXCA_SP,standardize.y = FALSE, unchanged = NULL)
OXCA_H_Z<-standardize(OXCA_H,standardize.y = FALSE, unchanged = NULL)
OXCA_V_Z<-standardize(OXCA_V,standardize.y = FALSE, unchanged = NULL)
OXCA_SA_Z<-standardize(OXCA_SA,standardize.y = TRUE, unchanged = NULL)
OXCA_HAB_Z<-standardize(OXCA_HAB,standardize.y = TRUE, unchanged = NULL)



# stats in table
OXCA_SP_tidy <-tidy(OXCA_SP_Z)
OXCA_H_tidy <-tidy(OXCA_H_Z)
OXCA_V_tidy <-tidy(OXCA_V_Z)
OXCA_SA_tidy <-tidy(OXCA_SA_Z)
OXCA_HAB_tidy <-tidy(OXCA_HAB_Z)


# all together
OXCA_ME<-rbind(OXCA_H_tidy, 
               OXCA_V_tidy, 
               OXCA_SA_tidy,
               OXCA_SP_tidy,
               OXCA_HAB_tidy)


# remove coefficients for  and intercept
OXCA_GO<-OXCA_ME %>%
  filter(!grepl('Intercept|Habitat|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))



OXCA_print<-OXCA_ME %>%
  filter(!grepl('Intercept|Site|Tran|Site', term)) %>% 
  filter(!grepl('Residual|Tran_ID', group))

write.csv(OXCA_print,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Tables/OXCA.csv",row.names = FALSE)



# table 
View(CHPU_print[c(3,8,4)])
View(LYDA_print[c(3,8,4)])
View(PACL_print[c(3,8,4)])
View(HASE_print[c(3,8,4)])
View(HYRU_print[c(3,8,4)])
View(OXCA_print[c(3,8,4)])






#### Improved summary (1/15/2021); sensu Burnham and Anderson 2002

ms_CHPU <- data.frame(model.sel(CHPU_SP,
                                CHPU_H,
                                CHPU_V,
                                CHPU_SA,
                                CHPU_HAB, extra = "R2"))

ms_LYDA <- data.frame(model.sel(LYDA_site_spp,
                                LYDA_site_H,
                                LYDA_site_V,
                                LYDA_site_SA,
                                LYDA_site_HAB, extra = "R2"))

ms_PACL <- data.frame(model.sel(PACL_site_spp,
                                PACL_site_H,
                                PACL_site_V,
                                PACL_site_SA,
                                PACL_site_HAB, extra = "R2"))

ms_HASE   <- data.frame(model.sel(HASE_site_spp,
                                  HASE_site_H,
                                  HASE_site_V,
                                  HASE_site_SA,
                                  HASE_site_HAB, extra = "R2"))

ms_HYRU <- data.frame(model.sel(HYRU_site_spp,
                                HYRU_site_H,
                                HYRU_site_V,
                                HYRU_site_SA,
                                HYRU_site_HAB, extra = "R2"))

ms_OXCA <- data.frame(model.sel(OXCA_site_spp,
                                OXCA_site_H,
                                OXCA_site_V,
                                OXCA_site_SA,
                                OXCA_site_HAB, extra = "R2"))

spp_ms<-bind_rows(ms_CHPU,ms_LYDA,ms_PACL,ms_HASE,ms_HYRU,ms_OXCA)

spp_ms_v2<-spp_ms[-c(1:11,13,15)]

spp_ms_v2<-spp_ms[-c(1:11,13,14)]











#### Multivariate methods ####
library(vegan)
library(broom)
library(ggvegan)
library(ggExtra)

#### MLR-type models ####
fish_transformed<-allspp_data_log[c(4:37)]

# do a dbRDA with MLR model strucutre: surface area, height, etc.
overall_site_spp_dbRDA<-capscale(fish_transformed ~ MLR_data_mod$Tran_ID + MLR_data_mod$Season + MLR_data_mod$Site + MLR_data_mod$Habitat_PC1 + MLR_data_mod$Habitat_PC2 + algae_sqrt$HORN_DENS + algae_sqrt$SAPA_DENS + algae_log$MAPY_DENS +algae_sqrt$STOS_DENS,dist = "euc")
overall_site_height_dbRDA<-capscale(fish_transformed ~ MLR_data_mod$Tran_ID + MLR_data_mod$Season + MLR_data_mod$Site + MLR_data_mod$Habitat_PC1 + MLR_data_mod$Habitat_PC2 + MLR_data_mod$Height, data = algae_log,dist = "euc")
overall_site_SA_dbRDA<-capscale(fish_transformed ~ MLR_data_mod$Tran_ID + MLR_data_mod$Season + MLR_data_mod$Site + MLR_data_mod$Habitat_PC1 + MLR_data_mod$Habitat_PC2 + MLR_data_mod$SA , data = algae_log,dist = "euc")
overall_site_volume_dbRDA<-capscale(fish_transformed ~ MLR_data_mod$Tran_ID + MLR_data_mod$Season + MLR_data_mod$Site + MLR_data_mod$Habitat_PC1 + MLR_data_mod$Habitat_PC2 + MLR_data_mod$Volume , data = algae_log,dist = "euc")
overall_site_site_dbRDA<-capscale(fish_transformed ~ MLR_data_mod$Tran_ID + MLR_data_mod$Season + MLR_data_mod$Site + MLR_data_mod$Habitat_PC1 + MLR_data_mod$Habitat_PC2, data = algae_log,dist = "euc")


summary(overall_site_spp_dbRDA)
summary(overall_site_height_dbRDA)
summary(overall_site_SA_dbRDA)
summary(overall_site_volume_dbRDA)
summary(overall_site_site_dbRDA)



anova(overall_site_spp_dbRDA)
anova(overall_site_height_dbRDA)
anova(overall_site_SA_dbRDA)
anova(overall_site_volume_dbRDA)
anova(overall_site_site_dbRDA)


plot(overall_site_spp_dbRDA)
plot(overall_site_height_dbRDA)
plot(overall_site_SA_dbRDA)
plot(overall_site_volume_dbRDA)
plot(overall_site_site_dbRDA)



# anova outputs for models
multi_spp_aov<-anova.cca(overall_site_spp_dbRDA, by="terms", permu=999) # distLM-like result; these do have a single rquared: can be extracted 
multi_height_aov<-anova.cca(overall_site_height_dbRDA, by="terms", permu=999) # distLM-like result
multi_SA_aov<-anova.cca(overall_site_SA_dbRDA, by="terms", permu=999) # distLM-like result
multi_V_aov<-anova.cca(overall_site_volume_dbRDA, by="terms", permu=999) # distLM-like result
multi_site_aov<-anova.cca(overall_site_site_dbRDA, by="terms", permu=999) # distLM-like result

# pull r-squared from models
overall_spp_r<-RsquareAdj(overall_site_spp_dbRDA)$r.squared
overall_h_r<-RsquareAdj(overall_site_height_dbRDA)$r.squared
overall_sa_r<-RsquareAdj(overall_site_SA_dbRDA)$r.squared
overall_v_r<-RsquareAdj(overall_site_volume_dbRDA)$r.squared
overall_site_r<-RsquareAdj(overall_site_site_dbRDA)$r.squared



extractAIC(overall_site_spp_dbRDA)
deviance(overall_site_spp_dbRDA)

overall_site_spp_dbRDA



multi_spp_aov_t<-tidy(multi_spp_aov)
multi_height_aov_t<-tidy(multi_height_aov)
multi_SA_aov_t<-tidy(multi_SA_aov)
multi_V_aov_t<-tidy(multi_V_aov)
multi_site_aov_t<-tidy(multi_site_aov)





##### 1/8/2021 --> Ordination plot ####



# bring in environmental variables
env_variables<-cbind(MLR_data[c(1:3,7:9,16,17)],algae_log[c(1:3,5)])


# envfit them
fish_mds <- metaMDS(comm = fish_transformed, distance = "euclidean", trace = FALSE, autotransform = FALSE)
stress_no<-fish_mds$stress

stressplot(fish_mds) 


ef <- envfit(fish_mds, env_variables, permu = 999)
#plot(fish_mds)
#plot(ef)
summary(ef)

fish_mds_points<-as.data.frame(fish_mds$points)

vectors<-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
vectors$species<-rownames(vectors)



#### MDS for seasonal variation in fish assemblage ####


env_variables<-cbind(MLR_data[c(7:11)],algae_log[c(1:3,5)])
fish_mds_points

community_MDS<-cbind(fish_mds_points,MLR_data[c(1:3,7:11)])


## For plotting means +/- SE on objects
se_fn <- function(x) sd(x) / sqrt(length(x)) # Create own function

community_MDS_summarize <- community_MDS %>% 
  group_by(Season,Site) %>%
  dplyr::summarize(mean_MDS1 = mean(MDS1),
                   mean_MDS2 = mean(MDS2),
                   se_MDS1 = se_fn(MDS1),
                   se_MDS2 = se_fn(MDS2))


community_MDS_summarize$Season <- factor(community_MDS_summarize$Season, levels = c("SUMMER","FALL","WINTER","SPRING","SUMMER-2"))
levels(community_MDS_summarize$Season) <- c("1","2","3","4","5")

my.cols <-c("#E41A1C","#377EB8","black","green4","#FF7F00","grey","#A65628")

summarized_MDS<-ggplot(data = community_MDS_summarize[order(community_MDS_summarize$Season),], 
                       aes(x=mean_MDS1,y=mean_MDS2, color = Site)) + 
  geom_point(aes(color = Site),size=5.5, color = "black") +  # underlying point so that every point has a black outline
  geom_point(aes(color = Site),size=5) +  # actual points
  geom_errorbarh(aes(color = Site,xmax = mean_MDS1 + se_MDS1, xmin = mean_MDS1 - se_MDS1)) + # horizontal error bars (SE) 
  geom_errorbar(aes(color = Site,ymax = mean_MDS2 + se_MDS2, ymin = mean_MDS2 - se_MDS2)) + # vertical error bars (SE) 
  theme_bw() + # theme
  theme(text = element_text(size = 14, color = "black"), # theme formatting 
        axis.text.x = element_text(color="black"), 
        axis.text.y = element_text(color="black"),
        legend.position = c(0.2,.15),
        legend.background = element_blank(),
        legend.title = element_blank()) +
  labs(color = "Site", x = "MDS1",y="MDS2") + # labels 
  scale_color_manual(values = my.cols) +
  removeGrid() +# no grid please
  guides(color=guide_legend(ncol=2))


  

  
  ### add vectors for overlay 
summarized_MDS_V2<-summarized_MDS +
  geom_segment(data=vectors,                           
               aes(x=0,xend=NMDS1*10,y=0,yend=NMDS2*10),inherit.aes = F,
               arrow = arrow(length = unit(0.2,"cm")), 
               colour="black") + 
    geom_text(data=vectors, 
              aes(x=NMDS1*10,y=NMDS2*10,label=species),inherit.aes = F, 
              size=4,
              color = "black")

  
  

ggsave("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Figures/spatiotemporal_MDS.pdf", 
       plot = summarized_MDS_V2, width = 7.5, height = 5.5)
  




#### Supplement - Correlation matrix for algae and structure variables #### 
#library(corrplot)

corr_table<-cbind(MLR_data_ext[c(10:12)],algae_log[c(5)],algae_sqrt[c(1:3)])

names(corr_table) <-c("Height",
                      "Volume",
                      "Surface area",
                      "M. pyrifera",
                      "S. horneri",
                      "S. palmeri",
                      "S. neglecta")

cor_matrix<-cor(corr_table)

col2 = colorRampPalette(c('red', 'white', 'blue'))  


pdf(file = "/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Results/Figures/variable_corrplot.pdf")
library(corrplot)
corrplot(cor_matrix, type = "upper", 
         tl.col = "black", 
         tl.srt = 45,
         col = col2(10))

dev.off()





### END ###
### Fish biomass ####

# Written by: G. Srednick - 3/25/2022


# Steps: 
# (1) convert fish total length to biomass per species
# (2) aggregate that to total biomass per transect across species


# packages
library(tidyverse)
library(rfishbase)
library(plyr)

names(allspp_data)


# ok so make a list of all of these spp in scientific name, then source the length weight relationships
spp_name_df<-read.csv("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/spp_table.csv",strip.white=TRUE)

spp_name<-spp_name_df$species

L2W_table<-length_weight(spp_name)

# manually filter by local measurements
#write.csv(L2W_table,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure/Analyses/submitted_file/Data_for_pub/L2W_table.csv", row.names = F)

# back in 
L2W_table_actual<-read.csv("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Analyses/submitted_file/Data_for_pub/L2W_table_in.csv",strip.white=TRUE)

L2W_table_actual<-L2W_table_actual[c(2,14,16)]


# the function is here

# x = length in cm 
# a = spp unique body shape and condition
# b = spp unique isometric growth

L2W<-function(x) {
  W<- (a(x)^b)
  return(W)
}

# then covert in dplyr with a function based on metrics provided ("a" and "b")


## Steps
# can take each species which will have colmans as sizes and rows as number per transect
# convert those colnames (TL in cm) to biomasses
# then sum to get biomass of each spp per transect



## a bunch of csvs; one for every species
# read them all in and convert them to long
# covert TL to biomass
# multiply biomass by count 
# sum for transect 

setwd("/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/species")


species_csvs = list.files(pattern="*.csv")

ALHO<-read.csv("ALHO.csv",check.names = FALSE,strip.white=TRUE)
ANDA<-read.csv("ANDA.csv",check.names = FALSE,strip.white=TRUE)
APGU<-read.csv("APGU.csv",check.names = FALSE,strip.white=TRUE)
ATAF<-read.csv("ATAF.csv",check.names = FALSE,strip.white=TRUE)
BAPO<-read.csv("BAPO.csv",check.names = FALSE,strip.white=TRUE)
BRFR<-read.csv("BRFR.csv",check.names = FALSE,strip.white=TRUE)
CAPR<-read.csv("CAPR.csv",check.names = FALSE,strip.white=TRUE)
CHPU<-read.csv("CHPU.csv",check.names = FALSE,strip.white=TRUE)
EMJA<-read.csv("EMJA.csv",check.names = FALSE,strip.white=TRUE)
GIEL<-read.csv("GIEL.csv",check.names = FALSE,strip.white=TRUE)
GINI<-read.csv("GINI.csv",check.names = FALSE,strip.white=TRUE)
HASE<-read.csv("HASE.csv",check.names = FALSE,strip.white=TRUE)
HEFR<-read.csv("HEFR.csv",check.names = FALSE,strip.white=TRUE)
HERO<-read.csv("HERO.csv",check.names = FALSE,strip.white=TRUE)
HYRU<-read.csv("HYRU.csv",check.names = FALSE,strip.white=TRUE)
LYDA<-read.csv("LYDA.csv",check.names = FALSE,strip.white=TRUE)
LYZE<-read.csv("LYZE.csv",check.names = FALSE,strip.white=TRUE)
MECA<-read.csv("MECA.csv",check.names = FALSE,strip.white=TRUE)
OXCA<-read.csv("OXCA.csv",check.names = FALSE,strip.white=TRUE)
PACL<-read.csv("PACL.csv",check.names = FALSE,strip.white=TRUE)
PAIN<-read.csv("PAIN.csv",check.names = FALSE,strip.white=TRUE)
PANE<-read.csv("PANE.csv",check.names = FALSE,strip.white=TRUE)
RHNI<-read.csv("RHNI.csv",check.names = FALSE,strip.white=TRUE)
RHTO<-read.csv("RHTO.csv",check.names = FALSE,strip.white=TRUE)
RHVA<-read.csv("RHVA.csv",check.names = FALSE,strip.white=TRUE)
SCGU<-read.csv("SCGU.csv",check.names = FALSE,strip.white=TRUE)
SEAT<-read.csv("SEAT.csv",check.names = FALSE,strip.white=TRUE)
SEAU<-read.csv("SEAU.csv",check.names = FALSE,strip.white=TRUE)
SELA<-read.csv("SELA.csv",check.names = FALSE,strip.white=TRUE)
SEPU<-read.csv("SEPU.csv",check.names = FALSE,strip.white=TRUE)
SESE<-read.csv("SESE.csv",check.names = FALSE,strip.white=TRUE)
STGI<-read.csv("STGI.csv",check.names = FALSE,strip.white=TRUE)
TRSC<-read.csv("TRSC.csv",check.names = FALSE,strip.white=TRUE)

all_species_lengths<-rbind.fill(
  ALHO,
      ANDA,
      APGU,
      ATAF,
      BAPO,
      BRFR,
      CAPR,
      CHPU,
      EMJA,
      GIEL,
      GINI,
      HASE,
      HEFR,
      HERO,
      HYRU,
      LYDA,
      LYZE,
      MECA,
      OXCA,
      PACL,
      PAIN,
      PANE,
      RHNI,
      RHTO,
      RHVA,
      SCGU,
      SEAT,
      SEAU,
      SELA,
      SEPU,
      SESE,
      STGI,
      TRSC)

names(all_species_lengths)

# a function to turn all of these to long
all_species_lengths_long<-all_species_lengths %>% 
  gather(size, abundance,-c(Season,Site,Transect,Strata,species,QUAD))

dim(all_species_lengths_long) # [1] 1085700      8

# merge this and bring in a and b
all_species_lengths_long_V2<-merge(all_species_lengths_long,L2W_table_actual, by = "species")

dim(all_species_lengths_long_V2) #[1] 1085700      10




# summarize fish to transect 
all_species_lengths_long_V3<-all_species_lengths_long_V2 %>% group_by(species) %>% 
  group_by(Site,Season,Transect,Strata,species,size,a,b) %>% 
  summarise_if(is.numeric, mean) %>% 
  select(1:10)

all_species_lengths_long_V3$QUAD<-NULL
dim(all_species_lengths_long_V3)

# merge back with common names
biomass_merge<-merge(all_species_lengths_long_V3,spp_name_df)
dim(biomass_merge)



#### Convert to proper spatial scale (transect)

# cryptic
cryptic_coverted<-biomass_merge %>% 
  filter(strata_counted == "crpytic") %>%
  mutate(meter_abundance = (as.numeric(abundance)*4))
  
dim(cryptic_coverted)

# consp
cons_converted<-biomass_merge %>% 
  filter(strata_counted == "cons") %>%
  mutate(meter_abundance = (as.numeric(abundance)/60))
 
dim(cons_converted)


# good to go
abundance_converted<-rbind(cryptic_coverted,cons_converted)

# now convert to transect scale
abundance_converted$tran_abundance <- abundance_converted$meter_abundance*60


# now convert sizes to biomass, then total biomass, then average to transect

# ok now for every species convert the size to biomass
abundance_biomass<-abundance_converted %>% 
  group_by(species) %>% 
  mutate(biomass = a*as.numeric(size)^b,
         total_biomass = biomass*meter_abundance)




# biomass
all_species_biomass<-abundance_biomass %>% 
  group_by(Site,Season,Transect,Strata,species) %>% 
  summarize_at(vars(total_biomass),sum, na.rm = T)


# sum across strata 
biomass_merge_strata<-all_species_biomass %>% 
  group_by(Site,Season,Transect,species) %>% 
  summarize_at(vars(total_biomass),sum, na.rm = T)


# per transect 
# have a metric that removes LYDA, as well as LYDA and CHPU, so 3 metrics; 2 of which I will use for the present analysis

biomass_data<-biomass_merge_strata %>% 
  group_by(Site,Season,Transect) %>% 
  summarize_at(vars(total_biomass),sum, na.rm = T)


biomass_small<-biomass_merge_strata %>% 
  filter(!species == "Lythrypnus dalli") %>% 
  group_by(Site,Season,Transect) %>% 
  summarize_at(vars(total_biomass),
               funs(biomass_small = sum), na.rm = T)
  
biomass_smaller<-biomass_merge_strata %>% 
  filter(!species %in% c("Lythrypnus dalli", "Chromis punctipinnis")) %>% 
  group_by(Site,Season,Transect) %>% 
  summarize_at(vars(total_biomass),
               funs(biomass_reduced = sum), na.rm = T)

biomass_actual<-biomass_data
names(biomass_actual) <- c("SITE", "SEASON","TRAN","total_biomass")


# test assumptions
qqnorm(log(biomass_actual$total_biomass+1)) # yes for total biomass
biomass_actual$log_total_biomass<-log(biomass_actual$total_biomass +1)

ggplot(biomass_actual,aes(x = Season, y=log_total_biomass, fill = Site)) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               width = 0.01,
               position=position_dodge(0.6)) +
  stat_summary(aes(fill = Site),
               geom = "point",
               pch=21,
               fun = "mean", 
               size = 6,
               position=position_dodge(0.6)) +
  stat_summary(aes(group = Site), 
               geom = "line",
               fun = "mean", 
               position=position_dodge(0.6)) +
  theme(text = element_text(size=15)) +
  labs(x = "Season", y = "total biomass (g)") +
  theme_bw() +
  removeGrid()



names(biomass_actual) <- c("Site", "Season","TRAN","total_biomass","log_total_biomass")



write_csv(biomass_actual,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data_for_pub/biomass.csv")



#### END #####
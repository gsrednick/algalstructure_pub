### Calculate Diversity for Algal structure paper ###
# Written by G. Srednick
# 3/28/2022


# Goal calculate Shannon Weiner Diversity index for every replicate transect
# conduct formal tests back in author script for new metric


library(vegan)


diversity_df<-allspp_data[c(4:37)]

H<-diversity(diversity_df)


diversity_data<-cbind(allspp_data[c(1:3)],H)

qqnorm(log(diversity_data$H+1)) # unnecessary


write_csv(diversity_data,"/Users/icarus2/Documents/Software/R/Tidy Workshop/Git/algalstructure_pub/Data/diversity.csv")



#### DONE #### 
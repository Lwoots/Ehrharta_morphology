#Data cleaning of Ehrharta morphology measurements

#----------------------------------------------------

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(here, tidyverse, visdat, ggfortify, corrplot)

dat <- read.csv(here("Raw/Unprocessed", "Morphology_20191003.csv"), sep = ";")

#Data exploration
head(dat)
summary(dat)
vis_dat(dat)
glimpse(dat)

pairs(dat[,4:21])

#dat_omit <- na.omit(dat)
#
#corr_dat <- cor(dat[,4:24])
#corrplot(corr_dat[,4:21])
#
#hist(dat$Plant_height)
#hist(dat$Veg_height)
#hist((dat$Plant_height-dat$Veg_height)/dat$Plant_height)


datforpca <- scale(dat[,c(4:24)])

hist(dat$Plant_height)

dat_pca <- prcomp(datforpca)

plot(dat_pca$x[,1], dat_pca$x[,2])

autoplot(dat_pca,
         data = dat,
         colour = 'Species',
         #loadings = TRUE,
         #loadings.label = T
         )

autoplot(dat_pca,
         data = dat,
         colour = 'Species',
         loadings = TRUE,
         loadings.label = T,
         size = 3
)

derivative_dat <- dat %>%
    select(Collection,
           Species,
           Plant_height,
           Leaf_width,
           Leaf_length,
           Internode_dist,
           Culm_diam,
           Spikelet_no,
           Pedicel_length_epu,
           Striation_no) %>%
    mutate(Clearance = dat$Plant_height/dat$Veg_height,
           Leaf_importance = dat$Leaf_length/dat$Sheath_length,
           Glume_sheathing = dat$Glume_length_outer/dat$Spikelet_length,
           Lemma_balance = dat$Lemma_outer_length/dat$Lemma_inner_length,
           Spikelet_shape = dat$Spikelet_length/dat$Spikelet_width,
           Leaf_packing = dat$Sheath_length/dat$Internode_dist)

dat_pca <- prcomp(derivative_dat[,3:16], scale = T)

plot(dat_pca$x[,1], dat_pca$x[,3])

autoplot(dat_pca,
         data = dat,
         colour = 'Species',
         #loadings = TRUE,
         #loadings.label = T,
         size = 5
)

autoplot(dat_pca,
         data = dat,
         colour = 'Species',
         loadings = TRUE,
         loadings.label = T,
         size = 5
)


rupes <- dat[36:84,]
vis_dat(rupes)
rupes <- na.omit(rupes)

ruppca <- scale(rupes[,c(4:21)])
ruppca <- na.omit(ruppca)

dat_pca <- prcomp(ruppca)

plot(dat_pca$x[,1], dat_pca$x[,3])

autoplot(dat_pca,
         data = rupes,
         colour = 'Species',
         #loadings = TRUE,
         #loadings.label = T,
         size = 5
)

autoplot(dat_pca,
         data = rupes,
         colour = 'Species',
         loadings = TRUE,
         loadings.label = T,
         size = 5
)
chc


derivative_rup <- rupes %>%
    select(Collection,
           Species,
           Plant_height,
           Leaf_width,
           Leaf_length,
           Internode_dist,
           Culm_diam,
           Spikelet_no,
           Pedicel_length_epu) %>%
    mutate(Clearance = rupes$Plant_height/rupes$Veg_height,
           Leaf_importance = rupes$Leaf_length/rupes$Sheath_length,
           Glume_sheathing = rupes$Glume_length_outer/rupes$Spikelet_length,
           Lemma_balance = rupes$Lemma_outer_length/rupes$Lemma_inner_length,
           Spikelet_shape = rupes$Spikelet_length/rupes$Spikelet_width,
           Leaf_packing = rupes$Sheath_length/rupes$Internode_dist)


dat_pca <- prcomp(derivative_rup[,3:15], scale. = T)

autoplot(dat_pca,
         data = rupes,
         colour = 'Species',
         loadings = TRUE,
         loadings.label = T,
         size = 5
)
jh

#Data cleaning of Ehrharta morphology measurements

#----------------------------------------------------

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(here, tidyverse, visdat, ggfortify)

dat <- read.csv(here("Raw/Unprocessed", "Morphology_20191003.csv"), sep = ";")
dat <- na.omit(dat)

head(dat)
summary(dat)

vis_dat(dat)

pairs(dat[,4:21])
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
         size = 5
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
           Spikelet_shape = dat$Spikelet_length/dat$Spikelet_width)

dat_pca <- prcomp(derivative_dat[,3:14], scale = T)

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


rupes <- dat[36:60,]
vis_dat(rupes)
rupes <- na.omit(rupes)

ruppca <- scale(rupes[,c(4:21)])
ruppca <- na.omit(ruppca)

dat_pca <- prcomp(ruppca)

plot(dat_pca$x[,1], dat_pca$x[,2])

autoplot(dat_pca,
         data = rupes,
         colour = 'Species',
         #loadings = TRUE,
         #loadings.label = T
)

autoplot(dat_pca,
         data = rupes,
         colour = 'Species',
         loadings = TRUE,
         loadings.label = T,
         size = 5
)
chc
testing branches

and again

#Data cleaning of Ehrharta morphology measurements

#----------------------------------------------------

if(!require(pacman)){install.packages("pacman", dependencies=TRUE); library(pacman)}
p_load(here, tidyverse, visdat, ggfortify, corrplot)

dat <- read.csv(here("Raw/Unprocessed", "Morphology_20200120.csv"), sep = ",")

#Data exploration
head(dat)
summary(dat)
vis_dat(dat)
glimpse(dat)

dat <- dat[,1:24] #Get rid of extra columns

pairs(dat[,4:21])
#vegetative pairs
pairs(dat[,4:12])
#Spikelet pairs
pairs(dat[,15:21])


ehrharta <- na.omit(dat[,c(1:24)])

#Subset data

rupestris_clade <- ehrharta[c(grep("rupestris", ehrharta$Species), grep("setacea", ehrharta$Species)),]
ramosa_clade <- ehrharta[c(grep("ramosa", ehrharta$Species), grep("rehmannii", ehrharta$Species)),]


#PCA with the complete data

pca_ehrharta <- prcomp(ehrharta[,4:24], scale. = T)

plot(pca_ehrharta$x[,1], pca_ehrharta$x[,2])

autoplot(pca_ehrharta,
         data = ehrharta,
         colour = 'Species',
         size = 3
)

autoplot(pca_ehrharta,
         data = ehrharta,
         colour = 'Species',
         size = 3,
         loadings = TRUE,
         loadings.label = T
         
)

#PCA rupestris clade

pca_rupestris <- prcomp(rupestris_clade[,c(4:21)], scale. = T) #Drop last 3 columns as they don't vary

plot(pca_rupestris$x[,1], pca_rupestris$x[,2])

autoplot(pca_rupestris,
         data = rupestris_clade,
         colour = 'Species',
         size = 3
)

autoplot(pca_rupestris,
         data = rupestris_clade,
         colour = 'Species',
         size = 3,
         loadings = TRUE,
         loadings.label = T
)

#PCA ramosa clade

pca_ramosa <- prcomp(ramosa_clade[,c(4:24)], scale. = T)

plot(pca_ramosa$x[,1], pca_ramosa$x[,2])

autoplot(pca_ramosa,
         data = ramosa_clade,
         colour = 'Species',
         size = 3
)

autoplot(pca_ramosa,
         data = ramosa_clade,
         colour = 'Species',
         size = 3,
         loadings = TRUE,
         loadings.label = T
)

#How do ratios change the PCA? ####

#Create dataset with derivative measures
derivative_ehrharta <- ehrharta %>%
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
    mutate(Clearance = ehrharta$Plant_height/ehrharta$Veg_height,
           Leaf_importance = ehrharta$Leaf_length/ehrharta$Sheath_length,
           Glume_sheathing = ehrharta$Glume_length_outer/ehrharta$Spikelet_length,
           Lemma_balance = ehrharta$Lemma_outer_length/ehrharta$Lemma_inner_length,
           Spikelet_shape = ehrharta$Spikelet_length/ehrharta$Spikelet_width,
           Leaf_packing = ehrharta$Sheath_length/ehrharta$Internode_dist)

#Subset data

rupestris_derivative <- derivative_ehrharta[c(grep("rupestris", derivative_ehrharta$Species), grep("setacea", derivative_ehrharta$Species)),]
ramosa_derivative <- derivative_ehrharta[c(grep("ramosa", derivative_ehrharta$Species), grep("rehmannii", derivative_ehrharta$Species)),]

#PCA full dataset

pca_derivative_ehrharta <- prcomp(derivative_ehrharta[,3:16], scale. = T)

plot(pca_derivative_ehrharta$x[,1], pca_derivative_ehrharta$x[,2])

autoplot(pca_derivative_ehrharta,
         data = derivative_ehrharta,
         colour = 'Species',
         size = 3
)

autoplot(pca_derivative_ehrharta,
         data = derivative_ehrharta,
         colour = 'Species',
         size = 3,
         loadings = TRUE,
         loadings.label = T
         
)

#PCA rupestris clade

pca_derivative_rupestris <- prcomp(rupestris_derivative[,c(3:9, 11:16)], scale. = T) #Drop unvarying column

plot(pca_derivative_rupestris$x[,1], pca_derivative_rupestris$x[,2])

autoplot(pca_derivative_rupestris,
         data = rupestris_derivative,
         colour = 'Species',
         size = 3
)

autoplot(pca_derivative_rupestris,
         data = rupestris_derivative,
         colour = 'Species',
         size = 3,
         loadings = TRUE,
         loadings.label = T
)

#PCA ramosa clade

pca_derivative_ramosa <- prcomp(ramosa_derivative[,c(3:16)], scale. = T)

plot(pca_derivative_ramosa$x[,1], pca_derivative_ramosa$x[,2])

autoplot(pca_derivative_ramosa,
         data = ramosa_derivative,
         colour = 'Species',
         size = 3
)

autoplot(pca_derivative_ramosa,
         data = ramosa_derivative,
         colour = 'Species',
         size = 3,
         loadings = TRUE,
         loadings.label = T
)



#Some data visualisation

ggplot(ehrharta, aes(Plant_height)) +
    theme_classic() +
    geom_histogram(alpha = 0.15) +
    geom_histogram(data = rupestris_clade, fill = "dark green", alpha = 0.2) +
    geom_histogram(data = ramosa_clade, fill = "red", alpha = 0.2)

ggplot(ehrharta, aes(Leaf_width)) +
    theme_classic() +
    geom_histogram(alpha = 0.15) +
    geom_histogram(data = rupestris_clade, fill = "dark green", alpha = 0.2) +
    geom_histogram(data = ramosa_clade, fill = "red", alpha = 0.2)

ggplot(ehrharta, aes(Internode_dist)) +
    theme_classic() +
    geom_histogram(alpha = 0.15) +
    geom_histogram(data = rupestris_clade, fill = "dark green", alpha = 0.2) +
    geom_histogram(data = ramosa_clade, fill = "red", alpha = 0.2)



datforpca <- scale(dat_omit[,c(4:24)])

hist(dat$Plant_height)

dat_pca <- prcomp(datforpca)

plot(dat_pca$x[,1], dat_pca$x[,2])

autoplot(dat_pca,
         data = dat_omit,
         colour = 'Species',
         #loadings = TRUE,
         #loadings.label = T
         )

autoplot(dat_pca,
         data = dat_omit,
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

derivative_dat <- na.omit(derivative_dat)
dat_pca <- prcomp(derivative_dat[,3:16], scale = T)

plot(dat_pca$x[,1], dat_pca$x[,3])

autoplot(dat_pca,
         data = derivative_dat,
         colour = 'Species',
         #loadings = TRUE,
         #loadings.label = T,
         size = 5,
         alpha(0.7)
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

ramclade <- dat %>% filter(Species == "E_rehmannii_subspicata" | 
                               Species == "E_rehmannii_filiformis" |
                               Species == "E_rehmannii_rehmannii" |
                               Species == "E_ramosa_ramosa" |
                               Species == "E_ramosa_aphylla")
na.omit(ramclade)
ramclade$Collection <- as.factor(ramclade$Collection)
ramclade_pca <- prcomp(ramclade[,4:15], scale. = T)
autoplot(ramclade_pca,
         data = ramclade,
         colour = 'Species',
         #loadings = TRUE,
         #loadings.label = T,
         size = 5
)


derivative_ram <- ramclade %>%
    select(Collection,
           Species,
           Plant_height,
           Leaf_width,
           Leaf_length,
           Internode_dist,
           Culm_diam,
           Spikelet_no,
           Pedicel_length_epu) %>%
    mutate(Clearance = ramclade$Plant_height/ramclade$Veg_height,
           Leaf_importance = ramclade$Leaf_length/ramclade$Sheath_length,
           Glume_sheathing = ramclade$Glume_length_outer/ramclade$Spikelet_length,
           Lemma_balance = ramclade$Lemma_outer_length/ramclade$Lemma_inner_length,
           Spikelet_shape = ramclade$Spikelet_length/ramclade$Spikelet_width,
           Leaf_packing = ramclade$Sheath_length/ramclade$Internode_dist)


dat_pca <- prcomp(derivative_ram[,4:15], scale. = T)

autoplot(dat_pca,
         data = ramclade,
         colour = 'Species',
         loadings = TRUE,
         loadings.label = T,
         size = 5
)

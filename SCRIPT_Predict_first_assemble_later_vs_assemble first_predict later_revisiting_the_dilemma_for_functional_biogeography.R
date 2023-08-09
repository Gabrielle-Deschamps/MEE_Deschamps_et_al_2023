rm(list=ls())
#This script computes all path analyses from :
#Deschamps G.; Giovanni P.; Brun P.; Thuiller, W. (2023) 
#Predict first - assemble later vs assemble first - predict later: revisiting the dilemma for functional biogeography.

## SETTINGS
setwd("...") # set the path where the function script is
source("FONCTIONS_SCRIPT_Predict_first_assemble_later_vs_assemble first_predict later_revisiting_the_dilemma_for_functional_biogeography")
# savings
path_to_load <- "..." # where raw data are saved
path_to_save <- "..." # where results have to be saved
path_to_plot <- "..." # where plots have to be plotted
dir.create(path_to_save)
dir.create(path_to_plot)

###########################
#### -- Data tables -- ####
###########################
# This part formats the data to be used in later analyses

#### -- Environmental data -- ####
load(paste0(path_to_load,"/environmental_variables.RData")) # load environmental variables, name is env_var
#This files contains following informations: "Site", "Temperature_Seasonality", "Mean_Temperature_of_Coldest_Quarter", "Aridity", "Precipitation_Seasonality", "Soil_ind.F_rdf", 
#"Soil_ind.R_rdf", "Slope".
#Further description of the sampling and variables is descripted in the Mat&Met section of the paper

#### -- botanical survey -- ####
load(paste0(path_to_load, "/plant_community.RData")) # load grassland community plot dataset, name is plant_recovery and plant_pa
#plant_recovery contains the cover-abundance classes of each species for each sampled plot.
#plant_pa contains the information of presence-absence for each species for each sampled plot.
#Further description of the sampling and variables is descripted in the Mat&Met section of the paper.

#### -- species traits -- ####
load(paste0(path_to_load, "/species_traits.RData")) # load species traits, name is traits.
#This files contains following informations: "Taxa","LNC","SLA","PLH".
#Further description of the sampling and variables is descripted in the Mat&Met section of the paper.

#### -- functional indices of the community -- ####
CWM_obs <- CWM(plant_recovery) 
Site_order <- data.frame("Site"=rownames(plant_recovery))
CWM_obs <- left_join(Site_order, CWM_obs, by="Site")
rownames(CWM_obs)=CWM_obs$Site
## built a data table with for each sample site the community weighted means of each of the three traits (LNC, SLA, PLH).
# This files contains following informations : "Site","LNC","SLA","PLH".
CM_obs <- CWM(plant_pa) ## built a data table with for each sample site the community means of each of the three traits (LNC, SLA, PLH).
Site_order <- data.frame("Site"=rownames(plant_pa))
CM_obs <- left_join(Site_order, CM_obs, by="Site")
rownames(CM_obs)=CM_obs$Site
# This files contains following informations : "Site","LNC","SLA","PLH".
FDis_obs <- FDis(plant_recovery)## built a data table with for each sample site the weighted functional dispersion of each of the three traits (LNC, SLA, PLH).
Site_order <- data.frame("Site"=rownames(plant_recovery))
FDis_obs <- left_join(Site_order, FDis_obs, by="Site")
rownames(FDis_obs)=FDis_obs$Site
# This files contains following informations : "Site","LNC","SLA","PLH".
uFDis_obs <- FDis(plant_pa)## built a data table with for each sample site the unweighted functional dispersion of each of the three traits (LNC, SLA, PLH).
Site_order <- data.frame("Site"=rownames(plant_pa))
uFDis_obs <- left_join(Site_order, uFDis_obs, by="Site")
rownames(uFDis_obs)=uFDis_obs$Site
# This files contains following informations : "Site","LNC","SLA","PLH".

#### -- groups for cross-validation procedures -- ####
load(paste0(path_to_load, "/group_cross_validation.RData")) # load groups for the two different four-fold cross-validation procedures, name is group_extra1, group_extra2, group_extra3, group_extra4, group_inter1, group_inter2, group_inter3, group_inter4.
#Further description of the sampling and variables is descripted in the Mat&Met section of the paper.

######################################
#### -- SPECIES-BASED APPROACH -- ####
######################################
####### -- Species Distribution Models for the Species-based approach-- ####
###### -- PREDICTION OF SPECIES ABUNDANCE ######
#### -- This computes path predicts the relative abundance of species in the context of interpolation -- ####
gam_data <- data.frame(plant_recovery,env_var) # assembly of species recovery data and environmental variables.

## - For interpolation group 1 - ##
species_models_recovery_grow1<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,]) # construction of an abundance model for the species Aceras anthropophorum at sites in interpolation groups 2, 3 and 4.
species_models_recovery_inter1 <- predict(species_models_recovery_grow1, gam_data[group_inter1,])$predicted # prediction of abundance of species Aceras anthropophorum at sites in interpolation group 1.
for (i in 2 : ncol(plant_recovery[group_inter1,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter1 <- cbind(species_models_recovery_inter1, predict(model, gam_data[group_inter1,])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter1) <- colnames(plant_recovery)
rownames(species_models_recovery_inter1) <- rownames(plant_recovery[group_inter1,])


## - For interpolation group 2 - ##
species_models_recovery_grow2<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,]) # construction of an abundance model for the species Aceras anthropophorum at sites in interpolation groups 2, 3 and 4.
species_models_recovery_inter2 <- predict(species_models_recovery_grow2, gam_data[group_inter2,])$predicted # prediction of abundance of species Aceras anthropophorum at sites in interpolation group 2.
for (i in 2 : ncol(plant_recovery[group_inter2,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter2 <- cbind(species_models_recovery_inter2, predict(model, gam_data[group_inter2,])$predicted) # prediction of abundance for each species at sites i interpolation group 2.
}
colnames(species_models_recovery_inter2) <- colnames(plant_recovery)
rownames(species_models_recovery_inter2) <- rownames(plant_recovery[group_inter2,])

## - For interpolation group 3 - ##
species_models_recovery_grow3<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,]) # construction of an abundance model for the species Aceras anthropophorum at sites in interpolation groups 2, 3 and 4.
species_models_recovery_inter3 <- predict(species_models_recovery_grow3, gam_data[group_inter3,])$predicted # prediction of abundance of species Aceras anthropophorum at sites in interpolation group 3.
for (i in 2 : ncol(plant_recovery[group_inter3,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter3 <- cbind(species_models_recovery_inter3, predict(model, gam_data[group_inter3,])$predicted) # prediction of abundance for each species at sites i interpolation group 3.
}
colnames(species_models_recovery_inter3) <- colnames(plant_recovery)
rownames(species_models_recovery_inter3) <- rownames(plant_recovery[group_inter3,])

## - For interpolation group 4 - ##
species_models_recovery_grow4<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,]) # construction of an abundance model for the species Aceras anthropophorum at sites in interpolation groups 2, 3 and 4.
species_models_recovery_inter4 <- predict(species_models_recovery_grow4, gam_data[group_inter4,])$predicted # prediction of abundance of species Aceras anthropophorum at sites in interpolation group 4.
for (i in 2 : ncol(plant_recovery[group_inter4,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter4 <- cbind(species_models_recovery_inter4, predict(model, gam_data[group_inter4,])$predicted) # prediction of abundance for each species at sites i interpolation group 4.
}
colnames(species_models_recovery_inter4) <- colnames(plant_recovery)
rownames(species_models_recovery_inter4) <- rownames(plant_recovery[group_inter4,])


## assembly in a single prediction table for the 4 groups :
species_models_recovery_inter <- data.frame(rbind(species_models_recovery_inter1, species_models_recovery_inter2, species_models_recovery_inter3, species_models_recovery_inter4))
colnames(species_models_recovery_inter)<- colnames(species_models_recovery_inter1)
species_models_recovery_inter <- cbind("Site"=rownames(species_models_recovery_inter), species_models_recovery_inter)
Site_order <- data.frame("Site"=rownames(plant_recovery))
species_models_recovery_inter <- left_join(Site_order, species_models_recovery_inter, by="Site")
rownames(species_models_recovery_inter)=species_models_recovery_inter$Site
species_models_recovery_inter <- species_models_recovery_inter[,-1]


#### -- This computes path predicts the relative abundance of species in the context of extrapolation -- ####
## - For extrapolation group 1 - ##
species_models_recovery_grow1<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,]) # construction of an abundance model for the species Aceras anthropophorum at sites in extrapolation groups 2, 3 and 4.
species_models_recovery_extra1 <- predict(species_models_recovery_grow1, gam_data[group_extra1,])$predicted # prediction of abundance of species Aceras anthropophorum at sites in extrapolation group 1.
for (i in 2 : ncol(plant_recovery[group_extra1,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra1 <- cbind(species_models_recovery_extra1, predict(model, gam_data[group_extra1,])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra1) <- colnames(plant_recovery)
rownames(species_models_recovery_extra1) <- rownames(plant_recovery[group_extra1,])


## - For extrapolation group 2 - ##
species_models_recovery_grow2<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,]) # construction of an abundance model for the species Aceras anthropophorum at sites in extrapolation groups 2, 3 and 4.
species_models_recovery_extra2 <- predict(species_models_recovery_grow2, gam_data[group_extra2,])$predicted # prediction of abundance of species Aceras anthropophorum at sites in extrapolation group 2.
for (i in 2 : ncol(plant_recovery[group_extra2,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra2 <- cbind(species_models_recovery_extra2, predict(model, gam_data[group_extra2,])$predicted) # prediction of abundance for each species at sites i extrapolation group 2.
}
colnames(species_models_recovery_extra2) <- colnames(plant_recovery)
rownames(species_models_recovery_extra2) <- rownames(plant_recovery[group_extra2,])

## - For extrapolation group 3 - ##
species_models_recovery_grow3<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,]) # construction of an abundance model for the species Aceras anthropophorum at sites in extrapolation groups 2, 3 and 4.
species_models_recovery_extra3 <- predict(species_models_recovery_grow3, gam_data[group_extra3,])$predicted # prediction of abundance of species Aceras anthropophorum at sites in extrapolation group 3.
for (i in 2 : ncol(plant_recovery[group_extra3,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra3 <- cbind(species_models_recovery_extra3, predict(model, gam_data[group_extra3,])$predicted) # prediction of abundance for each species at sites i extrapolation group 3.
}
colnames(species_models_recovery_extra3) <- colnames(plant_recovery)
rownames(species_models_recovery_extra3) <- rownames(plant_recovery[group_extra3,])

## - For extrapolation group 4 - ##
species_models_recovery_grow4<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,]) # construction of an abundance model for the species Aceras anthropophorum at sites in extrapolation groups 2, 3 and 4.
species_models_recovery_extra4 <- predict(species_models_recovery_grow4, gam_data[group_extra4,])$predicted # prediction of abundance of species Aceras anthropophorum at sites in extrapolation group 4.
for (i in 2 : ncol(plant_recovery[group_extra4,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra4 <- cbind(species_models_recovery_extra4, predict(model, gam_data[group_extra4,])$predicted) # prediction of abundance for each species at sites i extrapolation group 4.
}
colnames(species_models_recovery_extra4) <- colnames(plant_recovery)
rownames(species_models_recovery_extra4) <- rownames(plant_recovery[group_extra4,])

## assembly in a single prediction table for the 4 groups :
species_models_recovery_extra <- data.frame(rbind(species_models_recovery_extra1, species_models_recovery_extra2, species_models_recovery_extra3, species_models_recovery_extra4))
colnames(species_models_recovery_extra)<- colnames(species_models_recovery_extra1)
species_models_recovery_extra <- cbind("Site"=rownames(species_models_recovery_extra), species_models_recovery_extra)
Site_order <- data.frame("Site"=rownames(plant_recovery))
species_models_recovery_extra <- left_join(Site_order, species_models_recovery_extra, by="Site")
rownames(species_models_recovery_extra)=species_models_recovery_extra$Site
species_models_recovery_extra <- species_models_recovery_extra[,-1]

#### performance of species abundance models :
# for interpolation : 
R2_species_models_recovery_inter <- rep(0,831)
for (i in 1 : 831){
  R2_species_models_recovery_inter[i] <- cor(plant_recovery[,i],species_models_recovery_inter[,i])^2
}
# for extrapolation :
R2_species_models_recovery_extra <- rep(0,831)
for (i in 1 : 831){
  R2_species_models_recovery_extra[i] <- cor(plant_recovery[,i],species_models_recovery_extra[,i])^2
}

###### PREDICTION OF PRESENCE PROBABILITIES ######
#### -- This computes path predicts the presence probabilities of species in the context of interpolation -- ####
gam_data <- data.frame(plant_pa,env_var) # assembly of species presences/absences data and environmental variables.

## - For interpolation group 1 - ##
species_models_pa_grow1<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,]) # construction of a probability of occurrence model for the species Aceras anthropophorum at sites in interpolation groups 2, 3 and 4.
species_models_pa_inter1 <- predict(species_models_pa_grow1, gam_data[group_inter1,])$predicted # prediction of a probability of occurrence model of species Aceras anthropophorum at sites in interpolation group 1.
for (i in 2 : ncol(plant_pa[group_inter1,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])"))) #Probability of occurrence modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_pa_inter1 <- cbind(species_models_pa_inter1, predict(model, gam_data[group_inter1,])$predicted) # prediction of probabilities of occurrence for each species at sites i interpolation group 1.
}
colnames(species_models_pa_inter1) <- colnames(plant_pa)
rownames(species_models_pa_inter1) <- rownames(plant_pa[group_inter1,])


## - For interpolation group 2 - ##
species_models_pa_grow2<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,]) # construction of a probability of occurrence model for the species Aceras anthropophorum at sites in interpolation groups 2, 3 and 4.
species_models_pa_inter2 <- predict(species_models_pa_grow2, gam_data[group_inter2,])$predicted # prediction of probability of occurrence of species Aceras anthropophorum at sites in interpolation group 2.
for (i in 2 : ncol(plant_pa[group_inter2,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])"))) #probability of occurrence modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_pa_inter2 <- cbind(species_models_pa_inter2, predict(model, gam_data[group_inter2,])$predicted) # prediction of probabilities of occurrence for each species at sites i interpolation group 2.
}
colnames(species_models_pa_inter2) <- colnames(plant_pa)
rownames(species_models_pa_inter2) <- rownames(plant_pa[group_inter2,])

## - For interpolation group 3 - ##
species_models_pa_grow3<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,]) # construction of a probability of occurrence model for the species Aceras anthropophorum at sites in interpolation groups 2, 3 and 4.
species_models_pa_inter3 <- predict(species_models_pa_grow3, gam_data[group_inter3,])$predicted # prediction of probability of occurrence of species Aceras anthropophorum at sites in interpolation group 3.
for (i in 2 : ncol(plant_pa[group_inter3,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])"))) #Probability of occurrence modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_pa_inter3 <- cbind(species_models_pa_inter3, predict(model, gam_data[group_inter3,])$predicted) # prediction of probability of occurrence for each species at sites i interpolation group 3.
}
colnames(species_models_pa_inter3) <- colnames(plant_pa)
rownames(species_models_pa_inter3) <- rownames(plant_pa[group_inter3,])

## - For interpolation group 4 - ##
species_models_pa_grow4<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,]) # construction of a probability of occurrence model for the species Aceras anthropophorum at sites in interpolation groups 2, 3 and 4.
species_models_pa_inter4 <- predict(species_models_pa_grow4, gam_data[group_inter4,])$predicted # prediction of probability of occurrence of species Aceras anthropophorum at sites in interpolation group 4.
for (i in 2 : ncol(plant_pa[group_inter4,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])"))) #Probability of occurrence modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_pa_inter4 <- cbind(species_models_pa_inter4, predict(model, gam_data[group_inter4,])$predicted) # prediction of probability of occurrence for each species at sites i interpolation group 4.
}
colnames(species_models_pa_inter4) <- colnames(plant_pa)
rownames(species_models_pa_inter4) <- rownames(plant_pa[group_inter4,])

## assembly in a single prediction table for the 4 groups :
species_models_pa_inter <- data.frame(rbind(species_models_pa_inter1, species_models_pa_inter2, species_models_pa_inter3, species_models_pa_inter4))
colnames(species_models_pa_inter) <- colnames(species_models_pa_inter1)
species_models_pa_inter <- cbind("Site"=rownames(species_models_pa_inter), species_models_pa_inter)
Site_order <- data.frame("Site"=rownames(plant_pa))
species_models_pa_inter <- left_join(Site_order, species_models_pa_inter, by="Site")
rownames(species_models_pa_inter)=species_models_pa_inter$Site
species_models_pa_inter <- species_models_pa_inter[,-1]

#### -- This computes path predicts the relative abundance of species in the context of extrapolation -- ####
## - For extrapolation group 1 - ##
species_models_pa_grow1<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,]) # construction of a probability of occurrence model for the species Aceras anthropophorum at sites in extrapolation groups 2, 3 and 4.
species_models_pa_extra1 <- predict(species_models_pa_grow1, gam_data[group_extra1,])$predicted # prediction of probability of occurrence of species Aceras anthropophorum at sites in extrapolation group 1.
for (i in 2 : ncol(plant_pa[group_extra1,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])"))) #Probability of occurrence modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_pa_extra1 <- cbind(species_models_pa_extra1, predict(model, gam_data[group_extra1,])$predicted) # prediction of probability of occurrence for each species at sites i extrapolation group 1.
}
colnames(species_models_pa_extra1) <- colnames(plant_pa)
rownames(species_models_pa_extra1) <- rownames(plant_pa[group_extra1,])

## - For extrapolation group 2 - ##
species_models_pa_grow2<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,]) # construction of a probability of occurrence model for the species Aceras anthropophorum at sites in extrapolation groups 2, 3 and 4.
species_models_pa_extra2 <- predict(species_models_pa_grow2, gam_data[group_extra2,])$predicted # prediction of probability of occurrence of species Aceras anthropophorum at sites in extrapolation group 2.
for (i in 2 : ncol(plant_pa[group_extra2,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])"))) #Probability of occurrence modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_pa_extra2 <- cbind(species_models_pa_extra2, predict(model, gam_data[group_extra2,])$predicted) # prediction of probability of occurrence for each species at sites i extrapolation group 2.
}
colnames(species_models_pa_extra2) <- colnames(plant_pa)
rownames(species_models_pa_extra2) <- rownames(plant_pa[group_extra2,])

## - For extrapolation group 3 - ##
species_models_pa_grow3<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,]) # construction of a probability of occurrence model for the species Aceras anthropophorum at sites in extrapolation groups 2, 3 and 4.
species_models_pa_extra3 <- predict(species_models_pa_grow3, gam_data[group_extra3,])$predicted # prediction of probability of occurrence of species Aceras anthropophorum at sites in extrapolation group 3.
for (i in 2 : ncol(plant_pa[group_extra3,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])"))) #Probability of occurrence modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_pa_extra3 <- cbind(species_models_pa_extra3, predict(model, gam_data[group_extra3,])$predicted) # prediction of probability of occurrence for each species at sites i extrapolation group 3.
}
colnames(species_models_pa_extra3) <- colnames(plant_pa)
rownames(species_models_pa_extra3) <- rownames(plant_pa[group_extra3,])

## - For extrapolation group 4 - ##
species_models_pa_grow4<- rfsrc(Aceras.anthropophorum~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,]) # construction of a probability of occurrence model for the species Aceras anthropophorum at sites in extrapolation groups 2, 3 and 4.
species_models_pa_extra4 <- predict(species_models_pa_grow4, gam_data[group_extra4,])$predicted # prediction of probability of occurrence of species Aceras anthropophorum at sites in extrapolation group 4.
for (i in 2 : ncol(plant_pa[group_extra4,])){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])"))) #Probability of occurrence modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_pa_extra4 <- cbind(species_models_pa_extra4, predict(model, gam_data[group_extra4,])$predicted) # prediction of probability of occurrence for each species at sites i extrapolation group 4.
}
colnames(species_models_pa_extra4) <- colnames(plant_pa)
rownames(species_models_pa_extra4) <- rownames(plant_pa[group_extra4,])

## assembly in a single prediction table for the 4 groups :
species_models_pa_extra <- data.frame(rbind(species_models_pa_extra1, species_models_pa_extra2, species_models_pa_extra3, species_models_pa_extra4))
colnames(species_models_pa_extra)<- colnames(species_models_pa_extra1)
species_models_pa_extra <- cbind("Site"=rownames(species_models_pa_extra), species_models_pa_extra)
Site_order <- data.frame("Site"=rownames(plant_pa))
species_models_pa_extra <- left_join(Site_order, species_models_pa_extra, by="Site")
rownames(species_models_pa_extra)=species_models_pa_extra$Site
species_models_pa_extra <- species_models_pa_extra[,-1]

#### -- performance of species presence/absence models ####
## for interpolation ##
# Transformation into presence/absence prediction #
threshold_species_models_inter <- rep(0,831)
for (i in 1 : 831){
  threshold_species_models_inter[i]=KappaRepet(plant_pa[,i], species_models_pa_inter[,i], TSS = TRUE)$CutOff
}
species_models_pa_inter_threshold <- species_models_pa_inter
for (i in 1 : 831){
  species_models_pa_inter_threshold[species_models_pa_inter_threshold[,i]<threshold_species_models_inter[i],i]=0
  species_models_pa_inter_threshold[species_models_pa_inter_threshold[,i]>0,i]=1
}
# TSS calculation
TSS_species_models_pa_inter <- rep(0,831)
for (i in 1 : 831){
  Misc <- table(data.frame("Observated"=as.numeric(plant_pa[,i]), "Predicted"=as.numeric(species_models_pa_inter_threshold[,i])))
  TSS_species_models_pa_inter[i]<-(Misc[2,2]/(Misc[2,2]+Misc[2,1])+(Misc[1,1]/(Misc[1,1]+Misc[1,2]))-1)
}

## for extrapolation ##
# Transformation into presence/absence prediction #
threshold_species_models_extra <- rep(0,831)
for (i in 1 : 831){
  threshold_species_models_extra[i]=KappaRepet(plant_pa[,i], species_models_pa_extra[,i], TSS = TRUE)$CutOff
}
species_models_pa_extra_threshold <- species_models_pa_extra
for (i in 1 : 831){
  species_models_pa_extra_threshold[species_models_pa_extra_threshold[,i]<threshold_species_models_extra[i],i]=0
  species_models_pa_extra_threshold[species_models_pa_extra_threshold[,i]>0,i]=1
}
# TSS calculation
TSS_species_models_pa_extra <- rep(0,831)
for (i in 1 : 831){
  Misc <- table(data.frame("Observated"=as.numeric(plant_pa[,i]), "Predicted"=as.numeric(species_models_pa_extra_threshold[,i])))
  TSS_species_models_pa_extra[i]<-(Misc[2,2]/(Misc[2,2]+Misc[2,1])+(Misc[1,1]/(Misc[1,1]+Misc[1,2]))-1)
}

##### -- transformation of species predictions into community functional metrics predictions -- ####
## CWM calculated with abundance prediction for interpolation :
CWM_species_based_inter1 <- CWM(species_models_recovery_inter1)
CWM_species_based_inter2 <- CWM(species_models_recovery_inter2)
CWM_species_based_inter3 <- CWM(species_models_recovery_inter3)
CWM_species_based_inter4 <- CWM(species_models_recovery_inter4)
CWM_species_based_inter <- CWM(species_models_recovery_inter)
Site_order <- data.frame("Site"=rownames(plant_recovery))
CWM_species_based_inter <- left_join(Site_order, CWM_species_based_inter, by="Site")
rownames(CWM_species_based_inter)=CWM_species_based_inter$Site
CWM_species_based_inter <- CWM_species_based_inter[,-1]

## CWM calculated with abundance prediction for extrapolation :
CWM_species_based_extra1 <- CWM(species_models_recovery_extra1)
CWM_species_based_extra2 <- CWM(species_models_recovery_extra2)
CWM_species_based_extra3 <- CWM(species_models_recovery_extra3)
CWM_species_based_extra4 <- CWM(species_models_recovery_extra4)
CWM_species_based_extra <- CWM(species_models_recovery_extra)
Site_order <- data.frame("Site"=rownames(plant_recovery))
CWM_species_based_extra <- left_join(Site_order, CWM_species_based_extra, by="Site")
rownames(CWM_species_based_extra)=CWM_species_based_extra$Site
CWM_species_based_extra <- CWM_species_based_extra[,-1]

## CM calculated with probability of occurrence for interpolation :
CM_species_based_inter1 <- CWM(species_models_pa_inter1)
CM_species_based_inter2 <- CWM(species_models_pa_inter2)
CM_species_based_inter3 <- CWM(species_models_pa_inter3)
CM_species_based_inter4 <- CWM(species_models_pa_inter4)
CM_species_based_inter <- CWM(species_models_pa_inter)
Site_order <- data.frame("Site"=rownames(plant_recovery))
CM_species_based_inter <- left_join(Site_order, CM_species_based_inter, by="Site")
rownames(CM_species_based_inter)=CM_species_based_inter$Site
CM_species_based_inter <- CM_species_based_inter[,-1]

## CM calculated with probability of occurrence for extrapolation :
CM_species_based_extra1 <- CWM(species_models_pa_extra1)
CM_species_based_extra2 <- CWM(species_models_pa_extra2)
CM_species_based_extra3 <- CWM(species_models_pa_extra3)
CM_species_based_extra4 <- CWM(species_models_pa_extra4)
CM_species_based_extra <- CWM(species_models_pa_extra)
Site_order <- data.frame("Site"=rownames(plant_recovery))
CM_species_based_extra <- left_join(Site_order, CM_species_based_extra, by="Site")
rownames(CM_species_based_extra)=CM_species_based_extra$Site
CM_species_based_extra <- CM_species_based_extra[,-1]

## FDis calculated with abundance prediction for interpolation :
FDis_species_based_inter1 <- FDis(species_models_recovery_inter1)
FDis_species_based_inter2 <- FDis(species_models_recovery_inter2)
FDis_species_based_inter3 <- FDis(species_models_recovery_inter3)
FDis_species_based_inter4 <- FDis(species_models_recovery_inter4)
FDis_species_based_inter <- FDis(species_models_recovery_inter)
Site_order <- data.frame("Site"=rownames(plant_recovery))
FDis_species_based_inter <- left_join(Site_order, FDis_species_based_inter, by="Site")
rownames(FDis_species_based_inter)=FDis_species_based_inter$Site
FDis_species_based_inter <- FDis_species_based_inter[,-1]

## FDis calculated with abundance prediction for extrapolation :
FDis_species_based_extra1 <- FDis(species_models_recovery_extra1)
FDis_species_based_extra2 <- FDis(species_models_recovery_extra2)
FDis_species_based_extra3 <- FDis(species_models_recovery_extra3)
FDis_species_based_extra4 <- FDis(species_models_recovery_extra4)
FDis_species_based_extra <- FDis(species_models_recovery_extra)
Site_order <- data.frame("Site"=rownames(plant_recovery))
FDis_species_based_extra <- left_join(Site_order, FDis_species_based_extra, by="Site")
rownames(FDis_species_based_extra)=FDis_species_based_extra$Site
FDis_species_based_extra <- FDis_species_based_extra[,-1]

## uFDis calculated with probability of occurrence for interpolation :
uFDis_species_based_inter1 <- FDis(species_models_pa_inter1)
uFDis_species_based_inter2 <- FDis(species_models_pa_inter2)
uFDis_species_based_inter3 <- FDis(species_models_pa_inter3)
uFDis_species_based_inter4 <- FDis(species_models_pa_inter4)
uFDis_species_based_inter <- FDis(species_models_pa_inter)
Site_order <- data.frame("Site"=rownames(plant_recovery))
uFDis_species_based_inter <- left_join(Site_order, uFDis_species_based_inter, by="Site")
rownames(uFDis_species_based_inter)=uFDis_species_based_inter$Site
uFDis_species_based_inter <- uFDis_species_based_inter[,-1]

## uFDis calculated with probability of occurrence for extrapolation :
uFDis_species_based_extra1 <- FDis(species_models_pa_extra1)
uFDis_species_based_extra2 <- FDis(species_models_pa_extra2)
uFDis_species_based_extra3 <- FDis(species_models_pa_extra3)
uFDis_species_based_extra4 <- FDis(species_models_pa_extra4)
uFDis_species_based_extra <- FDis(species_models_pa_extra)
Site_order <- data.frame("Site"=rownames(plant_recovery))
uFDis_species_based_extra <- left_join(Site_order, uFDis_species_based_extra, by="Site")
rownames(uFDis_species_based_extra)=uFDis_species_based_extra$Site
uFDis_species_based_extra <- uFDis_species_based_extra[,-1]


####################################
#### -- TRAIT-BASED APPROACH -- ####
####################################
####### -- Traits Distribution Models for the Trait-based approach-- 
###### PREDICTION OF COMMUNITY WEIGHTED MEAN  ######
#### -- This computes path predicts the community weighted mean in the context of interpolation -- ####
gam_data <- left_join(CWM_obs,env_var, by="Site")
## - For interpolation group 1 - ##
models_CWM_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
models_CWM_LNC_inter1 <- predict(models_CWM_LNC_grow1, gam_data[group_inter1,])$predicted
models_CWM_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
models_CWM_SLA_inter1 <- predict(models_CWM_SLA_grow1, gam_data[group_inter1,])$predicted
models_CWM_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
models_CWM_PLH_inter1 <- predict(models_CWM_PLH_grow1, gam_data[group_inter1,])$predicted

CWM_trait_based_inter1<-data.frame("LNC"=models_CWM_LNC_inter1, "SLA"=models_CWM_SLA_inter1, "PLH"=models_CWM_PLH_inter1)
rownames(CWM_trait_based_inter1) <- CWM_obs$Site[group_inter1]

## - For interpolation group 2 - ##
models_CWM_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
models_CWM_LNC_inter2 <- predict(models_CWM_LNC_grow2, gam_data[group_inter2,])$predicted
models_CWM_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
models_CWM_SLA_inter2 <- predict(models_CWM_SLA_grow2, gam_data[group_inter2,])$predicted
models_CWM_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
models_CWM_PLH_inter2 <- predict(models_CWM_PLH_grow2, gam_data[group_inter2,])$predicted

CWM_trait_based_inter2<-data.frame("LNC"=models_CWM_LNC_inter2, "SLA"=models_CWM_SLA_inter2, "PLH"=models_CWM_PLH_inter2)
rownames(CWM_trait_based_inter2) <- CWM_obs$Site[group_inter2]

## - For interpolation group 3 - ##
models_CWM_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
models_CWM_LNC_inter3 <- predict(models_CWM_LNC_grow3, gam_data[group_inter3,])$predicted
models_CWM_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
models_CWM_SLA_inter3 <- predict(models_CWM_SLA_grow3, gam_data[group_inter3,])$predicted
models_CWM_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
models_CWM_PLH_inter3 <- predict(models_CWM_PLH_grow3, gam_data[group_inter3,])$predicted

CWM_trait_based_inter3<-data.frame("LNC"=models_CWM_LNC_inter3, "SLA"=models_CWM_SLA_inter3, "PLH"=models_CWM_PLH_inter3)
rownames(CWM_trait_based_inter3) <- CWM_obs$Site[group_inter3]

## - For interpolation group 4 - ##
models_CWM_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
models_CWM_LNC_inter4 <- predict(models_CWM_LNC_grow4, gam_data[group_inter4,])$predicted
models_CWM_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
models_CWM_SLA_inter4 <- predict(models_CWM_SLA_grow4, gam_data[group_inter4,])$predicted
models_CWM_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
models_CWM_PLH_inter4 <- predict(models_CWM_PLH_grow4, gam_data[group_inter4,])$predicted

CWM_trait_based_inter4<-data.frame("LNC"=models_CWM_LNC_inter4, "SLA"=models_CWM_SLA_inter4, "PLH"=models_CWM_PLH_inter4)
rownames(CWM_trait_based_inter4) <- CWM_obs$Site[group_inter4]

## assembly in a single prediction table for the 4 groups :
CWM_trait_based_inter <- rbind(CWM_trait_based_inter1, CWM_trait_based_inter2, CWM_trait_based_inter3, CWM_trait_based_inter4)
CWM_trait_based_inter <- cbind("Site"=rownames(CWM_trait_based_inter), CWM_trait_based_inter)
Site_order <- data.frame("Site"=rownames(plant_recovery))
CWM_trait_based_inter <- left_join(Site_order, CWM_trait_based_inter, by="Site")
rownames(CWM_trait_based_inter) <- CWM_trait_based_inter$Site
CWM_trait_based_inter <- CWM_trait_based_inter[,-1]

#### -- This computes path predicts the community weighted mean in the context of extrapolation -- ####
gam_data <- left_join(CWM_obs,env_var, by="Site")
## - For extrapolation group 1 - ##
models_CWM_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
models_CWM_LNC_extra1 <- predict(models_CWM_LNC_grow1, gam_data[group_extra1,])$predicted
models_CWM_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
models_CWM_SLA_extra1 <- predict(models_CWM_SLA_grow1, gam_data[group_extra1,])$predicted
models_CWM_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
models_CWM_PLH_extra1 <- predict(models_CWM_PLH_grow1, gam_data[group_extra1,])$predicted

CWM_trait_based_extra1<-data.frame("LNC"=models_CWM_LNC_extra1, "SLA"=models_CWM_SLA_extra1, "PLH"=models_CWM_PLH_extra1)
rownames(CWM_trait_based_extra1) <- CWM_obs$Site[group_extra1]

## - For extrapolation group 2 - ##
models_CWM_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
models_CWM_LNC_extra2 <- predict(models_CWM_LNC_grow2, gam_data[group_extra2,])$predicted
models_CWM_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
models_CWM_SLA_extra2 <- predict(models_CWM_SLA_grow2, gam_data[group_extra2,])$predicted
models_CWM_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
models_CWM_PLH_extra2 <- predict(models_CWM_PLH_grow2, gam_data[group_extra2,])$predicted

CWM_trait_based_extra2<-data.frame("LNC"=models_CWM_LNC_extra2, "SLA"=models_CWM_SLA_extra2, "PLH"=models_CWM_PLH_extra2)
rownames(CWM_trait_based_extra2) <- CWM_obs$Site[group_extra2]

## - For extrapolation group 3 - ##
models_CWM_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
models_CWM_LNC_extra3 <- predict(models_CWM_LNC_grow3, gam_data[group_extra3,])$predicted
models_CWM_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
models_CWM_SLA_extra3 <- predict(models_CWM_SLA_grow3, gam_data[group_extra3,])$predicted
models_CWM_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
models_CWM_PLH_extra3 <- predict(models_CWM_PLH_grow3, gam_data[group_extra3,])$predicted

CWM_trait_based_extra3<-data.frame("LNC"=models_CWM_LNC_extra3, "SLA"=models_CWM_SLA_extra3, "PLH"=models_CWM_PLH_extra3)
rownames(CWM_trait_based_extra3) <- CWM_obs$Site[group_extra3]

## - For extrapolation group 4 - ##
models_CWM_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
models_CWM_LNC_extra4 <- predict(models_CWM_LNC_grow4, gam_data[group_extra4,])$predicted
models_CWM_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
models_CWM_SLA_extra4 <- predict(models_CWM_SLA_grow4, gam_data[group_extra4,])$predicted
models_CWM_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
models_CWM_PLH_extra4 <- predict(models_CWM_PLH_grow4, gam_data[group_extra4,])$predicted

CWM_trait_based_extra4<-data.frame("LNC"=models_CWM_LNC_extra4, "SLA"=models_CWM_SLA_extra4, "PLH"=models_CWM_PLH_extra4)
rownames(CWM_trait_based_extra4) <- CWM_obs$Site[group_extra4]

## assembly in a single prediction table for the 4 groups :
CWM_trait_based_extra <- rbind(CWM_trait_based_extra1, CWM_trait_based_extra2, CWM_trait_based_extra3, CWM_trait_based_extra4)
CWM_trait_based_extra <- cbind("Site"=rownames(CWM_trait_based_extra), CWM_trait_based_extra)
Site_order <- data.frame("Site"=rownames(plant_recovery))
CWM_trait_based_extra <- left_join(Site_order, CWM_trait_based_extra, by="Site")
rownames(CWM_trait_based_extra) <- CWM_trait_based_extra$Site
CWM_trait_based_extra <- CWM_trait_based_extra[,-1]

###### PREDICTION OF COMMUNITY MEAN  ######
#### -- This computes path predicts the community mean in the context of interpolation -- ####
## - For interpolation group 1 - ##
models_CM_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
models_CM_LNC_inter1 <- predict(models_CM_LNC_grow1, gam_data[group_inter1,])$predicted
models_CM_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
models_CM_SLA_inter1 <- predict(models_CM_SLA_grow1, gam_data[group_inter1,])$predicted
models_CM_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
models_CM_PLH_inter1 <- predict(models_CM_PLH_grow1, gam_data[group_inter1,])$predicted

CM_trait_based_inter1<-data.frame("LNC"=models_CM_LNC_inter1, "SLA"=models_CM_SLA_inter1, "PLH"=models_CM_PLH_inter1)
rownames(CM_trait_based_inter1) <- CM_obs$Site[group_inter1]

## - For interpolation group 2 - ##
models_CM_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
models_CM_LNC_inter2 <- predict(models_CM_LNC_grow2, gam_data[group_inter2,])$predicted
models_CM_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
models_CM_SLA_inter2 <- predict(models_CM_SLA_grow2, gam_data[group_inter2,])$predicted
models_CM_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
models_CM_PLH_inter2 <- predict(models_CM_PLH_grow2, gam_data[group_inter2,])$predicted

CM_trait_based_inter2<-data.frame("LNC"=models_CM_LNC_inter2, "SLA"=models_CM_SLA_inter2, "PLH"=models_CM_PLH_inter2)
rownames(CM_trait_based_inter2) <- CM_obs$Site[group_inter2]

## - For interpolation group 3 - ##
models_CM_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
models_CM_LNC_inter3 <- predict(models_CM_LNC_grow3, gam_data[group_inter3,])$predicted
models_CM_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
models_CM_SLA_inter3 <- predict(models_CM_SLA_grow3, gam_data[group_inter3,])$predicted
models_CM_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
models_CM_PLH_inter3 <- predict(models_CM_PLH_grow3, gam_data[group_inter3,])$predicted

CM_trait_based_inter3<-data.frame("LNC"=models_CM_LNC_inter3, "SLA"=models_CM_SLA_inter3, "PLH"=models_CM_PLH_inter3)
rownames(CM_trait_based_inter3) <- CM_obs$Site[group_inter3]

## - For interpolation group 4 - ##
models_CM_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
models_CM_LNC_inter4 <- predict(models_CM_LNC_grow4, gam_data[group_inter4,])$predicted
models_CM_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
models_CM_SLA_inter4 <- predict(models_CM_SLA_grow4, gam_data[group_inter4,])$predicted
models_CM_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
models_CM_PLH_inter4 <- predict(models_CM_PLH_grow4, gam_data[group_inter4,])$predicted

CM_trait_based_inter4<-data.frame("LNC"=models_CM_LNC_inter4, "SLA"=models_CM_SLA_inter4, "PLH"=models_CM_PLH_inter4)
rownames(CM_trait_based_inter4) <- CM_obs$Site[group_inter4]

## assembly in a single prediction table for the 4 groups :
CM_trait_based_inter <- rbind(CM_trait_based_inter1, CM_trait_based_inter2, CM_trait_based_inter3, CM_trait_based_inter4)
CM_trait_based_inter <- cbind("Site"=rownames(CM_trait_based_inter), CM_trait_based_inter)
Site_order <- data.frame("Site"=rownames(plant_pa))
CM_trait_based_inter <- left_join(Site_order, CM_trait_based_inter, by="Site")
rownames(CM_trait_based_inter) <- CM_trait_based_inter$Site
CM_trait_based_inter <- CM_trait_based_inter[,-1]
#### -- This computes path predicts the community mean in the context of extrapolation -- ####
gam_data <- left_join(CM_obs,env_var, by="Site")
## - For extrapolation group 1 - ##
models_CM_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
models_CM_LNC_extra1 <- predict(models_CM_LNC_grow1, gam_data[group_extra1,])$predicted
models_CM_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
models_CM_SLA_extra1 <- predict(models_CM_SLA_grow1, gam_data[group_extra1,])$predicted
models_CM_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
models_CM_PLH_extra1 <- predict(models_CM_PLH_grow1, gam_data[group_extra1,])$predicted

CM_trait_based_extra1<-data.frame("LNC"=models_CM_LNC_extra1, "SLA"=models_CM_SLA_extra1, "PLH"=models_CM_PLH_extra1)
rownames(CM_trait_based_extra1) <- CM_obs$Site[group_extra1]

## - For extrapolation group 2 - ##
models_CM_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
models_CM_LNC_extra2 <- predict(models_CM_LNC_grow2, gam_data[group_extra2,])$predicted
models_CM_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
models_CM_SLA_extra2 <- predict(models_CM_SLA_grow2, gam_data[group_extra2,])$predicted
models_CM_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
models_CM_PLH_extra2 <- predict(models_CM_PLH_grow2, gam_data[group_extra2,])$predicted

CM_trait_based_extra2<-data.frame("LNC"=models_CM_LNC_extra2, "SLA"=models_CM_SLA_extra2, "PLH"=models_CM_PLH_extra2)
rownames(CM_trait_based_extra2) <- CM_obs$Site[group_extra2]

## - For extrapolation group 3 - ##
models_CM_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
models_CM_LNC_extra3 <- predict(models_CM_LNC_grow3, gam_data[group_extra3,])$predicted
models_CM_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
models_CM_SLA_extra3 <- predict(models_CM_SLA_grow3, gam_data[group_extra3,])$predicted
models_CM_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
models_CM_PLH_extra3 <- predict(models_CM_PLH_grow3, gam_data[group_extra3,])$predicted

CM_trait_based_extra3<-data.frame("LNC"=models_CM_LNC_extra3, "SLA"=models_CM_SLA_extra3, "PLH"=models_CM_PLH_extra3)
rownames(CM_trait_based_extra3) <- CM_obs$Site[group_extra3]

## - For extrapolation group 4 - ##
models_CM_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
models_CM_LNC_extra4 <- predict(models_CM_LNC_grow4, gam_data[group_extra4,])$predicted
models_CM_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
models_CM_SLA_extra4 <- predict(models_CM_SLA_grow4, gam_data[group_extra4,])$predicted
models_CM_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
models_CM_PLH_extra4 <- predict(models_CM_PLH_grow4, gam_data[group_extra4,])$predicted

CM_trait_based_extra4<-data.frame("LNC"=models_CM_LNC_extra4, "SLA"=models_CM_SLA_extra4, "PLH"=models_CM_PLH_extra4)
rownames(CM_trait_based_extra4) <- CM_obs$Site[group_extra4]

## assembly in a single prediction table for the 4 groups :
CM_trait_based_extra <- rbind(CM_trait_based_extra1, CM_trait_based_extra2, CM_trait_based_extra3, CM_trait_based_extra4)
CM_trait_based_extra <- cbind("Site"=rownames(CM_trait_based_extra), CM_trait_based_extra)
Site_order <- data.frame("Site"=rownames(plant_pa))
CM_trait_based_extra <- left_join(Site_order, CM_trait_based_extra, by="Site")
rownames(CM_trait_based_extra) <- CM_trait_based_extra$Site
CM_trait_based_extra <- CM_trait_based_extra[,-1]

###### PREDICTION OF FUNCTIONAL DISPERSION  ######
#### -- This computes path predicts the functional dispersion in the context of interpolation -- ####
gam_data <- left_join(FDis_obs,env_var, by="Site")
## - For interpolation group 1 - ##
models_FDis_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
models_FDis_LNC_inter1 <- predict(models_FDis_LNC_grow1, gam_data[group_inter1,])$predicted
models_FDis_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
models_FDis_SLA_inter1 <- predict(models_FDis_SLA_grow1, gam_data[group_inter1,])$predicted
models_FDis_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
models_FDis_PLH_inter1 <- predict(models_FDis_PLH_grow1, gam_data[group_inter1,])$predicted

FDis_trait_based_inter1<-data.frame("LNC"=models_FDis_LNC_inter1, "SLA"=models_FDis_SLA_inter1, "PLH"=models_FDis_PLH_inter1)
rownames(FDis_trait_based_inter1) <- FDis_obs$Site[group_inter1]

## - For interpolation group 2 - ##
models_FDis_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
models_FDis_LNC_inter2 <- predict(models_FDis_LNC_grow2, gam_data[group_inter2,])$predicted
models_FDis_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
models_FDis_SLA_inter2 <- predict(models_FDis_SLA_grow2, gam_data[group_inter2,])$predicted
models_FDis_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
models_FDis_PLH_inter2 <- predict(models_FDis_PLH_grow2, gam_data[group_inter2,])$predicted

FDis_trait_based_inter2<-data.frame("LNC"=models_FDis_LNC_inter2, "SLA"=models_FDis_SLA_inter2, "PLH"=models_FDis_PLH_inter2)
rownames(FDis_trait_based_inter2) <- FDis_obs$Site[group_inter2]

## - For interpolation group 3 - ##
models_FDis_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
models_FDis_LNC_inter3 <- predict(models_FDis_LNC_grow3, gam_data[group_inter3,])$predicted
models_FDis_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
models_FDis_SLA_inter3 <- predict(models_FDis_SLA_grow3, gam_data[group_inter3,])$predicted
models_FDis_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
models_FDis_PLH_inter3 <- predict(models_FDis_PLH_grow3, gam_data[group_inter3,])$predicted

FDis_trait_based_inter3<-data.frame("LNC"=models_FDis_LNC_inter3, "SLA"=models_FDis_SLA_inter3, "PLH"=models_FDis_PLH_inter3)
rownames(FDis_trait_based_inter3) <- FDis_obs$Site[group_inter3]

## - For interpolation group 4 - ##
models_FDis_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
models_FDis_LNC_inter4 <- predict(models_FDis_LNC_grow4, gam_data[group_inter4,])$predicted
models_FDis_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
models_FDis_SLA_inter4 <- predict(models_FDis_SLA_grow4, gam_data[group_inter4,])$predicted
models_FDis_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
models_FDis_PLH_inter4 <- predict(models_FDis_PLH_grow4, gam_data[group_inter4,])$predicted

FDis_trait_based_inter4<-data.frame("LNC"=models_FDis_LNC_inter4, "SLA"=models_FDis_SLA_inter4, "PLH"=models_FDis_PLH_inter4)
rownames(FDis_trait_based_inter4) <- FDis_obs$Site[group_inter4]

## assembly in a single prediction table for the 4 groups :
FDis_trait_based_inter <- rbind(FDis_trait_based_inter1, FDis_trait_based_inter2, FDis_trait_based_inter3, FDis_trait_based_inter4)
FDis_trait_based_inter <- cbind("Site"=rownames(FDis_trait_based_inter), FDis_trait_based_inter)
Site_order <- data.frame("Site"=rownames(plant_recovery))
FDis_trait_based_inter <- left_join(Site_order, FDis_trait_based_inter, by="Site")
rownames(FDis_trait_based_inter) <- FDis_trait_based_inter$Site
FDis_trait_based_inter <- FDis_trait_based_inter[,-1]

#### -- This computes path predicts the functional dispersion in the context of extrapolation -- ####
gam_data <- left_join(FDis_obs,env_var, by="Site")
## - For extrapolation group 1 - ##
models_FDis_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
models_FDis_LNC_extra1 <- predict(models_FDis_LNC_grow1, gam_data[group_extra1,])$predicted
models_FDis_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
models_FDis_SLA_extra1 <- predict(models_FDis_SLA_grow1, gam_data[group_extra1,])$predicted
models_FDis_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
models_FDis_PLH_extra1 <- predict(models_FDis_PLH_grow1, gam_data[group_extra1,])$predicted

FDis_trait_based_extra1<-data.frame("LNC"=models_FDis_LNC_extra1, "SLA"=models_FDis_SLA_extra1, "PLH"=models_FDis_PLH_extra1)
rownames(FDis_trait_based_extra1) <- FDis_obs$Site[group_extra1]

## - For extrapolation group 2 - ##
models_FDis_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
models_FDis_LNC_extra2 <- predict(models_FDis_LNC_grow2, gam_data[group_extra2,])$predicted
models_FDis_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
models_FDis_SLA_extra2 <- predict(models_FDis_SLA_grow2, gam_data[group_extra2,])$predicted
models_FDis_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
models_FDis_PLH_extra2 <- predict(models_FDis_PLH_grow2, gam_data[group_extra2,])$predicted

FDis_trait_based_extra2<-data.frame("LNC"=models_FDis_LNC_extra2, "SLA"=models_FDis_SLA_extra2, "PLH"=models_FDis_PLH_extra2)
rownames(FDis_trait_based_extra2) <- FDis_obs$Site[group_extra2]

## - For extrapolation group 3 - ##
models_FDis_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
models_FDis_LNC_extra3 <- predict(models_FDis_LNC_grow3, gam_data[group_extra3,])$predicted
models_FDis_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
models_FDis_SLA_extra3 <- predict(models_FDis_SLA_grow3, gam_data[group_extra3,])$predicted
models_FDis_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
models_FDis_PLH_extra3 <- predict(models_FDis_PLH_grow3, gam_data[group_extra3,])$predicted

FDis_trait_based_extra3<-data.frame("LNC"=models_FDis_LNC_extra3, "SLA"=models_FDis_SLA_extra3, "PLH"=models_FDis_PLH_extra3)
rownames(FDis_trait_based_extra3) <- FDis_obs$Site[group_extra3]

## - For extrapolation group 4 - ##
models_FDis_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
models_FDis_LNC_extra4 <- predict(models_FDis_LNC_grow4, gam_data[group_extra4,])$predicted
models_FDis_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
models_FDis_SLA_extra4 <- predict(models_FDis_SLA_grow4, gam_data[group_extra4,])$predicted
models_FDis_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
models_FDis_PLH_extra4 <- predict(models_FDis_PLH_grow4, gam_data[group_extra4,])$predicted

FDis_trait_based_extra4<-data.frame("LNC"=models_FDis_LNC_extra4, "SLA"=models_FDis_SLA_extra4, "PLH"=models_FDis_PLH_extra4)
rownames(FDis_trait_based_extra4) <- FDis_obs$Site[group_extra4]

## assembly in a single prediction table for the 4 groups :
FDis_trait_based_extra <- rbind(FDis_trait_based_extra1, FDis_trait_based_extra2, FDis_trait_based_extra3, FDis_trait_based_extra4)
FDis_trait_based_extra <- cbind("Site"=rownames(FDis_trait_based_extra), FDis_trait_based_extra)
Site_order <- data.frame("Site"=rownames(plant_recovery))
FDis_trait_based_extra <- left_join(Site_order, FDis_trait_based_extra, by="Site")
rownames(FDis_trait_based_extra) <- FDis_trait_based_extra$Site
FDis_trait_based_extra <- FDis_trait_based_extra[,-1]

###### PREDICTION OF UNWEIGHTED FUNCTIONAL DISPERSION  ######
#### -- This computes path predicts the unweighted functional dispersion in the context of interpolation -- ####
gam_data <- left_join(uFDis_obs,env_var, by="Site")
## - For interpolation group 1 - ##
models_uFDis_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
models_uFDis_LNC_inter1 <- predict(models_uFDis_LNC_grow1, gam_data[group_inter1,])$predicted
models_uFDis_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
models_uFDis_SLA_inter1 <- predict(models_uFDis_SLA_grow1, gam_data[group_inter1,])$predicted
models_uFDis_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
models_uFDis_PLH_inter1 <- predict(models_uFDis_PLH_grow1, gam_data[group_inter1,])$predicted

uFDis_trait_based_inter1<-data.frame("LNC"=models_uFDis_LNC_inter1, "SLA"=models_uFDis_SLA_inter1, "PLH"=models_uFDis_PLH_inter1)
rownames(uFDis_trait_based_inter1) <- uFDis_obs$Site[group_inter1]

## - For interpolation group 2 - ##
models_uFDis_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
models_uFDis_LNC_inter2 <- predict(models_uFDis_LNC_grow2, gam_data[group_inter2,])$predicted
models_uFDis_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
models_uFDis_SLA_inter2 <- predict(models_uFDis_SLA_grow2, gam_data[group_inter2,])$predicted
models_uFDis_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
models_uFDis_PLH_inter2 <- predict(models_uFDis_PLH_grow2, gam_data[group_inter2,])$predicted

uFDis_trait_based_inter2<-data.frame("LNC"=models_uFDis_LNC_inter2, "SLA"=models_uFDis_SLA_inter2, "PLH"=models_uFDis_PLH_inter2)
rownames(uFDis_trait_based_inter2) <- uFDis_obs$Site[group_inter2]

## - For interpolation group 3 - ##
models_uFDis_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
models_uFDis_LNC_inter3 <- predict(models_uFDis_LNC_grow3, gam_data[group_inter3,])$predicted
models_uFDis_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
models_uFDis_SLA_inter3 <- predict(models_uFDis_SLA_grow3, gam_data[group_inter3,])$predicted
models_uFDis_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
models_uFDis_PLH_inter3 <- predict(models_uFDis_PLH_grow3, gam_data[group_inter3,])$predicted

uFDis_trait_based_inter3<-data.frame("LNC"=models_uFDis_LNC_inter3, "SLA"=models_uFDis_SLA_inter3, "PLH"=models_uFDis_PLH_inter3)
rownames(uFDis_trait_based_inter3) <- uFDis_obs$Site[group_inter3]

## - For interpolation group 4 - ##
models_uFDis_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
models_uFDis_LNC_inter4 <- predict(models_uFDis_LNC_grow4, gam_data[group_inter4,])$predicted
models_uFDis_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
models_uFDis_SLA_inter4 <- predict(models_uFDis_SLA_grow4, gam_data[group_inter4,])$predicted
models_uFDis_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
models_uFDis_PLH_inter4 <- predict(models_uFDis_PLH_grow4, gam_data[group_inter4,])$predicted

uFDis_trait_based_inter4<-data.frame("LNC"=models_uFDis_LNC_inter4, "SLA"=models_uFDis_SLA_inter4, "PLH"=models_uFDis_PLH_inter4)
rownames(uFDis_trait_based_inter4) <- uFDis_obs$Site[group_inter4]

## assembly in a single prediction table for the 4 groups :
uFDis_trait_based_inter <- rbind(uFDis_trait_based_inter1, uFDis_trait_based_inter2, uFDis_trait_based_inter3, uFDis_trait_based_inter4)
uFDis_trait_based_inter <- cbind("Site"=rownames(uFDis_trait_based_inter), uFDis_trait_based_inter)
Site_order <- data.frame("Site"=rownames(plant_pa))
uFDis_trait_based_inter <- left_join(Site_order, uFDis_trait_based_inter, by="Site")
rownames(uFDis_trait_based_inter) <- uFDis_trait_based_inter$Site
uFDis_trait_based_inter <- uFDis_trait_based_inter[,-1]

#### -- This computes path predicts the community mean in the context of extrapolation -- ####
## - For extrapolation group 1 - ##
models_uFDis_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
models_uFDis_LNC_extra1 <- predict(models_uFDis_LNC_grow1, gam_data[group_extra1,])$predicted
models_uFDis_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
models_uFDis_SLA_extra1 <- predict(models_uFDis_SLA_grow1, gam_data[group_extra1,])$predicted
models_uFDis_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
models_uFDis_PLH_extra1 <- predict(models_uFDis_PLH_grow1, gam_data[group_extra1,])$predicted

uFDis_trait_based_extra1<-data.frame("LNC"=models_uFDis_LNC_extra1, "SLA"=models_uFDis_SLA_extra1, "PLH"=models_uFDis_PLH_extra1)
rownames(uFDis_trait_based_extra1) <- uFDis_obs$Site[group_extra1]

## - For extrapolation group 2 - ##
models_uFDis_LNC_grow2 <- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
models_uFDis_LNC_extra2 <- predict(models_uFDis_LNC_grow2, gam_data[group_extra2,])$predicted
models_uFDis_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
models_uFDis_SLA_extra2 <- predict(models_uFDis_SLA_grow2, gam_data[group_extra2,])$predicted
models_uFDis_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
models_uFDis_PLH_extra2 <- predict(models_uFDis_PLH_grow2, gam_data[group_extra2,])$predicted

uFDis_trait_based_extra2<-data.frame("LNC"=models_uFDis_LNC_extra2, "SLA"=models_uFDis_SLA_extra2, "PLH"=models_uFDis_PLH_extra2)
rownames(uFDis_trait_based_extra2) <- uFDis_obs$Site[group_extra2]

## - For extrapolation group 3 - ##
models_uFDis_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
models_uFDis_LNC_extra3 <- predict(models_uFDis_LNC_grow3, gam_data[group_extra3,])$predicted
models_uFDis_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
models_uFDis_SLA_extra3 <- predict(models_uFDis_SLA_grow3, gam_data[group_extra3,])$predicted
models_uFDis_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
models_uFDis_PLH_extra3 <- predict(models_uFDis_PLH_grow3, gam_data[group_extra3,])$predicted

uFDis_trait_based_extra3<-data.frame("LNC"=models_uFDis_LNC_extra3, "SLA"=models_uFDis_SLA_extra3, "PLH"=models_uFDis_PLH_extra3)
rownames(uFDis_trait_based_extra3) <- uFDis_obs$Site[group_extra3]

## - For extrapolation group 4 - ##
models_uFDis_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
models_uFDis_LNC_extra4 <- predict(models_uFDis_LNC_grow4, gam_data[group_extra4,])$predicted
models_uFDis_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
models_uFDis_SLA_extra4 <- predict(models_uFDis_SLA_grow4, gam_data[group_extra4,])$predicted
models_uFDis_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
models_uFDis_PLH_extra4 <- predict(models_uFDis_PLH_grow4, gam_data[group_extra4,])$predicted

uFDis_trait_based_extra4<-data.frame("LNC"=models_uFDis_LNC_extra4, "SLA"=models_uFDis_SLA_extra4, "PLH"=models_uFDis_PLH_extra4)
rownames(uFDis_trait_based_extra4) <- uFDis_obs$Site[group_extra4]

## assembly in a single prediction table for the 4 groups :
uFDis_trait_based_extra <- rbind(uFDis_trait_based_extra1, uFDis_trait_based_extra2, uFDis_trait_based_extra3, uFDis_trait_based_extra4)
uFDis_trait_based_extra <- cbind("Site"=rownames(uFDis_trait_based_extra), uFDis_trait_based_extra)
Site_order <- data.frame("Site"=rownames(plant_pa))
uFDis_trait_based_extra <- left_join(Site_order, uFDis_trait_based_extra, by="Site")
rownames(uFDis_trait_based_extra) <- uFDis_trait_based_extra$Site
uFDis_trait_based_extra <- uFDis_trait_based_extra[,-1]

################################################################
#### -- CONVERGENCE BETWEEN THE TWO DIFFERENT APPROACHES -- ####
################################################################
#### -- Figure 2 : Functional indices of plant height predicted by trait-based approach plotted against those predicted by species-based approach ####
Figure_2 <- data.frame("Traits"=rep(c(rep("LNC", 4465),rep("SLA", 4465), rep("PLH", 4465)),24),
                         "Assemble-First Approach"=c(CWM_trait_based_inter$LNC,
                                                     28,18.0,
                                                   CWM_trait_based_inter$SLA,
                                                   28,14,
                                                   CWM_trait_based_inter$PLH,
                                                   45,5,
                                                   CM_trait_based_inter$LNC,
                                                   28,18.0,
                                                   CM_trait_based_inter$SLA,
                                                   28,14,
                                                   CM_trait_based_inter$PLH,
                                                   45,5,
                                                   FDis_trait_based_inter$LNC,
                                                   8,3,
                                                   FDis_trait_based_inter$SLA,
                                                   6.5,2,
                                                   FDis_trait_based_inter$PLH,
                                                   20,2,
                                                   uFDis_trait_based_inter$LNC,
                                                   8,3,
                                                   uFDis_trait_based_inter$SLA,
                                                   6.5,2,
                                                   uFDis_trait_based_inter$PLH,
                                                   20,2,
                                                   CWM_trait_based_extra$LNC,
                                                   28,18.0,
                                                   CWM_trait_based_extra$SLA,
                                                   28,14,
                                                   CWM_trait_based_extra$PLH,
                                                   45,5,
                                                   CM_trait_based_extra$LNC,
                                                   28,18.0,
                                                   CM_trait_based_extra$SLA,
                                                   28,14,
                                                   CM_trait_based_extra$PLH,
                                                   45,5,
                                                   FDis_trait_based_extra$LNC,
                                                   8,3,
                                                   FDis_trait_based_extra$SLA,
                                                   6.5,2,
                                                   FDis_trait_based_extra$PLH,
                                                   20,2,
                                                   uFDis_trait_based_extra$LNC,
                                                   8,3,
                                                   uFDis_trait_based_extra$SLA,
                                                   6.5,2,
                                                   uFDis_trait_based_extra$PLH,
                                                   20,2),
                         "Predict-First Approach"=c(CWM_species_based_inter$LNC,
                                                    28,19,
                                                    CWM_species_based_inter$SLA,
                                                    26,14,
                                                    CWM_species_based_inter$PLH,
                                                    45,5,
                                                    CM_species_based_inter$LNC,
                                                    28,19,
                                                    CM_species_based_inter$SLA,
                                                    26,14,
                                                    CM_species_based_inter$PLH,
                                                    45,5,
                                                    FDis_species_based_inter$LNC,
                                                    8,3.5,
                                                    FDis_species_based_inter$SLA,
                                                    7.5,2.5,
                                                    FDis_species_based_inter$PLH,
                                                    20,2,
                                                    uFDis_species_based_inter$LNC,
                                                    8,3.5,
                                                    uFDis_species_based_inter$SLA,
                                                    7.5,2.5,
                                                    uFDis_species_based_inter$PLH,
                                                    20,2,
                                                    CWM_species_based_extra$LNC,
                                                    28,19,
                                                    CWM_species_based_extra$SLA,
                                                    26,14,
                                                    CWM_species_based_extra$PLH,
                                                    45,5,
                                                    CM_species_based_extra$LNC,
                                                    28,19,
                                                    CM_species_based_extra$SLA,
                                                    26,14,
                                                    CM_species_based_extra$PLH,
                                                    45,5,
                                                    FDis_species_based_extra$LNC,
                                                    8,3.5,
                                                    FDis_species_based_extra$SLA,
                                                    7.5,2.5,
                                                    FDis_species_based_extra$PLH,
                                                    20,2,
                                                    uFDis_species_based_extra$LNC,
                                                    8,3.5,
                                                    uFDis_species_based_extra$SLA,
                                                    7.5,2.5,
                                                    uFDis_species_based_extra$PLH,
                                                    20,2),
                         "Indice"=rep(c(rep("CWM/CM", 4465*6), rep("FDis/uFDis", 4465*6)),2),
                         "Data"=rep(c(rep("Abundance data", 4465*3), rep("Occurrence data", 4465*3)),4),
                         "Crossvalidation"=c(rep("Interpolation",4465*12), rep("Extrapolation",4465*12)),
                         "alpha"=rep(c(rep(1,4463),c(0,0)),12)
)

Figure_2$Type <- fct_relevel(Figure_2$Data, c("Abundance data", "Occurrence data")) 
Figure_2$Crossvalidation <- factor(Figure_2$Crossvalidation, levels = c("Interpolation",  "Extrapolation"))
#### -- R2 between the predictions of the two approaches ####
cor(CWM_trait_based_inter$LNC, CWM_species_based_inter$LNC)^2
cor(CWM_trait_based_inter$SLA, CWM_species_based_inter$SLA)^2
cor(CWM_trait_based_inter$PLH, CWM_species_based_inter$PLH)^2

cor(CWM_trait_based_extra$LNC, CWM_species_based_extra$LNC)^2
cor(CWM_trait_based_extra$SLA, CWM_species_based_extra$SLA)^2
cor(CWM_trait_based_extra$PLH, CWM_species_based_extra$PLH)^2

cor(CM_trait_based_inter$LNC, CM_species_based_inter$LNC)^2
cor(CM_trait_based_inter$SLA, CM_species_based_inter$SLA)^2
cor(CM_trait_based_inter$PLH, CM_species_based_inter$PLH)^2

cor(CM_trait_based_extra$LNC, CM_species_based_extra$LNC)^2
cor(CM_trait_based_extra$SLA, CM_species_based_extra$SLA)^2
cor(CM_trait_based_extra$PLH, CM_species_based_extra$PLH)^2

cor(FDis_trait_based_inter$LNC, FDis_species_based_inter$LNC)^2
cor(FDis_trait_based_inter$SLA, FDis_species_based_inter$SLA)^2
cor(FDis_trait_based_inter$PLH, FDis_species_based_inter$PLH)^2

cor(FDis_trait_based_extra$LNC, FDis_species_based_extra$LNC)^2
cor(FDis_trait_based_extra$SLA, FDis_species_based_extra$SLA)^2
cor(FDis_trait_based_extra$PLH, FDis_species_based_extra$PLH)^2

cor(uFDis_trait_based_inter$LNC, uFDis_species_based_inter$LNC)^2
cor(uFDis_trait_based_inter$SLA, uFDis_species_based_inter$SLA)^2
cor(uFDis_trait_based_inter$PLH, uFDis_species_based_inter$PLH)^2

cor(uFDis_trait_based_extra$LNC, uFDis_species_based_extra$LNC)^2
cor(uFDis_trait_based_extra$SLA, uFDis_species_based_extra$SLA)^2
cor(uFDis_trait_based_extra$PLH, uFDis_species_based_extra$PLH)^2

#### -- RMSE between the predictions of the two approaches ####
rmse(CWM_trait_based_inter$LNC, CWM_species_based_inter$LNC)
rmse(CWM_trait_based_inter$SLA, CWM_species_based_inter$SLA)
rmse(CWM_trait_based_inter$PLH, CWM_species_based_inter$PLH)

rmse(CWM_trait_based_extra$LNC, CWM_species_based_extra$LNC)
rmse(CWM_trait_based_extra$SLA, CWM_species_based_extra$SLA)
rmse(CWM_trait_based_extra$PLH, CWM_species_based_extra$PLH)

rmse(CM_trait_based_inter$LNC, CM_species_based_inter$LNC)
rmse(CM_trait_based_inter$SLA, CM_species_based_inter$SLA)
rmse(CM_trait_based_inter$PLH, CM_species_based_inter$PLH)

rmse(CM_trait_based_extra$LNC, CM_species_based_extra$LNC)
rmse(CM_trait_based_extra$SLA, CM_species_based_extra$SLA)
rmse(CM_trait_based_extra$PLH, CM_species_based_extra$PLH)

rmse(FDis_trait_based_inter$LNC, FDis_species_based_inter$LNC)
rmse(FDis_trait_based_inter$SLA, FDis_species_based_inter$SLA)
rmse(FDis_trait_based_inter$PLH, FDis_species_based_inter$PLH)

rmse(FDis_trait_based_extra$LNC, FDis_species_based_extra$LNC)
rmse(FDis_trait_based_extra$SLA, FDis_species_based_extra$SLA)
rmse(FDis_trait_based_extra$PLH, FDis_species_based_extra$PLH)

rmse(uFDis_trait_based_inter$LNC, uFDis_species_based_inter$LNC)
rmse(uFDis_trait_based_inter$SLA, uFDis_species_based_inter$SLA)
rmse(uFDis_trait_based_inter$PLH, uFDis_species_based_inter$PLH)

rmse(uFDis_trait_based_extra$LNC, uFDis_species_based_extra$LNC)
rmse(uFDis_trait_based_extra$SLA, uFDis_species_based_extra$SLA)
rmse(uFDis_trait_based_extra$PLH, uFDis_species_based_extra$PLH)


data_lab_PLH <- data.frame(x=c(14, 14, 14, 14, 6,6,6,6), 
                           y=c(42,42,42, 42, 18.5,18.5,18.5,18.5), 
                           lab=c("R2 = 0.97",
                                 "R2 = 0.98",
                                 "R2 = 0.90",
                                 "R2 = 0.94",
                                 "R2 = 0.78",
                                 "R2 = 0.92",
                                 "R2 = 0.63",
                                 "R2 = 0.84"),
                           Indice=c("CWM/CM","CWM/CM","CWM/CM","CWM/CM","FDis/uFDis","FDis/uFDis","FDis/uFDis","FDis/uFDis"),
                           Type = rep(c("Abundance data", "Occurrence data"),4),
                           Crossvalidation = factor(rep(c(rep("Interpolation", 2), rep("Extrapolation",2)),2), levels = c("Interpolation","Extrapolation")))

data_lab_SLA <- data.frame(x=c(16.5, 16.5, 16.5, 16.5, 3.60,3.6,3.6,3.6), 
                           y=c(26.85,26.85,26.85, 26.85, 6.15,6.15,6.15,6.15), 
                           lab=c("R2 = 0.95",
                                 "R2 = 0.97",
                                 "R2 = 0.75",
                                 "R2 = 0.91",
                                 "R2 = 0.53",
                                 "R2 = 0.60",
                                 "R2 = 0.17",
                                 "R2 = 0.25"),
                           Indice=c("CWM/CM","CWM/CM","CWM/CM","CWM/CM","FDis/uFDis","FDis/uFDis","FDis/uFDis","FDis/uFDis"),
                           Type = rep(c("Abundance data", "Occurrence data"),4),
                           Crossvalidation = factor(rep(c(rep("Interpolation", 2), rep("Extrapolation",2)),2), levels = c("Interpolation","Extrapolation")))

data_lab_LNC <- data.frame(x=c(21, 21, 21, 21, 4.5,4.5,4.5,4.5), 
                           y=c(27.25,27.25,27.25, 27.25, 7.6,7.6,7.6,7.6), 
                           lab=c("R2 = 0.71",
                                 "R2 = 0.81",
                                 "R2 = 0.08",
                                 "R2 = 0.55",
                                 "R2 = 0.60",
                                 "R2 = 0.76",
                                 "R2 = 0.14",
                                 "R2 = 0.50"),
                           Indice=c("CWM/CM","CWM/CM","CWM/CM","CWM/CM","FDis/uFDis","FDis/uFDis","FDis/uFDis","FDis/uFDis"),
                           Type = rep(c("Abundance data", "Occurrence data"),4),
                           Crossvalidation = factor(rep(c(rep("Interpolation", 2), rep("Extrapolation",2)),2), levels = c("Interpolation","Extrapolation")))



Figure_2 %>%
  filter(Traits %in% "PLH") %>%
  filter(alpha == 1) %>%
  ggplot() +
  aes(x = Predict.First.Approach, y = Assemble.First.Approach) +
  geom_hex()+
  scale_fill_gradient(low = "grey", high = "black")+
  xlab("Predict-First Approach")+
  ylab("Assemble-First Approach")+
  geom_point(data = Figure_2[which(Figure_2$Traits=="PLH"& Figure_2$alpha==0),], aes(x = Predict.First.Approach, y = Assemble.First.Approach), alpha=0)+
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    strip.text.x = element_text(size = 18),
    strip.text.y = element_text(size = 18)) + 
  geom_abline(slope=1, intercept = 0)+
  facet_nested(Type ~ Indice + Crossvalidation, scales = "free", independent = "y")+
  geom_label(data=data_lab_PLH, mapping=aes(x=x, y=y, label=lab), size=3)


##############################################################################
#### -- COMPARAISON BETWEEN OBSERVED AND PREDICTED FUNCTIONAL INDICES --  ####
##############################################################################
figure_3 <- data.frame("rmse"=c(rmse(CWM_obs$LNC[group_inter1], CWM_trait_based_inter$LNC[group_inter1]), # CWM trait based approach interpolation
                              rmse(CWM_obs$SLA[group_inter1], CWM_trait_based_inter$SLA[group_inter1]),
                              rmse(CWM_obs$PLH[group_inter1], CWM_trait_based_inter$PLH[group_inter1]),
                              rmse(CWM_obs$LNC[group_inter2], CWM_trait_based_inter$LNC[group_inter2]),
                              rmse(CWM_obs$SLA[group_inter2], CWM_trait_based_inter$SLA[group_inter2]),
                              rmse(CWM_obs$PLH[group_inter2], CWM_trait_based_inter$PLH[group_inter2]),
                              rmse(CWM_obs$LNC[group_inter3], CWM_trait_based_inter$LNC[group_inter3]),
                              rmse(CWM_obs$SLA[group_inter3], CWM_trait_based_inter$SLA[group_inter3]),
                              rmse(CWM_obs$PLH[group_inter3], CWM_trait_based_inter$PLH[group_inter3]),
                              rmse(CWM_obs$LNC[group_inter4], CWM_trait_based_inter$LNC[group_inter4]),
                              rmse(CWM_obs$SLA[group_inter4], CWM_trait_based_inter$SLA[group_inter4]),
                              rmse(CWM_obs$PLH[group_inter4], CWM_trait_based_inter$PLH[group_inter4]),
                              rmse(CWM_obs$LNC[group_extra1], CWM_trait_based_extra$LNC[group_extra1]),# CWM trait based approach interpolation
                              rmse(CWM_obs$SLA[group_extra1], CWM_trait_based_extra$SLA[group_extra1]),
                              rmse(CWM_obs$PLH[group_extra1], CWM_trait_based_extra$PLH[group_extra1]),
                              rmse(CWM_obs$LNC[group_extra2], CWM_trait_based_extra$LNC[group_extra2]),
                              rmse(CWM_obs$SLA[group_extra2], CWM_trait_based_extra$SLA[group_extra2]),
                              rmse(CWM_obs$PLH[group_extra2], CWM_trait_based_extra$PLH[group_extra2]),
                              rmse(CWM_obs$LNC[group_extra3], CWM_trait_based_extra$LNC[group_extra3]),
                              rmse(CWM_obs$SLA[group_extra3], CWM_trait_based_extra$SLA[group_extra3]),
                              rmse(CWM_obs$PLH[group_extra3], CWM_trait_based_extra$PLH[group_extra3]),
                              rmse(CWM_obs$LNC[group_extra4], CWM_trait_based_extra$LNC[group_extra4]),
                              rmse(CWM_obs$SLA[group_extra4], CWM_trait_based_extra$SLA[group_extra4]),
                              rmse(CWM_obs$PLH[group_extra4], CWM_trait_based_extra$PLH[group_extra4]),
                              rmse(FDis_obs$LNC[group_inter1], FDis_trait_based_inter$LNC[group_inter1]),# FDis trait based approach interpolation
                              rmse(FDis_obs$SLA[group_inter1], FDis_trait_based_inter$SLA[group_inter1]),
                              rmse(FDis_obs$PLH[group_inter1], FDis_trait_based_inter$PLH[group_inter1]),
                              rmse(FDis_obs$LNC[group_inter2], FDis_trait_based_inter$LNC[group_inter2]),
                              rmse(FDis_obs$SLA[group_inter2], FDis_trait_based_inter$SLA[group_inter2]),
                              rmse(FDis_obs$PLH[group_inter2], FDis_trait_based_inter$PLH[group_inter2]),
                              rmse(FDis_obs$LNC[group_inter3], FDis_trait_based_inter$LNC[group_inter3]),
                              rmse(FDis_obs$SLA[group_inter3], FDis_trait_based_inter$SLA[group_inter3]),
                              rmse(FDis_obs$PLH[group_inter3], FDis_trait_based_inter$PLH[group_inter3]),
                              rmse(FDis_obs$LNC[group_inter4], FDis_trait_based_inter$LNC[group_inter4]),
                              rmse(FDis_obs$SLA[group_inter4], FDis_trait_based_inter$SLA[group_inter4]),
                              rmse(FDis_obs$PLH[group_inter4], FDis_trait_based_inter$PLH[group_inter4]),
                              rmse(FDis_obs$LNC[group_extra1], FDis_trait_based_extra$LNC[group_extra1]), # FDis trait based approach extrapolation
                              rmse(FDis_obs$SLA[group_extra1], FDis_trait_based_extra$SLA[group_extra1]),
                              rmse(FDis_obs$PLH[group_extra1], FDis_trait_based_extra$PLH[group_extra1]),
                              rmse(FDis_obs$LNC[group_extra2], FDis_trait_based_extra$LNC[group_extra2]),
                              rmse(FDis_obs$SLA[group_extra2], FDis_trait_based_extra$SLA[group_extra2]),
                              rmse(FDis_obs$PLH[group_extra2], FDis_trait_based_extra$PLH[group_extra2]),
                              rmse(FDis_obs$LNC[group_extra3], FDis_trait_based_extra$LNC[group_extra3]),
                              rmse(FDis_obs$SLA[group_extra3], FDis_trait_based_extra$SLA[group_extra3]),
                              rmse(FDis_obs$PLH[group_extra3], FDis_trait_based_extra$PLH[group_extra3]),
                              rmse(FDis_obs$LNC[group_extra4], FDis_trait_based_extra$LNC[group_extra4]),
                              rmse(FDis_obs$SLA[group_extra4], FDis_trait_based_extra$SLA[group_extra4]),
                              rmse(FDis_obs$PLH[group_extra4], FDis_trait_based_extra$PLH[group_extra4]),
                              rmse(CWM_obs$LNC[group_inter1], CWM_species_based_inter$LNC[group_inter1]), # CWM species based approach interpolation
                              rmse(CWM_obs$SLA[group_inter1], CWM_species_based_inter$SLA[group_inter1]),
                              rmse(CWM_obs$PLH[group_inter1], CWM_species_based_inter$PLH[group_inter1]),
                              rmse(CWM_obs$LNC[group_inter2], CWM_species_based_inter$LNC[group_inter2]),
                              rmse(CWM_obs$SLA[group_inter2], CWM_species_based_inter$SLA[group_inter2]),
                              rmse(CWM_obs$PLH[group_inter2], CWM_species_based_inter$PLH[group_inter2]),
                              rmse(CWM_obs$LNC[group_inter3], CWM_species_based_inter$LNC[group_inter3]),
                              rmse(CWM_obs$SLA[group_inter3], CWM_species_based_inter$SLA[group_inter3]),
                              rmse(CWM_obs$PLH[group_inter3], CWM_species_based_inter$PLH[group_inter3]),
                              rmse(CWM_obs$LNC[group_inter4], CWM_species_based_inter$LNC[group_inter4]),
                              rmse(CWM_obs$SLA[group_inter4], CWM_species_based_inter$SLA[group_inter4]),
                              rmse(CWM_obs$PLH[group_inter4], CWM_species_based_inter$PLH[group_inter4]),
                              rmse(CWM_obs$LNC[group_extra1], CWM_species_based_extra$LNC[group_extra1]), # CWM species based approach extrapolation
                              rmse(CWM_obs$SLA[group_extra1], CWM_species_based_extra$SLA[group_extra1]),
                              rmse(CWM_obs$PLH[group_extra1], CWM_species_based_extra$PLH[group_extra1]),
                              rmse(CWM_obs$LNC[group_extra2], CWM_species_based_extra$LNC[group_extra2]),
                              rmse(CWM_obs$SLA[group_extra2], CWM_species_based_extra$SLA[group_extra2]),
                              rmse(CWM_obs$PLH[group_extra2], CWM_species_based_extra$PLH[group_extra2]),
                              rmse(CWM_obs$LNC[group_extra3], CWM_species_based_extra$LNC[group_extra3]),
                              rmse(CWM_obs$SLA[group_extra3], CWM_species_based_extra$SLA[group_extra3]),
                              rmse(CWM_obs$PLH[group_extra3], CWM_species_based_extra$PLH[group_extra3]),
                              rmse(CWM_obs$LNC[group_extra4], CWM_species_based_extra$LNC[group_extra4]),
                              rmse(CWM_obs$SLA[group_extra4], CWM_species_based_extra$SLA[group_extra4]),
                              rmse(CWM_obs$PLH[group_extra4], CWM_species_based_extra$PLH[group_extra4]),
                              rmse(FDis_obs$LNC[group_inter1], FDis_species_based_inter$LNC[group_inter1]), # FDis species based approach interpolation
                              rmse(FDis_obs$SLA[group_inter1], FDis_species_based_inter$SLA[group_inter1]),
                              rmse(FDis_obs$PLH[group_inter1], FDis_species_based_inter$PLH[group_inter1]),
                              rmse(FDis_obs$LNC[group_inter2], FDis_species_based_inter$LNC[group_inter2]),
                              rmse(FDis_obs$SLA[group_inter2], FDis_species_based_inter$SLA[group_inter2]),
                              rmse(FDis_obs$PLH[group_inter2], FDis_species_based_inter$PLH[group_inter2]),
                              rmse(FDis_obs$LNC[group_inter3], FDis_species_based_inter$LNC[group_inter3]),
                              rmse(FDis_obs$SLA[group_inter3], FDis_species_based_inter$SLA[group_inter3]),
                              rmse(FDis_obs$PLH[group_inter3], FDis_species_based_inter$PLH[group_inter3]),
                              rmse(FDis_obs$LNC[group_inter4], FDis_species_based_inter$LNC[group_inter4]),
                              rmse(FDis_obs$SLA[group_inter4], FDis_species_based_inter$SLA[group_inter4]),
                              rmse(FDis_obs$PLH[group_inter4], FDis_species_based_inter$PLH[group_inter4]),
                              rmse(FDis_obs$LNC[group_extra1], FDis_species_based_extra$LNC[group_extra1]), # FDis species based approach extrapolation
                              rmse(FDis_obs$SLA[group_extra1], FDis_species_based_extra$SLA[group_extra1]),
                              rmse(FDis_obs$PLH[group_extra1], FDis_species_based_extra$PLH[group_extra1]),
                              rmse(FDis_obs$LNC[group_extra2], FDis_species_based_extra$LNC[group_extra2]),
                              rmse(FDis_obs$SLA[group_extra2], FDis_species_based_extra$SLA[group_extra2]),
                              rmse(FDis_obs$PLH[group_extra2], FDis_species_based_extra$PLH[group_extra2]),
                              rmse(FDis_obs$LNC[group_extra3], FDis_species_based_extra$LNC[group_extra3]),
                              rmse(FDis_obs$SLA[group_extra3], FDis_species_based_extra$SLA[group_extra3]),
                              rmse(FDis_obs$PLH[group_extra3], FDis_species_based_extra$PLH[group_extra3]),
                              rmse(FDis_obs$LNC[group_extra4], FDis_species_based_extra$LNC[group_extra4]),
                              rmse(FDis_obs$SLA[group_extra4], FDis_species_based_extra$SLA[group_extra4]),
                              rmse(FDis_obs$PLH[group_extra4], FDis_species_based_extra$PLH[group_extra4]),
                              
                              rmse(CM_obs$LNC[group_inter1], CM_trait_based_inter$LNC[group_inter1]), # CM trait based approach interpolation
                              rmse(CM_obs$SLA[group_inter1], CM_trait_based_inter$SLA[group_inter1]),
                              rmse(CM_obs$PLH[group_inter1], CM_trait_based_inter$PLH[group_inter1]),
                              rmse(CM_obs$LNC[group_inter2], CM_trait_based_inter$LNC[group_inter2]),
                              rmse(CM_obs$SLA[group_inter2], CM_trait_based_inter$SLA[group_inter2]),
                              rmse(CM_obs$PLH[group_inter2], CM_trait_based_inter$PLH[group_inter2]),
                              rmse(CM_obs$LNC[group_inter3], CM_trait_based_inter$LNC[group_inter3]),
                              rmse(CM_obs$SLA[group_inter3], CM_trait_based_inter$SLA[group_inter3]),
                              rmse(CM_obs$PLH[group_inter3], CM_trait_based_inter$PLH[group_inter3]),
                              rmse(CM_obs$LNC[group_inter4], CM_trait_based_inter$LNC[group_inter4]),
                              rmse(CM_obs$SLA[group_inter4], CM_trait_based_inter$SLA[group_inter4]),
                              rmse(CM_obs$PLH[group_inter4], CM_trait_based_inter$PLH[group_inter4]),
                              rmse(CM_obs$LNC[group_extra1], CM_trait_based_extra$LNC[group_extra1]),# CM trait based approach interpolation
                              rmse(CM_obs$SLA[group_extra1], CM_trait_based_extra$SLA[group_extra1]),
                              rmse(CM_obs$PLH[group_extra1], CM_trait_based_extra$PLH[group_extra1]),
                              rmse(CM_obs$LNC[group_extra2], CM_trait_based_extra$LNC[group_extra2]),
                              rmse(CM_obs$SLA[group_extra2], CM_trait_based_extra$SLA[group_extra2]),
                              rmse(CM_obs$PLH[group_extra2], CM_trait_based_extra$PLH[group_extra2]),
                              rmse(CM_obs$LNC[group_extra3], CM_trait_based_extra$LNC[group_extra3]),
                              rmse(CM_obs$SLA[group_extra3], CM_trait_based_extra$SLA[group_extra3]),
                              rmse(CM_obs$PLH[group_extra3], CM_trait_based_extra$PLH[group_extra3]),
                              rmse(CM_obs$LNC[group_extra4], CM_trait_based_extra$LNC[group_extra4]),
                              rmse(CM_obs$SLA[group_extra4], CM_trait_based_extra$SLA[group_extra4]),
                              rmse(CM_obs$PLH[group_extra4], CM_trait_based_extra$PLH[group_extra4]),
                              rmse(uFDis_obs$LNC[group_inter1], uFDis_trait_based_inter$LNC[group_inter1]),# uFDis trait based approach interpolation
                              rmse(uFDis_obs$SLA[group_inter1], uFDis_trait_based_inter$SLA[group_inter1]),
                              rmse(uFDis_obs$PLH[group_inter1], uFDis_trait_based_inter$PLH[group_inter1]),
                              rmse(uFDis_obs$LNC[group_inter2], uFDis_trait_based_inter$LNC[group_inter2]),
                              rmse(uFDis_obs$SLA[group_inter2], uFDis_trait_based_inter$SLA[group_inter2]),
                              rmse(uFDis_obs$PLH[group_inter2], uFDis_trait_based_inter$PLH[group_inter2]),
                              rmse(uFDis_obs$LNC[group_inter3], uFDis_trait_based_inter$LNC[group_inter3]),
                              rmse(uFDis_obs$SLA[group_inter3], uFDis_trait_based_inter$SLA[group_inter3]),
                              rmse(uFDis_obs$PLH[group_inter3], uFDis_trait_based_inter$PLH[group_inter3]),
                              rmse(uFDis_obs$LNC[group_inter4], uFDis_trait_based_inter$LNC[group_inter4]),
                              rmse(uFDis_obs$SLA[group_inter4], uFDis_trait_based_inter$SLA[group_inter4]),
                              rmse(uFDis_obs$PLH[group_inter4], uFDis_trait_based_inter$PLH[group_inter4]),
                              rmse(uFDis_obs$LNC[group_extra1], uFDis_trait_based_extra$LNC[group_extra1]), # uFDis trait based approach extrapolation
                              rmse(uFDis_obs$SLA[group_extra1], uFDis_trait_based_extra$SLA[group_extra1]),
                              rmse(uFDis_obs$PLH[group_extra1], uFDis_trait_based_extra$PLH[group_extra1]),
                              rmse(uFDis_obs$LNC[group_extra2], uFDis_trait_based_extra$LNC[group_extra2]),
                              rmse(uFDis_obs$SLA[group_extra2], uFDis_trait_based_extra$SLA[group_extra2]),
                              rmse(uFDis_obs$PLH[group_extra2], uFDis_trait_based_extra$PLH[group_extra2]),
                              rmse(uFDis_obs$LNC[group_extra3], uFDis_trait_based_extra$LNC[group_extra3]),
                              rmse(uFDis_obs$SLA[group_extra3], uFDis_trait_based_extra$SLA[group_extra3]),
                              rmse(uFDis_obs$PLH[group_extra3], uFDis_trait_based_extra$PLH[group_extra3]),
                              rmse(uFDis_obs$LNC[group_extra4], uFDis_trait_based_extra$LNC[group_extra4]),
                              rmse(uFDis_obs$SLA[group_extra4], uFDis_trait_based_extra$SLA[group_extra4]),
                              rmse(uFDis_obs$PLH[group_extra4], uFDis_trait_based_extra$PLH[group_extra4]),
                              rmse(CM_obs$LNC[group_inter1], CM_species_based_inter$LNC[group_inter1]), # CM species based approach interpolation
                              rmse(CM_obs$SLA[group_inter1], CM_species_based_inter$SLA[group_inter1]),
                              rmse(CM_obs$PLH[group_inter1], CM_species_based_inter$PLH[group_inter1]),
                              rmse(CM_obs$LNC[group_inter2], CM_species_based_inter$LNC[group_inter2]),
                              rmse(CM_obs$SLA[group_inter2], CM_species_based_inter$SLA[group_inter2]),
                              rmse(CM_obs$PLH[group_inter2], CM_species_based_inter$PLH[group_inter2]),
                              rmse(CM_obs$LNC[group_inter3], CM_species_based_inter$LNC[group_inter3]),
                              rmse(CM_obs$SLA[group_inter3], CM_species_based_inter$SLA[group_inter3]),
                              rmse(CM_obs$PLH[group_inter3], CM_species_based_inter$PLH[group_inter3]),
                              rmse(CM_obs$LNC[group_inter4], CM_species_based_inter$LNC[group_inter4]),
                              rmse(CM_obs$SLA[group_inter4], CM_species_based_inter$SLA[group_inter4]),
                              rmse(CM_obs$PLH[group_inter4], CM_species_based_inter$PLH[group_inter4]),
                              rmse(CM_obs$LNC[group_extra1], CM_species_based_extra$LNC[group_extra1]), # CM species based approach extrapolation
                              rmse(CM_obs$SLA[group_extra1], CM_species_based_extra$SLA[group_extra1]),
                              rmse(CM_obs$PLH[group_extra1], CM_species_based_extra$PLH[group_extra1]),
                              rmse(CM_obs$LNC[group_extra2], CM_species_based_extra$LNC[group_extra2]),
                              rmse(CM_obs$SLA[group_extra2], CM_species_based_extra$SLA[group_extra2]),
                              rmse(CM_obs$PLH[group_extra2], CM_species_based_extra$PLH[group_extra2]),
                              rmse(CM_obs$LNC[group_extra3], CM_species_based_extra$LNC[group_extra3]),
                              rmse(CM_obs$SLA[group_extra3], CM_species_based_extra$SLA[group_extra3]),
                              rmse(CM_obs$PLH[group_extra3], CM_species_based_extra$PLH[group_extra3]),
                              rmse(CM_obs$LNC[group_extra4], CM_species_based_extra$LNC[group_extra4]),
                              rmse(CM_obs$SLA[group_extra4], CM_species_based_extra$SLA[group_extra4]),
                              rmse(CM_obs$PLH[group_extra4], CM_species_based_extra$PLH[group_extra4]),
                              rmse(uFDis_obs$LNC[group_inter1], uFDis_species_based_inter$LNC[group_inter1]), # uFDis species based approach interpolation
                              rmse(uFDis_obs$SLA[group_inter1], uFDis_species_based_inter$SLA[group_inter1]),
                              rmse(uFDis_obs$PLH[group_inter1], uFDis_species_based_inter$PLH[group_inter1]),
                              rmse(uFDis_obs$LNC[group_inter2], uFDis_species_based_inter$LNC[group_inter2]),
                              rmse(uFDis_obs$SLA[group_inter2], uFDis_species_based_inter$SLA[group_inter2]),
                              rmse(uFDis_obs$PLH[group_inter2], uFDis_species_based_inter$PLH[group_inter2]),
                              rmse(uFDis_obs$LNC[group_inter3], uFDis_species_based_inter$LNC[group_inter3]),
                              rmse(uFDis_obs$SLA[group_inter3], uFDis_species_based_inter$SLA[group_inter3]),
                              rmse(uFDis_obs$PLH[group_inter3], uFDis_species_based_inter$PLH[group_inter3]),
                              rmse(uFDis_obs$LNC[group_inter4], uFDis_species_based_inter$LNC[group_inter4]),
                              rmse(uFDis_obs$SLA[group_inter4], uFDis_species_based_inter$SLA[group_inter4]),
                              rmse(uFDis_obs$PLH[group_inter4], uFDis_species_based_inter$PLH[group_inter4]),
                              rmse(uFDis_obs$LNC[group_extra1], uFDis_species_based_extra$LNC[group_extra1]), # uFDis species based approach extrapolation
                              rmse(uFDis_obs$SLA[group_extra1], uFDis_species_based_extra$SLA[group_extra1]),
                              rmse(uFDis_obs$PLH[group_extra1], uFDis_species_based_extra$PLH[group_extra1]),
                              rmse(uFDis_obs$LNC[group_extra2], uFDis_species_based_extra$LNC[group_extra2]),
                              rmse(uFDis_obs$SLA[group_extra2], uFDis_species_based_extra$SLA[group_extra2]),
                              rmse(uFDis_obs$PLH[group_extra2], uFDis_species_based_extra$PLH[group_extra2]),
                              rmse(uFDis_obs$LNC[group_extra3], uFDis_species_based_extra$LNC[group_extra3]),
                              rmse(uFDis_obs$SLA[group_extra3], uFDis_species_based_extra$SLA[group_extra3]),
                              rmse(uFDis_obs$PLH[group_extra3], uFDis_species_based_extra$PLH[group_extra3]),
                              rmse(uFDis_obs$LNC[group_extra4], uFDis_species_based_extra$LNC[group_extra4]),
                              rmse(uFDis_obs$SLA[group_extra4], uFDis_species_based_extra$SLA[group_extra4]),
                              rmse(uFDis_obs$PLH[group_extra4], uFDis_species_based_extra$PLH[group_extra4])),
                       "R2"=c(cor(CWM_obs$LNC[group_inter1], CWM_trait_based_inter$LNC[group_inter1])^2, # CWM trait based approach interpolation
                                cor(CWM_obs$SLA[group_inter1], CWM_trait_based_inter$SLA[group_inter1])^2,
                                cor(CWM_obs$PLH[group_inter1], CWM_trait_based_inter$PLH[group_inter1])^2,
                                cor(CWM_obs$LNC[group_inter2], CWM_trait_based_inter$LNC[group_inter2])^2,
                                cor(CWM_obs$SLA[group_inter2], CWM_trait_based_inter$SLA[group_inter2])^2,
                                cor(CWM_obs$PLH[group_inter2], CWM_trait_based_inter$PLH[group_inter2])^2,
                                cor(CWM_obs$LNC[group_inter3], CWM_trait_based_inter$LNC[group_inter3])^2,
                                cor(CWM_obs$SLA[group_inter3], CWM_trait_based_inter$SLA[group_inter3])^2,
                                cor(CWM_obs$PLH[group_inter3], CWM_trait_based_inter$PLH[group_inter3])^2,
                                cor(CWM_obs$LNC[group_inter4], CWM_trait_based_inter$LNC[group_inter4])^2,
                                cor(CWM_obs$SLA[group_inter4], CWM_trait_based_inter$SLA[group_inter4])^2,
                                cor(CWM_obs$PLH[group_inter4], CWM_trait_based_inter$PLH[group_inter4])^2,
                                cor(CWM_obs$LNC[group_extra1], CWM_trait_based_extra$LNC[group_extra1])^2,# CWM trait based approach interpolation
                                cor(CWM_obs$SLA[group_extra1], CWM_trait_based_extra$SLA[group_extra1])^2,
                                cor(CWM_obs$PLH[group_extra1], CWM_trait_based_extra$PLH[group_extra1])^2,
                                cor(CWM_obs$LNC[group_extra2], CWM_trait_based_extra$LNC[group_extra2])^2,
                                cor(CWM_obs$SLA[group_extra2], CWM_trait_based_extra$SLA[group_extra2])^2,
                                cor(CWM_obs$PLH[group_extra2], CWM_trait_based_extra$PLH[group_extra2])^2,
                                cor(CWM_obs$LNC[group_extra3], CWM_trait_based_extra$LNC[group_extra3])^2,
                                cor(CWM_obs$SLA[group_extra3], CWM_trait_based_extra$SLA[group_extra3])^2,
                                cor(CWM_obs$PLH[group_extra3], CWM_trait_based_extra$PLH[group_extra3])^2,
                                cor(CWM_obs$LNC[group_extra4], CWM_trait_based_extra$LNC[group_extra4])^2,
                                cor(CWM_obs$SLA[group_extra4], CWM_trait_based_extra$SLA[group_extra4])^2,
                                cor(CWM_obs$PLH[group_extra4], CWM_trait_based_extra$PLH[group_extra4])^2,
                                cor(FDis_obs$LNC[group_inter1], FDis_trait_based_inter$LNC[group_inter1])^2,# FDis trait based approach interpolation
                                cor(FDis_obs$SLA[group_inter1], FDis_trait_based_inter$SLA[group_inter1])^2,
                                cor(FDis_obs$PLH[group_inter1], FDis_trait_based_inter$PLH[group_inter1])^2,
                                cor(FDis_obs$LNC[group_inter2], FDis_trait_based_inter$LNC[group_inter2])^2,
                                cor(FDis_obs$SLA[group_inter2], FDis_trait_based_inter$SLA[group_inter2])^2,
                                cor(FDis_obs$PLH[group_inter2], FDis_trait_based_inter$PLH[group_inter2])^2,
                                cor(FDis_obs$LNC[group_inter3], FDis_trait_based_inter$LNC[group_inter3])^2,
                                cor(FDis_obs$SLA[group_inter3], FDis_trait_based_inter$SLA[group_inter3])^2,
                                cor(FDis_obs$PLH[group_inter3], FDis_trait_based_inter$PLH[group_inter3])^2,
                                cor(FDis_obs$LNC[group_inter4], FDis_trait_based_inter$LNC[group_inter4])^2,
                                cor(FDis_obs$SLA[group_inter4], FDis_trait_based_inter$SLA[group_inter4])^2,
                                cor(FDis_obs$PLH[group_inter4], FDis_trait_based_inter$PLH[group_inter4])^2,
                                cor(FDis_obs$LNC[group_extra1], FDis_trait_based_extra$LNC[group_extra1])^2, # FDis trait based approach extrapolation
                                cor(FDis_obs$SLA[group_extra1], FDis_trait_based_extra$SLA[group_extra1])^2,
                                cor(FDis_obs$PLH[group_extra1], FDis_trait_based_extra$PLH[group_extra1])^2,
                                cor(FDis_obs$LNC[group_extra2], FDis_trait_based_extra$LNC[group_extra2])^2,
                                cor(FDis_obs$SLA[group_extra2], FDis_trait_based_extra$SLA[group_extra2])^2,
                                cor(FDis_obs$PLH[group_extra2], FDis_trait_based_extra$PLH[group_extra2])^2,
                                cor(FDis_obs$LNC[group_extra3], FDis_trait_based_extra$LNC[group_extra3])^2,
                                cor(FDis_obs$SLA[group_extra3], FDis_trait_based_extra$SLA[group_extra3])^2,
                                cor(FDis_obs$PLH[group_extra3], FDis_trait_based_extra$PLH[group_extra3])^2,
                                cor(FDis_obs$LNC[group_extra4], FDis_trait_based_extra$LNC[group_extra4])^2,
                                cor(FDis_obs$SLA[group_extra4], FDis_trait_based_extra$SLA[group_extra4])^2,
                                cor(FDis_obs$PLH[group_extra4], FDis_trait_based_extra$PLH[group_extra4])^2,
                                cor(CWM_obs$LNC[group_inter1], CWM_species_based_inter$LNC[group_inter1])^2, # CWM species based approach interpolation
                                cor(CWM_obs$SLA[group_inter1], CWM_species_based_inter$SLA[group_inter1])^2,
                                cor(CWM_obs$PLH[group_inter1], CWM_species_based_inter$PLH[group_inter1])^2,
                                cor(CWM_obs$LNC[group_inter2], CWM_species_based_inter$LNC[group_inter2])^2,
                                cor(CWM_obs$SLA[group_inter2], CWM_species_based_inter$SLA[group_inter2])^2,
                                cor(CWM_obs$PLH[group_inter2], CWM_species_based_inter$PLH[group_inter2])^2,
                                cor(CWM_obs$LNC[group_inter3], CWM_species_based_inter$LNC[group_inter3])^2,
                                cor(CWM_obs$SLA[group_inter3], CWM_species_based_inter$SLA[group_inter3])^2,
                                cor(CWM_obs$PLH[group_inter3], CWM_species_based_inter$PLH[group_inter3])^2,
                                cor(CWM_obs$LNC[group_inter4], CWM_species_based_inter$LNC[group_inter4])^2,
                                cor(CWM_obs$SLA[group_inter4], CWM_species_based_inter$SLA[group_inter4])^2,
                                cor(CWM_obs$PLH[group_inter4], CWM_species_based_inter$PLH[group_inter4])^2,
                                cor(CWM_obs$LNC[group_extra1], CWM_species_based_extra$LNC[group_extra1])^2, # CWM species based approach extrapolation
                                cor(CWM_obs$SLA[group_extra1], CWM_species_based_extra$SLA[group_extra1])^2,
                                cor(CWM_obs$PLH[group_extra1], CWM_species_based_extra$PLH[group_extra1])^2,
                                cor(CWM_obs$LNC[group_extra2], CWM_species_based_extra$LNC[group_extra2])^2,
                                cor(CWM_obs$SLA[group_extra2], CWM_species_based_extra$SLA[group_extra2])^2,
                                cor(CWM_obs$PLH[group_extra2], CWM_species_based_extra$PLH[group_extra2])^2,
                                cor(CWM_obs$LNC[group_extra3], CWM_species_based_extra$LNC[group_extra3])^2,
                                cor(CWM_obs$SLA[group_extra3], CWM_species_based_extra$SLA[group_extra3])^2,
                                cor(CWM_obs$PLH[group_extra3], CWM_species_based_extra$PLH[group_extra3])^2,
                                cor(CWM_obs$LNC[group_extra4], CWM_species_based_extra$LNC[group_extra4])^2,
                                cor(CWM_obs$SLA[group_extra4], CWM_species_based_extra$SLA[group_extra4])^2,
                                cor(CWM_obs$PLH[group_extra4], CWM_species_based_extra$PLH[group_extra4])^2,
                                cor(FDis_obs$LNC[group_inter1], FDis_species_based_inter$LNC[group_inter1])^2, # FDis species based approach interpolation
                                cor(FDis_obs$SLA[group_inter1], FDis_species_based_inter$SLA[group_inter1])^2,
                                cor(FDis_obs$PLH[group_inter1], FDis_species_based_inter$PLH[group_inter1])^2,
                                cor(FDis_obs$LNC[group_inter2], FDis_species_based_inter$LNC[group_inter2])^2,
                                cor(FDis_obs$SLA[group_inter2], FDis_species_based_inter$SLA[group_inter2])^2,
                                cor(FDis_obs$PLH[group_inter2], FDis_species_based_inter$PLH[group_inter2])^2,
                                cor(FDis_obs$LNC[group_inter3], FDis_species_based_inter$LNC[group_inter3])^2,
                                cor(FDis_obs$SLA[group_inter3], FDis_species_based_inter$SLA[group_inter3])^2,
                                cor(FDis_obs$PLH[group_inter3], FDis_species_based_inter$PLH[group_inter3])^2,
                                cor(FDis_obs$LNC[group_inter4], FDis_species_based_inter$LNC[group_inter4])^2,
                                cor(FDis_obs$SLA[group_inter4], FDis_species_based_inter$SLA[group_inter4])^2,
                                cor(FDis_obs$PLH[group_inter4], FDis_species_based_inter$PLH[group_inter4])^2,
                                cor(FDis_obs$LNC[group_extra1], FDis_species_based_extra$LNC[group_extra1])^2, # FDis species based approach extrapolation
                                cor(FDis_obs$SLA[group_extra1], FDis_species_based_extra$SLA[group_extra1])^2,
                                cor(FDis_obs$PLH[group_extra1], FDis_species_based_extra$PLH[group_extra1])^2,
                                cor(FDis_obs$LNC[group_extra2], FDis_species_based_extra$LNC[group_extra2])^2,
                                cor(FDis_obs$SLA[group_extra2], FDis_species_based_extra$SLA[group_extra2])^2,
                                cor(FDis_obs$PLH[group_extra2], FDis_species_based_extra$PLH[group_extra2])^2,
                                cor(FDis_obs$LNC[group_extra3], FDis_species_based_extra$LNC[group_extra3])^2,
                                cor(FDis_obs$SLA[group_extra3], FDis_species_based_extra$SLA[group_extra3])^2,
                                cor(FDis_obs$PLH[group_extra3], FDis_species_based_extra$PLH[group_extra3])^2,
                                cor(FDis_obs$LNC[group_extra4], FDis_species_based_extra$LNC[group_extra4])^2,
                                cor(FDis_obs$SLA[group_extra4], FDis_species_based_extra$SLA[group_extra4])^2,
                                cor(FDis_obs$PLH[group_extra4], FDis_species_based_extra$PLH[group_extra4])^2,
                                
                                cor(CM_obs$LNC[group_inter1], CM_trait_based_inter$LNC[group_inter1])^2, # CM trait based approach interpolation
                                cor(CM_obs$SLA[group_inter1], CM_trait_based_inter$SLA[group_inter1])^2,
                                cor(CM_obs$PLH[group_inter1], CM_trait_based_inter$PLH[group_inter1])^2,
                                cor(CM_obs$LNC[group_inter2], CM_trait_based_inter$LNC[group_inter2])^2,
                                cor(CM_obs$SLA[group_inter2], CM_trait_based_inter$SLA[group_inter2])^2,
                                cor(CM_obs$PLH[group_inter2], CM_trait_based_inter$PLH[group_inter2])^2,
                                cor(CM_obs$LNC[group_inter3], CM_trait_based_inter$LNC[group_inter3])^2,
                                cor(CM_obs$SLA[group_inter3], CM_trait_based_inter$SLA[group_inter3])^2,
                                cor(CM_obs$PLH[group_inter3], CM_trait_based_inter$PLH[group_inter3])^2,
                                cor(CM_obs$LNC[group_inter4], CM_trait_based_inter$LNC[group_inter4])^2,
                                cor(CM_obs$SLA[group_inter4], CM_trait_based_inter$SLA[group_inter4])^2,
                                cor(CM_obs$PLH[group_inter4], CM_trait_based_inter$PLH[group_inter4])^2,
                                cor(CM_obs$LNC[group_extra1], CM_trait_based_extra$LNC[group_extra1])^2,# CM trait based approach interpolation
                                cor(CM_obs$SLA[group_extra1], CM_trait_based_extra$SLA[group_extra1])^2,
                                cor(CM_obs$PLH[group_extra1], CM_trait_based_extra$PLH[group_extra1])^2,
                                cor(CM_obs$LNC[group_extra2], CM_trait_based_extra$LNC[group_extra2])^2,
                                cor(CM_obs$SLA[group_extra2], CM_trait_based_extra$SLA[group_extra2])^2,
                                cor(CM_obs$PLH[group_extra2], CM_trait_based_extra$PLH[group_extra2])^2,
                                cor(CM_obs$LNC[group_extra3], CM_trait_based_extra$LNC[group_extra3])^2,
                                cor(CM_obs$SLA[group_extra3], CM_trait_based_extra$SLA[group_extra3])^2,
                                cor(CM_obs$PLH[group_extra3], CM_trait_based_extra$PLH[group_extra3])^2,
                                cor(CM_obs$LNC[group_extra4], CM_trait_based_extra$LNC[group_extra4])^2,
                                cor(CM_obs$SLA[group_extra4], CM_trait_based_extra$SLA[group_extra4])^2,
                                cor(CM_obs$PLH[group_extra4], CM_trait_based_extra$PLH[group_extra4])^2,
                                cor(uFDis_obs$LNC[group_inter1], uFDis_trait_based_inter$LNC[group_inter1])^2,# uFDis trait based approach interpolation
                                cor(uFDis_obs$SLA[group_inter1], uFDis_trait_based_inter$SLA[group_inter1])^2,
                                cor(uFDis_obs$PLH[group_inter1], uFDis_trait_based_inter$PLH[group_inter1])^2,
                                cor(uFDis_obs$LNC[group_inter2], uFDis_trait_based_inter$LNC[group_inter2])^2,
                                cor(uFDis_obs$SLA[group_inter2], uFDis_trait_based_inter$SLA[group_inter2])^2,
                                cor(uFDis_obs$PLH[group_inter2], uFDis_trait_based_inter$PLH[group_inter2])^2,
                                cor(uFDis_obs$LNC[group_inter3], uFDis_trait_based_inter$LNC[group_inter3])^2,
                                cor(uFDis_obs$SLA[group_inter3], uFDis_trait_based_inter$SLA[group_inter3])^2,
                                cor(uFDis_obs$PLH[group_inter3], uFDis_trait_based_inter$PLH[group_inter3])^2,
                                cor(uFDis_obs$LNC[group_inter4], uFDis_trait_based_inter$LNC[group_inter4])^2,
                                cor(uFDis_obs$SLA[group_inter4], uFDis_trait_based_inter$SLA[group_inter4])^2,
                                cor(uFDis_obs$PLH[group_inter4], uFDis_trait_based_inter$PLH[group_inter4])^2,
                                cor(uFDis_obs$LNC[group_extra1], uFDis_trait_based_extra$LNC[group_extra1])^2, # uFDis trait based approach extrapolation
                                cor(uFDis_obs$SLA[group_extra1], uFDis_trait_based_extra$SLA[group_extra1])^2,
                                cor(uFDis_obs$PLH[group_extra1], uFDis_trait_based_extra$PLH[group_extra1])^2,
                                cor(uFDis_obs$LNC[group_extra2], uFDis_trait_based_extra$LNC[group_extra2])^2,
                                cor(uFDis_obs$SLA[group_extra2], uFDis_trait_based_extra$SLA[group_extra2])^2,
                                cor(uFDis_obs$PLH[group_extra2], uFDis_trait_based_extra$PLH[group_extra2])^2,
                                cor(uFDis_obs$LNC[group_extra3], uFDis_trait_based_extra$LNC[group_extra3])^2,
                                cor(uFDis_obs$SLA[group_extra3], uFDis_trait_based_extra$SLA[group_extra3])^2,
                                cor(uFDis_obs$PLH[group_extra3], uFDis_trait_based_extra$PLH[group_extra3])^2,
                                cor(uFDis_obs$LNC[group_extra4], uFDis_trait_based_extra$LNC[group_extra4])^2,
                                cor(uFDis_obs$SLA[group_extra4], uFDis_trait_based_extra$SLA[group_extra4])^2,
                                cor(uFDis_obs$PLH[group_extra4], uFDis_trait_based_extra$PLH[group_extra4])^2,
                                cor(CM_obs$LNC[group_inter1], CM_species_based_inter$LNC[group_inter1])^2, # CM species based approach interpolation
                                cor(CM_obs$SLA[group_inter1], CM_species_based_inter$SLA[group_inter1])^2,
                                cor(CM_obs$PLH[group_inter1], CM_species_based_inter$PLH[group_inter1])^2,
                                cor(CM_obs$LNC[group_inter2], CM_species_based_inter$LNC[group_inter2])^2,
                                cor(CM_obs$SLA[group_inter2], CM_species_based_inter$SLA[group_inter2])^2,
                                cor(CM_obs$PLH[group_inter2], CM_species_based_inter$PLH[group_inter2])^2,
                                cor(CM_obs$LNC[group_inter3], CM_species_based_inter$LNC[group_inter3])^2,
                                cor(CM_obs$SLA[group_inter3], CM_species_based_inter$SLA[group_inter3])^2,
                                cor(CM_obs$PLH[group_inter3], CM_species_based_inter$PLH[group_inter3])^2,
                                cor(CM_obs$LNC[group_inter4], CM_species_based_inter$LNC[group_inter4])^2,
                                cor(CM_obs$SLA[group_inter4], CM_species_based_inter$SLA[group_inter4])^2,
                                cor(CM_obs$PLH[group_inter4], CM_species_based_inter$PLH[group_inter4])^2,
                                cor(CM_obs$LNC[group_extra1], CM_species_based_extra$LNC[group_extra1])^2, # CM species based approach extrapolation
                                cor(CM_obs$SLA[group_extra1], CM_species_based_extra$SLA[group_extra1])^2,
                                cor(CM_obs$PLH[group_extra1], CM_species_based_extra$PLH[group_extra1])^2,
                                cor(CM_obs$LNC[group_extra2], CM_species_based_extra$LNC[group_extra2])^2,
                                cor(CM_obs$SLA[group_extra2], CM_species_based_extra$SLA[group_extra2])^2,
                                cor(CM_obs$PLH[group_extra2], CM_species_based_extra$PLH[group_extra2])^2,
                                cor(CM_obs$LNC[group_extra3], CM_species_based_extra$LNC[group_extra3])^2,
                                cor(CM_obs$SLA[group_extra3], CM_species_based_extra$SLA[group_extra3])^2,
                                cor(CM_obs$PLH[group_extra3], CM_species_based_extra$PLH[group_extra3])^2,
                                cor(CM_obs$LNC[group_extra4], CM_species_based_extra$LNC[group_extra4])^2,
                                cor(CM_obs$SLA[group_extra4], CM_species_based_extra$SLA[group_extra4])^2,
                                cor(CM_obs$PLH[group_extra4], CM_species_based_extra$PLH[group_extra4])^2,
                                cor(uFDis_obs$LNC[group_inter1], uFDis_species_based_inter$LNC[group_inter1])^2, # uFDis species based approach interpolation
                                cor(uFDis_obs$SLA[group_inter1], uFDis_species_based_inter$SLA[group_inter1])^2,
                                cor(uFDis_obs$PLH[group_inter1], uFDis_species_based_inter$PLH[group_inter1])^2,
                                cor(uFDis_obs$LNC[group_inter2], uFDis_species_based_inter$LNC[group_inter2])^2,
                                cor(uFDis_obs$SLA[group_inter2], uFDis_species_based_inter$SLA[group_inter2])^2,
                                cor(uFDis_obs$PLH[group_inter2], uFDis_species_based_inter$PLH[group_inter2])^2,
                                cor(uFDis_obs$LNC[group_inter3], uFDis_species_based_inter$LNC[group_inter3])^2,
                                cor(uFDis_obs$SLA[group_inter3], uFDis_species_based_inter$SLA[group_inter3])^2,
                                cor(uFDis_obs$PLH[group_inter3], uFDis_species_based_inter$PLH[group_inter3])^2,
                                cor(uFDis_obs$LNC[group_inter4], uFDis_species_based_inter$LNC[group_inter4])^2,
                                cor(uFDis_obs$SLA[group_inter4], uFDis_species_based_inter$SLA[group_inter4])^2,
                                cor(uFDis_obs$PLH[group_inter4], uFDis_species_based_inter$PLH[group_inter4])^2,
                                cor(uFDis_obs$LNC[group_extra1], uFDis_species_based_extra$LNC[group_extra1])^2, # uFDis species based approach extrapolation
                                cor(uFDis_obs$SLA[group_extra1], uFDis_species_based_extra$SLA[group_extra1])^2,
                                cor(uFDis_obs$PLH[group_extra1], uFDis_species_based_extra$PLH[group_extra1])^2,
                                cor(uFDis_obs$LNC[group_extra2], uFDis_species_based_extra$LNC[group_extra2])^2,
                                cor(uFDis_obs$SLA[group_extra2], uFDis_species_based_extra$SLA[group_extra2])^2,
                                cor(uFDis_obs$PLH[group_extra2], uFDis_species_based_extra$PLH[group_extra2])^2,
                                cor(uFDis_obs$LNC[group_extra3], uFDis_species_based_extra$LNC[group_extra3])^2,
                                cor(uFDis_obs$SLA[group_extra3], uFDis_species_based_extra$SLA[group_extra3])^2,
                                cor(uFDis_obs$PLH[group_extra3], uFDis_species_based_extra$PLH[group_extra3])^2,
                                cor(uFDis_obs$LNC[group_extra4], uFDis_species_based_extra$LNC[group_extra4])^2,
                                cor(uFDis_obs$SLA[group_extra4], uFDis_species_based_extra$SLA[group_extra4])^2,
                                cor(uFDis_obs$PLH[group_extra4], uFDis_species_based_extra$PLH[group_extra4])^2),
                       
                       "Traits"=rep(c(rep(c("LNC", "SLA", "PLH"), 4*8)),2),
                       "Indices"=c(c(rep(c(rep("CWM", 24), rep("FDis", 24)),2)),c(rep(c(rep("CM", 24), rep("uFDis", 24)),2))),
                       "Approaches"=rep(c(rep("Assemble-First Approach", 12*4), rep("Predict-First Approach", 12*4)),2),
                       "Crossvalidation"=rep(c(rep(c(rep("Interpolation", 12), rep("Extrapolation", 12)),4)),2),
                       "Type"=c(rep("Abundance data", 4*8*3), rep("Occurrence data", 4*8*3)),
                       "repet"=rep(c(rep(1,3), rep(2,3), rep(3,3), rep(4,3)), 16)
)

figure_3$Traits <- fct_relevel(figure_3$Traits, c("PLH", "SLA", "LNC"))
figure_3$Crossvalidation <- fct_relevel(figure_3$Crossvalidation, c("Interpolation", "Extrapolation"))
figure_3$Type <- fct_relevel(figure_3$Type, c("Abundance data","Occurrence data"))
figure_3$Approaches <- fct_relevel(figure_3$Approaches, c("Predict-First Approach","Assemble-First Approach"))
figure_3$Indices <- fct_relevel(figure_3$Indices, c("CWM","FDis", "CM", 'uFDis'))

### 10 null tests with permutation of the trait-species matrix rows
##### 
shuffle_traits <- traits[sample(1:nrow(traits)),]
shuffle_traits$Taxa <- traits$Taxa
shuffle_CWM_obs <- CWM_shuffle(plant_recovery)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_CWM_obs <- left_join(Site_order, shuffle_CWM_obs, by="Site")
rownames(shuffle_CWM_obs)=shuffle_CWM_obs$Site
shuffle_FDis_obs <- FDis_shuffle(plant_recovery)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_FDis_obs <- left_join(Site_order, shuffle_FDis_obs, by="Site")
rownames(shuffle_FDis_obs)=shuffle_FDis_obs$Site

shuffle_CM_obs <- CWM_shuffle(plant_pa)
Site_order <- data.frame("Site"=rownames(plant_pa))
shuffle_CM_obs <- left_join(Site_order, shuffle_CM_obs, by="Site")
rownames(shuffle_CM_obs)=shuffle_CM_obs$Site
shuffle_uFDis_obs <- FDis_shuffle(plant_pa)
Site_order <- data.frame("Site"=rownames(plant_pa))
shuffle_uFDis_obs <- left_join(Site_order, shuffle_uFDis_obs, by="Site")
rownames(shuffle_uFDis_obs)=shuffle_uFDis_obs$Site

shuffle_CWM_species_based_inter <- CWM_shuffle(species_models_recovery_inter)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_CWM_species_based_inter <- left_join(Site_order, shuffle_CWM_species_based_inter, by="Site")
rownames(shuffle_CWM_species_based_inter)=shuffle_CWM_species_based_inter$Site
shuffle_CWM_species_based_inter <- shuffle_CWM_species_based_inter[,-1]

shuffle_CWM_species_based_extra <- CWM_shuffle(species_models_recovery_extra)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_CWM_species_based_extra <- left_join(Site_order, shuffle_CWM_species_based_extra, by="Site")
rownames(shuffle_CWM_species_based_extra)=shuffle_CWM_species_based_extra$Site
shuffle_CWM_species_based_extra <- shuffle_CWM_species_based_extra[,-1]

shuffle_CM_species_based_inter <- CWM_shuffle(species_models_pa_inter)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_CM_species_based_inter <- left_join(Site_order, shuffle_CM_species_based_inter, by="Site")
rownames(shuffle_CM_species_based_inter)=shuffle_CM_species_based_inter$Site
shuffle_CM_species_based_inter <- shuffle_CM_species_based_inter[,-1]

shuffle_CM_species_based_extra <- CWM_shuffle(species_models_pa_extra)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_CM_species_based_extra <- left_join(Site_order, shuffle_CM_species_based_extra, by="Site")
rownames(shuffle_CM_species_based_extra)=shuffle_CM_species_based_extra$Site
shuffle_CM_species_based_extra <- shuffle_CM_species_based_extra[,-1]

shuffle_FDis_species_based_inter <- FDis_shuffle(species_models_recovery_inter)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_FDis_species_based_inter <- left_join(Site_order, shuffle_FDis_species_based_inter, by="Site")
rownames(shuffle_FDis_species_based_inter)=shuffle_FDis_species_based_inter$Site
shuffle_FDis_species_based_inter <- shuffle_FDis_species_based_inter[,-1]

shuffle_FDis_species_based_extra <- FDis_shuffle(species_models_recovery_extra)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_FDis_species_based_extra <- left_join(Site_order, shuffle_FDis_species_based_extra, by="Site")
rownames(shuffle_FDis_species_based_extra)=shuffle_FDis_species_based_extra$Site
shuffle_FDis_species_based_extra <- shuffle_FDis_species_based_extra[,-1]

shuffle_uFDis_species_based_inter <- FDis_shuffle(species_models_pa_inter)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_uFDis_species_based_inter <- left_join(Site_order, shuffle_uFDis_species_based_inter, by="Site")
rownames(shuffle_uFDis_species_based_inter)=shuffle_uFDis_species_based_inter$Site
shuffle_uFDis_species_based_inter <- shuffle_uFDis_species_based_inter[,-1]

shuffle_uFDis_species_based_extra <- FDis_shuffle(species_models_pa_extra)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_uFDis_species_based_extra <- left_join(Site_order, shuffle_uFDis_species_based_extra, by="Site")
rownames(shuffle_uFDis_species_based_extra)=shuffle_uFDis_species_based_extra$Site
shuffle_uFDis_species_based_extra <- shuffle_uFDis_species_based_extra[,-1]

##################trait_based 
gam_data <- left_join(shuffle_CWM_obs,env_var, by="Site")
## - For interpolation group 1 - ##
shuffle_models_CWM_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
shuffle_models_CWM_LNC_inter1 <- predict(shuffle_models_CWM_LNC_grow1, gam_data[group_inter1,])$predicted
shuffle_models_CWM_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
shuffle_models_CWM_SLA_inter1 <- predict(shuffle_models_CWM_SLA_grow1, gam_data[group_inter1,])$predicted
shuffle_models_CWM_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
shuffle_models_CWM_PLH_inter1 <- predict(shuffle_models_CWM_PLH_grow1, gam_data[group_inter1,])$predicted

shuffle_CWM_trait_based_inter1<-data.frame("LNC"=shuffle_models_CWM_LNC_inter1, "SLA"=shuffle_models_CWM_SLA_inter1, "PLH"=shuffle_models_CWM_PLH_inter1)
rownames(shuffle_CWM_trait_based_inter1) <- shuffle_CWM_obs$Site[group_inter1]

## - For interpolation group 2 - ##
shuffle_models_CWM_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
shuffle_models_CWM_LNC_inter2 <- predict(shuffle_models_CWM_LNC_grow2, gam_data[group_inter2,])$predicted
shuffle_models_CWM_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
shuffle_models_CWM_SLA_inter2 <- predict(shuffle_models_CWM_SLA_grow2, gam_data[group_inter2,])$predicted
shuffle_models_CWM_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
shuffle_models_CWM_PLH_inter2 <- predict(shuffle_models_CWM_PLH_grow2, gam_data[group_inter2,])$predicted

shuffle_CWM_trait_based_inter2<-data.frame("LNC"=shuffle_models_CWM_LNC_inter2, "SLA"=shuffle_models_CWM_SLA_inter2, "PLH"=shuffle_models_CWM_PLH_inter2)
rownames(shuffle_CWM_trait_based_inter2) <- shuffle_CWM_obs$Site[group_inter2]

## - For interpolation group 3 - ##
shuffle_models_CWM_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
shuffle_models_CWM_LNC_inter3 <- predict(shuffle_models_CWM_LNC_grow3, gam_data[group_inter3,])$predicted
shuffle_models_CWM_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
shuffle_models_CWM_SLA_inter3 <- predict(shuffle_models_CWM_SLA_grow3, gam_data[group_inter3,])$predicted
shuffle_models_CWM_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
shuffle_models_CWM_PLH_inter3 <- predict(shuffle_models_CWM_PLH_grow3, gam_data[group_inter3,])$predicted

shuffle_CWM_trait_based_inter3<-data.frame("LNC"=shuffle_models_CWM_LNC_inter3, "SLA"=shuffle_models_CWM_SLA_inter3, "PLH"=shuffle_models_CWM_PLH_inter3)
rownames(shuffle_CWM_trait_based_inter3) <- shuffle_CWM_obs$Site[group_inter3]

## - For interpolation group 4 - ##
shuffle_models_CWM_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
shuffle_models_CWM_LNC_inter4 <- predict(shuffle_models_CWM_LNC_grow4, gam_data[group_inter4,])$predicted
shuffle_models_CWM_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
shuffle_models_CWM_SLA_inter4 <- predict(shuffle_models_CWM_SLA_grow4, gam_data[group_inter4,])$predicted
shuffle_models_CWM_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
shuffle_models_CWM_PLH_inter4 <- predict(shuffle_models_CWM_PLH_grow4, gam_data[group_inter4,])$predicted

shuffle_CWM_trait_based_inter4<-data.frame("LNC"=shuffle_models_CWM_LNC_inter4, "SLA"=shuffle_models_CWM_SLA_inter4, "PLH"=shuffle_models_CWM_PLH_inter4)
rownames(shuffle_CWM_trait_based_inter4) <- shuffle_CWM_obs$Site[group_inter4]

## assembly in a single prediction table for the 4 groups :
shuffle_CWM_trait_based_inter <- rbind(shuffle_CWM_trait_based_inter1, shuffle_CWM_trait_based_inter2, shuffle_CWM_trait_based_inter3, shuffle_CWM_trait_based_inter4)
shuffle_CWM_trait_based_inter <- cbind("Site"=rownames(shuffle_CWM_trait_based_inter), shuffle_CWM_trait_based_inter)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_CWM_trait_based_inter <- left_join(Site_order, shuffle_CWM_trait_based_inter, by="Site")
rownames(shuffle_CWM_trait_based_inter) <- shuffle_CWM_trait_based_inter$Site
shuffle_CWM_trait_based_inter <- shuffle_CWM_trait_based_inter[,-1]

#### -- This computes path predicts the community weighted mean in the context of extrapolation -- ####
gam_data <- left_join(shuffle_CWM_obs,env_var, by="Site")
## - For extrapolation group 1 - ##
shuffle_models_CWM_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
shuffle_models_CWM_LNC_extra1 <- predict(shuffle_models_CWM_LNC_grow1, gam_data[group_extra1,])$predicted
shuffle_models_CWM_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
shuffle_models_CWM_SLA_extra1 <- predict(shuffle_models_CWM_SLA_grow1, gam_data[group_extra1,])$predicted
shuffle_models_CWM_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
shuffle_models_CWM_PLH_extra1 <- predict(shuffle_models_CWM_PLH_grow1, gam_data[group_extra1,])$predicted

shuffle_CWM_trait_based_extra1<-data.frame("LNC"=shuffle_models_CWM_LNC_extra1, "SLA"=shuffle_models_CWM_SLA_extra1, "PLH"=shuffle_models_CWM_PLH_extra1)
rownames(shuffle_CWM_trait_based_extra1) <- shuffle_CWM_obs$Site[group_extra1]

## - For extrapolation group 2 - ##
shuffle_models_CWM_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
shuffle_models_CWM_LNC_extra2 <- predict(shuffle_models_CWM_LNC_grow2, gam_data[group_extra2,])$predicted
shuffle_models_CWM_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
shuffle_models_CWM_SLA_extra2 <- predict(shuffle_models_CWM_SLA_grow2, gam_data[group_extra2,])$predicted
shuffle_models_CWM_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
shuffle_models_CWM_PLH_extra2 <- predict(shuffle_models_CWM_PLH_grow2, gam_data[group_extra2,])$predicted

shuffle_CWM_trait_based_extra2<-data.frame("LNC"=shuffle_models_CWM_LNC_extra2, "SLA"=shuffle_models_CWM_SLA_extra2, "PLH"=shuffle_models_CWM_PLH_extra2)
rownames(shuffle_CWM_trait_based_extra2) <- shuffle_CWM_obs$Site[group_extra2]

## - For extrapolation group 3 - ##
shuffle_models_CWM_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
shuffle_models_CWM_LNC_extra3 <- predict(shuffle_models_CWM_LNC_grow3, gam_data[group_extra3,])$predicted
shuffle_models_CWM_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
shuffle_models_CWM_SLA_extra3 <- predict(shuffle_models_CWM_SLA_grow3, gam_data[group_extra3,])$predicted
shuffle_models_CWM_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
shuffle_models_CWM_PLH_extra3 <- predict(shuffle_models_CWM_PLH_grow3, gam_data[group_extra3,])$predicted

shuffle_CWM_trait_based_extra3<-data.frame("LNC"=shuffle_models_CWM_LNC_extra3, "SLA"=shuffle_models_CWM_SLA_extra3, "PLH"=shuffle_models_CWM_PLH_extra3)
rownames(shuffle_CWM_trait_based_extra3) <- shuffle_CWM_obs$Site[group_extra3]

## - For extrapolation group 4 - ##
shuffle_models_CWM_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
shuffle_models_CWM_LNC_extra4 <- predict(shuffle_models_CWM_LNC_grow4, gam_data[group_extra4,])$predicted
shuffle_models_CWM_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
shuffle_models_CWM_SLA_extra4 <- predict(shuffle_models_CWM_SLA_grow4, gam_data[group_extra4,])$predicted
shuffle_models_CWM_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
shuffle_models_CWM_PLH_extra4 <- predict(shuffle_models_CWM_PLH_grow4, gam_data[group_extra4,])$predicted

shuffle_CWM_trait_based_extra4<-data.frame("LNC"=shuffle_models_CWM_LNC_extra4, "SLA"=shuffle_models_CWM_SLA_extra4, "PLH"=shuffle_models_CWM_PLH_extra4)
rownames(shuffle_CWM_trait_based_extra4) <- shuffle_CWM_obs$Site[group_extra4]

## assembly in a single prediction table for the 4 groups :
shuffle_CWM_trait_based_extra <- rbind(shuffle_CWM_trait_based_extra1, shuffle_CWM_trait_based_extra2, shuffle_CWM_trait_based_extra3, shuffle_CWM_trait_based_extra4)
shuffle_CWM_trait_based_extra <- cbind("Site"=rownames(shuffle_CWM_trait_based_extra), shuffle_CWM_trait_based_extra)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_CWM_trait_based_extra <- left_join(Site_order, shuffle_CWM_trait_based_extra, by="Site")
rownames(shuffle_CWM_trait_based_extra) <- shuffle_CWM_trait_based_extra$Site
shuffle_CWM_trait_based_extra <- shuffle_CWM_trait_based_extra[,-1]

###### PREDICTION OF COMMUNITY MEAN  ######
#### -- This computes path predicts the community mean in the context of interpolation -- ####
gam_data <- left_join(shuffle_CM_obs,env_var, by="Site")
## - For interpolation group 1 - ##
shuffle_models_CM_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
shuffle_models_CM_LNC_inter1 <- predict(shuffle_models_CM_LNC_grow1, gam_data[group_inter1,])$predicted
shuffle_models_CM_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
shuffle_models_CM_SLA_inter1 <- predict(shuffle_models_CM_SLA_grow1, gam_data[group_inter1,])$predicted
shuffle_models_CM_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
shuffle_models_CM_PLH_inter1 <- predict(shuffle_models_CM_PLH_grow1, gam_data[group_inter1,])$predicted

shuffle_CM_trait_based_inter1<-data.frame("LNC"=shuffle_models_CM_LNC_inter1, "SLA"=shuffle_models_CM_SLA_inter1, "PLH"=shuffle_models_CM_PLH_inter1)
rownames(shuffle_CM_trait_based_inter1) <- shuffle_CM_obs$Site[group_inter1]

## - For interpolation group 2 - ##
shuffle_models_CM_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
shuffle_models_CM_LNC_inter2 <- predict(shuffle_models_CM_LNC_grow2, gam_data[group_inter2,])$predicted
shuffle_models_CM_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
shuffle_models_CM_SLA_inter2 <- predict(shuffle_models_CM_SLA_grow2, gam_data[group_inter2,])$predicted
shuffle_models_CM_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
shuffle_models_CM_PLH_inter2 <- predict(shuffle_models_CM_PLH_grow2, gam_data[group_inter2,])$predicted

shuffle_CM_trait_based_inter2<-data.frame("LNC"=shuffle_models_CM_LNC_inter2, "SLA"=shuffle_models_CM_SLA_inter2, "PLH"=shuffle_models_CM_PLH_inter2)
rownames(shuffle_CM_trait_based_inter2) <- shuffle_CM_obs$Site[group_inter2]

## - For interpolation group 3 - ##
shuffle_models_CM_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
shuffle_models_CM_LNC_inter3 <- predict(shuffle_models_CM_LNC_grow3, gam_data[group_inter3,])$predicted
shuffle_models_CM_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
shuffle_models_CM_SLA_inter3 <- predict(shuffle_models_CM_SLA_grow3, gam_data[group_inter3,])$predicted
shuffle_models_CM_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
shuffle_models_CM_PLH_inter3 <- predict(shuffle_models_CM_PLH_grow3, gam_data[group_inter3,])$predicted

shuffle_CM_trait_based_inter3<-data.frame("LNC"=shuffle_models_CM_LNC_inter3, "SLA"=shuffle_models_CM_SLA_inter3, "PLH"=shuffle_models_CM_PLH_inter3)
rownames(shuffle_CM_trait_based_inter3) <- shuffle_CM_obs$Site[group_inter3]

## - For interpolation group 4 - ##
shuffle_models_CM_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
shuffle_models_CM_LNC_inter4 <- predict(shuffle_models_CM_LNC_grow4, gam_data[group_inter4,])$predicted
shuffle_models_CM_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
shuffle_models_CM_SLA_inter4 <- predict(shuffle_models_CM_SLA_grow4, gam_data[group_inter4,])$predicted
shuffle_models_CM_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
shuffle_models_CM_PLH_inter4 <- predict(shuffle_models_CM_PLH_grow4, gam_data[group_inter4,])$predicted

shuffle_CM_trait_based_inter4<-data.frame("LNC"=shuffle_models_CM_LNC_inter4, "SLA"=shuffle_models_CM_SLA_inter4, "PLH"=shuffle_models_CM_PLH_inter4)
rownames(shuffle_CM_trait_based_inter4) <- shuffle_CM_obs$Site[group_inter4]

## assembly in a single prediction table for the 4 groups :
shuffle_CM_trait_based_inter <- rbind(shuffle_CM_trait_based_inter1, shuffle_CM_trait_based_inter2, shuffle_CM_trait_based_inter3, shuffle_CM_trait_based_inter4)
shuffle_CM_trait_based_inter <- cbind("Site"=rownames(shuffle_CM_trait_based_inter), shuffle_CM_trait_based_inter)
Site_order <- data.frame("Site"=rownames(plant_pa))
shuffle_CM_trait_based_inter <- left_join(Site_order, shuffle_CM_trait_based_inter, by="Site")
rownames(shuffle_CM_trait_based_inter) <- shuffle_CM_trait_based_inter$Site
shuffle_CM_trait_based_inter <- shuffle_CM_trait_based_inter[,-1]

#### -- This computes path predicts the community mean in the context of extrapolation -- ####
gam_data <- left_join(shuffle_CM_obs,env_var, by="Site")
## - For extrapolation group 1 - ##
shuffle_models_CM_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
shuffle_models_CM_LNC_extra1 <- predict(shuffle_models_CM_LNC_grow1, gam_data[group_extra1,])$predicted
shuffle_models_CM_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
shuffle_models_CM_SLA_extra1 <- predict(shuffle_models_CM_SLA_grow1, gam_data[group_extra1,])$predicted
shuffle_models_CM_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
shuffle_models_CM_PLH_extra1 <- predict(shuffle_models_CM_PLH_grow1, gam_data[group_extra1,])$predicted

shuffle_CM_trait_based_extra1<-data.frame("LNC"=shuffle_models_CM_LNC_extra1, "SLA"=shuffle_models_CM_SLA_extra1, "PLH"=shuffle_models_CM_PLH_extra1)
rownames(shuffle_CM_trait_based_extra1) <- shuffle_CM_obs$Site[group_extra1]

## - For extrapolation group 2 - ##
shuffle_models_CM_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
shuffle_models_CM_LNC_extra2 <- predict(shuffle_models_CM_LNC_grow2, gam_data[group_extra2,])$predicted
shuffle_models_CM_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
shuffle_models_CM_SLA_extra2 <- predict(shuffle_models_CM_SLA_grow2, gam_data[group_extra2,])$predicted
shuffle_models_CM_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
shuffle_models_CM_PLH_extra2 <- predict(shuffle_models_CM_PLH_grow2, gam_data[group_extra2,])$predicted

shuffle_CM_trait_based_extra2<-data.frame("LNC"=shuffle_models_CM_LNC_extra2, "SLA"=shuffle_models_CM_SLA_extra2, "PLH"=shuffle_models_CM_PLH_extra2)
rownames(shuffle_CM_trait_based_extra2) <- shuffle_CM_obs$Site[group_extra2]

## - For extrapolation group 3 - ##
shuffle_models_CM_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
shuffle_models_CM_LNC_extra3 <- predict(shuffle_models_CM_LNC_grow3, gam_data[group_extra3,])$predicted
shuffle_models_CM_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
shuffle_models_CM_SLA_extra3 <- predict(shuffle_models_CM_SLA_grow3, gam_data[group_extra3,])$predicted
shuffle_models_CM_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
shuffle_models_CM_PLH_extra3 <- predict(shuffle_models_CM_PLH_grow3, gam_data[group_extra3,])$predicted

shuffle_CM_trait_based_extra3<-data.frame("LNC"=shuffle_models_CM_LNC_extra3, "SLA"=shuffle_models_CM_SLA_extra3, "PLH"=shuffle_models_CM_PLH_extra3)
rownames(shuffle_CM_trait_based_extra3) <- shuffle_CM_obs$Site[group_extra3]

## - For extrapolation group 4 - ##
shuffle_models_CM_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
shuffle_models_CM_LNC_extra4 <- predict(shuffle_models_CM_LNC_grow4, gam_data[group_extra4,])$predicted
shuffle_models_CM_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
shuffle_models_CM_SLA_extra4 <- predict(shuffle_models_CM_SLA_grow4, gam_data[group_extra4,])$predicted
shuffle_models_CM_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
shuffle_models_CM_PLH_extra4 <- predict(shuffle_models_CM_PLH_grow4, gam_data[group_extra4,])$predicted

shuffle_CM_trait_based_extra4<-data.frame("LNC"=shuffle_models_CM_LNC_extra4, "SLA"=shuffle_models_CM_SLA_extra4, "PLH"=shuffle_models_CM_PLH_extra4)
rownames(shuffle_CM_trait_based_extra4) <- shuffle_CM_obs$Site[group_extra4]

## assembly in a single prediction table for the 4 groups :
shuffle_CM_trait_based_extra <- rbind(shuffle_CM_trait_based_extra1, shuffle_CM_trait_based_extra2, shuffle_CM_trait_based_extra3, shuffle_CM_trait_based_extra4)
shuffle_CM_trait_based_extra <- cbind("Site"=rownames(shuffle_CM_trait_based_extra), shuffle_CM_trait_based_extra)
Site_order <- data.frame("Site"=rownames(plant_pa))
shuffle_CM_trait_based_extra <- left_join(Site_order, shuffle_CM_trait_based_extra, by="Site")
rownames(shuffle_CM_trait_based_extra) <- shuffle_CM_trait_based_extra$Site
shuffle_CM_trait_based_extra <- shuffle_CM_trait_based_extra[,-1]

###### PREDICTION OF FUNCTIONAL DISPERSION  ######
#### -- This computes path predicts the functional dispersion in the context of interpolation -- ####
gam_data <- left_join(shuffle_FDis_obs,env_var, by="Site")
## - For interpolation group 1 - ##
shuffle_models_FDis_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
shuffle_models_FDis_LNC_inter1 <- predict(shuffle_models_FDis_LNC_grow1, gam_data[group_inter1,])$predicted
shuffle_models_FDis_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
shuffle_models_FDis_SLA_inter1 <- predict(shuffle_models_FDis_SLA_grow1, gam_data[group_inter1,])$predicted
shuffle_models_FDis_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
shuffle_models_FDis_PLH_inter1 <- predict(shuffle_models_FDis_PLH_grow1, gam_data[group_inter1,])$predicted

shuffle_FDis_trait_based_inter1<-data.frame("LNC"=shuffle_models_FDis_LNC_inter1, "SLA"=shuffle_models_FDis_SLA_inter1, "PLH"=shuffle_models_FDis_PLH_inter1)
rownames(shuffle_FDis_trait_based_inter1) <- shuffle_FDis_obs$Site[group_inter1]

## - For interpolation group 2 - ##
shuffle_models_FDis_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
shuffle_models_FDis_LNC_inter2 <- predict(shuffle_models_FDis_LNC_grow2, gam_data[group_inter2,])$predicted
shuffle_models_FDis_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
shuffle_models_FDis_SLA_inter2 <- predict(shuffle_models_FDis_SLA_grow2, gam_data[group_inter2,])$predicted
shuffle_models_FDis_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
shuffle_models_FDis_PLH_inter2 <- predict(shuffle_models_FDis_PLH_grow2, gam_data[group_inter2,])$predicted

shuffle_FDis_trait_based_inter2<-data.frame("LNC"=shuffle_models_FDis_LNC_inter2, "SLA"=shuffle_models_FDis_SLA_inter2, "PLH"=shuffle_models_FDis_PLH_inter2)
rownames(shuffle_FDis_trait_based_inter2) <- shuffle_FDis_obs$Site[group_inter2]

## - For interpolation group 3 - ##
shuffle_models_FDis_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
shuffle_models_FDis_LNC_inter3 <- predict(shuffle_models_FDis_LNC_grow3, gam_data[group_inter3,])$predicted
shuffle_models_FDis_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
shuffle_models_FDis_SLA_inter3 <- predict(shuffle_models_FDis_SLA_grow3, gam_data[group_inter3,])$predicted
shuffle_models_FDis_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
shuffle_models_FDis_PLH_inter3 <- predict(shuffle_models_FDis_PLH_grow3, gam_data[group_inter3,])$predicted

shuffle_FDis_trait_based_inter3<-data.frame("LNC"=shuffle_models_FDis_LNC_inter3, "SLA"=shuffle_models_FDis_SLA_inter3, "PLH"=shuffle_models_FDis_PLH_inter3)
rownames(shuffle_FDis_trait_based_inter3) <- shuffle_FDis_obs$Site[group_inter3]

## - For interpolation group 4 - ##
shuffle_models_FDis_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
shuffle_models_FDis_LNC_inter4 <- predict(shuffle_models_FDis_LNC_grow4, gam_data[group_inter4,])$predicted
shuffle_models_FDis_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
shuffle_models_FDis_SLA_inter4 <- predict(shuffle_models_FDis_SLA_grow4, gam_data[group_inter4,])$predicted
shuffle_models_FDis_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
shuffle_models_FDis_PLH_inter4 <- predict(shuffle_models_FDis_PLH_grow4, gam_data[group_inter4,])$predicted

shuffle_FDis_trait_based_inter4<-data.frame("LNC"=shuffle_models_FDis_LNC_inter4, "SLA"=shuffle_models_FDis_SLA_inter4, "PLH"=shuffle_models_FDis_PLH_inter4)
rownames(shuffle_FDis_trait_based_inter4) <- shuffle_FDis_obs$Site[group_inter4]

## assembly in a single prediction table for the 4 groups :
shuffle_FDis_trait_based_inter <- rbind(shuffle_FDis_trait_based_inter1, shuffle_FDis_trait_based_inter2, shuffle_FDis_trait_based_inter3, shuffle_FDis_trait_based_inter4)
shuffle_FDis_trait_based_inter <- cbind("Site"=rownames(shuffle_FDis_trait_based_inter), shuffle_FDis_trait_based_inter)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_FDis_trait_based_inter <- left_join(Site_order, shuffle_FDis_trait_based_inter, by="Site")
rownames(shuffle_FDis_trait_based_inter) <- shuffle_FDis_trait_based_inter$Site
shuffle_FDis_trait_based_inter <- shuffle_FDis_trait_based_inter[,-1]

#### -- This computes path predicts the functional dispersion in the context of extrapolation -- ####
gam_data <- left_join(shuffle_FDis_obs,env_var, by="Site")
## - For extrapolation group 1 - ##
shuffle_models_FDis_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
shuffle_models_FDis_LNC_extra1 <- predict(shuffle_models_FDis_LNC_grow1, gam_data[group_extra1,])$predicted
shuffle_models_FDis_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
shuffle_models_FDis_SLA_extra1 <- predict(shuffle_models_FDis_SLA_grow1, gam_data[group_extra1,])$predicted
shuffle_models_FDis_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
shuffle_models_FDis_PLH_extra1 <- predict(shuffle_models_FDis_PLH_grow1, gam_data[group_extra1,])$predicted

shuffle_FDis_trait_based_extra1<-data.frame("LNC"=shuffle_models_FDis_LNC_extra1, "SLA"=shuffle_models_FDis_SLA_extra1, "PLH"=shuffle_models_FDis_PLH_extra1)
rownames(shuffle_FDis_trait_based_extra1) <- shuffle_FDis_obs$Site[group_extra1]

## - For extrapolation group 2 - ##
shuffle_models_FDis_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
shuffle_models_FDis_LNC_extra2 <- predict(shuffle_models_FDis_LNC_grow2, gam_data[group_extra2,])$predicted
shuffle_models_FDis_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
shuffle_models_FDis_SLA_extra2 <- predict(shuffle_models_FDis_SLA_grow2, gam_data[group_extra2,])$predicted
shuffle_models_FDis_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
shuffle_models_FDis_PLH_extra2 <- predict(shuffle_models_FDis_PLH_grow2, gam_data[group_extra2,])$predicted

shuffle_FDis_trait_based_extra2<-data.frame("LNC"=shuffle_models_FDis_LNC_extra2, "SLA"=shuffle_models_FDis_SLA_extra2, "PLH"=shuffle_models_FDis_PLH_extra2)
rownames(shuffle_FDis_trait_based_extra2) <- shuffle_FDis_obs$Site[group_extra2]

## - For extrapolation group 3 - ##
shuffle_models_FDis_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
shuffle_models_FDis_LNC_extra3 <- predict(shuffle_models_FDis_LNC_grow3, gam_data[group_extra3,])$predicted
shuffle_models_FDis_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
shuffle_models_FDis_SLA_extra3 <- predict(shuffle_models_FDis_SLA_grow3, gam_data[group_extra3,])$predicted
shuffle_models_FDis_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
shuffle_models_FDis_PLH_extra3 <- predict(shuffle_models_FDis_PLH_grow3, gam_data[group_extra3,])$predicted

shuffle_FDis_trait_based_extra3<-data.frame("LNC"=shuffle_models_FDis_LNC_extra3, "SLA"=shuffle_models_FDis_SLA_extra3, "PLH"=shuffle_models_FDis_PLH_extra3)
rownames(shuffle_FDis_trait_based_extra3) <- shuffle_FDis_obs$Site[group_extra3]

## - For extrapolation group 4 - ##
shuffle_models_FDis_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
shuffle_models_FDis_LNC_extra4 <- predict(shuffle_models_FDis_LNC_grow4, gam_data[group_extra4,])$predicted
shuffle_models_FDis_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
shuffle_models_FDis_SLA_extra4 <- predict(shuffle_models_FDis_SLA_grow4, gam_data[group_extra4,])$predicted
shuffle_models_FDis_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
shuffle_models_FDis_PLH_extra4 <- predict(shuffle_models_FDis_PLH_grow4, gam_data[group_extra4,])$predicted

shuffle_FDis_trait_based_extra4<-data.frame("LNC"=shuffle_models_FDis_LNC_extra4, "SLA"=shuffle_models_FDis_SLA_extra4, "PLH"=shuffle_models_FDis_PLH_extra4)
rownames(shuffle_FDis_trait_based_extra4) <- shuffle_FDis_obs$Site[group_extra4]

## assembly in a single prediction table for the 4 groups :
shuffle_FDis_trait_based_extra <- rbind(shuffle_FDis_trait_based_extra1, shuffle_FDis_trait_based_extra2, shuffle_FDis_trait_based_extra3, shuffle_FDis_trait_based_extra4)
shuffle_FDis_trait_based_extra <- cbind("Site"=rownames(shuffle_FDis_trait_based_extra), shuffle_FDis_trait_based_extra)
Site_order <- data.frame("Site"=rownames(plant_recovery))
shuffle_FDis_trait_based_extra <- left_join(Site_order, shuffle_FDis_trait_based_extra, by="Site")
rownames(shuffle_FDis_trait_based_extra) <- shuffle_FDis_trait_based_extra$Site
shuffle_FDis_trait_based_extra <- shuffle_FDis_trait_based_extra[,-1]

###### PREDICTION OF UNWEIGHTED FUNCTIONAL DISPERSION  ######
#### -- This computes path predicts the unweighted functional dispersion in the context of interpolation -- ####
gam_data <- left_join(shuffle_uFDis_obs,env_var, by="Site")
## - For interpolation group 1 - ##
shuffle_models_uFDis_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
shuffle_models_uFDis_LNC_inter1 <- predict(shuffle_models_uFDis_LNC_grow1, gam_data[group_inter1,])$predicted
shuffle_models_uFDis_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
shuffle_models_uFDis_SLA_inter1 <- predict(shuffle_models_uFDis_SLA_grow1, gam_data[group_inter1,])$predicted
shuffle_models_uFDis_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
shuffle_models_uFDis_PLH_inter1 <- predict(shuffle_models_uFDis_PLH_grow1, gam_data[group_inter1,])$predicted

shuffle_uFDis_trait_based_inter1<-data.frame("LNC"=shuffle_models_uFDis_LNC_inter1, "SLA"=shuffle_models_uFDis_SLA_inter1, "PLH"=shuffle_models_uFDis_PLH_inter1)
rownames(shuffle_uFDis_trait_based_inter1) <- shuffle_uFDis_obs$Site[group_inter1]

## - For interpolation group 2 - ##
shuffle_models_uFDis_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
shuffle_models_uFDis_LNC_inter2 <- predict(shuffle_models_uFDis_LNC_grow2, gam_data[group_inter2,])$predicted
shuffle_models_uFDis_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
shuffle_models_uFDis_SLA_inter2 <- predict(shuffle_models_uFDis_SLA_grow2, gam_data[group_inter2,])$predicted
shuffle_models_uFDis_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
shuffle_models_uFDis_PLH_inter2 <- predict(shuffle_models_uFDis_PLH_grow2, gam_data[group_inter2,])$predicted

shuffle_uFDis_trait_based_inter2<-data.frame("LNC"=shuffle_models_uFDis_LNC_inter2, "SLA"=shuffle_models_uFDis_SLA_inter2, "PLH"=shuffle_models_uFDis_PLH_inter2)
rownames(shuffle_uFDis_trait_based_inter2) <- shuffle_uFDis_obs$Site[group_inter2]

## - For interpolation group 3 - ##
shuffle_models_uFDis_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
shuffle_models_uFDis_LNC_inter3 <- predict(shuffle_models_uFDis_LNC_grow3, gam_data[group_inter3,])$predicted
shuffle_models_uFDis_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
shuffle_models_uFDis_SLA_inter3 <- predict(shuffle_models_uFDis_SLA_grow3, gam_data[group_inter3,])$predicted
shuffle_models_uFDis_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
shuffle_models_uFDis_PLH_inter3 <- predict(shuffle_models_uFDis_PLH_grow3, gam_data[group_inter3,])$predicted

shuffle_uFDis_trait_based_inter3<-data.frame("LNC"=shuffle_models_uFDis_LNC_inter3, "SLA"=shuffle_models_uFDis_SLA_inter3, "PLH"=shuffle_models_uFDis_PLH_inter3)
rownames(shuffle_uFDis_trait_based_inter3) <- shuffle_uFDis_obs$Site[group_inter3]

## - For interpolation group 4 - ##
shuffle_models_uFDis_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
shuffle_models_uFDis_LNC_inter4 <- predict(shuffle_models_uFDis_LNC_grow4, gam_data[group_inter4,])$predicted
shuffle_models_uFDis_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
shuffle_models_uFDis_SLA_inter4 <- predict(shuffle_models_uFDis_SLA_grow4, gam_data[group_inter4,])$predicted
shuffle_models_uFDis_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
shuffle_models_uFDis_PLH_inter4 <- predict(shuffle_models_uFDis_PLH_grow4, gam_data[group_inter4,])$predicted

shuffle_uFDis_trait_based_inter4<-data.frame("LNC"=shuffle_models_uFDis_LNC_inter4, "SLA"=shuffle_models_uFDis_SLA_inter4, "PLH"=shuffle_models_uFDis_PLH_inter4)
rownames(shuffle_uFDis_trait_based_inter4) <- shuffle_uFDis_obs$Site[group_inter4]

## assembly in a single prediction table for the 4 groups :
shuffle_uFDis_trait_based_inter <- rbind(shuffle_uFDis_trait_based_inter1, shuffle_uFDis_trait_based_inter2, shuffle_uFDis_trait_based_inter3, shuffle_uFDis_trait_based_inter4)
shuffle_uFDis_trait_based_inter <- cbind("Site"=rownames(shuffle_uFDis_trait_based_inter), shuffle_uFDis_trait_based_inter)
Site_order <- data.frame("Site"=rownames(plant_pa))
shuffle_uFDis_trait_based_inter <- left_join(Site_order, shuffle_uFDis_trait_based_inter, by="Site")
rownames(shuffle_uFDis_trait_based_inter) <- shuffle_uFDis_trait_based_inter$Site
shuffle_uFDis_trait_based_inter <- shuffle_uFDis_trait_based_inter[,-1]

#### -- This computes path predicts the community mean in the context of extrapolation -- ####
## - For extrapolation group 1 - ##
shuffle_models_uFDis_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
shuffle_models_uFDis_LNC_extra1 <- predict(shuffle_models_uFDis_LNC_grow1, gam_data[group_extra1,])$predicted
shuffle_models_uFDis_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
shuffle_models_uFDis_SLA_extra1 <- predict(shuffle_models_uFDis_SLA_grow1, gam_data[group_extra1,])$predicted
shuffle_models_uFDis_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
shuffle_models_uFDis_PLH_extra1 <- predict(shuffle_models_uFDis_PLH_grow1, gam_data[group_extra1,])$predicted

shuffle_uFDis_trait_based_extra1<-data.frame("LNC"=shuffle_models_uFDis_LNC_extra1, "SLA"=shuffle_models_uFDis_SLA_extra1, "PLH"=shuffle_models_uFDis_PLH_extra1)
rownames(shuffle_uFDis_trait_based_extra1) <- shuffle_uFDis_obs$Site[group_extra1]

## - For extrapolation group 2 - ##
shuffle_models_uFDis_LNC_grow2 <- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
shuffle_models_uFDis_LNC_extra2 <- predict(shuffle_models_uFDis_LNC_grow2, gam_data[group_extra2,])$predicted
shuffle_models_uFDis_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
shuffle_models_uFDis_SLA_extra2 <- predict(shuffle_models_uFDis_SLA_grow2, gam_data[group_extra2,])$predicted
shuffle_models_uFDis_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
shuffle_models_uFDis_PLH_extra2 <- predict(shuffle_models_uFDis_PLH_grow2, gam_data[group_extra2,])$predicted

shuffle_uFDis_trait_based_extra2<-data.frame("LNC"=shuffle_models_uFDis_LNC_extra2, "SLA"=shuffle_models_uFDis_SLA_extra2, "PLH"=shuffle_models_uFDis_PLH_extra2)
rownames(shuffle_uFDis_trait_based_extra2) <- shuffle_uFDis_obs$Site[group_extra2]

## - For extrapolation group 3 - ##
shuffle_models_uFDis_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
shuffle_models_uFDis_LNC_extra3 <- predict(shuffle_models_uFDis_LNC_grow3, gam_data[group_extra3,])$predicted
shuffle_models_uFDis_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
shuffle_models_uFDis_SLA_extra3 <- predict(shuffle_models_uFDis_SLA_grow3, gam_data[group_extra3,])$predicted
shuffle_models_uFDis_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
shuffle_models_uFDis_PLH_extra3 <- predict(shuffle_models_uFDis_PLH_grow3, gam_data[group_extra3,])$predicted

shuffle_uFDis_trait_based_extra3<-data.frame("LNC"=shuffle_models_uFDis_LNC_extra3, "SLA"=shuffle_models_uFDis_SLA_extra3, "PLH"=shuffle_models_uFDis_PLH_extra3)
rownames(shuffle_uFDis_trait_based_extra3) <- shuffle_uFDis_obs$Site[group_extra3]

## - For extrapolation group 4 - ##
shuffle_models_uFDis_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
shuffle_models_uFDis_LNC_extra4 <- predict(shuffle_models_uFDis_LNC_grow4, gam_data[group_extra4,])$predicted
shuffle_models_uFDis_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
shuffle_models_uFDis_SLA_extra4 <- predict(shuffle_models_uFDis_SLA_grow4, gam_data[group_extra4,])$predicted
shuffle_models_uFDis_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
shuffle_models_uFDis_PLH_extra4 <- predict(shuffle_models_uFDis_PLH_grow4, gam_data[group_extra4,])$predicted

shuffle_uFDis_trait_based_extra4<-data.frame("LNC"=shuffle_models_uFDis_LNC_extra4, "SLA"=shuffle_models_uFDis_SLA_extra4, "PLH"=shuffle_models_uFDis_PLH_extra4)
rownames(shuffle_uFDis_trait_based_extra4) <- shuffle_uFDis_obs$Site[group_extra4]

## assembly in a single prediction table for the 4 groups :
shuffle_uFDis_trait_based_extra <- rbind(shuffle_uFDis_trait_based_extra1, shuffle_uFDis_trait_based_extra2, shuffle_uFDis_trait_based_extra3, shuffle_uFDis_trait_based_extra4)
shuffle_uFDis_trait_based_extra <- cbind("Site"=rownames(shuffle_uFDis_trait_based_extra), shuffle_uFDis_trait_based_extra)
Site_order <- data.frame("Site"=rownames(plant_pa))
shuffle_uFDis_trait_based_extra <- left_join(Site_order, shuffle_uFDis_trait_based_extra, by="Site")
rownames(shuffle_uFDis_trait_based_extra) <- shuffle_uFDis_trait_based_extra$Site
shuffle_uFDis_trait_based_extra <- shuffle_uFDis_trait_based_extra[,-1]


shuffle_figure3 <- data.frame("RMSE"=c(rmse(shuffle_CWM_obs$LNC[group_inter1], shuffle_CWM_trait_based_inter$LNC[group_inter1]), # CWM trait based approach interpolation
                                       rmse(shuffle_CWM_obs$SLA[group_inter1], shuffle_CWM_trait_based_inter$SLA[group_inter1]),
                                       rmse(shuffle_CWM_obs$PLH[group_inter1], shuffle_CWM_trait_based_inter$PLH[group_inter1]),
                                       rmse(shuffle_CWM_obs$LNC[group_inter2], shuffle_CWM_trait_based_inter$LNC[group_inter2]),
                                       rmse(shuffle_CWM_obs$SLA[group_inter2], shuffle_CWM_trait_based_inter$SLA[group_inter2]),
                                       rmse(shuffle_CWM_obs$PLH[group_inter2], shuffle_CWM_trait_based_inter$PLH[group_inter2]),
                                       rmse(shuffle_CWM_obs$LNC[group_inter3], shuffle_CWM_trait_based_inter$LNC[group_inter3]),
                                       rmse(shuffle_CWM_obs$SLA[group_inter3], shuffle_CWM_trait_based_inter$SLA[group_inter3]),
                                       rmse(shuffle_CWM_obs$PLH[group_inter3], shuffle_CWM_trait_based_inter$PLH[group_inter3]),
                                       rmse(shuffle_CWM_obs$LNC[group_inter4], shuffle_CWM_trait_based_inter$LNC[group_inter4]),
                                       rmse(shuffle_CWM_obs$SLA[group_inter4], shuffle_CWM_trait_based_inter$SLA[group_inter4]),
                                       rmse(shuffle_CWM_obs$PLH[group_inter4], shuffle_CWM_trait_based_inter$PLH[group_inter4]),
                                       rmse(shuffle_CWM_obs$LNC[group_extra1], shuffle_CWM_trait_based_extra$LNC[group_extra1]),# CWM trait based approach interpolation
                                       rmse(shuffle_CWM_obs$SLA[group_extra1], shuffle_CWM_trait_based_extra$SLA[group_extra1]),
                                       rmse(shuffle_CWM_obs$PLH[group_extra1], shuffle_CWM_trait_based_extra$PLH[group_extra1]),
                                       rmse(shuffle_CWM_obs$LNC[group_extra2], shuffle_CWM_trait_based_extra$LNC[group_extra2]),
                                       rmse(shuffle_CWM_obs$SLA[group_extra2], shuffle_CWM_trait_based_extra$SLA[group_extra2]),
                                       rmse(shuffle_CWM_obs$PLH[group_extra2], shuffle_CWM_trait_based_extra$PLH[group_extra2]),
                                       rmse(shuffle_CWM_obs$LNC[group_extra3], shuffle_CWM_trait_based_extra$LNC[group_extra3]),
                                       rmse(shuffle_CWM_obs$SLA[group_extra3], shuffle_CWM_trait_based_extra$SLA[group_extra3]),
                                       rmse(shuffle_CWM_obs$PLH[group_extra3], shuffle_CWM_trait_based_extra$PLH[group_extra3]),
                                       rmse(shuffle_CWM_obs$LNC[group_extra4], shuffle_CWM_trait_based_extra$LNC[group_extra4]),
                                       rmse(shuffle_CWM_obs$SLA[group_extra4], shuffle_CWM_trait_based_extra$SLA[group_extra4]),
                                       rmse(shuffle_CWM_obs$PLH[group_extra4], shuffle_CWM_trait_based_extra$PLH[group_extra4]),
                                       rmse(shuffle_FDis_obs$LNC[group_inter1], shuffle_FDis_trait_based_inter$LNC[group_inter1]),# FDis trait based approach interpolation
                                       rmse(shuffle_FDis_obs$SLA[group_inter1], shuffle_FDis_trait_based_inter$SLA[group_inter1]),
                                       rmse(shuffle_FDis_obs$PLH[group_inter1], shuffle_FDis_trait_based_inter$PLH[group_inter1]),
                                       rmse(shuffle_FDis_obs$LNC[group_inter2], shuffle_FDis_trait_based_inter$LNC[group_inter2]),
                                       rmse(shuffle_FDis_obs$SLA[group_inter2], shuffle_FDis_trait_based_inter$SLA[group_inter2]),
                                       rmse(shuffle_FDis_obs$PLH[group_inter2], shuffle_FDis_trait_based_inter$PLH[group_inter2]),
                                       rmse(shuffle_FDis_obs$LNC[group_inter3], shuffle_FDis_trait_based_inter$LNC[group_inter3]),
                                       rmse(shuffle_FDis_obs$SLA[group_inter3], shuffle_FDis_trait_based_inter$SLA[group_inter3]),
                                       rmse(shuffle_FDis_obs$PLH[group_inter3], shuffle_FDis_trait_based_inter$PLH[group_inter3]),
                                       rmse(shuffle_FDis_obs$LNC[group_inter4], shuffle_FDis_trait_based_inter$LNC[group_inter4]),
                                       rmse(shuffle_FDis_obs$SLA[group_inter4], shuffle_FDis_trait_based_inter$SLA[group_inter4]),
                                       rmse(shuffle_FDis_obs$PLH[group_inter4], shuffle_FDis_trait_based_inter$PLH[group_inter4]),
                                       rmse(shuffle_FDis_obs$LNC[group_extra1], shuffle_FDis_trait_based_extra$LNC[group_extra1]), # FDis trait based approach extrapolation
                                       rmse(shuffle_FDis_obs$SLA[group_extra1], shuffle_FDis_trait_based_extra$SLA[group_extra1]),
                                       rmse(shuffle_FDis_obs$PLH[group_extra1], shuffle_FDis_trait_based_extra$PLH[group_extra1]),
                                       rmse(shuffle_FDis_obs$LNC[group_extra2], shuffle_FDis_trait_based_extra$LNC[group_extra2]),
                                       rmse(shuffle_FDis_obs$SLA[group_extra2], shuffle_FDis_trait_based_extra$SLA[group_extra2]),
                                       rmse(shuffle_FDis_obs$PLH[group_extra2], shuffle_FDis_trait_based_extra$PLH[group_extra2]),
                                       rmse(shuffle_FDis_obs$LNC[group_extra3], shuffle_FDis_trait_based_extra$LNC[group_extra3]),
                                       rmse(shuffle_FDis_obs$SLA[group_extra3], shuffle_FDis_trait_based_extra$SLA[group_extra3]),
                                       rmse(shuffle_FDis_obs$PLH[group_extra3], shuffle_FDis_trait_based_extra$PLH[group_extra3]),
                                       rmse(shuffle_FDis_obs$LNC[group_extra4], shuffle_FDis_trait_based_extra$LNC[group_extra4]),
                                       rmse(shuffle_FDis_obs$SLA[group_extra4], shuffle_FDis_trait_based_extra$SLA[group_extra4]),
                                       rmse(shuffle_FDis_obs$PLH[group_extra4], shuffle_FDis_trait_based_extra$PLH[group_extra4]),
                                       rmse(shuffle_CWM_obs$LNC[group_inter1], shuffle_CWM_species_based_inter$LNC[group_inter1]), # CWM species based approach interpolation
                                       rmse(shuffle_CWM_obs$SLA[group_inter1], shuffle_CWM_species_based_inter$SLA[group_inter1]),
                                       rmse(shuffle_CWM_obs$PLH[group_inter1], shuffle_CWM_species_based_inter$PLH[group_inter1]),
                                       rmse(shuffle_CWM_obs$LNC[group_inter2], shuffle_CWM_species_based_inter$LNC[group_inter2]),
                                       rmse(shuffle_CWM_obs$SLA[group_inter2], shuffle_CWM_species_based_inter$SLA[group_inter2]),
                                       rmse(shuffle_CWM_obs$PLH[group_inter2], shuffle_CWM_species_based_inter$PLH[group_inter2]),
                                       rmse(shuffle_CWM_obs$LNC[group_inter3], shuffle_CWM_species_based_inter$LNC[group_inter3]),
                                       rmse(shuffle_CWM_obs$SLA[group_inter3], shuffle_CWM_species_based_inter$SLA[group_inter3]),
                                       rmse(shuffle_CWM_obs$PLH[group_inter3], shuffle_CWM_species_based_inter$PLH[group_inter3]),
                                       rmse(shuffle_CWM_obs$LNC[group_inter4], shuffle_CWM_species_based_inter$LNC[group_inter4]),
                                       rmse(shuffle_CWM_obs$SLA[group_inter4], shuffle_CWM_species_based_inter$SLA[group_inter4]),
                                       rmse(shuffle_CWM_obs$PLH[group_inter4], shuffle_CWM_species_based_inter$PLH[group_inter4]),
                                       rmse(shuffle_CWM_obs$LNC[group_extra1], shuffle_CWM_species_based_extra$LNC[group_extra1]), # CWM species based approach extrapolation
                                       rmse(shuffle_CWM_obs$SLA[group_extra1], shuffle_CWM_species_based_extra$SLA[group_extra1]),
                                       rmse(shuffle_CWM_obs$PLH[group_extra1], shuffle_CWM_species_based_extra$PLH[group_extra1]),
                                       rmse(shuffle_CWM_obs$LNC[group_extra2], shuffle_CWM_species_based_extra$LNC[group_extra2]),
                                       rmse(shuffle_CWM_obs$SLA[group_extra2], shuffle_CWM_species_based_extra$SLA[group_extra2]),
                                       rmse(shuffle_CWM_obs$PLH[group_extra2], shuffle_CWM_species_based_extra$PLH[group_extra2]),
                                       rmse(shuffle_CWM_obs$LNC[group_extra3], shuffle_CWM_species_based_extra$LNC[group_extra3]),
                                       rmse(shuffle_CWM_obs$SLA[group_extra3], shuffle_CWM_species_based_extra$SLA[group_extra3]),
                                       rmse(shuffle_CWM_obs$PLH[group_extra3], shuffle_CWM_species_based_extra$PLH[group_extra3]),
                                       rmse(shuffle_CWM_obs$LNC[group_extra4], shuffle_CWM_species_based_extra$LNC[group_extra4]),
                                       rmse(shuffle_CWM_obs$SLA[group_extra4], shuffle_CWM_species_based_extra$SLA[group_extra4]),
                                       rmse(shuffle_CWM_obs$PLH[group_extra4], shuffle_CWM_species_based_extra$PLH[group_extra4]),
                                       rmse(shuffle_FDis_obs$LNC[group_inter1], shuffle_FDis_species_based_inter$LNC[group_inter1]), # FDis species based approach interpolation
                                       rmse(shuffle_FDis_obs$SLA[group_inter1], shuffle_FDis_species_based_inter$SLA[group_inter1]),
                                       rmse(shuffle_FDis_obs$PLH[group_inter1], shuffle_FDis_species_based_inter$PLH[group_inter1]),
                                       rmse(shuffle_FDis_obs$LNC[group_inter2], shuffle_FDis_species_based_inter$LNC[group_inter2]),
                                       rmse(shuffle_FDis_obs$SLA[group_inter2], shuffle_FDis_species_based_inter$SLA[group_inter2]),
                                       rmse(shuffle_FDis_obs$PLH[group_inter2], shuffle_FDis_species_based_inter$PLH[group_inter2]),
                                       rmse(shuffle_FDis_obs$LNC[group_inter3], shuffle_FDis_species_based_inter$LNC[group_inter3]),
                                       rmse(shuffle_FDis_obs$SLA[group_inter3], shuffle_FDis_species_based_inter$SLA[group_inter3]),
                                       rmse(shuffle_FDis_obs$PLH[group_inter3], shuffle_FDis_species_based_inter$PLH[group_inter3]),
                                       rmse(shuffle_FDis_obs$LNC[group_inter4], shuffle_FDis_species_based_inter$LNC[group_inter4]),
                                       rmse(shuffle_FDis_obs$SLA[group_inter4], shuffle_FDis_species_based_inter$SLA[group_inter4]),
                                       rmse(shuffle_FDis_obs$PLH[group_inter4], shuffle_FDis_species_based_inter$PLH[group_inter4]),
                                       rmse(shuffle_FDis_obs$LNC[group_extra1], shuffle_FDis_species_based_extra$LNC[group_extra1]), # FDis species based approach extrapolation
                                       rmse(shuffle_FDis_obs$SLA[group_extra1], shuffle_FDis_species_based_extra$SLA[group_extra1]),
                                       rmse(shuffle_FDis_obs$PLH[group_extra1], shuffle_FDis_species_based_extra$PLH[group_extra1]),
                                       rmse(shuffle_FDis_obs$LNC[group_extra2], shuffle_FDis_species_based_extra$LNC[group_extra2]),
                                       rmse(shuffle_FDis_obs$SLA[group_extra2], shuffle_FDis_species_based_extra$SLA[group_extra2]),
                                       rmse(shuffle_FDis_obs$PLH[group_extra2], shuffle_FDis_species_based_extra$PLH[group_extra2]),
                                       rmse(shuffle_FDis_obs$LNC[group_extra3], shuffle_FDis_species_based_extra$LNC[group_extra3]),
                                       rmse(shuffle_FDis_obs$SLA[group_extra3], shuffle_FDis_species_based_extra$SLA[group_extra3]),
                                       rmse(shuffle_FDis_obs$PLH[group_extra3], shuffle_FDis_species_based_extra$PLH[group_extra3]),
                                       rmse(shuffle_FDis_obs$LNC[group_extra4], shuffle_FDis_species_based_extra$LNC[group_extra4]),
                                       rmse(shuffle_FDis_obs$SLA[group_extra4], shuffle_FDis_species_based_extra$SLA[group_extra4]),
                                       rmse(shuffle_FDis_obs$PLH[group_extra4], shuffle_FDis_species_based_extra$PLH[group_extra4]),
                                       
                                       rmse(shuffle_CM_obs$LNC[group_inter1], shuffle_CM_trait_based_inter$LNC[group_inter1]), # CM trait based approach interpolation
                                       rmse(shuffle_CM_obs$SLA[group_inter1], shuffle_CM_trait_based_inter$SLA[group_inter1]),
                                       rmse(shuffle_CM_obs$PLH[group_inter1], shuffle_CM_trait_based_inter$PLH[group_inter1]),
                                       rmse(shuffle_CM_obs$LNC[group_inter2], shuffle_CM_trait_based_inter$LNC[group_inter2]),
                                       rmse(shuffle_CM_obs$SLA[group_inter2], shuffle_CM_trait_based_inter$SLA[group_inter2]),
                                       rmse(shuffle_CM_obs$PLH[group_inter2], shuffle_CM_trait_based_inter$PLH[group_inter2]),
                                       rmse(shuffle_CM_obs$LNC[group_inter3], shuffle_CM_trait_based_inter$LNC[group_inter3]),
                                       rmse(shuffle_CM_obs$SLA[group_inter3], shuffle_CM_trait_based_inter$SLA[group_inter3]),
                                       rmse(shuffle_CM_obs$PLH[group_inter3], shuffle_CM_trait_based_inter$PLH[group_inter3]),
                                       rmse(shuffle_CM_obs$LNC[group_inter4], shuffle_CM_trait_based_inter$LNC[group_inter4]),
                                       rmse(shuffle_CM_obs$SLA[group_inter4], shuffle_CM_trait_based_inter$SLA[group_inter4]),
                                       rmse(shuffle_CM_obs$PLH[group_inter4], shuffle_CM_trait_based_inter$PLH[group_inter4]),
                                       rmse(shuffle_CM_obs$LNC[group_extra1], shuffle_CM_trait_based_extra$LNC[group_extra1]),# CM trait based approach interpolation
                                       rmse(shuffle_CM_obs$SLA[group_extra1], shuffle_CM_trait_based_extra$SLA[group_extra1]),
                                       rmse(shuffle_CM_obs$PLH[group_extra1], shuffle_CM_trait_based_extra$PLH[group_extra1]),
                                       rmse(shuffle_CM_obs$LNC[group_extra2], shuffle_CM_trait_based_extra$LNC[group_extra2]),
                                       rmse(shuffle_CM_obs$SLA[group_extra2], shuffle_CM_trait_based_extra$SLA[group_extra2]),
                                       rmse(shuffle_CM_obs$PLH[group_extra2], shuffle_CM_trait_based_extra$PLH[group_extra2]),
                                       rmse(shuffle_CM_obs$LNC[group_extra3], shuffle_CM_trait_based_extra$LNC[group_extra3]),
                                       rmse(shuffle_CM_obs$SLA[group_extra3], shuffle_CM_trait_based_extra$SLA[group_extra3]),
                                       rmse(shuffle_CM_obs$PLH[group_extra3], shuffle_CM_trait_based_extra$PLH[group_extra3]),
                                       rmse(shuffle_CM_obs$LNC[group_extra4], shuffle_CM_trait_based_extra$LNC[group_extra4]),
                                       rmse(shuffle_CM_obs$SLA[group_extra4], shuffle_CM_trait_based_extra$SLA[group_extra4]),
                                       rmse(shuffle_CM_obs$PLH[group_extra4], shuffle_CM_trait_based_extra$PLH[group_extra4]),
                                       rmse(shuffle_uFDis_obs$LNC[group_inter1], shuffle_uFDis_trait_based_inter$LNC[group_inter1]),# uuFDis trait based approach interpolation
                                       rmse(shuffle_uFDis_obs$SLA[group_inter1], shuffle_uFDis_trait_based_inter$SLA[group_inter1]),
                                       rmse(shuffle_uFDis_obs$PLH[group_inter1], shuffle_uFDis_trait_based_inter$PLH[group_inter1]),
                                       rmse(shuffle_uFDis_obs$LNC[group_inter2], shuffle_uFDis_trait_based_inter$LNC[group_inter2]),
                                       rmse(shuffle_uFDis_obs$SLA[group_inter2], shuffle_uFDis_trait_based_inter$SLA[group_inter2]),
                                       rmse(shuffle_uFDis_obs$PLH[group_inter2], shuffle_uFDis_trait_based_inter$PLH[group_inter2]),
                                       rmse(shuffle_uFDis_obs$LNC[group_inter3], shuffle_uFDis_trait_based_inter$LNC[group_inter3]),
                                       rmse(shuffle_uFDis_obs$SLA[group_inter3], shuffle_uFDis_trait_based_inter$SLA[group_inter3]),
                                       rmse(shuffle_uFDis_obs$PLH[group_inter3], shuffle_uFDis_trait_based_inter$PLH[group_inter3]),
                                       rmse(shuffle_uFDis_obs$LNC[group_inter4], shuffle_uFDis_trait_based_inter$LNC[group_inter4]),
                                       rmse(shuffle_uFDis_obs$SLA[group_inter4], shuffle_uFDis_trait_based_inter$SLA[group_inter4]),
                                       rmse(shuffle_uFDis_obs$PLH[group_inter4], shuffle_uFDis_trait_based_inter$PLH[group_inter4]),
                                       rmse(shuffle_uFDis_obs$LNC[group_extra1], shuffle_uFDis_trait_based_extra$LNC[group_extra1]), # uuFDis trait based approach extrapolation
                                       rmse(shuffle_uFDis_obs$SLA[group_extra1], shuffle_uFDis_trait_based_extra$SLA[group_extra1]),
                                       rmse(shuffle_uFDis_obs$PLH[group_extra1], shuffle_uFDis_trait_based_extra$PLH[group_extra1]),
                                       rmse(shuffle_uFDis_obs$LNC[group_extra2], shuffle_uFDis_trait_based_extra$LNC[group_extra2]),
                                       rmse(shuffle_uFDis_obs$SLA[group_extra2], shuffle_uFDis_trait_based_extra$SLA[group_extra2]),
                                       rmse(shuffle_uFDis_obs$PLH[group_extra2], shuffle_uFDis_trait_based_extra$PLH[group_extra2]),
                                       rmse(shuffle_uFDis_obs$LNC[group_extra3], shuffle_uFDis_trait_based_extra$LNC[group_extra3]),
                                       rmse(shuffle_uFDis_obs$SLA[group_extra3], shuffle_uFDis_trait_based_extra$SLA[group_extra3]),
                                       rmse(shuffle_uFDis_obs$PLH[group_extra3], shuffle_uFDis_trait_based_extra$PLH[group_extra3]),
                                       rmse(shuffle_uFDis_obs$LNC[group_extra4], shuffle_uFDis_trait_based_extra$LNC[group_extra4]),
                                       rmse(shuffle_uFDis_obs$SLA[group_extra4], shuffle_uFDis_trait_based_extra$SLA[group_extra4]),
                                       rmse(shuffle_uFDis_obs$PLH[group_extra4], shuffle_uFDis_trait_based_extra$PLH[group_extra4]),
                                       rmse(shuffle_CM_obs$LNC[group_inter1], shuffle_CM_species_based_inter$LNC[group_inter1]), # CM species based approach interpolation
                                       rmse(shuffle_CM_obs$SLA[group_inter1], shuffle_CM_species_based_inter$SLA[group_inter1]),
                                       rmse(shuffle_CM_obs$PLH[group_inter1], shuffle_CM_species_based_inter$PLH[group_inter1]),
                                       rmse(shuffle_CM_obs$LNC[group_inter2], shuffle_CM_species_based_inter$LNC[group_inter2]),
                                       rmse(shuffle_CM_obs$SLA[group_inter2], shuffle_CM_species_based_inter$SLA[group_inter2]),
                                       rmse(shuffle_CM_obs$PLH[group_inter2], shuffle_CM_species_based_inter$PLH[group_inter2]),
                                       rmse(shuffle_CM_obs$LNC[group_inter3], shuffle_CM_species_based_inter$LNC[group_inter3]),
                                       rmse(shuffle_CM_obs$SLA[group_inter3], shuffle_CM_species_based_inter$SLA[group_inter3]),
                                       rmse(shuffle_CM_obs$PLH[group_inter3], shuffle_CM_species_based_inter$PLH[group_inter3]),
                                       rmse(shuffle_CM_obs$LNC[group_inter4], shuffle_CM_species_based_inter$LNC[group_inter4]),
                                       rmse(shuffle_CM_obs$SLA[group_inter4], shuffle_CM_species_based_inter$SLA[group_inter4]),
                                       rmse(shuffle_CM_obs$PLH[group_inter4], shuffle_CM_species_based_inter$PLH[group_inter4]),
                                       rmse(shuffle_CM_obs$LNC[group_extra1], shuffle_CM_species_based_extra$LNC[group_extra1]), # CM species based approach extrapolation
                                       rmse(shuffle_CM_obs$SLA[group_extra1], shuffle_CM_species_based_extra$SLA[group_extra1]),
                                       rmse(shuffle_CM_obs$PLH[group_extra1], shuffle_CM_species_based_extra$PLH[group_extra1]),
                                       rmse(shuffle_CM_obs$LNC[group_extra2], shuffle_CM_species_based_extra$LNC[group_extra2]),
                                       rmse(shuffle_CM_obs$SLA[group_extra2], shuffle_CM_species_based_extra$SLA[group_extra2]),
                                       rmse(shuffle_CM_obs$PLH[group_extra2], shuffle_CM_species_based_extra$PLH[group_extra2]),
                                       rmse(shuffle_CM_obs$LNC[group_extra3], shuffle_CM_species_based_extra$LNC[group_extra3]),
                                       rmse(shuffle_CM_obs$SLA[group_extra3], shuffle_CM_species_based_extra$SLA[group_extra3]),
                                       rmse(shuffle_CM_obs$PLH[group_extra3], shuffle_CM_species_based_extra$PLH[group_extra3]),
                                       rmse(shuffle_CM_obs$LNC[group_extra4], shuffle_CM_species_based_extra$LNC[group_extra4]),
                                       rmse(shuffle_CM_obs$SLA[group_extra4], shuffle_CM_species_based_extra$SLA[group_extra4]),
                                       rmse(shuffle_CM_obs$PLH[group_extra4], shuffle_CM_species_based_extra$PLH[group_extra4]),
                                       rmse(shuffle_uFDis_obs$LNC[group_inter1], shuffle_uFDis_species_based_inter$LNC[group_inter1]), # uuFDis species based approach interpolation
                                       rmse(shuffle_uFDis_obs$SLA[group_inter1], shuffle_uFDis_species_based_inter$SLA[group_inter1]),
                                       rmse(shuffle_uFDis_obs$PLH[group_inter1], shuffle_uFDis_species_based_inter$PLH[group_inter1]),
                                       rmse(shuffle_uFDis_obs$LNC[group_inter2], shuffle_uFDis_species_based_inter$LNC[group_inter2]),
                                       rmse(shuffle_uFDis_obs$SLA[group_inter2], shuffle_uFDis_species_based_inter$SLA[group_inter2]),
                                       rmse(shuffle_uFDis_obs$PLH[group_inter2], shuffle_uFDis_species_based_inter$PLH[group_inter2]),
                                       rmse(shuffle_uFDis_obs$LNC[group_inter3], shuffle_uFDis_species_based_inter$LNC[group_inter3]),
                                       rmse(shuffle_uFDis_obs$SLA[group_inter3], shuffle_uFDis_species_based_inter$SLA[group_inter3]),
                                       rmse(shuffle_uFDis_obs$PLH[group_inter3], shuffle_uFDis_species_based_inter$PLH[group_inter3]),
                                       rmse(shuffle_uFDis_obs$LNC[group_inter4], shuffle_uFDis_species_based_inter$LNC[group_inter4]),
                                       rmse(shuffle_uFDis_obs$SLA[group_inter4], shuffle_uFDis_species_based_inter$SLA[group_inter4]),
                                       rmse(shuffle_uFDis_obs$PLH[group_inter4], shuffle_uFDis_species_based_inter$PLH[group_inter4]),
                                       rmse(shuffle_uFDis_obs$LNC[group_extra1], shuffle_uFDis_species_based_extra$LNC[group_extra1]), # uuFDis species based approach extrapolation
                                       rmse(shuffle_uFDis_obs$SLA[group_extra1], shuffle_uFDis_species_based_extra$SLA[group_extra1]),
                                       rmse(shuffle_uFDis_obs$PLH[group_extra1], shuffle_uFDis_species_based_extra$PLH[group_extra1]),
                                       rmse(shuffle_uFDis_obs$LNC[group_extra2], shuffle_uFDis_species_based_extra$LNC[group_extra2]),
                                       rmse(shuffle_uFDis_obs$SLA[group_extra2], shuffle_uFDis_species_based_extra$SLA[group_extra2]),
                                       rmse(shuffle_uFDis_obs$PLH[group_extra2], shuffle_uFDis_species_based_extra$PLH[group_extra2]),
                                       rmse(shuffle_uFDis_obs$LNC[group_extra3], shuffle_uFDis_species_based_extra$LNC[group_extra3]),
                                       rmse(shuffle_uFDis_obs$SLA[group_extra3], shuffle_uFDis_species_based_extra$SLA[group_extra3]),
                                       rmse(shuffle_uFDis_obs$PLH[group_extra3], shuffle_uFDis_species_based_extra$PLH[group_extra3]),
                                       rmse(shuffle_uFDis_obs$LNC[group_extra4], shuffle_uFDis_species_based_extra$LNC[group_extra4]),
                                       rmse(shuffle_uFDis_obs$SLA[group_extra4], shuffle_uFDis_species_based_extra$SLA[group_extra4]),
                                       rmse(shuffle_uFDis_obs$PLH[group_extra4], shuffle_uFDis_species_based_extra$PLH[group_extra4])),
                              "R2"=c(cor(shuffle_CWM_obs$LNC[group_inter1], shuffle_CWM_trait_based_inter$LNC[group_inter1])^2, # CWM trait based approach interpolation
                                     cor(shuffle_CWM_obs$SLA[group_inter1], shuffle_CWM_trait_based_inter$SLA[group_inter1])^2,
                                     cor(shuffle_CWM_obs$PLH[group_inter1], shuffle_CWM_trait_based_inter$PLH[group_inter1])^2,
                                     cor(shuffle_CWM_obs$LNC[group_inter2], shuffle_CWM_trait_based_inter$LNC[group_inter2])^2,
                                     cor(shuffle_CWM_obs$SLA[group_inter2], shuffle_CWM_trait_based_inter$SLA[group_inter2])^2,
                                     cor(shuffle_CWM_obs$PLH[group_inter2], shuffle_CWM_trait_based_inter$PLH[group_inter2])^2,
                                     cor(shuffle_CWM_obs$LNC[group_inter3], shuffle_CWM_trait_based_inter$LNC[group_inter3])^2,
                                     cor(shuffle_CWM_obs$SLA[group_inter3], shuffle_CWM_trait_based_inter$SLA[group_inter3])^2,
                                     cor(shuffle_CWM_obs$PLH[group_inter3], shuffle_CWM_trait_based_inter$PLH[group_inter3])^2,
                                     cor(shuffle_CWM_obs$LNC[group_inter4], shuffle_CWM_trait_based_inter$LNC[group_inter4])^2,
                                     cor(shuffle_CWM_obs$SLA[group_inter4], shuffle_CWM_trait_based_inter$SLA[group_inter4])^2,
                                     cor(shuffle_CWM_obs$PLH[group_inter4], shuffle_CWM_trait_based_inter$PLH[group_inter4])^2,
                                     cor(shuffle_CWM_obs$LNC[group_extra1], shuffle_CWM_trait_based_extra$LNC[group_extra1])^2,# CWM trait based approach interpolation
                                     cor(shuffle_CWM_obs$SLA[group_extra1], shuffle_CWM_trait_based_extra$SLA[group_extra1])^2,
                                     cor(shuffle_CWM_obs$PLH[group_extra1], shuffle_CWM_trait_based_extra$PLH[group_extra1])^2,
                                     cor(shuffle_CWM_obs$LNC[group_extra2], shuffle_CWM_trait_based_extra$LNC[group_extra2])^2,
                                     cor(shuffle_CWM_obs$SLA[group_extra2], shuffle_CWM_trait_based_extra$SLA[group_extra2])^2,
                                     cor(shuffle_CWM_obs$PLH[group_extra2], shuffle_CWM_trait_based_extra$PLH[group_extra2])^2,
                                     cor(shuffle_CWM_obs$LNC[group_extra3], shuffle_CWM_trait_based_extra$LNC[group_extra3])^2,
                                     cor(shuffle_CWM_obs$SLA[group_extra3], shuffle_CWM_trait_based_extra$SLA[group_extra3])^2,
                                     cor(shuffle_CWM_obs$PLH[group_extra3], shuffle_CWM_trait_based_extra$PLH[group_extra3])^2,
                                     cor(shuffle_CWM_obs$LNC[group_extra4], shuffle_CWM_trait_based_extra$LNC[group_extra4])^2,
                                     cor(shuffle_CWM_obs$SLA[group_extra4], shuffle_CWM_trait_based_extra$SLA[group_extra4])^2,
                                     cor(shuffle_CWM_obs$PLH[group_extra4], shuffle_CWM_trait_based_extra$PLH[group_extra4])^2,
                                     cor(shuffle_FDis_obs$LNC[group_inter1], shuffle_FDis_trait_based_inter$LNC[group_inter1])^2,# FDis trait based approach interpolation
                                     cor(shuffle_FDis_obs$SLA[group_inter1], shuffle_FDis_trait_based_inter$SLA[group_inter1])^2,
                                     cor(shuffle_FDis_obs$PLH[group_inter1], shuffle_FDis_trait_based_inter$PLH[group_inter1])^2,
                                     cor(shuffle_FDis_obs$LNC[group_inter2], shuffle_FDis_trait_based_inter$LNC[group_inter2])^2,
                                     cor(shuffle_FDis_obs$SLA[group_inter2], shuffle_FDis_trait_based_inter$SLA[group_inter2])^2,
                                     cor(shuffle_FDis_obs$PLH[group_inter2], shuffle_FDis_trait_based_inter$PLH[group_inter2])^2,
                                     cor(shuffle_FDis_obs$LNC[group_inter3], shuffle_FDis_trait_based_inter$LNC[group_inter3])^2,
                                     cor(shuffle_FDis_obs$SLA[group_inter3], shuffle_FDis_trait_based_inter$SLA[group_inter3])^2,
                                     cor(shuffle_FDis_obs$PLH[group_inter3], shuffle_FDis_trait_based_inter$PLH[group_inter3])^2,
                                     cor(shuffle_FDis_obs$LNC[group_inter4], shuffle_FDis_trait_based_inter$LNC[group_inter4])^2,
                                     cor(shuffle_FDis_obs$SLA[group_inter4], shuffle_FDis_trait_based_inter$SLA[group_inter4])^2,
                                     cor(shuffle_FDis_obs$PLH[group_inter4], shuffle_FDis_trait_based_inter$PLH[group_inter4])^2,
                                     cor(shuffle_FDis_obs$LNC[group_extra1], shuffle_FDis_trait_based_extra$LNC[group_extra1])^2, # FDis trait based approach extrapolation
                                     cor(shuffle_FDis_obs$SLA[group_extra1], shuffle_FDis_trait_based_extra$SLA[group_extra1])^2,
                                     cor(shuffle_FDis_obs$PLH[group_extra1], shuffle_FDis_trait_based_extra$PLH[group_extra1])^2,
                                     cor(shuffle_FDis_obs$LNC[group_extra2], shuffle_FDis_trait_based_extra$LNC[group_extra2])^2,
                                     cor(shuffle_FDis_obs$SLA[group_extra2], shuffle_FDis_trait_based_extra$SLA[group_extra2])^2,
                                     cor(shuffle_FDis_obs$PLH[group_extra2], shuffle_FDis_trait_based_extra$PLH[group_extra2])^2,
                                     cor(shuffle_FDis_obs$LNC[group_extra3], shuffle_FDis_trait_based_extra$LNC[group_extra3])^2,
                                     cor(shuffle_FDis_obs$SLA[group_extra3], shuffle_FDis_trait_based_extra$SLA[group_extra3])^2,
                                     cor(shuffle_FDis_obs$PLH[group_extra3], shuffle_FDis_trait_based_extra$PLH[group_extra3])^2,
                                     cor(shuffle_FDis_obs$LNC[group_extra4], shuffle_FDis_trait_based_extra$LNC[group_extra4])^2,
                                     cor(shuffle_FDis_obs$SLA[group_extra4], shuffle_FDis_trait_based_extra$SLA[group_extra4])^2,
                                     cor(shuffle_FDis_obs$PLH[group_extra4], shuffle_FDis_trait_based_extra$PLH[group_extra4])^2,
                                     cor(shuffle_CWM_obs$LNC[group_inter1], shuffle_CWM_species_based_inter$LNC[group_inter1])^2, # CWM species based approach interpolation
                                     cor(shuffle_CWM_obs$SLA[group_inter1], shuffle_CWM_species_based_inter$SLA[group_inter1])^2,
                                     cor(shuffle_CWM_obs$PLH[group_inter1], shuffle_CWM_species_based_inter$PLH[group_inter1])^2,
                                     cor(shuffle_CWM_obs$LNC[group_inter2], shuffle_CWM_species_based_inter$LNC[group_inter2])^2,
                                     cor(shuffle_CWM_obs$SLA[group_inter2], shuffle_CWM_species_based_inter$SLA[group_inter2])^2,
                                     cor(shuffle_CWM_obs$PLH[group_inter2], shuffle_CWM_species_based_inter$PLH[group_inter2])^2,
                                     cor(shuffle_CWM_obs$LNC[group_inter3], shuffle_CWM_species_based_inter$LNC[group_inter3])^2,
                                     cor(shuffle_CWM_obs$SLA[group_inter3], shuffle_CWM_species_based_inter$SLA[group_inter3])^2,
                                     cor(shuffle_CWM_obs$PLH[group_inter3], shuffle_CWM_species_based_inter$PLH[group_inter3])^2,
                                     cor(shuffle_CWM_obs$LNC[group_inter4], shuffle_CWM_species_based_inter$LNC[group_inter4])^2,
                                     cor(shuffle_CWM_obs$SLA[group_inter4], shuffle_CWM_species_based_inter$SLA[group_inter4])^2,
                                     cor(shuffle_CWM_obs$PLH[group_inter4], shuffle_CWM_species_based_inter$PLH[group_inter4])^2,
                                     cor(shuffle_CWM_obs$LNC[group_extra1], shuffle_CWM_species_based_extra$LNC[group_extra1])^2, # CWM species based approach extrapolation
                                     cor(shuffle_CWM_obs$SLA[group_extra1], shuffle_CWM_species_based_extra$SLA[group_extra1])^2,
                                     cor(shuffle_CWM_obs$PLH[group_extra1], shuffle_CWM_species_based_extra$PLH[group_extra1])^2,
                                     cor(shuffle_CWM_obs$LNC[group_extra2], shuffle_CWM_species_based_extra$LNC[group_extra2])^2,
                                     cor(shuffle_CWM_obs$SLA[group_extra2], shuffle_CWM_species_based_extra$SLA[group_extra2])^2,
                                     cor(shuffle_CWM_obs$PLH[group_extra2], shuffle_CWM_species_based_extra$PLH[group_extra2])^2,
                                     cor(shuffle_CWM_obs$LNC[group_extra3], shuffle_CWM_species_based_extra$LNC[group_extra3])^2,
                                     cor(shuffle_CWM_obs$SLA[group_extra3], shuffle_CWM_species_based_extra$SLA[group_extra3])^2,
                                     cor(shuffle_CWM_obs$PLH[group_extra3], shuffle_CWM_species_based_extra$PLH[group_extra3])^2,
                                     cor(shuffle_CWM_obs$LNC[group_extra4], shuffle_CWM_species_based_extra$LNC[group_extra4])^2,
                                     cor(shuffle_CWM_obs$SLA[group_extra4], shuffle_CWM_species_based_extra$SLA[group_extra4])^2,
                                     cor(shuffle_CWM_obs$PLH[group_extra4], shuffle_CWM_species_based_extra$PLH[group_extra4])^2,
                                     cor(shuffle_FDis_obs$LNC[group_inter1], shuffle_FDis_species_based_inter$LNC[group_inter1])^2, # FDis species based approach interpolation
                                     cor(shuffle_FDis_obs$SLA[group_inter1], shuffle_FDis_species_based_inter$SLA[group_inter1])^2,
                                     cor(shuffle_FDis_obs$PLH[group_inter1], shuffle_FDis_species_based_inter$PLH[group_inter1])^2,
                                     cor(shuffle_FDis_obs$LNC[group_inter2], shuffle_FDis_species_based_inter$LNC[group_inter2])^2,
                                     cor(shuffle_FDis_obs$SLA[group_inter2], shuffle_FDis_species_based_inter$SLA[group_inter2])^2,
                                     cor(shuffle_FDis_obs$PLH[group_inter2], shuffle_FDis_species_based_inter$PLH[group_inter2])^2,
                                     cor(shuffle_FDis_obs$LNC[group_inter3], shuffle_FDis_species_based_inter$LNC[group_inter3])^2,
                                     cor(shuffle_FDis_obs$SLA[group_inter3], shuffle_FDis_species_based_inter$SLA[group_inter3])^2,
                                     cor(shuffle_FDis_obs$PLH[group_inter3], shuffle_FDis_species_based_inter$PLH[group_inter3])^2,
                                     cor(shuffle_FDis_obs$LNC[group_inter4], shuffle_FDis_species_based_inter$LNC[group_inter4])^2,
                                     cor(shuffle_FDis_obs$SLA[group_inter4], shuffle_FDis_species_based_inter$SLA[group_inter4])^2,
                                     cor(shuffle_FDis_obs$PLH[group_inter4], shuffle_FDis_species_based_inter$PLH[group_inter4])^2,
                                     cor(shuffle_FDis_obs$LNC[group_extra1], shuffle_FDis_species_based_extra$LNC[group_extra1])^2, # FDis species based approach extrapolation
                                     cor(shuffle_FDis_obs$SLA[group_extra1], shuffle_FDis_species_based_extra$SLA[group_extra1])^2,
                                     cor(shuffle_FDis_obs$PLH[group_extra1], shuffle_FDis_species_based_extra$PLH[group_extra1])^2,
                                     cor(shuffle_FDis_obs$LNC[group_extra2], shuffle_FDis_species_based_extra$LNC[group_extra2])^2,
                                     cor(shuffle_FDis_obs$SLA[group_extra2], shuffle_FDis_species_based_extra$SLA[group_extra2])^2,
                                     cor(shuffle_FDis_obs$PLH[group_extra2], shuffle_FDis_species_based_extra$PLH[group_extra2])^2,
                                     cor(shuffle_FDis_obs$LNC[group_extra3], shuffle_FDis_species_based_extra$LNC[group_extra3])^2,
                                     cor(shuffle_FDis_obs$SLA[group_extra3], shuffle_FDis_species_based_extra$SLA[group_extra3])^2,
                                     cor(shuffle_FDis_obs$PLH[group_extra3], shuffle_FDis_species_based_extra$PLH[group_extra3])^2,
                                     cor(shuffle_FDis_obs$LNC[group_extra4], shuffle_FDis_species_based_extra$LNC[group_extra4])^2,
                                     cor(shuffle_FDis_obs$SLA[group_extra4], shuffle_FDis_species_based_extra$SLA[group_extra4])^2,
                                     cor(shuffle_FDis_obs$PLH[group_extra4], shuffle_FDis_species_based_extra$PLH[group_extra4])^2,
                                     
                                     cor(shuffle_CM_obs$LNC[group_inter1], shuffle_CM_trait_based_inter$LNC[group_inter1])^2, # CM trait based approach interpolation
                                     cor(shuffle_CM_obs$SLA[group_inter1], shuffle_CM_trait_based_inter$SLA[group_inter1])^2,
                                     cor(shuffle_CM_obs$PLH[group_inter1], shuffle_CM_trait_based_inter$PLH[group_inter1])^2,
                                     cor(shuffle_CM_obs$LNC[group_inter2], shuffle_CM_trait_based_inter$LNC[group_inter2])^2,
                                     cor(shuffle_CM_obs$SLA[group_inter2], shuffle_CM_trait_based_inter$SLA[group_inter2])^2,
                                     cor(shuffle_CM_obs$PLH[group_inter2], shuffle_CM_trait_based_inter$PLH[group_inter2])^2,
                                     cor(shuffle_CM_obs$LNC[group_inter3], shuffle_CM_trait_based_inter$LNC[group_inter3])^2,
                                     cor(shuffle_CM_obs$SLA[group_inter3], shuffle_CM_trait_based_inter$SLA[group_inter3])^2,
                                     cor(shuffle_CM_obs$PLH[group_inter3], shuffle_CM_trait_based_inter$PLH[group_inter3])^2,
                                     cor(shuffle_CM_obs$LNC[group_inter4], shuffle_CM_trait_based_inter$LNC[group_inter4])^2,
                                     cor(shuffle_CM_obs$SLA[group_inter4], shuffle_CM_trait_based_inter$SLA[group_inter4])^2,
                                     cor(shuffle_CM_obs$PLH[group_inter4], shuffle_CM_trait_based_inter$PLH[group_inter4])^2,
                                     cor(shuffle_CM_obs$LNC[group_extra1], shuffle_CM_trait_based_extra$LNC[group_extra1])^2,# CM trait based approach interpolation
                                     cor(shuffle_CM_obs$SLA[group_extra1], shuffle_CM_trait_based_extra$SLA[group_extra1])^2,
                                     cor(shuffle_CM_obs$PLH[group_extra1], shuffle_CM_trait_based_extra$PLH[group_extra1])^2,
                                     cor(shuffle_CM_obs$LNC[group_extra2], shuffle_CM_trait_based_extra$LNC[group_extra2])^2,
                                     cor(shuffle_CM_obs$SLA[group_extra2], shuffle_CM_trait_based_extra$SLA[group_extra2])^2,
                                     cor(shuffle_CM_obs$PLH[group_extra2], shuffle_CM_trait_based_extra$PLH[group_extra2])^2,
                                     cor(shuffle_CM_obs$LNC[group_extra3], shuffle_CM_trait_based_extra$LNC[group_extra3])^2,
                                     cor(shuffle_CM_obs$SLA[group_extra3], shuffle_CM_trait_based_extra$SLA[group_extra3])^2,
                                     cor(shuffle_CM_obs$PLH[group_extra3], shuffle_CM_trait_based_extra$PLH[group_extra3])^2,
                                     cor(shuffle_CM_obs$LNC[group_extra4], shuffle_CM_trait_based_extra$LNC[group_extra4])^2,
                                     cor(shuffle_CM_obs$SLA[group_extra4], shuffle_CM_trait_based_extra$SLA[group_extra4])^2,
                                     cor(shuffle_CM_obs$PLH[group_extra4], shuffle_CM_trait_based_extra$PLH[group_extra4])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_inter1], shuffle_uFDis_trait_based_inter$LNC[group_inter1])^2,# uuFDis trait based approach interpolation
                                     cor(shuffle_uFDis_obs$SLA[group_inter1], shuffle_uFDis_trait_based_inter$SLA[group_inter1])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_inter1], shuffle_uFDis_trait_based_inter$PLH[group_inter1])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_inter2], shuffle_uFDis_trait_based_inter$LNC[group_inter2])^2,
                                     cor(shuffle_uFDis_obs$SLA[group_inter2], shuffle_uFDis_trait_based_inter$SLA[group_inter2])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_inter2], shuffle_uFDis_trait_based_inter$PLH[group_inter2])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_inter3], shuffle_uFDis_trait_based_inter$LNC[group_inter3])^2,
                                     cor(shuffle_uFDis_obs$SLA[group_inter3], shuffle_uFDis_trait_based_inter$SLA[group_inter3])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_inter3], shuffle_uFDis_trait_based_inter$PLH[group_inter3])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_inter4], shuffle_uFDis_trait_based_inter$LNC[group_inter4])^2,
                                     cor(shuffle_uFDis_obs$SLA[group_inter4], shuffle_uFDis_trait_based_inter$SLA[group_inter4])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_inter4], shuffle_uFDis_trait_based_inter$PLH[group_inter4])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_extra1], shuffle_uFDis_trait_based_extra$LNC[group_extra1])^2, # uuFDis trait based approach extrapolation
                                     cor(shuffle_uFDis_obs$SLA[group_extra1], shuffle_uFDis_trait_based_extra$SLA[group_extra1])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_extra1], shuffle_uFDis_trait_based_extra$PLH[group_extra1])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_extra2], shuffle_uFDis_trait_based_extra$LNC[group_extra2])^2,
                                     cor(shuffle_uFDis_obs$SLA[group_extra2], shuffle_uFDis_trait_based_extra$SLA[group_extra2])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_extra2], shuffle_uFDis_trait_based_extra$PLH[group_extra2])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_extra3], shuffle_uFDis_trait_based_extra$LNC[group_extra3])^2,
                                     cor(shuffle_uFDis_obs$SLA[group_extra3], shuffle_uFDis_trait_based_extra$SLA[group_extra3])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_extra3], shuffle_uFDis_trait_based_extra$PLH[group_extra3])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_extra4], shuffle_uFDis_trait_based_extra$LNC[group_extra4])^2,
                                     cor(shuffle_uFDis_obs$SLA[group_extra4], shuffle_uFDis_trait_based_extra$SLA[group_extra4])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_extra4], shuffle_uFDis_trait_based_extra$PLH[group_extra4])^2,
                                     cor(shuffle_CM_obs$LNC[group_inter1], shuffle_CM_species_based_inter$LNC[group_inter1])^2, # CM species based approach interpolation
                                     cor(shuffle_CM_obs$SLA[group_inter1], shuffle_CM_species_based_inter$SLA[group_inter1])^2,
                                     cor(shuffle_CM_obs$PLH[group_inter1], shuffle_CM_species_based_inter$PLH[group_inter1])^2,
                                     cor(shuffle_CM_obs$LNC[group_inter2], shuffle_CM_species_based_inter$LNC[group_inter2])^2,
                                     cor(shuffle_CM_obs$SLA[group_inter2], shuffle_CM_species_based_inter$SLA[group_inter2])^2,
                                     cor(shuffle_CM_obs$PLH[group_inter2], shuffle_CM_species_based_inter$PLH[group_inter2])^2,
                                     cor(shuffle_CM_obs$LNC[group_inter3], shuffle_CM_species_based_inter$LNC[group_inter3])^2,
                                     cor(shuffle_CM_obs$SLA[group_inter3], shuffle_CM_species_based_inter$SLA[group_inter3])^2,
                                     cor(shuffle_CM_obs$PLH[group_inter3], shuffle_CM_species_based_inter$PLH[group_inter3])^2,
                                     cor(shuffle_CM_obs$LNC[group_inter4], shuffle_CM_species_based_inter$LNC[group_inter4])^2,
                                     cor(shuffle_CM_obs$SLA[group_inter4], shuffle_CM_species_based_inter$SLA[group_inter4])^2,
                                     cor(shuffle_CM_obs$PLH[group_inter4], shuffle_CM_species_based_inter$PLH[group_inter4])^2,
                                     cor(shuffle_CM_obs$LNC[group_extra1], shuffle_CM_species_based_extra$LNC[group_extra1])^2, # CM species based approach extrapolation
                                     cor(shuffle_CM_obs$SLA[group_extra1], shuffle_CM_species_based_extra$SLA[group_extra1])^2,
                                     cor(shuffle_CM_obs$PLH[group_extra1], shuffle_CM_species_based_extra$PLH[group_extra1])^2,
                                     cor(shuffle_CM_obs$LNC[group_extra2], shuffle_CM_species_based_extra$LNC[group_extra2])^2,
                                     cor(shuffle_CM_obs$SLA[group_extra2], shuffle_CM_species_based_extra$SLA[group_extra2])^2,
                                     cor(shuffle_CM_obs$PLH[group_extra2], shuffle_CM_species_based_extra$PLH[group_extra2])^2,
                                     cor(shuffle_CM_obs$LNC[group_extra3], shuffle_CM_species_based_extra$LNC[group_extra3])^2,
                                     cor(shuffle_CM_obs$SLA[group_extra3], shuffle_CM_species_based_extra$SLA[group_extra3])^2,
                                     cor(shuffle_CM_obs$PLH[group_extra3], shuffle_CM_species_based_extra$PLH[group_extra3])^2,
                                     cor(shuffle_CM_obs$LNC[group_extra4], shuffle_CM_species_based_extra$LNC[group_extra4])^2,
                                     cor(shuffle_CM_obs$SLA[group_extra4], shuffle_CM_species_based_extra$SLA[group_extra4])^2,
                                     cor(shuffle_CM_obs$PLH[group_extra4], shuffle_CM_species_based_extra$PLH[group_extra4])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_inter1], shuffle_uFDis_species_based_inter$LNC[group_inter1])^2, # uuFDis species based approach interpolation
                                     cor(shuffle_uFDis_obs$SLA[group_inter1], shuffle_uFDis_species_based_inter$SLA[group_inter1])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_inter1], shuffle_uFDis_species_based_inter$PLH[group_inter1])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_inter2], shuffle_uFDis_species_based_inter$LNC[group_inter2])^2,
                                     cor(shuffle_uFDis_obs$SLA[group_inter2], shuffle_uFDis_species_based_inter$SLA[group_inter2])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_inter2], shuffle_uFDis_species_based_inter$PLH[group_inter2])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_inter3], shuffle_uFDis_species_based_inter$LNC[group_inter3])^2,
                                     cor(shuffle_uFDis_obs$SLA[group_inter3], shuffle_uFDis_species_based_inter$SLA[group_inter3])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_inter3], shuffle_uFDis_species_based_inter$PLH[group_inter3])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_inter4], shuffle_uFDis_species_based_inter$LNC[group_inter4])^2,
                                     cor(shuffle_uFDis_obs$SLA[group_inter4], shuffle_uFDis_species_based_inter$SLA[group_inter4])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_inter4], shuffle_uFDis_species_based_inter$PLH[group_inter4])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_extra1], shuffle_uFDis_species_based_extra$LNC[group_extra1])^2, # uuFDis species based approach extrapolation
                                     cor(shuffle_uFDis_obs$SLA[group_extra1], shuffle_uFDis_species_based_extra$SLA[group_extra1])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_extra1], shuffle_uFDis_species_based_extra$PLH[group_extra1])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_extra2], shuffle_uFDis_species_based_extra$LNC[group_extra2])^2,
                                     cor(shuffle_uFDis_obs$SLA[group_extra2], shuffle_uFDis_species_based_extra$SLA[group_extra2])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_extra2], shuffle_uFDis_species_based_extra$PLH[group_extra2])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_extra3], shuffle_uFDis_species_based_extra$LNC[group_extra3])^2,
                                     cor(shuffle_uFDis_obs$SLA[group_extra3], shuffle_uFDis_species_based_extra$SLA[group_extra3])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_extra3], shuffle_uFDis_species_based_extra$PLH[group_extra3])^2,
                                     cor(shuffle_uFDis_obs$LNC[group_extra4], shuffle_uFDis_species_based_extra$LNC[group_extra4])^2,
                                     cor(shuffle_uFDis_obs$SLA[group_extra4], shuffle_uFDis_species_based_extra$SLA[group_extra4])^2,
                                     cor(shuffle_uFDis_obs$PLH[group_extra4], shuffle_uFDis_species_based_extra$PLH[group_extra4])^2),
                              "Traits"=rep(c(rep(c("LNC", "SLA", "PLH"), 4*8)),2),
                              "Indices"=c(c(rep(c(rep("CWM", 24), rep("FDis", 24)),2)),c(rep(c(rep("CM", 24), rep("uFDis", 24)),2))),
                              "Approaches"=rep(c(rep("Assemble-First Approach", 12*4), rep("Predict-First Approach", 12*4)),2),
                              "Crossvalidation"=rep(c(rep(c(rep("Interpolation", 12), rep("Extrapolation", 12)),4)),2),
                              "Type"=c(rep("Abundance data", 4*8*3), rep("Occurrence data", 4*8*3)),
                              "repet"=rep(c(rep(1,3), rep(2,3), rep(3,3), rep(4,3)), 16),
                              "alpha"=rep(0.1, 3*4*8*2))

shuffle_figure3$Traits <- fct_relevel(shuffle_figure3$Traits, c("PLH", "SLA", "LNC"))
shuffle_figure3$Crossvalidation <- fct_relevel(shuffle_figure3$Crossvalidation, c("Interpolation", "Extrapolation"))
shuffle_figure3$Type <- fct_relevel(shuffle_figure3$Type, c("Abundance data","Occurrence data"))
shuffle_figure3$Approaches <- fct_relevel(shuffle_figure3$Approaches, c("Predict-First Approach","Assemble-First Approach"))
shuffle_figure3$Indices <- fct_relevel(shuffle_figure3$Indices, c("CWM","FDis", "CM", 'uFDis'))


for (i in 1:9){
  shuffle_traits <- traits[sample(1:nrow(traits)),]
  shuffle_traits$Taxa <- traits$Taxa
  shuffle_CWM_obs <- CWM_shuffle(plant_recovery)
  while(length(which(is.na(shuffle_CWM_obs)))>0){
    shuffle_traits <- traits[sample(1:nrow(traits)),]
    shuffle_traits$Taxa <- traits$Taxa
    shuffle_CWM_obs <- CWM_shuffle(plant_recovery)
  }
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_CWM_obs <- left_join(Site_order, shuffle_CWM_obs, by="Site")
  rownames(shuffle_CWM_obs)=shuffle_CWM_obs$Site
  shuffle_FDis_obs <- FDis_shuffle(plant_recovery)
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_FDis_obs <- left_join(Site_order, shuffle_FDis_obs, by="Site")
  rownames(shuffle_FDis_obs)=shuffle_FDis_obs$Site
  
  shuffle_CM_obs <- CWM_shuffle(plant_pa)
  Site_order <- data.frame("Site"=rownames(plant_pa))
  shuffle_CM_obs <- left_join(Site_order, shuffle_CM_obs, by="Site")
  rownames(shuffle_CM_obs)=shuffle_CM_obs$Site
  shuffle_uFDis_obs <- FDis_shuffle(plant_pa)
  Site_order <- data.frame("Site"=rownames(plant_pa))
  shuffle_uFDis_obs <- left_join(Site_order, shuffle_uFDis_obs, by="Site")
  rownames(shuffle_uFDis_obs)=shuffle_uFDis_obs$Site
  
  shuffle_CWM_species_based_inter <- CWM_shuffle(species_models_recovery_inter)
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_CWM_species_based_inter <- left_join(Site_order, shuffle_CWM_species_based_inter, by="Site")
  rownames(shuffle_CWM_species_based_inter)=shuffle_CWM_species_based_inter$Site
  shuffle_CWM_species_based_inter <- shuffle_CWM_species_based_inter[,-1]
  
  shuffle_CWM_species_based_extra <- CWM_shuffle(species_models_recovery_extra)
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_CWM_species_based_extra <- left_join(Site_order, shuffle_CWM_species_based_extra, by="Site")
  rownames(shuffle_CWM_species_based_extra)=shuffle_CWM_species_based_extra$Site
  shuffle_CWM_species_based_extra <- shuffle_CWM_species_based_extra[,-1]
  
  shuffle_CM_species_based_inter <- CWM_shuffle(species_models_pa_inter)
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_CM_species_based_inter <- left_join(Site_order, shuffle_CM_species_based_inter, by="Site")
  rownames(shuffle_CM_species_based_inter)=shuffle_CM_species_based_inter$Site
  shuffle_CM_species_based_inter <- shuffle_CM_species_based_inter[,-1]
  
  shuffle_CM_species_based_extra <- CWM_shuffle(species_models_pa_extra)
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_CM_species_based_extra <- left_join(Site_order, shuffle_CM_species_based_extra, by="Site")
  rownames(shuffle_CM_species_based_extra)=shuffle_CM_species_based_extra$Site
  shuffle_CM_species_based_extra <- shuffle_CM_species_based_extra[,-1]
  
  shuffle_FDis_species_based_inter <- FDis_shuffle(species_models_recovery_inter)
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_FDis_species_based_inter <- left_join(Site_order, shuffle_FDis_species_based_inter, by="Site")
  rownames(shuffle_FDis_species_based_inter)=shuffle_FDis_species_based_inter$Site
  shuffle_FDis_species_based_inter <- shuffle_FDis_species_based_inter[,-1]
  
  shuffle_FDis_species_based_extra <- FDis_shuffle(species_models_recovery_extra)
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_FDis_species_based_extra <- left_join(Site_order, shuffle_FDis_species_based_extra, by="Site")
  rownames(shuffle_FDis_species_based_extra)=shuffle_FDis_species_based_extra$Site
  shuffle_FDis_species_based_extra <- shuffle_FDis_species_based_extra[,-1]
  
  shuffle_uFDis_species_based_inter <- FDis_shuffle(species_models_pa_inter)
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_uFDis_species_based_inter <- left_join(Site_order, shuffle_uFDis_species_based_inter, by="Site")
  rownames(shuffle_uFDis_species_based_inter)=shuffle_uFDis_species_based_inter$Site
  shuffle_uFDis_species_based_inter <- shuffle_uFDis_species_based_inter[,-1]
  
  shuffle_uFDis_species_based_extra <- FDis_shuffle(species_models_pa_extra)
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_uFDis_species_based_extra <- left_join(Site_order, shuffle_uFDis_species_based_extra, by="Site")
  rownames(shuffle_uFDis_species_based_extra)=shuffle_uFDis_species_based_extra$Site
  shuffle_uFDis_species_based_extra <- shuffle_uFDis_species_based_extra[,-1]
  
  ##################trait_based 
  gam_data <- left_join(shuffle_CWM_obs,env_var, by="Site")
  ## - For interpolation group 1 - ##
  shuffle_models_CWM_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
  shuffle_models_CWM_LNC_inter1 <- predict(shuffle_models_CWM_LNC_grow1, gam_data[group_inter1,])$predicted
  shuffle_models_CWM_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
  shuffle_models_CWM_SLA_inter1 <- predict(shuffle_models_CWM_SLA_grow1, gam_data[group_inter1,])$predicted
  shuffle_models_CWM_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
  shuffle_models_CWM_PLH_inter1 <- predict(shuffle_models_CWM_PLH_grow1, gam_data[group_inter1,])$predicted
  
  shuffle_CWM_trait_based_inter1<-data.frame("LNC"=shuffle_models_CWM_LNC_inter1, "SLA"=shuffle_models_CWM_SLA_inter1, "PLH"=shuffle_models_CWM_PLH_inter1)
  rownames(shuffle_CWM_trait_based_inter1) <- shuffle_CWM_obs$Site[group_inter1]
  
  ## - For interpolation group 2 - ##
  shuffle_models_CWM_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
  shuffle_models_CWM_LNC_inter2 <- predict(shuffle_models_CWM_LNC_grow2, gam_data[group_inter2,])$predicted
  shuffle_models_CWM_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
  shuffle_models_CWM_SLA_inter2 <- predict(shuffle_models_CWM_SLA_grow2, gam_data[group_inter2,])$predicted
  shuffle_models_CWM_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
  shuffle_models_CWM_PLH_inter2 <- predict(shuffle_models_CWM_PLH_grow2, gam_data[group_inter2,])$predicted
  
  shuffle_CWM_trait_based_inter2<-data.frame("LNC"=shuffle_models_CWM_LNC_inter2, "SLA"=shuffle_models_CWM_SLA_inter2, "PLH"=shuffle_models_CWM_PLH_inter2)
  rownames(shuffle_CWM_trait_based_inter2) <- shuffle_CWM_obs$Site[group_inter2]
  
  ## - For interpolation group 3 - ##
  shuffle_models_CWM_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
  shuffle_models_CWM_LNC_inter3 <- predict(shuffle_models_CWM_LNC_grow3, gam_data[group_inter3,])$predicted
  shuffle_models_CWM_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
  shuffle_models_CWM_SLA_inter3 <- predict(shuffle_models_CWM_SLA_grow3, gam_data[group_inter3,])$predicted
  shuffle_models_CWM_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
  shuffle_models_CWM_PLH_inter3 <- predict(shuffle_models_CWM_PLH_grow3, gam_data[group_inter3,])$predicted
  
  shuffle_CWM_trait_based_inter3<-data.frame("LNC"=shuffle_models_CWM_LNC_inter3, "SLA"=shuffle_models_CWM_SLA_inter3, "PLH"=shuffle_models_CWM_PLH_inter3)
  rownames(shuffle_CWM_trait_based_inter3) <- shuffle_CWM_obs$Site[group_inter3]
  
  ## - For interpolation group 4 - ##
  shuffle_models_CWM_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
  shuffle_models_CWM_LNC_inter4 <- predict(shuffle_models_CWM_LNC_grow4, gam_data[group_inter4,])$predicted
  shuffle_models_CWM_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
  shuffle_models_CWM_SLA_inter4 <- predict(shuffle_models_CWM_SLA_grow4, gam_data[group_inter4,])$predicted
  shuffle_models_CWM_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
  shuffle_models_CWM_PLH_inter4 <- predict(shuffle_models_CWM_PLH_grow4, gam_data[group_inter4,])$predicted
  
  shuffle_CWM_trait_based_inter4<-data.frame("LNC"=shuffle_models_CWM_LNC_inter4, "SLA"=shuffle_models_CWM_SLA_inter4, "PLH"=shuffle_models_CWM_PLH_inter4)
  rownames(shuffle_CWM_trait_based_inter4) <- shuffle_CWM_obs$Site[group_inter4]
  
  ## assembly in a single prediction table for the 4 groups :
  shuffle_CWM_trait_based_inter <- rbind(shuffle_CWM_trait_based_inter1, shuffle_CWM_trait_based_inter2, shuffle_CWM_trait_based_inter3, shuffle_CWM_trait_based_inter4)
  shuffle_CWM_trait_based_inter <- cbind("Site"=rownames(shuffle_CWM_trait_based_inter), shuffle_CWM_trait_based_inter)
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_CWM_trait_based_inter <- left_join(Site_order, shuffle_CWM_trait_based_inter, by="Site")
  rownames(shuffle_CWM_trait_based_inter) <- shuffle_CWM_trait_based_inter$Site
  shuffle_CWM_trait_based_inter <- shuffle_CWM_trait_based_inter[,-1]
  
  #### -- This computes path predicts the community weighted mean in the context of extrapolation -- ####
  gam_data <- left_join(shuffle_CWM_obs,env_var, by="Site")
  ## - For extrapolation group 1 - ##
  shuffle_models_CWM_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
  shuffle_models_CWM_LNC_extra1 <- predict(shuffle_models_CWM_LNC_grow1, gam_data[group_extra1,])$predicted
  shuffle_models_CWM_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
  shuffle_models_CWM_SLA_extra1 <- predict(shuffle_models_CWM_SLA_grow1, gam_data[group_extra1,])$predicted
  shuffle_models_CWM_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
  shuffle_models_CWM_PLH_extra1 <- predict(shuffle_models_CWM_PLH_grow1, gam_data[group_extra1,])$predicted
  
  shuffle_CWM_trait_based_extra1<-data.frame("LNC"=shuffle_models_CWM_LNC_extra1, "SLA"=shuffle_models_CWM_SLA_extra1, "PLH"=shuffle_models_CWM_PLH_extra1)
  rownames(shuffle_CWM_trait_based_extra1) <- shuffle_CWM_obs$Site[group_extra1]
  
  ## - For extrapolation group 2 - ##
  shuffle_models_CWM_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
  shuffle_models_CWM_LNC_extra2 <- predict(shuffle_models_CWM_LNC_grow2, gam_data[group_extra2,])$predicted
  shuffle_models_CWM_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
  shuffle_models_CWM_SLA_extra2 <- predict(shuffle_models_CWM_SLA_grow2, gam_data[group_extra2,])$predicted
  shuffle_models_CWM_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
  shuffle_models_CWM_PLH_extra2 <- predict(shuffle_models_CWM_PLH_grow2, gam_data[group_extra2,])$predicted
  
  shuffle_CWM_trait_based_extra2<-data.frame("LNC"=shuffle_models_CWM_LNC_extra2, "SLA"=shuffle_models_CWM_SLA_extra2, "PLH"=shuffle_models_CWM_PLH_extra2)
  rownames(shuffle_CWM_trait_based_extra2) <- shuffle_CWM_obs$Site[group_extra2]
  
  ## - For extrapolation group 3 - ##
  shuffle_models_CWM_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
  shuffle_models_CWM_LNC_extra3 <- predict(shuffle_models_CWM_LNC_grow3, gam_data[group_extra3,])$predicted
  shuffle_models_CWM_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
  shuffle_models_CWM_SLA_extra3 <- predict(shuffle_models_CWM_SLA_grow3, gam_data[group_extra3,])$predicted
  shuffle_models_CWM_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
  shuffle_models_CWM_PLH_extra3 <- predict(shuffle_models_CWM_PLH_grow3, gam_data[group_extra3,])$predicted
  
  shuffle_CWM_trait_based_extra3<-data.frame("LNC"=shuffle_models_CWM_LNC_extra3, "SLA"=shuffle_models_CWM_SLA_extra3, "PLH"=shuffle_models_CWM_PLH_extra3)
  rownames(shuffle_CWM_trait_based_extra3) <- shuffle_CWM_obs$Site[group_extra3]
  
  ## - For extrapolation group 4 - ##
  shuffle_models_CWM_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
  shuffle_models_CWM_LNC_extra4 <- predict(shuffle_models_CWM_LNC_grow4, gam_data[group_extra4,])$predicted
  shuffle_models_CWM_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
  shuffle_models_CWM_SLA_extra4 <- predict(shuffle_models_CWM_SLA_grow4, gam_data[group_extra4,])$predicted
  shuffle_models_CWM_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
  shuffle_models_CWM_PLH_extra4 <- predict(shuffle_models_CWM_PLH_grow4, gam_data[group_extra4,])$predicted
  
  shuffle_CWM_trait_based_extra4<-data.frame("LNC"=shuffle_models_CWM_LNC_extra4, "SLA"=shuffle_models_CWM_SLA_extra4, "PLH"=shuffle_models_CWM_PLH_extra4)
  rownames(shuffle_CWM_trait_based_extra4) <- shuffle_CWM_obs$Site[group_extra4]
  
  ## assembly in a single prediction table for the 4 groups :
  shuffle_CWM_trait_based_extra <- rbind(shuffle_CWM_trait_based_extra1, shuffle_CWM_trait_based_extra2, shuffle_CWM_trait_based_extra3, shuffle_CWM_trait_based_extra4)
  shuffle_CWM_trait_based_extra <- cbind("Site"=rownames(shuffle_CWM_trait_based_extra), shuffle_CWM_trait_based_extra)
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_CWM_trait_based_extra <- left_join(Site_order, shuffle_CWM_trait_based_extra, by="Site")
  rownames(shuffle_CWM_trait_based_extra) <- shuffle_CWM_trait_based_extra$Site
  shuffle_CWM_trait_based_extra <- shuffle_CWM_trait_based_extra[,-1]
  
  ###### PREDICTION OF COMMUNITY MEAN  ######
  #### -- This computes path predicts the community mean in the context of interpolation -- ####
  gam_data <- left_join(shuffle_CM_obs,env_var, by="Site")
  ## - For interpolation group 1 - ##
  shuffle_models_CM_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
  shuffle_models_CM_LNC_inter1 <- predict(shuffle_models_CM_LNC_grow1, gam_data[group_inter1,])$predicted
  shuffle_models_CM_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
  shuffle_models_CM_SLA_inter1 <- predict(shuffle_models_CM_SLA_grow1, gam_data[group_inter1,])$predicted
  shuffle_models_CM_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
  shuffle_models_CM_PLH_inter1 <- predict(shuffle_models_CM_PLH_grow1, gam_data[group_inter1,])$predicted
  
  shuffle_CM_trait_based_inter1<-data.frame("LNC"=shuffle_models_CM_LNC_inter1, "SLA"=shuffle_models_CM_SLA_inter1, "PLH"=shuffle_models_CM_PLH_inter1)
  rownames(shuffle_CM_trait_based_inter1) <- shuffle_CM_obs$Site[group_inter1]
  
  ## - For interpolation group 2 - ##
  shuffle_models_CM_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
  shuffle_models_CM_LNC_inter2 <- predict(shuffle_models_CM_LNC_grow2, gam_data[group_inter2,])$predicted
  shuffle_models_CM_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
  shuffle_models_CM_SLA_inter2 <- predict(shuffle_models_CM_SLA_grow2, gam_data[group_inter2,])$predicted
  shuffle_models_CM_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
  shuffle_models_CM_PLH_inter2 <- predict(shuffle_models_CM_PLH_grow2, gam_data[group_inter2,])$predicted
  
  shuffle_CM_trait_based_inter2<-data.frame("LNC"=shuffle_models_CM_LNC_inter2, "SLA"=shuffle_models_CM_SLA_inter2, "PLH"=shuffle_models_CM_PLH_inter2)
  rownames(shuffle_CM_trait_based_inter2) <- shuffle_CM_obs$Site[group_inter2]
  
  ## - For interpolation group 3 - ##
  shuffle_models_CM_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
  shuffle_models_CM_LNC_inter3 <- predict(shuffle_models_CM_LNC_grow3, gam_data[group_inter3,])$predicted
  shuffle_models_CM_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
  shuffle_models_CM_SLA_inter3 <- predict(shuffle_models_CM_SLA_grow3, gam_data[group_inter3,])$predicted
  shuffle_models_CM_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
  shuffle_models_CM_PLH_inter3 <- predict(shuffle_models_CM_PLH_grow3, gam_data[group_inter3,])$predicted
  
  shuffle_CM_trait_based_inter3<-data.frame("LNC"=shuffle_models_CM_LNC_inter3, "SLA"=shuffle_models_CM_SLA_inter3, "PLH"=shuffle_models_CM_PLH_inter3)
  rownames(shuffle_CM_trait_based_inter3) <- shuffle_CM_obs$Site[group_inter3]
  
  ## - For interpolation group 4 - ##
  shuffle_models_CM_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
  shuffle_models_CM_LNC_inter4 <- predict(shuffle_models_CM_LNC_grow4, gam_data[group_inter4,])$predicted
  shuffle_models_CM_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
  shuffle_models_CM_SLA_inter4 <- predict(shuffle_models_CM_SLA_grow4, gam_data[group_inter4,])$predicted
  shuffle_models_CM_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
  shuffle_models_CM_PLH_inter4 <- predict(shuffle_models_CM_PLH_grow4, gam_data[group_inter4,])$predicted
  
  shuffle_CM_trait_based_inter4<-data.frame("LNC"=shuffle_models_CM_LNC_inter4, "SLA"=shuffle_models_CM_SLA_inter4, "PLH"=shuffle_models_CM_PLH_inter4)
  rownames(shuffle_CM_trait_based_inter4) <- shuffle_CM_obs$Site[group_inter4]
  
  ## assembly in a single prediction table for the 4 groups :
  shuffle_CM_trait_based_inter <- rbind(shuffle_CM_trait_based_inter1, shuffle_CM_trait_based_inter2, shuffle_CM_trait_based_inter3, shuffle_CM_trait_based_inter4)
  shuffle_CM_trait_based_inter <- cbind("Site"=rownames(shuffle_CM_trait_based_inter), shuffle_CM_trait_based_inter)
  Site_order <- data.frame("Site"=rownames(plant_pa))
  shuffle_CM_trait_based_inter <- left_join(Site_order, shuffle_CM_trait_based_inter, by="Site")
  rownames(shuffle_CM_trait_based_inter) <- shuffle_CM_trait_based_inter$Site
  shuffle_CM_trait_based_inter <- shuffle_CM_trait_based_inter[,-1]
  
  #### -- This computes path predicts the community mean in the context of extrapolation -- ####
  gam_data <- left_join(shuffle_CM_obs,env_var, by="Site")
  ## - For extrapolation group 1 - ##
  shuffle_models_CM_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
  shuffle_models_CM_LNC_extra1 <- predict(shuffle_models_CM_LNC_grow1, gam_data[group_extra1,])$predicted
  shuffle_models_CM_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
  shuffle_models_CM_SLA_extra1 <- predict(shuffle_models_CM_SLA_grow1, gam_data[group_extra1,])$predicted
  shuffle_models_CM_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
  shuffle_models_CM_PLH_extra1 <- predict(shuffle_models_CM_PLH_grow1, gam_data[group_extra1,])$predicted
  
  shuffle_CM_trait_based_extra1<-data.frame("LNC"=shuffle_models_CM_LNC_extra1, "SLA"=shuffle_models_CM_SLA_extra1, "PLH"=shuffle_models_CM_PLH_extra1)
  rownames(shuffle_CM_trait_based_extra1) <- shuffle_CM_obs$Site[group_extra1]
  
  ## - For extrapolation group 2 - ##
  shuffle_models_CM_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
  shuffle_models_CM_LNC_extra2 <- predict(shuffle_models_CM_LNC_grow2, gam_data[group_extra2,])$predicted
  shuffle_models_CM_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
  shuffle_models_CM_SLA_extra2 <- predict(shuffle_models_CM_SLA_grow2, gam_data[group_extra2,])$predicted
  shuffle_models_CM_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
  shuffle_models_CM_PLH_extra2 <- predict(shuffle_models_CM_PLH_grow2, gam_data[group_extra2,])$predicted
  
  shuffle_CM_trait_based_extra2<-data.frame("LNC"=shuffle_models_CM_LNC_extra2, "SLA"=shuffle_models_CM_SLA_extra2, "PLH"=shuffle_models_CM_PLH_extra2)
  rownames(shuffle_CM_trait_based_extra2) <- shuffle_CM_obs$Site[group_extra2]
  
  ## - For extrapolation group 3 - ##
  shuffle_models_CM_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
  shuffle_models_CM_LNC_extra3 <- predict(shuffle_models_CM_LNC_grow3, gam_data[group_extra3,])$predicted
  shuffle_models_CM_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
  shuffle_models_CM_SLA_extra3 <- predict(shuffle_models_CM_SLA_grow3, gam_data[group_extra3,])$predicted
  shuffle_models_CM_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
  shuffle_models_CM_PLH_extra3 <- predict(shuffle_models_CM_PLH_grow3, gam_data[group_extra3,])$predicted
  
  shuffle_CM_trait_based_extra3<-data.frame("LNC"=shuffle_models_CM_LNC_extra3, "SLA"=shuffle_models_CM_SLA_extra3, "PLH"=shuffle_models_CM_PLH_extra3)
  rownames(shuffle_CM_trait_based_extra3) <- shuffle_CM_obs$Site[group_extra3]
  
  ## - For extrapolation group 4 - ##
  shuffle_models_CM_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
  shuffle_models_CM_LNC_extra4 <- predict(shuffle_models_CM_LNC_grow4, gam_data[group_extra4,])$predicted
  shuffle_models_CM_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
  shuffle_models_CM_SLA_extra4 <- predict(shuffle_models_CM_SLA_grow4, gam_data[group_extra4,])$predicted
  shuffle_models_CM_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
  shuffle_models_CM_PLH_extra4 <- predict(shuffle_models_CM_PLH_grow4, gam_data[group_extra4,])$predicted
  
  shuffle_CM_trait_based_extra4<-data.frame("LNC"=shuffle_models_CM_LNC_extra4, "SLA"=shuffle_models_CM_SLA_extra4, "PLH"=shuffle_models_CM_PLH_extra4)
  rownames(shuffle_CM_trait_based_extra4) <- shuffle_CM_obs$Site[group_extra4]
  
  ## assembly in a single prediction table for the 4 groups :
  shuffle_CM_trait_based_extra <- rbind(shuffle_CM_trait_based_extra1, shuffle_CM_trait_based_extra2, shuffle_CM_trait_based_extra3, shuffle_CM_trait_based_extra4)
  shuffle_CM_trait_based_extra <- cbind("Site"=rownames(shuffle_CM_trait_based_extra), shuffle_CM_trait_based_extra)
  Site_order <- data.frame("Site"=rownames(plant_pa))
  shuffle_CM_trait_based_extra <- left_join(Site_order, shuffle_CM_trait_based_extra, by="Site")
  rownames(shuffle_CM_trait_based_extra) <- shuffle_CM_trait_based_extra$Site
  shuffle_CM_trait_based_extra <- shuffle_CM_trait_based_extra[,-1]
  
  ###### PREDICTION OF FUNCTIONAL DISPERSION  ######
  #### -- This computes path predicts the functional dispersion in the context of interpolation -- ####
  gam_data <- left_join(shuffle_FDis_obs,env_var, by="Site")
  ## - For interpolation group 1 - ##
  shuffle_models_FDis_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
  shuffle_models_FDis_LNC_inter1 <- predict(shuffle_models_FDis_LNC_grow1, gam_data[group_inter1,])$predicted
  shuffle_models_FDis_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
  shuffle_models_FDis_SLA_inter1 <- predict(shuffle_models_FDis_SLA_grow1, gam_data[group_inter1,])$predicted
  shuffle_models_FDis_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
  shuffle_models_FDis_PLH_inter1 <- predict(shuffle_models_FDis_PLH_grow1, gam_data[group_inter1,])$predicted
  
  shuffle_FDis_trait_based_inter1<-data.frame("LNC"=shuffle_models_FDis_LNC_inter1, "SLA"=shuffle_models_FDis_SLA_inter1, "PLH"=shuffle_models_FDis_PLH_inter1)
  rownames(shuffle_FDis_trait_based_inter1) <- shuffle_FDis_obs$Site[group_inter1]
  
  ## - For interpolation group 2 - ##
  shuffle_models_FDis_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
  shuffle_models_FDis_LNC_inter2 <- predict(shuffle_models_FDis_LNC_grow2, gam_data[group_inter2,])$predicted
  shuffle_models_FDis_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
  shuffle_models_FDis_SLA_inter2 <- predict(shuffle_models_FDis_SLA_grow2, gam_data[group_inter2,])$predicted
  shuffle_models_FDis_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
  shuffle_models_FDis_PLH_inter2 <- predict(shuffle_models_FDis_PLH_grow2, gam_data[group_inter2,])$predicted
  
  shuffle_FDis_trait_based_inter2<-data.frame("LNC"=shuffle_models_FDis_LNC_inter2, "SLA"=shuffle_models_FDis_SLA_inter2, "PLH"=shuffle_models_FDis_PLH_inter2)
  rownames(shuffle_FDis_trait_based_inter2) <- shuffle_FDis_obs$Site[group_inter2]
  
  ## - For interpolation group 3 - ##
  shuffle_models_FDis_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
  shuffle_models_FDis_LNC_inter3 <- predict(shuffle_models_FDis_LNC_grow3, gam_data[group_inter3,])$predicted
  shuffle_models_FDis_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
  shuffle_models_FDis_SLA_inter3 <- predict(shuffle_models_FDis_SLA_grow3, gam_data[group_inter3,])$predicted
  shuffle_models_FDis_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
  shuffle_models_FDis_PLH_inter3 <- predict(shuffle_models_FDis_PLH_grow3, gam_data[group_inter3,])$predicted
  
  shuffle_FDis_trait_based_inter3<-data.frame("LNC"=shuffle_models_FDis_LNC_inter3, "SLA"=shuffle_models_FDis_SLA_inter3, "PLH"=shuffle_models_FDis_PLH_inter3)
  rownames(shuffle_FDis_trait_based_inter3) <- shuffle_FDis_obs$Site[group_inter3]
  
  ## - For interpolation group 4 - ##
  shuffle_models_FDis_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
  shuffle_models_FDis_LNC_inter4 <- predict(shuffle_models_FDis_LNC_grow4, gam_data[group_inter4,])$predicted
  shuffle_models_FDis_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
  shuffle_models_FDis_SLA_inter4 <- predict(shuffle_models_FDis_SLA_grow4, gam_data[group_inter4,])$predicted
  shuffle_models_FDis_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
  shuffle_models_FDis_PLH_inter4 <- predict(shuffle_models_FDis_PLH_grow4, gam_data[group_inter4,])$predicted
  
  shuffle_FDis_trait_based_inter4<-data.frame("LNC"=shuffle_models_FDis_LNC_inter4, "SLA"=shuffle_models_FDis_SLA_inter4, "PLH"=shuffle_models_FDis_PLH_inter4)
  rownames(shuffle_FDis_trait_based_inter4) <- shuffle_FDis_obs$Site[group_inter4]
  
  ## assembly in a single prediction table for the 4 groups :
  shuffle_FDis_trait_based_inter <- rbind(shuffle_FDis_trait_based_inter1, shuffle_FDis_trait_based_inter2, shuffle_FDis_trait_based_inter3, shuffle_FDis_trait_based_inter4)
  shuffle_FDis_trait_based_inter <- cbind("Site"=rownames(shuffle_FDis_trait_based_inter), shuffle_FDis_trait_based_inter)
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_FDis_trait_based_inter <- left_join(Site_order, shuffle_FDis_trait_based_inter, by="Site")
  rownames(shuffle_FDis_trait_based_inter) <- shuffle_FDis_trait_based_inter$Site
  shuffle_FDis_trait_based_inter <- shuffle_FDis_trait_based_inter[,-1]
  
  #### -- This computes path predicts the functional dispersion in the context of extrapolation -- ####
  gam_data <- left_join(shuffle_FDis_obs,env_var, by="Site")
  ## - For extrapolation group 1 - ##
  shuffle_models_FDis_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
  shuffle_models_FDis_LNC_extra1 <- predict(shuffle_models_FDis_LNC_grow1, gam_data[group_extra1,])$predicted
  shuffle_models_FDis_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
  shuffle_models_FDis_SLA_extra1 <- predict(shuffle_models_FDis_SLA_grow1, gam_data[group_extra1,])$predicted
  shuffle_models_FDis_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
  shuffle_models_FDis_PLH_extra1 <- predict(shuffle_models_FDis_PLH_grow1, gam_data[group_extra1,])$predicted
  
  shuffle_FDis_trait_based_extra1<-data.frame("LNC"=shuffle_models_FDis_LNC_extra1, "SLA"=shuffle_models_FDis_SLA_extra1, "PLH"=shuffle_models_FDis_PLH_extra1)
  rownames(shuffle_FDis_trait_based_extra1) <- shuffle_FDis_obs$Site[group_extra1]
  
  ## - For extrapolation group 2 - ##
  shuffle_models_FDis_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
  shuffle_models_FDis_LNC_extra2 <- predict(shuffle_models_FDis_LNC_grow2, gam_data[group_extra2,])$predicted
  shuffle_models_FDis_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
  shuffle_models_FDis_SLA_extra2 <- predict(shuffle_models_FDis_SLA_grow2, gam_data[group_extra2,])$predicted
  shuffle_models_FDis_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
  shuffle_models_FDis_PLH_extra2 <- predict(shuffle_models_FDis_PLH_grow2, gam_data[group_extra2,])$predicted
  
  shuffle_FDis_trait_based_extra2<-data.frame("LNC"=shuffle_models_FDis_LNC_extra2, "SLA"=shuffle_models_FDis_SLA_extra2, "PLH"=shuffle_models_FDis_PLH_extra2)
  rownames(shuffle_FDis_trait_based_extra2) <- shuffle_FDis_obs$Site[group_extra2]
  
  ## - For extrapolation group 3 - ##
  shuffle_models_FDis_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
  shuffle_models_FDis_LNC_extra3 <- predict(shuffle_models_FDis_LNC_grow3, gam_data[group_extra3,])$predicted
  shuffle_models_FDis_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
  shuffle_models_FDis_SLA_extra3 <- predict(shuffle_models_FDis_SLA_grow3, gam_data[group_extra3,])$predicted
  shuffle_models_FDis_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
  shuffle_models_FDis_PLH_extra3 <- predict(shuffle_models_FDis_PLH_grow3, gam_data[group_extra3,])$predicted
  
  shuffle_FDis_trait_based_extra3<-data.frame("LNC"=shuffle_models_FDis_LNC_extra3, "SLA"=shuffle_models_FDis_SLA_extra3, "PLH"=shuffle_models_FDis_PLH_extra3)
  rownames(shuffle_FDis_trait_based_extra3) <- shuffle_FDis_obs$Site[group_extra3]
  
  ## - For extrapolation group 4 - ##
  shuffle_models_FDis_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
  shuffle_models_FDis_LNC_extra4 <- predict(shuffle_models_FDis_LNC_grow4, gam_data[group_extra4,])$predicted
  shuffle_models_FDis_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
  shuffle_models_FDis_SLA_extra4 <- predict(shuffle_models_FDis_SLA_grow4, gam_data[group_extra4,])$predicted
  shuffle_models_FDis_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
  shuffle_models_FDis_PLH_extra4 <- predict(shuffle_models_FDis_PLH_grow4, gam_data[group_extra4,])$predicted
  
  shuffle_FDis_trait_based_extra4<-data.frame("LNC"=shuffle_models_FDis_LNC_extra4, "SLA"=shuffle_models_FDis_SLA_extra4, "PLH"=shuffle_models_FDis_PLH_extra4)
  rownames(shuffle_FDis_trait_based_extra4) <- shuffle_FDis_obs$Site[group_extra4]
  
  ## assembly in a single prediction table for the 4 groups :
  shuffle_FDis_trait_based_extra <- rbind(shuffle_FDis_trait_based_extra1, shuffle_FDis_trait_based_extra2, shuffle_FDis_trait_based_extra3, shuffle_FDis_trait_based_extra4)
  shuffle_FDis_trait_based_extra <- cbind("Site"=rownames(shuffle_FDis_trait_based_extra), shuffle_FDis_trait_based_extra)
  Site_order <- data.frame("Site"=rownames(plant_recovery))
  shuffle_FDis_trait_based_extra <- left_join(Site_order, shuffle_FDis_trait_based_extra, by="Site")
  rownames(shuffle_FDis_trait_based_extra) <- shuffle_FDis_trait_based_extra$Site
  shuffle_FDis_trait_based_extra <- shuffle_FDis_trait_based_extra[,-1]
  
  ###### PREDICTION OF UNWEIGHTED FUNCTIONAL DISPERSION  ######
  #### -- This computes path predicts the unweighted functional dispersion in the context of interpolation -- ####
  gam_data <- left_join(shuffle_uFDis_obs,env_var, by="Site")
  ## - For interpolation group 1 - ##
  shuffle_models_uFDis_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
  shuffle_models_uFDis_LNC_inter1 <- predict(shuffle_models_uFDis_LNC_grow1, gam_data[group_inter1,])$predicted
  shuffle_models_uFDis_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
  shuffle_models_uFDis_SLA_inter1 <- predict(shuffle_models_uFDis_SLA_grow1, gam_data[group_inter1,])$predicted
  shuffle_models_uFDis_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter1,])
  shuffle_models_uFDis_PLH_inter1 <- predict(shuffle_models_uFDis_PLH_grow1, gam_data[group_inter1,])$predicted
  
  shuffle_uFDis_trait_based_inter1<-data.frame("LNC"=shuffle_models_uFDis_LNC_inter1, "SLA"=shuffle_models_uFDis_SLA_inter1, "PLH"=shuffle_models_uFDis_PLH_inter1)
  rownames(shuffle_uFDis_trait_based_inter1) <- shuffle_uFDis_obs$Site[group_inter1]
  
  ## - For interpolation group 2 - ##
  shuffle_models_uFDis_LNC_grow2<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
  shuffle_models_uFDis_LNC_inter2 <- predict(shuffle_models_uFDis_LNC_grow2, gam_data[group_inter2,])$predicted
  shuffle_models_uFDis_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
  shuffle_models_uFDis_SLA_inter2 <- predict(shuffle_models_uFDis_SLA_grow2, gam_data[group_inter2,])$predicted
  shuffle_models_uFDis_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter2,])
  shuffle_models_uFDis_PLH_inter2 <- predict(shuffle_models_uFDis_PLH_grow2, gam_data[group_inter2,])$predicted
  
  shuffle_uFDis_trait_based_inter2<-data.frame("LNC"=shuffle_models_uFDis_LNC_inter2, "SLA"=shuffle_models_uFDis_SLA_inter2, "PLH"=shuffle_models_uFDis_PLH_inter2)
  rownames(shuffle_uFDis_trait_based_inter2) <- shuffle_uFDis_obs$Site[group_inter2]
  
  ## - For interpolation group 3 - ##
  shuffle_models_uFDis_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
  shuffle_models_uFDis_LNC_inter3 <- predict(shuffle_models_uFDis_LNC_grow3, gam_data[group_inter3,])$predicted
  shuffle_models_uFDis_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
  shuffle_models_uFDis_SLA_inter3 <- predict(shuffle_models_uFDis_SLA_grow3, gam_data[group_inter3,])$predicted
  shuffle_models_uFDis_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter3,])
  shuffle_models_uFDis_PLH_inter3 <- predict(shuffle_models_uFDis_PLH_grow3, gam_data[group_inter3,])$predicted
  
  shuffle_uFDis_trait_based_inter3<-data.frame("LNC"=shuffle_models_uFDis_LNC_inter3, "SLA"=shuffle_models_uFDis_SLA_inter3, "PLH"=shuffle_models_uFDis_PLH_inter3)
  rownames(shuffle_uFDis_trait_based_inter3) <- shuffle_uFDis_obs$Site[group_inter3]
  
  ## - For interpolation group 4 - ##
  shuffle_models_uFDis_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
  shuffle_models_uFDis_LNC_inter4 <- predict(shuffle_models_uFDis_LNC_grow4, gam_data[group_inter4,])$predicted
  shuffle_models_uFDis_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
  shuffle_models_uFDis_SLA_inter4 <- predict(shuffle_models_uFDis_SLA_grow4, gam_data[group_inter4,])$predicted
  shuffle_models_uFDis_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_inter4,])
  shuffle_models_uFDis_PLH_inter4 <- predict(shuffle_models_uFDis_PLH_grow4, gam_data[group_inter4,])$predicted
  
  shuffle_uFDis_trait_based_inter4<-data.frame("LNC"=shuffle_models_uFDis_LNC_inter4, "SLA"=shuffle_models_uFDis_SLA_inter4, "PLH"=shuffle_models_uFDis_PLH_inter4)
  rownames(shuffle_uFDis_trait_based_inter4) <- shuffle_uFDis_obs$Site[group_inter4]
  
  ## assembly in a single prediction table for the 4 groups :
  shuffle_uFDis_trait_based_inter <- rbind(shuffle_uFDis_trait_based_inter1, shuffle_uFDis_trait_based_inter2, shuffle_uFDis_trait_based_inter3, shuffle_uFDis_trait_based_inter4)
  shuffle_uFDis_trait_based_inter <- cbind("Site"=rownames(shuffle_uFDis_trait_based_inter), shuffle_uFDis_trait_based_inter)
  Site_order <- data.frame("Site"=rownames(plant_pa))
  shuffle_uFDis_trait_based_inter <- left_join(Site_order, shuffle_uFDis_trait_based_inter, by="Site")
  rownames(shuffle_uFDis_trait_based_inter) <- shuffle_uFDis_trait_based_inter$Site
  shuffle_uFDis_trait_based_inter <- shuffle_uFDis_trait_based_inter[,-1]
  
  #### -- This computes path predicts the community mean in the context of extrapolation -- ####
  ## - For extrapolation group 1 - ##
  shuffle_models_uFDis_LNC_grow1<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
  shuffle_models_uFDis_LNC_extra1 <- predict(shuffle_models_uFDis_LNC_grow1, gam_data[group_extra1,])$predicted
  shuffle_models_uFDis_SLA_grow1<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
  shuffle_models_uFDis_SLA_extra1 <- predict(shuffle_models_uFDis_SLA_grow1, gam_data[group_extra1,])$predicted
  shuffle_models_uFDis_PLH_grow1<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra1,])
  shuffle_models_uFDis_PLH_extra1 <- predict(shuffle_models_uFDis_PLH_grow1, gam_data[group_extra1,])$predicted
  
  shuffle_uFDis_trait_based_extra1<-data.frame("LNC"=shuffle_models_uFDis_LNC_extra1, "SLA"=shuffle_models_uFDis_SLA_extra1, "PLH"=shuffle_models_uFDis_PLH_extra1)
  rownames(shuffle_uFDis_trait_based_extra1) <- shuffle_uFDis_obs$Site[group_extra1]
  
  ## - For extrapolation group 2 - ##
  shuffle_models_uFDis_LNC_grow2 <- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
  shuffle_models_uFDis_LNC_extra2 <- predict(shuffle_models_uFDis_LNC_grow2, gam_data[group_extra2,])$predicted
  shuffle_models_uFDis_SLA_grow2<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
  shuffle_models_uFDis_SLA_extra2 <- predict(shuffle_models_uFDis_SLA_grow2, gam_data[group_extra2,])$predicted
  shuffle_models_uFDis_PLH_grow2<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra2,])
  shuffle_models_uFDis_PLH_extra2 <- predict(shuffle_models_uFDis_PLH_grow2, gam_data[group_extra2,])$predicted
  
  shuffle_uFDis_trait_based_extra2<-data.frame("LNC"=shuffle_models_uFDis_LNC_extra2, "SLA"=shuffle_models_uFDis_SLA_extra2, "PLH"=shuffle_models_uFDis_PLH_extra2)
  rownames(shuffle_uFDis_trait_based_extra2) <- shuffle_uFDis_obs$Site[group_extra2]
  
  ## - For extrapolation group 3 - ##
  shuffle_models_uFDis_LNC_grow3<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
  shuffle_models_uFDis_LNC_extra3 <- predict(shuffle_models_uFDis_LNC_grow3, gam_data[group_extra3,])$predicted
  shuffle_models_uFDis_SLA_grow3<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
  shuffle_models_uFDis_SLA_extra3 <- predict(shuffle_models_uFDis_SLA_grow3, gam_data[group_extra3,])$predicted
  shuffle_models_uFDis_PLH_grow3<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra3,])
  shuffle_models_uFDis_PLH_extra3 <- predict(shuffle_models_uFDis_PLH_grow3, gam_data[group_extra3,])$predicted
  
  shuffle_uFDis_trait_based_extra3<-data.frame("LNC"=shuffle_models_uFDis_LNC_extra3, "SLA"=shuffle_models_uFDis_SLA_extra3, "PLH"=shuffle_models_uFDis_PLH_extra3)
  rownames(shuffle_uFDis_trait_based_extra3) <- shuffle_uFDis_obs$Site[group_extra3]
  
  ## - For extrapolation group 4 - ##
  shuffle_models_uFDis_LNC_grow4<- rfsrc(LNC ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
  shuffle_models_uFDis_LNC_extra4 <- predict(shuffle_models_uFDis_LNC_grow4, gam_data[group_extra4,])$predicted
  shuffle_models_uFDis_SLA_grow4<- rfsrc(SLA ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
  shuffle_models_uFDis_SLA_extra4 <- predict(shuffle_models_uFDis_SLA_grow4, gam_data[group_extra4,])$predicted
  shuffle_models_uFDis_PLH_grow4<- rfsrc(PLH ~ Aridity + Temperature_Seasonality + Mean_Temperature_of_Coldest_Quarter + Slope + Precipitation_Seasonality + Soil_ind.F_rdf + Soil_ind.R_rdf, data=gam_data[-group_extra4,])
  shuffle_models_uFDis_PLH_extra4 <- predict(shuffle_models_uFDis_PLH_grow4, gam_data[group_extra4,])$predicted
  
  shuffle_uFDis_trait_based_extra4<-data.frame("LNC"=shuffle_models_uFDis_LNC_extra4, "SLA"=shuffle_models_uFDis_SLA_extra4, "PLH"=shuffle_models_uFDis_PLH_extra4)
  rownames(shuffle_uFDis_trait_based_extra4) <- shuffle_uFDis_obs$Site[group_extra4]
  
  ## assembly in a single prediction table for the 4 groups :
  shuffle_uFDis_trait_based_extra <- rbind(shuffle_uFDis_trait_based_extra1, shuffle_uFDis_trait_based_extra2, shuffle_uFDis_trait_based_extra3, shuffle_uFDis_trait_based_extra4)
  shuffle_uFDis_trait_based_extra <- cbind("Site"=rownames(shuffle_uFDis_trait_based_extra), shuffle_uFDis_trait_based_extra)
  Site_order <- data.frame("Site"=rownames(plant_pa))
  shuffle_uFDis_trait_based_extra <- left_join(Site_order, shuffle_uFDis_trait_based_extra, by="Site")
  rownames(shuffle_uFDis_trait_based_extra) <- shuffle_uFDis_trait_based_extra$Site
  shuffle_uFDis_trait_based_extra <- shuffle_uFDis_trait_based_extra[,-1]
  
  
  shuffle_figure3 <- rbind(shuffle_figure3, data.frame("RMSE"=c(rmse(shuffle_CWM_obs$LNC[group_inter1], shuffle_CWM_trait_based_inter$LNC[group_inter1]), # CWM trait based approach interpolation
                                                                rmse(shuffle_CWM_obs$SLA[group_inter1], shuffle_CWM_trait_based_inter$SLA[group_inter1]),
                                                                rmse(shuffle_CWM_obs$PLH[group_inter1], shuffle_CWM_trait_based_inter$PLH[group_inter1]),
                                                                rmse(shuffle_CWM_obs$LNC[group_inter2], shuffle_CWM_trait_based_inter$LNC[group_inter2]),
                                                                rmse(shuffle_CWM_obs$SLA[group_inter2], shuffle_CWM_trait_based_inter$SLA[group_inter2]),
                                                                rmse(shuffle_CWM_obs$PLH[group_inter2], shuffle_CWM_trait_based_inter$PLH[group_inter2]),
                                                                rmse(shuffle_CWM_obs$LNC[group_inter3], shuffle_CWM_trait_based_inter$LNC[group_inter3]),
                                                                rmse(shuffle_CWM_obs$SLA[group_inter3], shuffle_CWM_trait_based_inter$SLA[group_inter3]),
                                                                rmse(shuffle_CWM_obs$PLH[group_inter3], shuffle_CWM_trait_based_inter$PLH[group_inter3]),
                                                                rmse(shuffle_CWM_obs$LNC[group_inter4], shuffle_CWM_trait_based_inter$LNC[group_inter4]),
                                                                rmse(shuffle_CWM_obs$SLA[group_inter4], shuffle_CWM_trait_based_inter$SLA[group_inter4]),
                                                                rmse(shuffle_CWM_obs$PLH[group_inter4], shuffle_CWM_trait_based_inter$PLH[group_inter4]),
                                                                rmse(shuffle_CWM_obs$LNC[group_extra1], shuffle_CWM_trait_based_extra$LNC[group_extra1]),# CWM trait based approach interpolation
                                                                rmse(shuffle_CWM_obs$SLA[group_extra1], shuffle_CWM_trait_based_extra$SLA[group_extra1]),
                                                                rmse(shuffle_CWM_obs$PLH[group_extra1], shuffle_CWM_trait_based_extra$PLH[group_extra1]),
                                                                rmse(shuffle_CWM_obs$LNC[group_extra2], shuffle_CWM_trait_based_extra$LNC[group_extra2]),
                                                                rmse(shuffle_CWM_obs$SLA[group_extra2], shuffle_CWM_trait_based_extra$SLA[group_extra2]),
                                                                rmse(shuffle_CWM_obs$PLH[group_extra2], shuffle_CWM_trait_based_extra$PLH[group_extra2]),
                                                                rmse(shuffle_CWM_obs$LNC[group_extra3], shuffle_CWM_trait_based_extra$LNC[group_extra3]),
                                                                rmse(shuffle_CWM_obs$SLA[group_extra3], shuffle_CWM_trait_based_extra$SLA[group_extra3]),
                                                                rmse(shuffle_CWM_obs$PLH[group_extra3], shuffle_CWM_trait_based_extra$PLH[group_extra3]),
                                                                rmse(shuffle_CWM_obs$LNC[group_extra4], shuffle_CWM_trait_based_extra$LNC[group_extra4]),
                                                                rmse(shuffle_CWM_obs$SLA[group_extra4], shuffle_CWM_trait_based_extra$SLA[group_extra4]),
                                                                rmse(shuffle_CWM_obs$PLH[group_extra4], shuffle_CWM_trait_based_extra$PLH[group_extra4]),
                                                                rmse(shuffle_FDis_obs$LNC[group_inter1], shuffle_FDis_trait_based_inter$LNC[group_inter1]),# FDis trait based approach interpolation
                                                                rmse(shuffle_FDis_obs$SLA[group_inter1], shuffle_FDis_trait_based_inter$SLA[group_inter1]),
                                                                rmse(shuffle_FDis_obs$PLH[group_inter1], shuffle_FDis_trait_based_inter$PLH[group_inter1]),
                                                                rmse(shuffle_FDis_obs$LNC[group_inter2], shuffle_FDis_trait_based_inter$LNC[group_inter2]),
                                                                rmse(shuffle_FDis_obs$SLA[group_inter2], shuffle_FDis_trait_based_inter$SLA[group_inter2]),
                                                                rmse(shuffle_FDis_obs$PLH[group_inter2], shuffle_FDis_trait_based_inter$PLH[group_inter2]),
                                                                rmse(shuffle_FDis_obs$LNC[group_inter3], shuffle_FDis_trait_based_inter$LNC[group_inter3]),
                                                                rmse(shuffle_FDis_obs$SLA[group_inter3], shuffle_FDis_trait_based_inter$SLA[group_inter3]),
                                                                rmse(shuffle_FDis_obs$PLH[group_inter3], shuffle_FDis_trait_based_inter$PLH[group_inter3]),
                                                                rmse(shuffle_FDis_obs$LNC[group_inter4], shuffle_FDis_trait_based_inter$LNC[group_inter4]),
                                                                rmse(shuffle_FDis_obs$SLA[group_inter4], shuffle_FDis_trait_based_inter$SLA[group_inter4]),
                                                                rmse(shuffle_FDis_obs$PLH[group_inter4], shuffle_FDis_trait_based_inter$PLH[group_inter4]),
                                                                rmse(shuffle_FDis_obs$LNC[group_extra1], shuffle_FDis_trait_based_extra$LNC[group_extra1]), # FDis trait based approach extrapolation
                                                                rmse(shuffle_FDis_obs$SLA[group_extra1], shuffle_FDis_trait_based_extra$SLA[group_extra1]),
                                                                rmse(shuffle_FDis_obs$PLH[group_extra1], shuffle_FDis_trait_based_extra$PLH[group_extra1]),
                                                                rmse(shuffle_FDis_obs$LNC[group_extra2], shuffle_FDis_trait_based_extra$LNC[group_extra2]),
                                                                rmse(shuffle_FDis_obs$SLA[group_extra2], shuffle_FDis_trait_based_extra$SLA[group_extra2]),
                                                                rmse(shuffle_FDis_obs$PLH[group_extra2], shuffle_FDis_trait_based_extra$PLH[group_extra2]),
                                                                rmse(shuffle_FDis_obs$LNC[group_extra3], shuffle_FDis_trait_based_extra$LNC[group_extra3]),
                                                                rmse(shuffle_FDis_obs$SLA[group_extra3], shuffle_FDis_trait_based_extra$SLA[group_extra3]),
                                                                rmse(shuffle_FDis_obs$PLH[group_extra3], shuffle_FDis_trait_based_extra$PLH[group_extra3]),
                                                                rmse(shuffle_FDis_obs$LNC[group_extra4], shuffle_FDis_trait_based_extra$LNC[group_extra4]),
                                                                rmse(shuffle_FDis_obs$SLA[group_extra4], shuffle_FDis_trait_based_extra$SLA[group_extra4]),
                                                                rmse(shuffle_FDis_obs$PLH[group_extra4], shuffle_FDis_trait_based_extra$PLH[group_extra4]),
                                                                rmse(shuffle_CWM_obs$LNC[group_inter1], shuffle_CWM_species_based_inter$LNC[group_inter1]), # CWM species based approach interpolation
                                                                rmse(shuffle_CWM_obs$SLA[group_inter1], shuffle_CWM_species_based_inter$SLA[group_inter1]),
                                                                rmse(shuffle_CWM_obs$PLH[group_inter1], shuffle_CWM_species_based_inter$PLH[group_inter1]),
                                                                rmse(shuffle_CWM_obs$LNC[group_inter2], shuffle_CWM_species_based_inter$LNC[group_inter2]),
                                                                rmse(shuffle_CWM_obs$SLA[group_inter2], shuffle_CWM_species_based_inter$SLA[group_inter2]),
                                                                rmse(shuffle_CWM_obs$PLH[group_inter2], shuffle_CWM_species_based_inter$PLH[group_inter2]),
                                                                rmse(shuffle_CWM_obs$LNC[group_inter3], shuffle_CWM_species_based_inter$LNC[group_inter3]),
                                                                rmse(shuffle_CWM_obs$SLA[group_inter3], shuffle_CWM_species_based_inter$SLA[group_inter3]),
                                                                rmse(shuffle_CWM_obs$PLH[group_inter3], shuffle_CWM_species_based_inter$PLH[group_inter3]),
                                                                rmse(shuffle_CWM_obs$LNC[group_inter4], shuffle_CWM_species_based_inter$LNC[group_inter4]),
                                                                rmse(shuffle_CWM_obs$SLA[group_inter4], shuffle_CWM_species_based_inter$SLA[group_inter4]),
                                                                rmse(shuffle_CWM_obs$PLH[group_inter4], shuffle_CWM_species_based_inter$PLH[group_inter4]),
                                                                rmse(shuffle_CWM_obs$LNC[group_extra1], shuffle_CWM_species_based_extra$LNC[group_extra1]), # CWM species based approach extrapolation
                                                                rmse(shuffle_CWM_obs$SLA[group_extra1], shuffle_CWM_species_based_extra$SLA[group_extra1]),
                                                                rmse(shuffle_CWM_obs$PLH[group_extra1], shuffle_CWM_species_based_extra$PLH[group_extra1]),
                                                                rmse(shuffle_CWM_obs$LNC[group_extra2], shuffle_CWM_species_based_extra$LNC[group_extra2]),
                                                                rmse(shuffle_CWM_obs$SLA[group_extra2], shuffle_CWM_species_based_extra$SLA[group_extra2]),
                                                                rmse(shuffle_CWM_obs$PLH[group_extra2], shuffle_CWM_species_based_extra$PLH[group_extra2]),
                                                                rmse(shuffle_CWM_obs$LNC[group_extra3], shuffle_CWM_species_based_extra$LNC[group_extra3]),
                                                                rmse(shuffle_CWM_obs$SLA[group_extra3], shuffle_CWM_species_based_extra$SLA[group_extra3]),
                                                                rmse(shuffle_CWM_obs$PLH[group_extra3], shuffle_CWM_species_based_extra$PLH[group_extra3]),
                                                                rmse(shuffle_CWM_obs$LNC[group_extra4], shuffle_CWM_species_based_extra$LNC[group_extra4]),
                                                                rmse(shuffle_CWM_obs$SLA[group_extra4], shuffle_CWM_species_based_extra$SLA[group_extra4]),
                                                                rmse(shuffle_CWM_obs$PLH[group_extra4], shuffle_CWM_species_based_extra$PLH[group_extra4]),
                                                                rmse(shuffle_FDis_obs$LNC[group_inter1], shuffle_FDis_species_based_inter$LNC[group_inter1]), # FDis species based approach interpolation
                                                                rmse(shuffle_FDis_obs$SLA[group_inter1], shuffle_FDis_species_based_inter$SLA[group_inter1]),
                                                                rmse(shuffle_FDis_obs$PLH[group_inter1], shuffle_FDis_species_based_inter$PLH[group_inter1]),
                                                                rmse(shuffle_FDis_obs$LNC[group_inter2], shuffle_FDis_species_based_inter$LNC[group_inter2]),
                                                                rmse(shuffle_FDis_obs$SLA[group_inter2], shuffle_FDis_species_based_inter$SLA[group_inter2]),
                                                                rmse(shuffle_FDis_obs$PLH[group_inter2], shuffle_FDis_species_based_inter$PLH[group_inter2]),
                                                                rmse(shuffle_FDis_obs$LNC[group_inter3], shuffle_FDis_species_based_inter$LNC[group_inter3]),
                                                                rmse(shuffle_FDis_obs$SLA[group_inter3], shuffle_FDis_species_based_inter$SLA[group_inter3]),
                                                                rmse(shuffle_FDis_obs$PLH[group_inter3], shuffle_FDis_species_based_inter$PLH[group_inter3]),
                                                                rmse(shuffle_FDis_obs$LNC[group_inter4], shuffle_FDis_species_based_inter$LNC[group_inter4]),
                                                                rmse(shuffle_FDis_obs$SLA[group_inter4], shuffle_FDis_species_based_inter$SLA[group_inter4]),
                                                                rmse(shuffle_FDis_obs$PLH[group_inter4], shuffle_FDis_species_based_inter$PLH[group_inter4]),
                                                                rmse(shuffle_FDis_obs$LNC[group_extra1], shuffle_FDis_species_based_extra$LNC[group_extra1]), # FDis species based approach extrapolation
                                                                rmse(shuffle_FDis_obs$SLA[group_extra1], shuffle_FDis_species_based_extra$SLA[group_extra1]),
                                                                rmse(shuffle_FDis_obs$PLH[group_extra1], shuffle_FDis_species_based_extra$PLH[group_extra1]),
                                                                rmse(shuffle_FDis_obs$LNC[group_extra2], shuffle_FDis_species_based_extra$LNC[group_extra2]),
                                                                rmse(shuffle_FDis_obs$SLA[group_extra2], shuffle_FDis_species_based_extra$SLA[group_extra2]),
                                                                rmse(shuffle_FDis_obs$PLH[group_extra2], shuffle_FDis_species_based_extra$PLH[group_extra2]),
                                                                rmse(shuffle_FDis_obs$LNC[group_extra3], shuffle_FDis_species_based_extra$LNC[group_extra3]),
                                                                rmse(shuffle_FDis_obs$SLA[group_extra3], shuffle_FDis_species_based_extra$SLA[group_extra3]),
                                                                rmse(shuffle_FDis_obs$PLH[group_extra3], shuffle_FDis_species_based_extra$PLH[group_extra3]),
                                                                rmse(shuffle_FDis_obs$LNC[group_extra4], shuffle_FDis_species_based_extra$LNC[group_extra4]),
                                                                rmse(shuffle_FDis_obs$SLA[group_extra4], shuffle_FDis_species_based_extra$SLA[group_extra4]),
                                                                rmse(shuffle_FDis_obs$PLH[group_extra4], shuffle_FDis_species_based_extra$PLH[group_extra4]),
                                                                
                                                                rmse(shuffle_CM_obs$LNC[group_inter1], shuffle_CM_trait_based_inter$LNC[group_inter1]), # CM trait based approach interpolation
                                                                rmse(shuffle_CM_obs$SLA[group_inter1], shuffle_CM_trait_based_inter$SLA[group_inter1]),
                                                                rmse(shuffle_CM_obs$PLH[group_inter1], shuffle_CM_trait_based_inter$PLH[group_inter1]),
                                                                rmse(shuffle_CM_obs$LNC[group_inter2], shuffle_CM_trait_based_inter$LNC[group_inter2]),
                                                                rmse(shuffle_CM_obs$SLA[group_inter2], shuffle_CM_trait_based_inter$SLA[group_inter2]),
                                                                rmse(shuffle_CM_obs$PLH[group_inter2], shuffle_CM_trait_based_inter$PLH[group_inter2]),
                                                                rmse(shuffle_CM_obs$LNC[group_inter3], shuffle_CM_trait_based_inter$LNC[group_inter3]),
                                                                rmse(shuffle_CM_obs$SLA[group_inter3], shuffle_CM_trait_based_inter$SLA[group_inter3]),
                                                                rmse(shuffle_CM_obs$PLH[group_inter3], shuffle_CM_trait_based_inter$PLH[group_inter3]),
                                                                rmse(shuffle_CM_obs$LNC[group_inter4], shuffle_CM_trait_based_inter$LNC[group_inter4]),
                                                                rmse(shuffle_CM_obs$SLA[group_inter4], shuffle_CM_trait_based_inter$SLA[group_inter4]),
                                                                rmse(shuffle_CM_obs$PLH[group_inter4], shuffle_CM_trait_based_inter$PLH[group_inter4]),
                                                                rmse(shuffle_CM_obs$LNC[group_extra1], shuffle_CM_trait_based_extra$LNC[group_extra1]),# CM trait based approach interpolation
                                                                rmse(shuffle_CM_obs$SLA[group_extra1], shuffle_CM_trait_based_extra$SLA[group_extra1]),
                                                                rmse(shuffle_CM_obs$PLH[group_extra1], shuffle_CM_trait_based_extra$PLH[group_extra1]),
                                                                rmse(shuffle_CM_obs$LNC[group_extra2], shuffle_CM_trait_based_extra$LNC[group_extra2]),
                                                                rmse(shuffle_CM_obs$SLA[group_extra2], shuffle_CM_trait_based_extra$SLA[group_extra2]),
                                                                rmse(shuffle_CM_obs$PLH[group_extra2], shuffle_CM_trait_based_extra$PLH[group_extra2]),
                                                                rmse(shuffle_CM_obs$LNC[group_extra3], shuffle_CM_trait_based_extra$LNC[group_extra3]),
                                                                rmse(shuffle_CM_obs$SLA[group_extra3], shuffle_CM_trait_based_extra$SLA[group_extra3]),
                                                                rmse(shuffle_CM_obs$PLH[group_extra3], shuffle_CM_trait_based_extra$PLH[group_extra3]),
                                                                rmse(shuffle_CM_obs$LNC[group_extra4], shuffle_CM_trait_based_extra$LNC[group_extra4]),
                                                                rmse(shuffle_CM_obs$SLA[group_extra4], shuffle_CM_trait_based_extra$SLA[group_extra4]),
                                                                rmse(shuffle_CM_obs$PLH[group_extra4], shuffle_CM_trait_based_extra$PLH[group_extra4]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_inter1], shuffle_uFDis_trait_based_inter$LNC[group_inter1]),# uuFDis trait based approach interpolation
                                                                rmse(shuffle_uFDis_obs$SLA[group_inter1], shuffle_uFDis_trait_based_inter$SLA[group_inter1]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_inter1], shuffle_uFDis_trait_based_inter$PLH[group_inter1]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_inter2], shuffle_uFDis_trait_based_inter$LNC[group_inter2]),
                                                                rmse(shuffle_uFDis_obs$SLA[group_inter2], shuffle_uFDis_trait_based_inter$SLA[group_inter2]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_inter2], shuffle_uFDis_trait_based_inter$PLH[group_inter2]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_inter3], shuffle_uFDis_trait_based_inter$LNC[group_inter3]),
                                                                rmse(shuffle_uFDis_obs$SLA[group_inter3], shuffle_uFDis_trait_based_inter$SLA[group_inter3]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_inter3], shuffle_uFDis_trait_based_inter$PLH[group_inter3]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_inter4], shuffle_uFDis_trait_based_inter$LNC[group_inter4]),
                                                                rmse(shuffle_uFDis_obs$SLA[group_inter4], shuffle_uFDis_trait_based_inter$SLA[group_inter4]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_inter4], shuffle_uFDis_trait_based_inter$PLH[group_inter4]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_extra1], shuffle_uFDis_trait_based_extra$LNC[group_extra1]), # uuFDis trait based approach extrapolation
                                                                rmse(shuffle_uFDis_obs$SLA[group_extra1], shuffle_uFDis_trait_based_extra$SLA[group_extra1]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_extra1], shuffle_uFDis_trait_based_extra$PLH[group_extra1]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_extra2], shuffle_uFDis_trait_based_extra$LNC[group_extra2]),
                                                                rmse(shuffle_uFDis_obs$SLA[group_extra2], shuffle_uFDis_trait_based_extra$SLA[group_extra2]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_extra2], shuffle_uFDis_trait_based_extra$PLH[group_extra2]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_extra3], shuffle_uFDis_trait_based_extra$LNC[group_extra3]),
                                                                rmse(shuffle_uFDis_obs$SLA[group_extra3], shuffle_uFDis_trait_based_extra$SLA[group_extra3]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_extra3], shuffle_uFDis_trait_based_extra$PLH[group_extra3]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_extra4], shuffle_uFDis_trait_based_extra$LNC[group_extra4]),
                                                                rmse(shuffle_uFDis_obs$SLA[group_extra4], shuffle_uFDis_trait_based_extra$SLA[group_extra4]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_extra4], shuffle_uFDis_trait_based_extra$PLH[group_extra4]),
                                                                rmse(shuffle_CM_obs$LNC[group_inter1], shuffle_CM_species_based_inter$LNC[group_inter1]), # CM species based approach interpolation
                                                                rmse(shuffle_CM_obs$SLA[group_inter1], shuffle_CM_species_based_inter$SLA[group_inter1]),
                                                                rmse(shuffle_CM_obs$PLH[group_inter1], shuffle_CM_species_based_inter$PLH[group_inter1]),
                                                                rmse(shuffle_CM_obs$LNC[group_inter2], shuffle_CM_species_based_inter$LNC[group_inter2]),
                                                                rmse(shuffle_CM_obs$SLA[group_inter2], shuffle_CM_species_based_inter$SLA[group_inter2]),
                                                                rmse(shuffle_CM_obs$PLH[group_inter2], shuffle_CM_species_based_inter$PLH[group_inter2]),
                                                                rmse(shuffle_CM_obs$LNC[group_inter3], shuffle_CM_species_based_inter$LNC[group_inter3]),
                                                                rmse(shuffle_CM_obs$SLA[group_inter3], shuffle_CM_species_based_inter$SLA[group_inter3]),
                                                                rmse(shuffle_CM_obs$PLH[group_inter3], shuffle_CM_species_based_inter$PLH[group_inter3]),
                                                                rmse(shuffle_CM_obs$LNC[group_inter4], shuffle_CM_species_based_inter$LNC[group_inter4]),
                                                                rmse(shuffle_CM_obs$SLA[group_inter4], shuffle_CM_species_based_inter$SLA[group_inter4]),
                                                                rmse(shuffle_CM_obs$PLH[group_inter4], shuffle_CM_species_based_inter$PLH[group_inter4]),
                                                                rmse(shuffle_CM_obs$LNC[group_extra1], shuffle_CM_species_based_extra$LNC[group_extra1]), # CM species based approach extrapolation
                                                                rmse(shuffle_CM_obs$SLA[group_extra1], shuffle_CM_species_based_extra$SLA[group_extra1]),
                                                                rmse(shuffle_CM_obs$PLH[group_extra1], shuffle_CM_species_based_extra$PLH[group_extra1]),
                                                                rmse(shuffle_CM_obs$LNC[group_extra2], shuffle_CM_species_based_extra$LNC[group_extra2]),
                                                                rmse(shuffle_CM_obs$SLA[group_extra2], shuffle_CM_species_based_extra$SLA[group_extra2]),
                                                                rmse(shuffle_CM_obs$PLH[group_extra2], shuffle_CM_species_based_extra$PLH[group_extra2]),
                                                                rmse(shuffle_CM_obs$LNC[group_extra3], shuffle_CM_species_based_extra$LNC[group_extra3]),
                                                                rmse(shuffle_CM_obs$SLA[group_extra3], shuffle_CM_species_based_extra$SLA[group_extra3]),
                                                                rmse(shuffle_CM_obs$PLH[group_extra3], shuffle_CM_species_based_extra$PLH[group_extra3]),
                                                                rmse(shuffle_CM_obs$LNC[group_extra4], shuffle_CM_species_based_extra$LNC[group_extra4]),
                                                                rmse(shuffle_CM_obs$SLA[group_extra4], shuffle_CM_species_based_extra$SLA[group_extra4]),
                                                                rmse(shuffle_CM_obs$PLH[group_extra4], shuffle_CM_species_based_extra$PLH[group_extra4]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_inter1], shuffle_uFDis_species_based_inter$LNC[group_inter1]), # uuFDis species based approach interpolation
                                                                rmse(shuffle_uFDis_obs$SLA[group_inter1], shuffle_uFDis_species_based_inter$SLA[group_inter1]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_inter1], shuffle_uFDis_species_based_inter$PLH[group_inter1]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_inter2], shuffle_uFDis_species_based_inter$LNC[group_inter2]),
                                                                rmse(shuffle_uFDis_obs$SLA[group_inter2], shuffle_uFDis_species_based_inter$SLA[group_inter2]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_inter2], shuffle_uFDis_species_based_inter$PLH[group_inter2]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_inter3], shuffle_uFDis_species_based_inter$LNC[group_inter3]),
                                                                rmse(shuffle_uFDis_obs$SLA[group_inter3], shuffle_uFDis_species_based_inter$SLA[group_inter3]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_inter3], shuffle_uFDis_species_based_inter$PLH[group_inter3]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_inter4], shuffle_uFDis_species_based_inter$LNC[group_inter4]),
                                                                rmse(shuffle_uFDis_obs$SLA[group_inter4], shuffle_uFDis_species_based_inter$SLA[group_inter4]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_inter4], shuffle_uFDis_species_based_inter$PLH[group_inter4]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_extra1], shuffle_uFDis_species_based_extra$LNC[group_extra1]), # uuFDis species based approach extrapolation
                                                                rmse(shuffle_uFDis_obs$SLA[group_extra1], shuffle_uFDis_species_based_extra$SLA[group_extra1]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_extra1], shuffle_uFDis_species_based_extra$PLH[group_extra1]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_extra2], shuffle_uFDis_species_based_extra$LNC[group_extra2]),
                                                                rmse(shuffle_uFDis_obs$SLA[group_extra2], shuffle_uFDis_species_based_extra$SLA[group_extra2]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_extra2], shuffle_uFDis_species_based_extra$PLH[group_extra2]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_extra3], shuffle_uFDis_species_based_extra$LNC[group_extra3]),
                                                                rmse(shuffle_uFDis_obs$SLA[group_extra3], shuffle_uFDis_species_based_extra$SLA[group_extra3]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_extra3], shuffle_uFDis_species_based_extra$PLH[group_extra3]),
                                                                rmse(shuffle_uFDis_obs$LNC[group_extra4], shuffle_uFDis_species_based_extra$LNC[group_extra4]),
                                                                rmse(shuffle_uFDis_obs$SLA[group_extra4], shuffle_uFDis_species_based_extra$SLA[group_extra4]),
                                                                rmse(shuffle_uFDis_obs$PLH[group_extra4], shuffle_uFDis_species_based_extra$PLH[group_extra4])),
                                                       "R2"=c(cor(shuffle_CWM_obs$LNC[group_inter1], shuffle_CWM_trait_based_inter$LNC[group_inter1])^2, # CWM trait based approach interpolation
                                                              cor(shuffle_CWM_obs$SLA[group_inter1], shuffle_CWM_trait_based_inter$SLA[group_inter1])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_inter1], shuffle_CWM_trait_based_inter$PLH[group_inter1])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_inter2], shuffle_CWM_trait_based_inter$LNC[group_inter2])^2,
                                                              cor(shuffle_CWM_obs$SLA[group_inter2], shuffle_CWM_trait_based_inter$SLA[group_inter2])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_inter2], shuffle_CWM_trait_based_inter$PLH[group_inter2])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_inter3], shuffle_CWM_trait_based_inter$LNC[group_inter3])^2,
                                                              cor(shuffle_CWM_obs$SLA[group_inter3], shuffle_CWM_trait_based_inter$SLA[group_inter3])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_inter3], shuffle_CWM_trait_based_inter$PLH[group_inter3])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_inter4], shuffle_CWM_trait_based_inter$LNC[group_inter4])^2,
                                                              cor(shuffle_CWM_obs$SLA[group_inter4], shuffle_CWM_trait_based_inter$SLA[group_inter4])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_inter4], shuffle_CWM_trait_based_inter$PLH[group_inter4])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_extra1], shuffle_CWM_trait_based_extra$LNC[group_extra1])^2,# CWM trait based approach interpolation
                                                              cor(shuffle_CWM_obs$SLA[group_extra1], shuffle_CWM_trait_based_extra$SLA[group_extra1])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_extra1], shuffle_CWM_trait_based_extra$PLH[group_extra1])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_extra2], shuffle_CWM_trait_based_extra$LNC[group_extra2])^2,
                                                              cor(shuffle_CWM_obs$SLA[group_extra2], shuffle_CWM_trait_based_extra$SLA[group_extra2])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_extra2], shuffle_CWM_trait_based_extra$PLH[group_extra2])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_extra3], shuffle_CWM_trait_based_extra$LNC[group_extra3])^2,
                                                              cor(shuffle_CWM_obs$SLA[group_extra3], shuffle_CWM_trait_based_extra$SLA[group_extra3])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_extra3], shuffle_CWM_trait_based_extra$PLH[group_extra3])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_extra4], shuffle_CWM_trait_based_extra$LNC[group_extra4])^2,
                                                              cor(shuffle_CWM_obs$SLA[group_extra4], shuffle_CWM_trait_based_extra$SLA[group_extra4])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_extra4], shuffle_CWM_trait_based_extra$PLH[group_extra4])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_inter1], shuffle_FDis_trait_based_inter$LNC[group_inter1])^2,# FDis trait based approach interpolation
                                                              cor(shuffle_FDis_obs$SLA[group_inter1], shuffle_FDis_trait_based_inter$SLA[group_inter1])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_inter1], shuffle_FDis_trait_based_inter$PLH[group_inter1])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_inter2], shuffle_FDis_trait_based_inter$LNC[group_inter2])^2,
                                                              cor(shuffle_FDis_obs$SLA[group_inter2], shuffle_FDis_trait_based_inter$SLA[group_inter2])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_inter2], shuffle_FDis_trait_based_inter$PLH[group_inter2])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_inter3], shuffle_FDis_trait_based_inter$LNC[group_inter3])^2,
                                                              cor(shuffle_FDis_obs$SLA[group_inter3], shuffle_FDis_trait_based_inter$SLA[group_inter3])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_inter3], shuffle_FDis_trait_based_inter$PLH[group_inter3])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_inter4], shuffle_FDis_trait_based_inter$LNC[group_inter4])^2,
                                                              cor(shuffle_FDis_obs$SLA[group_inter4], shuffle_FDis_trait_based_inter$SLA[group_inter4])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_inter4], shuffle_FDis_trait_based_inter$PLH[group_inter4])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_extra1], shuffle_FDis_trait_based_extra$LNC[group_extra1])^2, # FDis trait based approach extrapolation
                                                              cor(shuffle_FDis_obs$SLA[group_extra1], shuffle_FDis_trait_based_extra$SLA[group_extra1])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_extra1], shuffle_FDis_trait_based_extra$PLH[group_extra1])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_extra2], shuffle_FDis_trait_based_extra$LNC[group_extra2])^2,
                                                              cor(shuffle_FDis_obs$SLA[group_extra2], shuffle_FDis_trait_based_extra$SLA[group_extra2])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_extra2], shuffle_FDis_trait_based_extra$PLH[group_extra2])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_extra3], shuffle_FDis_trait_based_extra$LNC[group_extra3])^2,
                                                              cor(shuffle_FDis_obs$SLA[group_extra3], shuffle_FDis_trait_based_extra$SLA[group_extra3])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_extra3], shuffle_FDis_trait_based_extra$PLH[group_extra3])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_extra4], shuffle_FDis_trait_based_extra$LNC[group_extra4])^2,
                                                              cor(shuffle_FDis_obs$SLA[group_extra4], shuffle_FDis_trait_based_extra$SLA[group_extra4])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_extra4], shuffle_FDis_trait_based_extra$PLH[group_extra4])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_inter1], shuffle_CWM_species_based_inter$LNC[group_inter1])^2, # CWM species based approach interpolation
                                                              cor(shuffle_CWM_obs$SLA[group_inter1], shuffle_CWM_species_based_inter$SLA[group_inter1])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_inter1], shuffle_CWM_species_based_inter$PLH[group_inter1])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_inter2], shuffle_CWM_species_based_inter$LNC[group_inter2])^2,
                                                              cor(shuffle_CWM_obs$SLA[group_inter2], shuffle_CWM_species_based_inter$SLA[group_inter2])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_inter2], shuffle_CWM_species_based_inter$PLH[group_inter2])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_inter3], shuffle_CWM_species_based_inter$LNC[group_inter3])^2,
                                                              cor(shuffle_CWM_obs$SLA[group_inter3], shuffle_CWM_species_based_inter$SLA[group_inter3])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_inter3], shuffle_CWM_species_based_inter$PLH[group_inter3])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_inter4], shuffle_CWM_species_based_inter$LNC[group_inter4])^2,
                                                              cor(shuffle_CWM_obs$SLA[group_inter4], shuffle_CWM_species_based_inter$SLA[group_inter4])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_inter4], shuffle_CWM_species_based_inter$PLH[group_inter4])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_extra1], shuffle_CWM_species_based_extra$LNC[group_extra1])^2, # CWM species based approach extrapolation
                                                              cor(shuffle_CWM_obs$SLA[group_extra1], shuffle_CWM_species_based_extra$SLA[group_extra1])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_extra1], shuffle_CWM_species_based_extra$PLH[group_extra1])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_extra2], shuffle_CWM_species_based_extra$LNC[group_extra2])^2,
                                                              cor(shuffle_CWM_obs$SLA[group_extra2], shuffle_CWM_species_based_extra$SLA[group_extra2])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_extra2], shuffle_CWM_species_based_extra$PLH[group_extra2])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_extra3], shuffle_CWM_species_based_extra$LNC[group_extra3])^2,
                                                              cor(shuffle_CWM_obs$SLA[group_extra3], shuffle_CWM_species_based_extra$SLA[group_extra3])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_extra3], shuffle_CWM_species_based_extra$PLH[group_extra3])^2,
                                                              cor(shuffle_CWM_obs$LNC[group_extra4], shuffle_CWM_species_based_extra$LNC[group_extra4])^2,
                                                              cor(shuffle_CWM_obs$SLA[group_extra4], shuffle_CWM_species_based_extra$SLA[group_extra4])^2,
                                                              cor(shuffle_CWM_obs$PLH[group_extra4], shuffle_CWM_species_based_extra$PLH[group_extra4])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_inter1], shuffle_FDis_species_based_inter$LNC[group_inter1])^2, # FDis species based approach interpolation
                                                              cor(shuffle_FDis_obs$SLA[group_inter1], shuffle_FDis_species_based_inter$SLA[group_inter1])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_inter1], shuffle_FDis_species_based_inter$PLH[group_inter1])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_inter2], shuffle_FDis_species_based_inter$LNC[group_inter2])^2,
                                                              cor(shuffle_FDis_obs$SLA[group_inter2], shuffle_FDis_species_based_inter$SLA[group_inter2])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_inter2], shuffle_FDis_species_based_inter$PLH[group_inter2])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_inter3], shuffle_FDis_species_based_inter$LNC[group_inter3])^2,
                                                              cor(shuffle_FDis_obs$SLA[group_inter3], shuffle_FDis_species_based_inter$SLA[group_inter3])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_inter3], shuffle_FDis_species_based_inter$PLH[group_inter3])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_inter4], shuffle_FDis_species_based_inter$LNC[group_inter4])^2,
                                                              cor(shuffle_FDis_obs$SLA[group_inter4], shuffle_FDis_species_based_inter$SLA[group_inter4])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_inter4], shuffle_FDis_species_based_inter$PLH[group_inter4])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_extra1], shuffle_FDis_species_based_extra$LNC[group_extra1])^2, # FDis species based approach extrapolation
                                                              cor(shuffle_FDis_obs$SLA[group_extra1], shuffle_FDis_species_based_extra$SLA[group_extra1])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_extra1], shuffle_FDis_species_based_extra$PLH[group_extra1])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_extra2], shuffle_FDis_species_based_extra$LNC[group_extra2])^2,
                                                              cor(shuffle_FDis_obs$SLA[group_extra2], shuffle_FDis_species_based_extra$SLA[group_extra2])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_extra2], shuffle_FDis_species_based_extra$PLH[group_extra2])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_extra3], shuffle_FDis_species_based_extra$LNC[group_extra3])^2,
                                                              cor(shuffle_FDis_obs$SLA[group_extra3], shuffle_FDis_species_based_extra$SLA[group_extra3])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_extra3], shuffle_FDis_species_based_extra$PLH[group_extra3])^2,
                                                              cor(shuffle_FDis_obs$LNC[group_extra4], shuffle_FDis_species_based_extra$LNC[group_extra4])^2,
                                                              cor(shuffle_FDis_obs$SLA[group_extra4], shuffle_FDis_species_based_extra$SLA[group_extra4])^2,
                                                              cor(shuffle_FDis_obs$PLH[group_extra4], shuffle_FDis_species_based_extra$PLH[group_extra4])^2,
                                                              
                                                              cor(shuffle_CM_obs$LNC[group_inter1], shuffle_CM_trait_based_inter$LNC[group_inter1])^2, # CM trait based approach interpolation
                                                              cor(shuffle_CM_obs$SLA[group_inter1], shuffle_CM_trait_based_inter$SLA[group_inter1])^2,
                                                              cor(shuffle_CM_obs$PLH[group_inter1], shuffle_CM_trait_based_inter$PLH[group_inter1])^2,
                                                              cor(shuffle_CM_obs$LNC[group_inter2], shuffle_CM_trait_based_inter$LNC[group_inter2])^2,
                                                              cor(shuffle_CM_obs$SLA[group_inter2], shuffle_CM_trait_based_inter$SLA[group_inter2])^2,
                                                              cor(shuffle_CM_obs$PLH[group_inter2], shuffle_CM_trait_based_inter$PLH[group_inter2])^2,
                                                              cor(shuffle_CM_obs$LNC[group_inter3], shuffle_CM_trait_based_inter$LNC[group_inter3])^2,
                                                              cor(shuffle_CM_obs$SLA[group_inter3], shuffle_CM_trait_based_inter$SLA[group_inter3])^2,
                                                              cor(shuffle_CM_obs$PLH[group_inter3], shuffle_CM_trait_based_inter$PLH[group_inter3])^2,
                                                              cor(shuffle_CM_obs$LNC[group_inter4], shuffle_CM_trait_based_inter$LNC[group_inter4])^2,
                                                              cor(shuffle_CM_obs$SLA[group_inter4], shuffle_CM_trait_based_inter$SLA[group_inter4])^2,
                                                              cor(shuffle_CM_obs$PLH[group_inter4], shuffle_CM_trait_based_inter$PLH[group_inter4])^2,
                                                              cor(shuffle_CM_obs$LNC[group_extra1], shuffle_CM_trait_based_extra$LNC[group_extra1])^2,# CM trait based approach interpolation
                                                              cor(shuffle_CM_obs$SLA[group_extra1], shuffle_CM_trait_based_extra$SLA[group_extra1])^2,
                                                              cor(shuffle_CM_obs$PLH[group_extra1], shuffle_CM_trait_based_extra$PLH[group_extra1])^2,
                                                              cor(shuffle_CM_obs$LNC[group_extra2], shuffle_CM_trait_based_extra$LNC[group_extra2])^2,
                                                              cor(shuffle_CM_obs$SLA[group_extra2], shuffle_CM_trait_based_extra$SLA[group_extra2])^2,
                                                              cor(shuffle_CM_obs$PLH[group_extra2], shuffle_CM_trait_based_extra$PLH[group_extra2])^2,
                                                              cor(shuffle_CM_obs$LNC[group_extra3], shuffle_CM_trait_based_extra$LNC[group_extra3])^2,
                                                              cor(shuffle_CM_obs$SLA[group_extra3], shuffle_CM_trait_based_extra$SLA[group_extra3])^2,
                                                              cor(shuffle_CM_obs$PLH[group_extra3], shuffle_CM_trait_based_extra$PLH[group_extra3])^2,
                                                              cor(shuffle_CM_obs$LNC[group_extra4], shuffle_CM_trait_based_extra$LNC[group_extra4])^2,
                                                              cor(shuffle_CM_obs$SLA[group_extra4], shuffle_CM_trait_based_extra$SLA[group_extra4])^2,
                                                              cor(shuffle_CM_obs$PLH[group_extra4], shuffle_CM_trait_based_extra$PLH[group_extra4])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_inter1], shuffle_uFDis_trait_based_inter$LNC[group_inter1])^2,# uuFDis trait based approach interpolation
                                                              cor(shuffle_uFDis_obs$SLA[group_inter1], shuffle_uFDis_trait_based_inter$SLA[group_inter1])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_inter1], shuffle_uFDis_trait_based_inter$PLH[group_inter1])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_inter2], shuffle_uFDis_trait_based_inter$LNC[group_inter2])^2,
                                                              cor(shuffle_uFDis_obs$SLA[group_inter2], shuffle_uFDis_trait_based_inter$SLA[group_inter2])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_inter2], shuffle_uFDis_trait_based_inter$PLH[group_inter2])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_inter3], shuffle_uFDis_trait_based_inter$LNC[group_inter3])^2,
                                                              cor(shuffle_uFDis_obs$SLA[group_inter3], shuffle_uFDis_trait_based_inter$SLA[group_inter3])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_inter3], shuffle_uFDis_trait_based_inter$PLH[group_inter3])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_inter4], shuffle_uFDis_trait_based_inter$LNC[group_inter4])^2,
                                                              cor(shuffle_uFDis_obs$SLA[group_inter4], shuffle_uFDis_trait_based_inter$SLA[group_inter4])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_inter4], shuffle_uFDis_trait_based_inter$PLH[group_inter4])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_extra1], shuffle_uFDis_trait_based_extra$LNC[group_extra1])^2, # uuFDis trait based approach extrapolation
                                                              cor(shuffle_uFDis_obs$SLA[group_extra1], shuffle_uFDis_trait_based_extra$SLA[group_extra1])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_extra1], shuffle_uFDis_trait_based_extra$PLH[group_extra1])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_extra2], shuffle_uFDis_trait_based_extra$LNC[group_extra2])^2,
                                                              cor(shuffle_uFDis_obs$SLA[group_extra2], shuffle_uFDis_trait_based_extra$SLA[group_extra2])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_extra2], shuffle_uFDis_trait_based_extra$PLH[group_extra2])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_extra3], shuffle_uFDis_trait_based_extra$LNC[group_extra3])^2,
                                                              cor(shuffle_uFDis_obs$SLA[group_extra3], shuffle_uFDis_trait_based_extra$SLA[group_extra3])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_extra3], shuffle_uFDis_trait_based_extra$PLH[group_extra3])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_extra4], shuffle_uFDis_trait_based_extra$LNC[group_extra4])^2,
                                                              cor(shuffle_uFDis_obs$SLA[group_extra4], shuffle_uFDis_trait_based_extra$SLA[group_extra4])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_extra4], shuffle_uFDis_trait_based_extra$PLH[group_extra4])^2,
                                                              cor(shuffle_CM_obs$LNC[group_inter1], shuffle_CM_species_based_inter$LNC[group_inter1])^2, # CM species based approach interpolation
                                                              cor(shuffle_CM_obs$SLA[group_inter1], shuffle_CM_species_based_inter$SLA[group_inter1])^2,
                                                              cor(shuffle_CM_obs$PLH[group_inter1], shuffle_CM_species_based_inter$PLH[group_inter1])^2,
                                                              cor(shuffle_CM_obs$LNC[group_inter2], shuffle_CM_species_based_inter$LNC[group_inter2])^2,
                                                              cor(shuffle_CM_obs$SLA[group_inter2], shuffle_CM_species_based_inter$SLA[group_inter2])^2,
                                                              cor(shuffle_CM_obs$PLH[group_inter2], shuffle_CM_species_based_inter$PLH[group_inter2])^2,
                                                              cor(shuffle_CM_obs$LNC[group_inter3], shuffle_CM_species_based_inter$LNC[group_inter3])^2,
                                                              cor(shuffle_CM_obs$SLA[group_inter3], shuffle_CM_species_based_inter$SLA[group_inter3])^2,
                                                              cor(shuffle_CM_obs$PLH[group_inter3], shuffle_CM_species_based_inter$PLH[group_inter3])^2,
                                                              cor(shuffle_CM_obs$LNC[group_inter4], shuffle_CM_species_based_inter$LNC[group_inter4])^2,
                                                              cor(shuffle_CM_obs$SLA[group_inter4], shuffle_CM_species_based_inter$SLA[group_inter4])^2,
                                                              cor(shuffle_CM_obs$PLH[group_inter4], shuffle_CM_species_based_inter$PLH[group_inter4])^2,
                                                              cor(shuffle_CM_obs$LNC[group_extra1], shuffle_CM_species_based_extra$LNC[group_extra1])^2, # CM species based approach extrapolation
                                                              cor(shuffle_CM_obs$SLA[group_extra1], shuffle_CM_species_based_extra$SLA[group_extra1])^2,
                                                              cor(shuffle_CM_obs$PLH[group_extra1], shuffle_CM_species_based_extra$PLH[group_extra1])^2,
                                                              cor(shuffle_CM_obs$LNC[group_extra2], shuffle_CM_species_based_extra$LNC[group_extra2])^2,
                                                              cor(shuffle_CM_obs$SLA[group_extra2], shuffle_CM_species_based_extra$SLA[group_extra2])^2,
                                                              cor(shuffle_CM_obs$PLH[group_extra2], shuffle_CM_species_based_extra$PLH[group_extra2])^2,
                                                              cor(shuffle_CM_obs$LNC[group_extra3], shuffle_CM_species_based_extra$LNC[group_extra3])^2,
                                                              cor(shuffle_CM_obs$SLA[group_extra3], shuffle_CM_species_based_extra$SLA[group_extra3])^2,
                                                              cor(shuffle_CM_obs$PLH[group_extra3], shuffle_CM_species_based_extra$PLH[group_extra3])^2,
                                                              cor(shuffle_CM_obs$LNC[group_extra4], shuffle_CM_species_based_extra$LNC[group_extra4])^2,
                                                              cor(shuffle_CM_obs$SLA[group_extra4], shuffle_CM_species_based_extra$SLA[group_extra4])^2,
                                                              cor(shuffle_CM_obs$PLH[group_extra4], shuffle_CM_species_based_extra$PLH[group_extra4])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_inter1], shuffle_uFDis_species_based_inter$LNC[group_inter1])^2, # uuFDis species based approach interpolation
                                                              cor(shuffle_uFDis_obs$SLA[group_inter1], shuffle_uFDis_species_based_inter$SLA[group_inter1])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_inter1], shuffle_uFDis_species_based_inter$PLH[group_inter1])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_inter2], shuffle_uFDis_species_based_inter$LNC[group_inter2])^2,
                                                              cor(shuffle_uFDis_obs$SLA[group_inter2], shuffle_uFDis_species_based_inter$SLA[group_inter2])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_inter2], shuffle_uFDis_species_based_inter$PLH[group_inter2])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_inter3], shuffle_uFDis_species_based_inter$LNC[group_inter3])^2,
                                                              cor(shuffle_uFDis_obs$SLA[group_inter3], shuffle_uFDis_species_based_inter$SLA[group_inter3])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_inter3], shuffle_uFDis_species_based_inter$PLH[group_inter3])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_inter4], shuffle_uFDis_species_based_inter$LNC[group_inter4])^2,
                                                              cor(shuffle_uFDis_obs$SLA[group_inter4], shuffle_uFDis_species_based_inter$SLA[group_inter4])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_inter4], shuffle_uFDis_species_based_inter$PLH[group_inter4])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_extra1], shuffle_uFDis_species_based_extra$LNC[group_extra1])^2, # uuFDis species based approach extrapolation
                                                              cor(shuffle_uFDis_obs$SLA[group_extra1], shuffle_uFDis_species_based_extra$SLA[group_extra1])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_extra1], shuffle_uFDis_species_based_extra$PLH[group_extra1])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_extra2], shuffle_uFDis_species_based_extra$LNC[group_extra2])^2,
                                                              cor(shuffle_uFDis_obs$SLA[group_extra2], shuffle_uFDis_species_based_extra$SLA[group_extra2])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_extra2], shuffle_uFDis_species_based_extra$PLH[group_extra2])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_extra3], shuffle_uFDis_species_based_extra$LNC[group_extra3])^2,
                                                              cor(shuffle_uFDis_obs$SLA[group_extra3], shuffle_uFDis_species_based_extra$SLA[group_extra3])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_extra3], shuffle_uFDis_species_based_extra$PLH[group_extra3])^2,
                                                              cor(shuffle_uFDis_obs$LNC[group_extra4], shuffle_uFDis_species_based_extra$LNC[group_extra4])^2,
                                                              cor(shuffle_uFDis_obs$SLA[group_extra4], shuffle_uFDis_species_based_extra$SLA[group_extra4])^2,
                                                              cor(shuffle_uFDis_obs$PLH[group_extra4], shuffle_uFDis_species_based_extra$PLH[group_extra4])^2),
                                                       "Traits"=rep(c(rep(c("LNC", "SLA", "PLH"), 4*8)),2),
                                                       "Indices"=c(c(rep(c(rep("CWM", 24), rep("FDis", 24)),2)),c(rep(c(rep("CM", 24), rep("uFDis", 24)),2))),
                                                       "Approaches"=rep(c(rep("Assemble-First Approach", 12*4), rep("Predict-First Approach", 12*4)),2),
                                                       "Crossvalidation"=rep(c(rep(c(rep("Interpolation", 12), rep("Extrapolation", 12)),4)),2),
                                                       "Type"=c(rep("Abundance data", 4*8*3), rep("Occurrence data", 4*8*3)),
                                                       "repet"=rep(c(rep(1,3), rep(2,3), rep(3,3), rep(4,3)), 16),
                                                       "alpha"=rep(0.1, 3*4*8*2))
  )
  
  shuffle_figure3$Traits <- fct_relevel(shuffle_figure3$Traits, c("PLH", "SLA", "LNC"))
  shuffle_figure3$Crossvalidation <- fct_relevel(shuffle_figure3$Crossvalidation, c("Interpolation", "Extrapolation"))
  shuffle_figure3$Type <- fct_relevel(shuffle_figure3$Type, c("Abundance data","Occurrence data"))
  shuffle_figure3$Approaches <- fct_relevel(shuffle_figure3$Approaches, c("Predict-First Approach","Assemble-First Approach"))
  shuffle_figure3$Indices <- fct_relevel(shuffle_figure3$Indices, c("CWM","FDis", "CM", 'uFDis'))
}

### Figure 3 :
fig_3 <-figure_3%>%
  ggplot() +
  aes(x = Crossvalidation, y = R2, fill = Approaches) +
  scale_fill_manual(values=c("grey", "white"))+
  geom_boxplot() +
  theme_bw() +
  facet_nested(Traits ~ Type+Indices, scales = "free_y")+
  theme(axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 10))
fig_3 +
  geom_boxplot(data=tot_shuffle, mapping=aes(x=Crossvalidation, y=R2, fill=Approaches, color=Approaches, alpha=1)) +
  theme_bw() +
  scale_fill_manual(values=c("grey", "white"))+
  scale_color_manual(values=c("grey", "grey"))+
  facet_nested(Traits ~ Type+Indices , scales = "free_y")


#######################################################################################
#### -- INFLUENCE OF SPECIES OVERPREDICTION ON THE CALCULATION OF FUNCTIONAL INDICES --  ####
#######################################################################################
### -- Interpolation for abundance data ####
# For each site, species are ranked according to their predicted relative cover
order_species_recovery_inter <- plant_recovery
for (i in 1 : 4463){ 
  order_species_recovery_inter[i,] <- order(-species_models_recovery_inter[i,])
  print(i)
}

# calculation of functional indices and their correlation to those observed for each species addition
plant_recovery_inter_fig4 <- data.frame(matrix(0, ncol=831, nrow=4463), row.names = rownames(plant_recovery))
colnames(plant_recovery_inter_fig4)=colnames(plant_recovery)
R2_CWM_inter_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(R2_CWM_inter_fig4)<- c("LNC","SLA","PLH")
R2_FDis_inter_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(R2_FDis_inter_fig4)<- c("LNC","SLA","PLH")
RMSE_CWM_inter_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(RMSE_CWM_inter_fig4)<- c("LNC","SLA","PLH")
RMSE_FDis_inter_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(RMSE_FDis_inter_fig4)<- c("LNC","SLA","PLH")
for (m in 1:831){ ## One new species is added each time
  for (i in 1 : 4463){ ## we add the new species in the botanical survey for each site
    plant_recovery_inter_fig4[i,as.numeric(order_species_recovery_inter[i,m])] = species_models_recovery_inter[i,as.numeric(order_species_recovery_inter[i,m])]
  }
  print(m)
  CWM_inter_fig4 <- CWM(plant_recovery_inter_fig4)
  CWM_inter_fig4 <- left_join(CWM_inter_fig4, CWM_obs, by="Site")
  R2_CWM_inter_fig4[m,] <- c(cor(CWM_inter_fig4$LNC.x, CWM_inter_fig4$LNC.y)^2, cor(CWM_inter_fig4$SLA.x, CWM_inter_fig4$SLA.y)^2, cor(CWM_inter_fig4$PLH.x, CWM_inter_fig4$PLH.y)^2)
  RMSE_CWM_inter_fig4[m,] <- c(rmse(CWM_inter_fig4$LNC.x, CWM_inter_fig4$LNC.y), rmse(CWM_inter_fig4$SLA.x, CWM_inter_fig4$SLA.y), rmse(CWM_inter_fig4$PLH.x, CWM_inter_fig4$PLH.y))
  
  FDis_inter_fig4 <- FDis(plant_recovery_inter_fig4)
  FDis_inter_fig4 <- left_join(FDis_inter_fig4, FDis_obs, by="Site")
  R2_FDis_inter_fig4[m,] <- c(cor(FDis_inter_fig4$LNC.x, FDis_inter_fig4$LNC.y)^2, cor(FDis_inter_fig4$SLA.x, FDis_inter_fig4$SLA.y)^2, cor(FDis_inter_fig4$PLH.x, FDis_inter_fig4$PLH.y)^2)
  RMSE_FDis_inter_fig4[m,] <- c(rmse(FDis_inter_fig4$LNC.x, FDis_inter_fig4$LNC.y), rmse(FDis_inter_fig4$SLA.x, FDis_inter_fig4$SLA.y), rmse(FDis_inter_fig4$PLH.x, FDis_inter_fig4$PLH.y))
}
R2_CWM_inter_fig4 <- R2_CWM_inter_fig4[-1,]
R2_FDis_inter_fig4 <- R2_FDis_inter_fig4[-1,]
RMSE_CWM_inter_fig4 <- RMSE_CWM_inter_fig4[-1,]
RMSE_FDis_inter_fig4 <- RMSE_FDis_inter_fig4[-1,]

### -- Extrapolation for abundance data ####
# For each site, species are ranked according to their predicted relative cover
order_species_recovery_extra <- plant_recovery
for (i in 1 : 4463){ 
  order_species_recovery_extra[i,] <- order(-species_models_recovery_extra[i,])
}

# calculation of functional indices and their correlation to those observed for each species addition
plant_recovery_extra_fig4 <- data.frame(matrix(0, ncol=831, nrow=4463), row.names = rownames(plant_recovery))
colnames(plant_recovery_extra_fig4)=colnames(plant_recovery)
R2_CWM_extra_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(R2_CWM_extra_fig4)<- c("LNC","SLA","PLH")
R2_FDis_extra_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(R2_FDis_extra_fig4)<- c("LNC","SLA","PLH")
RMSE_CWM_extra_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(RMSE_CWM_extra_fig4)<- c("LNC","SLA","PLH")
RMSE_FDis_extra_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(RMSE_FDis_extra_fig4)<- c("LNC","SLA","PLH")
for (m in 1:831){ ## One new species is added each time
  for (i in 1 : 4463){ ## we add the new species in the botanical survey for each site
    plant_recovery_extra_fig4[i,as.numeric(order_species_recovery_extra[i,m])] = species_models_recovery_extra[i,as.numeric(order_species_recovery_extra[i,m])]
  }
  CWM_extra_fig4 <- CWM(plant_recovery_extra_fig4)
  CWM_extra_fig4 <- left_join(CWM_extra_fig4, CWM_obs, by="Site")
  R2_CWM_extra_fig4[m,] <- c(cor(CWM_extra_fig4$LNC.x, CWM_extra_fig4$LNC.y)^2, cor(CWM_extra_fig4$SLA.x, CWM_extra_fig4$SLA.y)^2, cor(CWM_extra_fig4$PLH.x, CWM_extra_fig4$PLH.y)^2)
  RMSE_CWM_extra_fig4[m,] <- c(rmse(CWM_extra_fig4$LNC.x, CWM_extra_fig4$LNC.y), rmse(CWM_extra_fig4$SLA.x, CWM_extra_fig4$SLA.y), rmse(CWM_extra_fig4$PLH.x, CWM_extra_fig4$PLH.y))
  

  FDis_extra_fig4 <- FDis(plant_recovery_extra_fig4)
  FDis_extra_fig4 <- left_join(FDis_extra_fig4, FDis_obs, by="Site")
  R2_FDis_extra_fig4[m,] <- c(cor(FDis_extra_fig4$LNC.x, FDis_extra_fig4$LNC.y)^2, cor(FDis_extra_fig4$SLA.x, FDis_extra_fig4$SLA.y)^2, cor(FDis_extra_fig4$PLH.x, FDis_extra_fig4$PLH.y)^2)
  RMSE_FDis_extra_fig4[m,] <- c(rmse(FDis_extra_fig4$LNC.x, FDis_extra_fig4$LNC.y), rmse(FDis_extra_fig4$SLA.x, FDis_extra_fig4$SLA.y), rmse(FDis_extra_fig4$PLH.x, FDis_extra_fig4$PLH.y))
}
R2_CWM_extra_fig4 <- R2_CWM_extra_fig4[-1,]
R2_FDis_extra_fig4 <- R2_FDis_extra_fig4[-1,]
RMSE_CWM_extra_fig4 <- RMSE_CWM_extra_fig4[-1,]
RMSE_FDis_extra_fig4 <- RMSE_FDis_extra_fig4[-1,]
### -- Interpolation for occurrence data ####
# For each site, species are ranked according to their predicted probability of presence
order_species_pa_inter <- plant_pa
for (i in 1 : 4463){ 
  order_species_pa_inter[i,] <- order(-species_models_pa_inter[i,])
}

# calculation of functional indices and their correlation to those observed for each species addition
plant_pa_inter_fig4 <- data.frame(matrix(0, ncol=831, nrow=4463), row.names = rownames(plant_pa))
colnames(plant_pa_inter_fig4)=colnames(plant_pa)
R2_CM_inter_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(R2_CM_inter_fig4)<- c("LNC","SLA","PLH")
R2_uFDis_inter_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(R2_uFDis_inter_fig4)<- c("LNC","SLA","PLH")
RMSE_CM_inter_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(RMSE_CM_inter_fig4)<- c("LNC","SLA","PLH")
RMSE_uFDis_inter_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(RMSE_uFDis_inter_fig4)<- c("LNC","SLA","PLH")

for (m in 1:831){ ## One new species is added each time
  for (i in 1 : 4463){ ## we add the new species in the botanical survey for each site
    plant_pa_inter_fig4[i,as.numeric(order_species_pa_inter[i,m])] = species_models_pa_inter[i,as.numeric(order_species_pa_inter[i,m])]
  }
  print(m)
  CM_inter_fig4 <- CWM(plant_pa_inter_fig4)
  CM_inter_fig4 <- left_join(CM_inter_fig4, CM_obs, by="Site")
  R2_CM_inter_fig4[m,] <- c(cor(CM_inter_fig4$LNC.x, CM_inter_fig4$LNC.y)^2, cor(CM_inter_fig4$SLA.x, CM_inter_fig4$SLA.y)^2, cor(CM_inter_fig4$PLH.x, CM_inter_fig4$PLH.y)^2)
  RMSE_CM_inter_fig4[m,] <- c(rmse(CM_inter_fig4$LNC.x, CM_inter_fig4$LNC.y), rmse(CM_inter_fig4$SLA.x, CM_inter_fig4$SLA.y), rmse(CM_inter_fig4$PLH.x, CM_inter_fig4$PLH.y))
  
  uFDis_inter_fig4 <- FDis(plant_pa_inter_fig4)
  uFDis_inter_fig4 <- left_join(uFDis_inter_fig4, uFDis_obs, by="Site")
  R2_uFDis_inter_fig4[m,] <- c(cor(uFDis_inter_fig4$LNC.x, uFDis_inter_fig4$LNC.y)^2, cor(uFDis_inter_fig4$SLA.x, uFDis_inter_fig4$SLA.y)^2, cor(uFDis_inter_fig4$PLH.x, uFDis_inter_fig4$PLH.y)^2)
  RMSE_uFDis_inter_fig4[m,] <- c(rmse(uFDis_inter_fig4$LNC.x, uFDis_inter_fig4$LNC.y), rmse(uFDis_inter_fig4$SLA.x, uFDis_inter_fig4$SLA.y), rmse(uFDis_inter_fig4$PLH.x, uFDis_inter_fig4$PLH.y))
}
R2_CM_inter_fig4 <- R2_CM_inter_fig4[-1,]
R2_uFDis_inter_fig4 <- R2_uFDis_inter_fig4[-1,]

RMSE_CM_inter_fig4 <- RMSE_CM_inter_fig4[-1,]
RMSE_uFDis_inter_fig4 <- RMSE_uFDis_inter_fig4[-1,]

### -- Extrapolation for occurrence data ####
# For each site, species are ranked according to their predicted probability of presence
order_species_pa_extra <- plant_pa
for (i in 1 : 4463){ 
  order_species_pa_extra[i,] <- order(-species_models_pa_extra[i,])
}

# calculation of functional indices and their correlation to those observed for each species addition
plant_pa_extra_fig4 <- data.frame(matrix(0, ncol=831, nrow=4463), row.names = rownames(plant_pa))
colnames(plant_pa_extra_fig4)=colnames(plant_pa)
R2_CM_extra_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(R2_CM_extra_fig4)<- c("LNC","SLA","PLH")
R2_uFDis_extra_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(R2_uFDis_extra_fig4)<- c("LNC","SLA","PLH")
RMSE_CM_extra_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(RMSE_CM_extra_fig4)<- c("LNC","SLA","PLH")
RMSE_uFDis_extra_fig4 <- data.frame(matrix(0, nrow=831, ncol=3), row.names=1:831)
colnames(RMSE_uFDis_extra_fig4)<- c("LNC","SLA","PLH")
for (m in 1:831){ ## One new species is added each time
  for (i in 1 : 4463){ ## we add the new species in the botanical survey for each site
    plant_pa_extra_fig4[i,as.numeric(order_species_pa_extra[i,m])] = species_models_pa_extra[i,as.numeric(order_species_pa_extra[i,m])]
  }
  CM_extra_fig4 <- CWM(plant_pa_extra_fig4)
  CM_extra_fig4 <- left_join(CM_extra_fig4, CM_obs, by="Site")
  R2_CM_extra_fig4[m,] <- c(cor(CM_extra_fig4$LNC.x, CM_extra_fig4$LNC.y)^2, cor(CM_extra_fig4$SLA.x, CM_extra_fig4$SLA.y)^2, cor(CM_extra_fig4$PLH.x, CM_extra_fig4$PLH.y)^2)
  RMSE_CM_extra_fig4[m,] <- c(rmse(CM_extra_fig4$LNC.x, CM_extra_fig4$LNC.y), rmse(CM_extra_fig4$SLA.x, CM_extra_fig4$SLA.y), rmse(CM_extra_fig4$PLH.x, CM_extra_fig4$PLH.y))
  
  uFDis_extra_fig4 <- FDis(plant_pa_extra_fig4)
  uFDis_extra_fig4 <- left_join(uFDis_extra_fig4, uFDis_obs, by="Site")
  R2_uFDis_extra_fig4[m,] <- c(cor(uFDis_extra_fig4$LNC.x, uFDis_extra_fig4$LNC.y)^2, cor(uFDis_extra_fig4$SLA.x, uFDis_extra_fig4$SLA.y)^2, cor(uFDis_extra_fig4$PLH.x, uFDis_extra_fig4$PLH.y)^2)
  RMSE_uFDis_extra_fig4[m,] <- c(rmse(uFDis_extra_fig4$LNC.x, uFDis_extra_fig4$LNC.y), rmse(uFDis_extra_fig4$SLA.x, uFDis_extra_fig4$SLA.y), rmse(uFDis_extra_fig4$PLH.x, uFDis_extra_fig4$PLH.y))
}
R2_CM_extra_fig4 <- R2_CM_extra_fig4[-1,]
R2_uFDis_extra_fig4 <- R2_uFDis_extra_fig4[-1,]
RMSE_CM_extra_fig4 <- RMSE_CM_extra_fig4[-1,]
RMSE_uFDis_extra_fig4 <- RMSE_uFDis_extra_fig4[-1,]

### -- Figure 4 ####
figure_4 <- data.frame("Number of species" =rep(2:831,12*2),
                       "R2"=c(R2_CWM_inter_fig4$LNC, R2_CWM_inter_fig4$SLA, R2_CWM_inter_fig4$PLH, R2_FDis_inter_fig4$LNC, R2_FDis_inter_fig4$SLA, R2_FDis_inter_fig4$PLH,
                              R2_CWM_extra_fig4$LNC, R2_CWM_extra_fig4$SLA, R2_CWM_extra_fig4$PLH, R2_FDis_extra_fig4$LNC, R2_FDis_extra_fig4$SLA, R2_FDis_extra_fig4$PLH,
                              R2_CM_inter_fig4$LNC, R2_CM_inter_fig4$SLA, R2_CM_inter_fig4$PLH, R2_uFDis_inter_fig4$LNC, R2_uFDis_inter_fig4$SLA, R2_uFDis_inter_fig4$PLH,
                              R2_CM_extra_fig4$LNC, R2_CM_extra_fig4$SLA, R2_CM_extra_fig4$PLH, R2_uFDis_extra_fig4$LNC, R2_uFDis_extra_fig4$SLA, R2_uFDis_extra_fig4$PLH),
                       "Trait"=rep(c(rep(c(rep("LNC", 830), rep("SLA", 830), rep("PLH", 830)),2)),4),
                       "Indices"=rep(c(rep("CWM/CM", 830*3), rep("FDis/uFDis", 830*3)),4), 
                       "Data"=c(rep("Abundance Data", 830*6*2), rep("Occurrence Data", 830*6*2)),
                       "Crossvalidation"=rep(c(rep("Interpolation", 830*6), rep("Extrapolation", 830*6)), 2), 
                       "Trait_based models"=c(rep(cor(CWM_obs$LNC, CWM_trait_based_inter$LNC)^2,830),
                                                  rep(cor(CWM_obs$SLA, CWM_trait_based_inter$SLA)^2,830),
                                                  rep(cor(CWM_obs$PLH, CWM_trait_based_inter$PLH)^2,830),
                                                  rep(cor(FDis_obs$LNC, FDis_trait_based_inter$LNC)^2,830),
                                                  rep(cor(FDis_obs$SLA, FDis_trait_based_inter$SLA)^2,830),
                                                  rep(cor(FDis_obs$PLH, FDis_trait_based_inter$PLH)^2,830),
                                                  rep(cor(CWM_obs$LNC, CWM_trait_based_extra$LNC)^2,830),
                                                  rep(cor(CWM_obs$SLA, CWM_trait_based_extra$SLA)^2,830),
                                                  rep(cor(CWM_obs$PLH, CWM_trait_based_extra$PLH)^2,830),
                                                  rep(cor(FDis_obs$LNC, FDis_trait_based_extra$LNC)^2,830),
                                                  rep(cor(FDis_obs$SLA, FDis_trait_based_extra$SLA)^2,830),
                                                  rep(cor(FDis_obs$PLH, FDis_trait_based_extra$PLH)^2,830),
                                                  rep(cor(CM_obs$LNC, CM_trait_based_inter$LNC)^2,830),
                                                  rep(cor(CM_obs$SLA, CM_trait_based_inter$SLA)^2,830),
                                                  rep(cor(CM_obs$PLH, CM_trait_based_inter$PLH)^2,830),
                                                  rep(cor(uFDis_obs$LNC, uFDis_trait_based_inter$LNC)^2,830),
                                                  rep(cor(uFDis_obs$SLA, uFDis_trait_based_inter$SLA)^2,830),
                                                  rep(cor(uFDis_obs$PLH, uFDis_trait_based_inter$PLH)^2,830),
                                                  rep(cor(CM_obs$LNC, CM_trait_based_extra$LNC)^2,830),
                                                  rep(cor(CM_obs$SLA, CM_trait_based_extra$SLA)^2,830),
                                                  rep(cor(CM_obs$PLH, CM_trait_based_extra$PLH)^2,830),
                                                  rep(cor(uFDis_obs$LNC, uFDis_trait_based_extra$LNC)^2,830),
                                                  rep(cor(uFDis_obs$SLA, uFDis_trait_based_extra$SLA)^2,830),
                                                  rep(cor(uFDis_obs$PLH, uFDis_trait_based_extra$PLH)^2,830)))
figure_4$Trait <- fct_relevel(figure_4$Trait, c("LNC", "SLA", "PLH"))
figure_4$Crossvalidation <- fct_relevel(figure_4$Crossvalidation, c("Interpolation", "Extrapolation"))


figure_4 %>%
  filter(Trait %in% "PLH") %>%
  ggplot() +
  aes(x = Number.of.species, y = R2) +
  scale_x_continuous(trans='log2')+
  coord_trans(x="log2")+
  geom_point(shape = "circle", size = 1, col="black")+
  geom_hline(aes(yintercept = Trait_based.models), col="darkgrey", linetype = "dashed", size=2)+
  theme_bw() +  
  theme(axis.text.x = element_text(size=12),
                      axis.text.y = element_text(size=15),
                      axis.title.x = element_text(size=20),
                      axis.title.y = element_text(size=20),
                      strip.text.x = element_text(size = 16),
                      strip.text.y = element_text(size = 18),
                      legend.title = element_text(size = 16),
                      legend.text = element_text(size = 10))+
  xlab("Number of species")+
  facet_nested(Data ~ Indices + Crossvalidation)

##########################################################################################
#### -- DIFFERENCES BETWEEN PREDICTED AND OBSERVED SPECIES ABUNDANCE AT EACH SITE --  ####
##########################################################################################
### Density of the coverage for each plant observation/prediction according to the distance between CWM and plant height
## distances between the traits of each species and the CWM:
pourcentage_obs <- plant_recovery
for (s in 1 : 4463){
  tot = rowSums(pourcentage_obs)[s]
  pourcentage_obs[s,]=pourcentage_obs[s,]*100/tot
}
pourcentage_pred <- species_models_recovery_inter
for (s in 1 : 4463){
  tot = rowSums(pourcentage_pred)[s]
  pourcentage_pred[s,]=pourcentage_pred[s,]*100/tot
}

# for LNC : 
site_distance_LNC <- pourcentage_obs
ordonne_site <- data.frame("Site"=rownames(site_distance_LNC))
CWM_obs <- left_join(ordonne_site, CWM_obs, by="Site")

for (s in 1:nrow(site_distance_LNC)){
  for(e in 1:ncol(site_distance_LNC)){
    if (!is.na(traits$LNC[e])){
      site_distance_LNC[s,e] <- CWM_obs$LNC[s]-traits$LNC[e]
    }else{
      site_distance_LNC[s,e] <-0
    }
  }
}

distance_LNC <- data.frame("distance LNC"=rep(0, 831*4463), "Frequence_obs"=rep(0, 831*4463), "Frequence_pred" = rep(0, 831*4463))

lim_inf=1
lim_sup = 831

for (s in 1 : 4463){
  distance_LNC[lim_inf:lim_sup,] <- data.frame(site_distance_LNC[s,], as.matrix(pourcentage_obs)[s,], as.matrix(pourcentage_pred)[s,])
  lim_inf = lim_sup+1
  lim_sup = lim_inf+830
}
distance.LNC_obs <- rep(distance_LNC$distance.LNC[1], round(distance_LNC$Frequence_obs[1]))
for (i in 1 : dim(distance_LNC)[1]){
  distance.LNC_obs <- c(distance.LNC_obs, rep(distance_LNC$distance.LNC[i], round(distance_LNC$Frequence_obs[i])))
}
distance.LNC_pred <- rep(distance_LNC$distance.LNC[1], round(distance_LNC$Frequence_pred[1]))
for (i in 1 : dim(distance_LNC)[1]){
  distance.LNC_pred <- c(distance.LNC_pred, rep(distance_LNC$distance.LNC[i], round(distance_LNC$Frequence_pred[i])))
}

# for SLA : 
site_distance_SLA <- pourcentage_obs
ordonne_site <- data.frame("Site"=rownames(site_distance_SLA))
CWM_obs <- left_join(ordonne_site, CWM_obs, by="Site")

for (s in 1:nrow(site_distance_SLA)){
  for(e in 1:ncol(site_distance_SLA)){
    if (!is.na(traits$SLA[e])){
      site_distance_SLA[s,e] <- CWM_obs$SLA[s]-traits$SLA[e]
    }else{
      site_distance_SLA[s,e] <-0
    }
  }
}

distance_SLA <- data.frame("distance SLA"=rep(0, 831*4463), "Frequence_obs"=rep(0, 831*4463), "Frequence_pred" = rep(0, 831*4463))

lim_inf=1
lim_sup = 831

for (s in 1 : 4463){
  distance_SLA[lim_inf:lim_sup,] <- data.frame(site_distance_SLA[s,], as.matrix(pourcentage_obs)[s,], as.matrix(pourcentage_pred)[s,])
  lim_inf = lim_sup+1
  lim_sup = lim_inf+830
}
distance.SLA_obs <- rep(distance_SLA$distance.SLA[1], round(distance_SLA$Frequence_obs[1]))
for (i in 1 : dim(distance_SLA)[1]){
  distance.SLA_obs <- c(distance.SLA_obs, rep(distance_SLA$distance.SLA[i], round(distance_SLA$Frequence_obs[i])))
}
distance.SLA_pred <- rep(distance_SLA$distance.SLA[1], round(distance_SLA$Frequence_pred[1]))
for (i in 1 : dim(distance_SLA)[1]){
  distance.SLA_pred <- c(distance.SLA_pred, rep(distance_SLA$distance.SLA[i], round(distance_SLA$Frequence_pred[i])))
}
# for PLH : 
site_distance_PLH <- pourcentage_obs
ordonne_site <- data.frame("Site"=rownames(site_distance_PLH))
CWM_obs <- left_join(ordonne_site, CWM_obs, by="Site")

for (s in 1:nrow(site_distance_PLH)){
  for(e in 1:ncol(site_distance_PLH)){
    if (!is.na(traits$PLH[e])){
      site_distance_PLH[s,e] <- CWM_obs$PLH[s]-traits$PLH[e]
    }else{
      site_distance_PLH[s,e] <-NA
    }
  }
}

distance_PLH <- data.frame("distance PLH"=rep(0, 831*4463), "Frequence_obs"=rep(0, 831*4463), "Frequence_pred" = rep(0, 831*4463))

lim_inf=1
lim_sup = 831

for (s in 1 : 4463){
  distance_PLH[lim_inf:lim_sup,] <- data.frame(site_distance_PLH[s,], as.matrix(pourcentage_obs)[s,], as.matrix(pourcentage_pred)[s,])
  lim_inf = lim_sup+1
  lim_sup = lim_inf+830
}
distance.PLH_obs <- rep(distance_PLH$distance.PLH[1], round(distance_PLH$Frequence_obs[1]))
for (i in 1 : dim(distance_PLH)[1]){
  distance.PLH_obs <- c(distance.PLH_obs, rep(distance_PLH$distance.PLH[i], round(distance_PLH$Frequence_obs[i])))
}
distance.PLH_pred <- rep(distance_PLH$distance.PLH[1], round(distance_PLH$Frequence_pred[1]))
for (i in 1 : dim(distance_PLH)[1]){
  distance.PLH_pred <- c(distance.PLH_pred, rep(distance_PLH$distance.PLH[i], round(distance_PLH$Frequence_pred[i])))
}

### Figure 5 ####
ggplot() +
  geom_histogram(data=distance.PLH_obs, aes(x = Distance, y = ..density..), position="identity", alpha=0.5, bins=100, fill="grey")+
  geom_density(data=distance.PLH_obs, aes(x = Distance, y = ..density..), alpha=0.6, color="grey")+
  geom_histogram(data=distance.PLH_pred, aes(x = Distance, y = -..density..), position="identity", alpha=0.5, bins=100, fill="#F8766D")+
  geom_density(data=distance.PLH_pred, aes(x = Distance, y = -..density..), alpha=0.6, color="#F8766D")+
  scale_colour_identity()+
  geom_segment(aes(x=-12.4007955, xend=-12.4007955, y=0, yend=Inf, color="grey"))+
  geom_segment(aes(x=11.6661972, xend=11.6661972, y=0, yend=Inf, color="grey"))+
  geom_segment(aes(x=45.6335, xend=45.6335, y=0, yend=Inf, color="grey"), linetype="dashed")+
  geom_segment(aes(x=-124.2556, xend=-124.2556, y=0, yend=Inf, color="grey"), linetype="dashed")+
  geom_segment(aes(x=-16.220440, xend=-16.220440, y=0, yend=-Inf, color="#F8766D"))+
  geom_segment(aes(x=15.221220, xend=15.221220, y=0, yend=-Inf, color="#F8766D"))+
  geom_segment(aes(x=-131.5561, xend=-131.5561, y=0, yend=-Inf, color="#F8766D"), linetype="dashed")+
  geom_segment(aes(x=60.98018, xend=60.98018, y=0, yend=-Inf, color="#F8766D"), linetype="dashed")+
  theme(
    axis.text.x = element_text(size = 12,face="bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(vjust= 1.8, size = 16),
    axis.title.x = element_text(vjust= -0.5, size = 16),
    axis.title = element_text(face = "bold"))+
  xlab("Distance between CWM and plant height")+
  ylab("Density of the coverage \n for each plant observation/prediction")+
  theme_minimal()

### 1000 predictions of null model ####
# Interpolation for abundance data ####
mean_recovery <- colSums(plant_recovery)/4463 # mean abundances of each species in all the dataset
CWM_LNC_null_model_recovery_inter <- data.frame(matrix(0, nrow=1000, ncol=4463))
colnames(CWM_LNC_null_model_recovery_inter) <- rownames(plant_recovery)
CWM_SLA_null_model_recovery_inter <- data.frame(matrix(0, nrow=1000, ncol=4463))
colnames(CWM_SLA_null_model_recovery_inter) <- rownames(plant_recovery)
CWM_PLH_null_model_recovery_inter <- data.frame(matrix(0, nrow=1000, ncol=4463))
colnames(CWM_PLH_null_model_recovery_inter) <- rownames(plant_recovery)

for (j in 1:1000){
  null_model_recovery_inter <- species_models_recovery_inter
  for (i in 1 : 4463){
    random_sampling <- sample(c(1:831), size = 831, replace = FALSE, prob = mean_recovery)
    sampled_species=c(rep(0,831))
    species_predict <- as.numeric(species_models_recovery_inter[i,])
    sampled_species[random_sampling]=species_models_recovery_inter[i,order(-species_predict)]
    null_model_recovery_inter[i,]<- sampled_species
  }
  print(j)
  CWM_null_model_recovery_inter <- CWM(null_model_recovery_inter)
  CWM_LNC_null_model_recovery_inter[j,] <- CWM_null_model_recovery_inter$LNC
  CWM_SLA_null_model_recovery_inter[j,] <- CWM_null_model_recovery_inter$SLA
  CWM_PLH_null_model_recovery_inter[j,] <- CWM_null_model_recovery_inter$PLH
  
  FDis_null_model_recovery_inter <- FDis(null_model_recovery_inter)
  FDis_LNC_null_model_recovery_inter[j,] <- FDis_null_model_recovery_inter$LNC
  FDis_SLA_null_model_recovery_inter[j,] <- FDis_null_model_recovery_inter$SLA
  FDis_PLH_null_model_recovery_inter[j,] <- FDis_null_model_recovery_inter$PLH
}

R2_CWM_null_model_recovery_inter<- data.frame("R2_LNC"=rep(0,1000), "RMSE_LNC"=rep(0,1000),
                                              "R2_SLA"=rep(0,1000), "RMSE_SLA"=rep(0,1000),
                                              "R2_PLH"=rep(0,1000), "RMSE_PLH"=rep(0,1000))
for (i in 1:1000){
  R2_CWM_null_model_recovery_inter[i,] <- c(cor(as.numeric(CWM_LNC_null_model_recovery_inter[i,]), CWM_obs$LNC)^2, rmse(as.numeric(CWM_LNC_null_model_recovery_inter[i,]), CWM_obs$LNC),
                                            cor(as.numeric(CWM_SLA_null_model_recovery_inter[i,]), CWM_obs$SLA)^2, rmse(as.numeric(CWM_SLA_null_model_recovery_inter[i,]), CWM_obs$SLA),
                                            cor(as.numeric(CWM_PLH_null_model_recovery_inter[i,]), CWM_obs$PLH)^2, rmse(as.numeric(CWM_PLH_null_model_recovery_inter[i,]), CWM_obs$PLH))
  print(i)
}

R2_FDis_null_model_recovery_inter <- data.frame("R2_LNC"=rep(0,1000), "RMSE_LNC"=rep(0,1000),
                                                "R2_SLA"=rep(0,1000), "RMSE_SLA"=rep(0,1000),
                                                "R2_PLH"=rep(0,1000), "RMSE_PLH"=rep(0,1000))
for (i in 1:1000){
  R2_FDis_null_model_recovery_inter[i,] <- c(cor(as.numeric(CWM_LNC_null_model_recovery_inter[i,]), CWM_obs$LNC)^2, rmse(as.numeric(CWM_LNC_null_model_recovery_inter[i,]), CWM_obs$LNC),
                                             cor(as.numeric(CWM_SLA_null_model_recovery_inter[i,]), CWM_obs$SLA)^2, rmse(as.numeric(CWM_SLA_null_model_recovery_inter[i,]), CWM_obs$SLA),
                                             cor(as.numeric(CWM_PLH_null_model_recovery_inter[i,]), CWM_obs$PLH)^2, rmse(as.numeric(CWM_PLH_null_model_recovery_inter[i,]), CWM_obs$PLH))
  print(i)
}


# Extrapolation for abundance data ####
CWM_LNC_null_model_recovery_extra <- data.frame(matrix(0, nrow=1000, ncol=4463))
colnames(CWM_LNC_null_model_recovery_extra) <- rownames(plant_recovery)
CWM_SLA_null_model_recovery_extra <- data.frame(matrix(0, nrow=1000, ncol=4463))
colnames(CWM_SLA_null_model_recovery_extra) <- rownames(plant_recovery)
CWM_PLH_null_model_recovery_extra <- data.frame(matrix(0, nrow=1000, ncol=4463))
colnames(CWM_PLH_null_model_recovery_extra) <- rownames(plant_recovery)

for (j in 1:1000){
  null_model_recovery_extra <- species_models_recovery_extra
  for (i in 1 : 4463){
    random_sampling <- sample(c(1:831), size = 831, replace = FALSE, prob = mean_recovery)
    sampled_species=c(rep(0,831))
    sampled_species[random_sampling]=species_models_recovery_extra[i,order(-species_models_recovery_extra[i,])]
    null_model_recovery_extra[i,]<- sampled_species
  }
  CWM_null_model_recovery_extra <- CWM(null_model_recovery_extra)
  CWM_LNC_null_model_recovery_extra[j,] <- CWM_null_model_recovery_extra$LNC
  CWM_SLA_null_model_recovery_extra[j,] <- CWM_null_model_recovery_extra$SLA
  CWM_PLH_null_model_recovery_extra[j,] <- CWM_null_model_recovery_extra$PLH
  FDis_null_model_recovery_extra <- FDis(null_model_recovery_extra)
  FDis_LNC_null_model_recovery_extra[j,] <- FDis_null_model_recovery_extra$LNC
  FDis_SLA_null_model_recovery_extra[j,] <- FDis_null_model_recovery_extra$SLA
  FDis_PLH_null_model_recovery_extra[j,] <- FDis_null_model_recovery_extra$PLH
}

R2_CWM_null_model_recovery_extra <- data.frame("R2_PLH"=rep(0,1000), "RMSE_PLH"=rep(0,1000))
for (i in 1:1000){
  R2_CWM_null_model_recovery_extra[i,] <- c(cor(CWM_PLH_null_model_recovery_extra[i,], CWM_obs$PLH[i])^2, rmse(CWM_PLH_null_model_recovery_extra[i,], CWM_obs$PLH[i]))
  print(i)
}

R2_FDis_null_model_recovery_extra <- data.frame("R2_LNC"=rep(0,1000), "RMSE_LNC"=rep(0,1000),
                                                "R2_SLA"=rep(0,1000), "RMSE_SLA"=rep(0,1000),
                                                "R2_PLH"=rep(0,1000), "RMSE_PLH"=rep(0,1000))
for (i in 1:1000){
  R2_FDis_null_model_recovery_extra[i,] <- c(cor(CWM_LNC_null_model_recovery_extra[i,], CWM_obs$LNC[i])^2, rmse(CWM_LNC_null_model_recovery_extra[i,], CWM_obs$LNC[i]),
                                            cor(CWM_SLA_null_model_recovery_extra[i,], CWM_obs$SLA[i])^2, rmse(CWM_SLA_null_model_recovery_extra[i,], CWM_obs$SLA[i]),
                                            cor(CWM_PLH_null_model_recovery_extra[i,], CWM_obs$PLH[i])^2, rmse(CWM_PLH_null_model_recovery_extra[i,], CWM_obs$PLH[i]))
  print(i)
}

### Figure 5 for the plot i ####
colnames(plant_recovery)=str_replace(colnames(plant_recovery),"\\."," ")
Obs_plot_i<- c()
Species_i <- which(plant_recovery[i,]>0)
for (s in 1:length(Species_i)){
  sp=Species_i[s]
  Obs_plot_i<-c(Obs_plot_i,rep(traits$PLH[which(traits$Taxa==colnames(plant_recovery)[sp])],plant_recovery[i,sp]*100))
}

Inter_plot_i<- c()
Species_i <- which(species_models_recovery_inter[i,]>0)
for (s in 1:length(Species_i)){
  sp=Species_i[s]
  Inter_plot_i<-c(Inter_plot_i,rep(traits$PLH[which(traits$Taxa==colnames(species_models_recovery_inter)[sp])],species_models_recovery_inter[i,sp]*100))
}

data.hist <- data.frame("Traits"=c(Obs_plot_i, Inter_plot_i), 
                        "Data"=c(rep("Observations", length(Obs_plot_i)), rep("Interpolations", length(Inter_plot_i))))
data.violin <- data.frame("Site"=CWM_PLH_null_model_recovery_inter[,i], "y"=rep(0,1000))

ggplot()+
  geom_histogram(data = data.hist[which(data.hist$Data=="Observations"),], aes(x = Traits, color="Observations", y=100*(..count..)/sum(..count..)),bins = 1000, position="dodge", alpha=0.2, show.legend = FALSE) +
  geom_histogram(data = data.hist[which(data.hist$Data=="Interpolations"),], aes(x = Traits, color="Interpolations", y=-100*(..count..)/sum(..count..)),bins = 1000, position="dodge", alpha=0.2, show.legend=FALSE) +
  scale_color_manual(values=c("#F8766D", "grey"))+
  scale_fill_hue(direction = 1) +
  ylab("Species coverage")+
  xlab("Plant Heigth")+
  annotate(geom = "rect", xmin = CWM_obs$PLH[i]-FDis_obs$PLH[i], xmax = CWM_obs$PLH[i]+FDis_obs$PLH[i], 
           ymin = 0, ymax = Inf, 
           fill = "grey", alpha = 0.2) +  
  annotate(geom = "rect", xmin = CWM_species_based_inter$PLH[i]-FDis_species_based_inter$PLH[i], xmax = CWM_species_based_inter$PLH[i]+FDis_species_based_inter$PLH[i], 
           ymin = -Inf, ymax = 0, 
           fill = "#F8766D", alpha = 0.2) +
  geom_abline(intercept=0, slope=0, color="#E6E6E6")+
  ggtitle(label="PLH")+
  theme_minimal()+
  scale_y_continuous(breaks=c(-90, -80, -70, -60,-50, -40, -30, -20, -10, 0, 10, 20, 30,40, 50,60, 70, 80, 90), labels = scales::percent_format(scale = 1))+
  geom_violin(data=data.violin,aes(x=Site, y=y))+
  geom_linerange(aes(x=CWM_obs$PLH[i], y=NULL, ymax=Inf, ymin=0), linetype="dashed", color="grey")+
  geom_linerange(aes(x=CWM_trait_based_inter$PLH[i], y=NULL, ymax=0, ymin=-Inf), linetype="dashed", color="#F8766D")+
  geom_rug(data = traits, mapping=aes(x=PLH), color="black")+
  annotate(geom="text", x=125, y=5, label="Observations")+
  annotate(geom="text", x=125, y=-5, label="Interpolations")+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))
  
##################################################################################
########### Mechanistically simulated virtual plant community data################
##################################################################################

######### To simulate our virtual plant community we used a pool of N=2000 species
x <- 10 ^ runif(2000, min = -3, max = 0) ## each species is associated with a prevalence randomly selected from a logistic rule from 1 to 10^-3

# This function returns a list of 20 numbers that represent the species composition of a simulated plant community.
com <- function(sigma, env){
  pool <- data.frame("species"=c(1:2000), "prob"=x) ## pool of 2000 species with a prevalence x
  com <- c()
  while(length(com)<20){ ## random selection of species stops when there are 20 individuals in the community
    sp = round(rnorm(1, env, sigma)) ## random selection of individuals from the species pool based on their trait from a normal distribution parameterized on the values associated with the environment
    if(sp>=1 & sp<=2000){ 
      des <- runif(1,0,1) ## random value between 0 and 1
      if(des<=pool$prob[which(pool$species==sp)]){ ## if the value obtained is lower than the prevalence associated to the species of the individual we incorporate it in the community
        com=c(com, sp)
      }
    }
  }
  return(com)
}

## In order to integrate communities of different composition in the same environments, we simulated 10 communities for each environment
env <- rep(501:1500, 10) # env (from 501 to 1500) represents the optimal trait value of the species on the environment.
env <- env[order(env)]
sigma <-rep(c(rep((1:20)*10,1:20)),5) ## sigma (from 1 to 200) a standard deviation value which represents the ability of the site to accept species with a trait value different from the optimal value env. The standard deviation is correlated with the optimal value env
sigma <- sigma[order(sigma)]
sigma <-sigma[1:1000]
sigma <- rep(sigma, 10)

ggenv_sigma <- data.frame("Optimal value y"=env, "Standard deviation sd"=sigma)
ggplot(ggenv_sigma, aes(x=Optimal.value.y, y=Standard.deviation.sd))+
  geom_point()


### Simulation of 10000 communities of 20 individus :
communautes <- as.data.frame(matrix(0, nrow=10000, ncol=20))

for (i in 1:10000){
  communautes[i,]<- com(sigma[i], env[i])
  print(i)
}

## transformation into site/species abundances data :
site_sp <- data.frame(matrix(0, ncol=2000, nrow=10000))
for (i in 1:10000){
  n=1
  m=unique(as.numeric(communautes[i,]))
  for (j in m[order(m)]){
    site_sp[i,j]<- table(as.numeric(communautes[i,]))[[n]]
    n=n+1
  }
  print(i)
}

site_sp<- cbind(env, site_sp)
com_simul <- site_sp[which(site_sp$env%in%c(501:1500)),]
colnames(com_simul)<-c("Env",1:2000)
com_simul<- com_simul[,-which(colSums(com_simul)==0)]


####### Final simulated data
### we degraded the environmental value by adding a random error (from a normal distribution centered on 0 with a standard deviation of 10)
env_simul_obs <- data.frame("Site"=1:10000,"Env"=com_simul$Env+rnorm(1,0,10))

com_simul <- com_simul[,-1]
traits_simul <- data.frame("Taxa"= paste0("species_",1:dim(com_simul)[2]), "Trait"=as.numeric(colnames(com_simul)))
colnames(com_simul) <- paste0("species_",1:dim(com_simul)[2])
rownames(com_simul) <- 1:10000

##### in order to limit the overfitting of the models, we have randomly selected 4 replicates of 300 of these environments to build communities.
replicate_1 <- sample(c(1:10000), 300)
replicate_2 <- sample(c(1:10000)[-c(replicate_1)], 300)
replicate_3 <- sample(c(1:10000)[-c(replicate_1, replicate_2)], 300)
replicate_4 <- sample(c(1:10000)[-c(replicate_1, replicate_2, replicate_3)], 300)

com_replicate_1 <- com_simul[replicate_1,]
com_replicate_1 <- com_replicate_1[,which(colSums(com_replicate_1)>0)]
env_replicate_1 <- env_simul_obs[replicate_1,]

com_replicate_2 <- com_simul[replicate_2,]
com_replicate_2 <- com_replicate_2[,which(colSums(com_replicate_2)>0)]
env_replicate_2 <- env_simul_obs[replicate_2,]

com_replicate_3 <- com_simul[replicate_3,]
com_replicate_3 <- com_replicate_3[,which(colSums(com_replicate_3)>0)]
env_replicate_3 <- env_simul_obs[replicate_3,]

com_replicate_4 <- com_simul[replicate_4,]
com_replicate_4 <- com_replicate_4[,which(colSums(com_replicate_4)>0)]
env_replicate_4 <- env_simul_obs[replicate_4,]

### The final simulated dataset is thus composed of 4*300 sites
com_tot <- full_join(com_replicate_1, com_replicate_2)
com_tot <- full_join(com_tot, com_replicate_3)
com_tot <- full_join(com_tot, com_replicate_4)
rownames(com_tot)<- c(rownames(com_replicate_1),rownames(com_replicate_2),rownames(com_replicate_3), rownames(com_replicate_4) )
com_tot <- cbind("Site"=as.numeric(rownames(com_tot)), com_tot)
com_tot <- com_tot[order(com_tot$Site), -1]
for (i in 1: 1200){
  for (j in 1 :dim(com_tot)[2]){
    if (is.na(com_tot[i,j])){
      com_tot[i,j]=0
    }
  }
}


CWM_simul <- CWM_simul(com_tot)
FDis_simul <- FDis_simul(com_tot)

CWM_simul_env <- left_join(CWM_simul, env_simul_obs, by="Site")
CWM_simul_env <- CWM_simul_env[order(CWM_simul_env$Site),-3]
FDis_simul_env <- left_join(FDis_simul, env_simul_obs, by="Site")
FDis_simul_env <- FDis_simul_env[order(FDis_simul_env$Site),-3]



########################################
####### Cross-validation groups ########
########################################

#### 4 interpolations groups for each replicated
replicate_1_inter_1 <- env_replicate_1$Site[sample(1:300, 75)]
replicate_1_inter_2 <- env_replicate_1$Site[-which(env_replicate_1$Site %in% replicate_1_inter_1)][sample(1:225, 75)]
replicate_1_inter_3 <- env_replicate_1$Site[-which(env_replicate_1$Site %in% c(replicate_1_inter_1, replicate_1_inter_2))][sample(1:150, 75)]
replicate_1_inter_4 <- env_replicate_1$Site[-which(env_replicate_1$Site %in% c(replicate_1_inter_1, replicate_1_inter_2, replicate_1_inter_3))]

replicate_2_inter_1 <- env_replicate_2$Site[sample(1:300, 75)]
replicate_2_inter_2 <- env_replicate_2$Site[-which(env_replicate_2$Site %in% replicate_2_inter_1)][sample(1:225, 75)]
replicate_2_inter_3 <- env_replicate_2$Site[-which(env_replicate_2$Site %in% c(replicate_2_inter_1, replicate_2_inter_2))][sample(1:150, 75)]
replicate_2_inter_4 <- env_replicate_2$Site[-which(env_replicate_2$Site %in% c(replicate_2_inter_1, replicate_2_inter_2, replicate_2_inter_3))]

replicate_3_inter_1 <- env_replicate_3$Site[sample(1:300, 75)]
replicate_3_inter_2 <- env_replicate_3$Site[-which(env_replicate_3$Site %in% replicate_3_inter_1)][sample(1:225, 75)]
replicate_3_inter_3 <- env_replicate_3$Site[-which(env_replicate_3$Site %in% c(replicate_3_inter_1, replicate_3_inter_2))][sample(1:150, 75)]
replicate_3_inter_4 <- env_replicate_3$Site[-which(env_replicate_3$Site %in% c(replicate_3_inter_1, replicate_3_inter_2, replicate_3_inter_3))]

replicate_4_inter_1 <- env_replicate_4$Site[sample(1:300, 75)]
replicate_4_inter_2 <- env_replicate_4$Site[-which(env_replicate_4$Site %in% replicate_4_inter_1)][sample(1:225, 75)]
replicate_4_inter_3 <- env_replicate_4$Site[-which(env_replicate_4$Site %in% c(replicate_4_inter_1, replicate_4_inter_2))][sample(1:150, 75)]
replicate_4_inter_4 <- env_replicate_4$Site[-which(env_replicate_4$Site %in% c(replicate_4_inter_1, replicate_4_inter_2, replicate_4_inter_3))]


#### For each of the four independent replicates in the simulated dataset, we created the four folds by randomly drawing 
#### sites according to a normal distribution centered on 4 different environmental value points of the sites
# replicate 1
norm_replicate_1_extra_1 <- rnorm(200, 37, 25)
norm_replicate_1_extra_2 <- rnorm(200, 112, 25)
norm_replicate_1_extra_3 <- rnorm(200, 187, 25)
norm_replicate_1_extra_4 <- rnorm(200, 262, 25)

replicate_1_extra_1 <- unique(round(norm_replicate_1_extra_1))[which(unique(round(norm_replicate_1_extra_1))>0)]
replicate_1_extra_2 <- unique(round(norm_replicate_1_extra_2))[which(unique(round(norm_replicate_1_extra_2))>0)]
replicate_1_extra_3 <- unique(round(norm_replicate_1_extra_3))[which(unique(round(norm_replicate_1_extra_3))<300)]
replicate_1_extra_4 <- unique(round(norm_replicate_1_extra_4))[which(unique(round(norm_replicate_1_extra_4))<300)]
double_1_extra_1 <- which(replicate_1_extra_1%in%replicate_1_extra_1[which(replicate_1_extra_1>=75)][which(replicate_1_extra_1[which(replicate_1_extra_1>=75)]%in% replicate_1_extra_2|replicate_1_extra_1[which(replicate_1_extra_1>=75)]%in% replicate_1_extra_3|replicate_1_extra_1[which(replicate_1_extra_1>75)]%in% replicate_1_extra_4)])
if(length(double_1_extra_1>0)){
  replicate_1_extra_1 <- replicate_1_extra_1[-double_1_extra_1]
}
double_1_extra_2 <- which(replicate_1_extra_2 %in% replicate_1_extra_2[which(replicate_1_extra_2>=150|replicate_1_extra_2<75)][which(replicate_1_extra_2[which(replicate_1_extra_2>=150|replicate_1_extra_2<75)]%in% replicate_1_extra_1|replicate_1_extra_2[which(replicate_1_extra_2>=150|replicate_1_extra_2<75)]%in% replicate_1_extra_3|replicate_1_extra_2[which(replicate_1_extra_2>=150|replicate_1_extra_2<75)]%in% replicate_1_extra_4)])
if(length(double_1_extra_2>0)){
  replicate_1_extra_2 <- replicate_1_extra_2[-double_1_extra_2]
}
double_1_extra_3 <- which(replicate_1_extra_3 %in% replicate_1_extra_3[which(replicate_1_extra_3>=225|replicate_1_extra_3<150)][which(replicate_1_extra_3[which(replicate_1_extra_3>=225|replicate_1_extra_3<150)]%in% replicate_1_extra_1|replicate_1_extra_3[which(replicate_1_extra_3>=225|replicate_1_extra_3<150)]%in% replicate_1_extra_2|replicate_1_extra_3[which(replicate_1_extra_3>=225|replicate_1_extra_3<150)]%in% replicate_1_extra_4)])
if(length(double_1_extra_3>0)){
  replicate_1_extra_3 <- replicate_1_extra_3[-double_1_extra_3]
}
double_1_extra_4 <- which(replicate_1_extra_4 %in% replicate_1_extra_4[which(replicate_1_extra_4<225)][which(replicate_1_extra_4[which(replicate_1_extra_4<225)]%in% replicate_1_extra_1|replicate_1_extra_4[which(replicate_1_extra_4<225)]%in% replicate_1_extra_2|replicate_1_extra_4[which(replicate_1_extra_4<225)]%in% replicate_1_extra_3)])
if(length(double_1_extra_4>0)){
  replicate_1_extra_4 <- replicate_1_extra_4[-double_1_extra_4]
}
replicate_1_extra_1 <- env_replicate_1$Site[order(env_replicate_1$Env)][replicate_1_extra_1]
replicate_1_extra_2 <- env_replicate_1$Site[order(env_replicate_1$Env)][replicate_1_extra_2]
replicate_1_extra_3 <- env_replicate_1$Site[order(env_replicate_1$Env)][replicate_1_extra_3]
replicate_1_extra_4 <- env_replicate_1$Site[order(env_replicate_1$Env)][replicate_1_extra_4]

# Replicate 2
norm_replicate_2_extra_1 <- rnorm(200, 37, 25)
norm_replicate_2_extra_2 <- rnorm(200, 112, 25)
norm_replicate_2_extra_3 <- rnorm(200, 187, 25)
norm_replicate_2_extra_4 <- rnorm(200, 262, 25)

replicate_2_extra_1 <- unique(round(norm_replicate_2_extra_1))[which(unique(round(norm_replicate_2_extra_1))>0)]
replicate_2_extra_2 <- unique(round(norm_replicate_2_extra_2))[which(unique(round(norm_replicate_2_extra_2))>0)]
replicate_2_extra_3 <- unique(round(norm_replicate_2_extra_3))[which(unique(round(norm_replicate_2_extra_3))<300)]
replicate_2_extra_4 <- unique(round(norm_replicate_2_extra_4))[which(unique(round(norm_replicate_2_extra_4))<300)]
double_2_extra_1 <- which(replicate_2_extra_1%in%replicate_2_extra_1[which(replicate_2_extra_1>=75)][which(replicate_2_extra_1[which(replicate_2_extra_1>=75)]%in% replicate_2_extra_2|replicate_2_extra_1[which(replicate_2_extra_1>=75)]%in% replicate_2_extra_3|replicate_2_extra_1[which(replicate_2_extra_1>75)]%in% replicate_2_extra_4)])
if(length(double_2_extra_1>0)){
  replicate_2_extra_1 <- replicate_2_extra_1[-double_2_extra_1]
}
double_2_extra_2 <- which(replicate_2_extra_2 %in% replicate_2_extra_2[which(replicate_2_extra_2>=150|replicate_2_extra_2<75)][which(replicate_2_extra_2[which(replicate_2_extra_2>=150|replicate_2_extra_2<75)]%in% replicate_2_extra_1|replicate_2_extra_2[which(replicate_2_extra_2>=150|replicate_2_extra_2<75)]%in% replicate_2_extra_3|replicate_2_extra_2[which(replicate_2_extra_2>=150|replicate_2_extra_2<75)]%in% replicate_2_extra_4)])
if(length(double_2_extra_2>0)){
  replicate_2_extra_2 <- replicate_2_extra_2[-double_2_extra_2]
}
double_2_extra_3 <- which(replicate_2_extra_3 %in% replicate_2_extra_3[which(replicate_2_extra_3>=225|replicate_2_extra_3<150)][which(replicate_2_extra_3[which(replicate_2_extra_3>=225|replicate_2_extra_3<150)]%in% replicate_2_extra_1|replicate_2_extra_3[which(replicate_2_extra_3>=225|replicate_2_extra_3<150)]%in% replicate_2_extra_2|replicate_2_extra_3[which(replicate_2_extra_3>=225|replicate_2_extra_3<150)]%in% replicate_2_extra_4)])
if(length(double_2_extra_3>0)){
  replicate_2_extra_3 <- replicate_2_extra_3[-double_2_extra_3]
}
double_2_extra_4 <- which(replicate_2_extra_4 %in% replicate_2_extra_4[which(replicate_2_extra_4<225)][which(replicate_2_extra_4[which(replicate_2_extra_4<225)]%in% replicate_2_extra_1|replicate_2_extra_4[which(replicate_2_extra_4<225)]%in% replicate_2_extra_2|replicate_2_extra_4[which(replicate_2_extra_4<225)]%in% replicate_2_extra_3)])
if(length(double_2_extra_4>0)){
  replicate_2_extra_4 <- replicate_2_extra_4[-double_2_extra_4]
}
replicate_2_extra_1 <- env_replicate_2$Site[order(env_replicate_2$Env)][replicate_2_extra_1]
replicate_2_extra_2 <- env_replicate_2$Site[order(env_replicate_2$Env)][replicate_2_extra_2]
replicate_2_extra_3 <- env_replicate_2$Site[order(env_replicate_2$Env)][replicate_2_extra_3]
replicate_2_extra_4 <- env_replicate_2$Site[order(env_replicate_2$Env)][replicate_2_extra_4]

# Replicate 3
norm_replicate_3_extra_1 <- rnorm(200, 37, 25)
norm_replicate_3_extra_2 <- rnorm(200, 112, 25)
norm_replicate_3_extra_3 <- rnorm(200, 187, 25)
norm_replicate_3_extra_4 <- rnorm(200, 262, 25)

replicate_3_extra_1 <- unique(round(norm_replicate_3_extra_1))[which(unique(round(norm_replicate_3_extra_1))>0)]
replicate_3_extra_2 <- unique(round(norm_replicate_3_extra_2))[which(unique(round(norm_replicate_3_extra_2))>0)]
replicate_3_extra_3 <- unique(round(norm_replicate_3_extra_3))[which(unique(round(norm_replicate_3_extra_3))<300)]
replicate_3_extra_4 <- unique(round(norm_replicate_3_extra_4))[which(unique(round(norm_replicate_3_extra_4))<300)]
double_3_extra_1 <- which(replicate_3_extra_1%in%replicate_3_extra_1[which(replicate_3_extra_1>=75)][which(replicate_3_extra_1[which(replicate_3_extra_1>=75)]%in% replicate_3_extra_2|replicate_3_extra_1[which(replicate_3_extra_1>=75)]%in% replicate_3_extra_3|replicate_3_extra_1[which(replicate_3_extra_1>75)]%in% replicate_3_extra_4)])
if(length(double_3_extra_1>0)){
  replicate_3_extra_1 <- replicate_3_extra_1[-double_3_extra_1]
}
double_3_extra_2 <- which(replicate_3_extra_2 %in% replicate_3_extra_2[which(replicate_3_extra_2>=150|replicate_3_extra_2<75)][which(replicate_3_extra_2[which(replicate_3_extra_2>=150|replicate_3_extra_2<75)]%in% replicate_3_extra_1|replicate_3_extra_2[which(replicate_3_extra_2>=150|replicate_3_extra_2<75)]%in% replicate_3_extra_3|replicate_3_extra_2[which(replicate_3_extra_2>=150|replicate_3_extra_2<75)]%in% replicate_3_extra_4)])
if(length(double_3_extra_2>0)){
  replicate_3_extra_2 <- replicate_3_extra_2[-double_3_extra_2]
}
double_3_extra_3 <- which(replicate_3_extra_3 %in% replicate_3_extra_3[which(replicate_3_extra_3>=225|replicate_3_extra_3<150)][which(replicate_3_extra_3[which(replicate_3_extra_3>=225|replicate_3_extra_3<150)]%in% replicate_3_extra_1|replicate_3_extra_3[which(replicate_3_extra_3>=225|replicate_3_extra_3<150)]%in% replicate_3_extra_2|replicate_3_extra_3[which(replicate_3_extra_3>=225|replicate_3_extra_3<150)]%in% replicate_3_extra_4)])
if(length(double_3_extra_3>0)){
  replicate_3_extra_3 <- replicate_3_extra_3[-double_3_extra_3]
}
double_3_extra_4 <- which(replicate_3_extra_4 %in% replicate_3_extra_4[which(replicate_3_extra_4<225)][which(replicate_3_extra_4[which(replicate_3_extra_4<225)]%in% replicate_3_extra_1|replicate_3_extra_4[which(replicate_3_extra_4<225)]%in% replicate_3_extra_2|replicate_3_extra_4[which(replicate_3_extra_4<225)]%in% replicate_3_extra_3)])
if(length(double_3_extra_4>0)){
  replicate_3_extra_4 <- replicate_3_extra_4[-double_3_extra_4]
}
replicate_3_extra_1 <- env_replicate_3$Site[order(env_replicate_3$Env)][replicate_3_extra_1]
replicate_3_extra_2 <- env_replicate_3$Site[order(env_replicate_3$Env)][replicate_3_extra_2]
replicate_3_extra_3 <- env_replicate_3$Site[order(env_replicate_3$Env)][replicate_3_extra_3]
replicate_3_extra_4 <- env_replicate_3$Site[order(env_replicate_3$Env)][replicate_3_extra_4]

# Replicate 4
norm_replicate_4_extra_1 <- rnorm(200, 37, 25)
norm_replicate_4_extra_2 <- rnorm(200, 112, 25)
norm_replicate_4_extra_3 <- rnorm(200, 187, 25)
norm_replicate_4_extra_4 <- rnorm(200, 262, 25)

replicate_4_extra_1 <- unique(round(norm_replicate_4_extra_1))[which(unique(round(norm_replicate_4_extra_1))>0)]
replicate_4_extra_2 <- unique(round(norm_replicate_4_extra_2))[which(unique(round(norm_replicate_4_extra_2))>0)]
replicate_4_extra_3 <- unique(round(norm_replicate_4_extra_3))[which(unique(round(norm_replicate_4_extra_3))<300)]
replicate_4_extra_4 <- unique(round(norm_replicate_4_extra_4))[which(unique(round(norm_replicate_4_extra_4))<300)]
double_4_extra_1 <- which(replicate_4_extra_1%in%replicate_4_extra_1[which(replicate_4_extra_1>=75)][which(replicate_4_extra_1[which(replicate_4_extra_1>=75)]%in% replicate_4_extra_2|replicate_4_extra_1[which(replicate_4_extra_1>=75)]%in% replicate_4_extra_3|replicate_4_extra_1[which(replicate_4_extra_1>75)]%in% replicate_4_extra_4)])
if(length(double_4_extra_1>0)){
  replicate_4_extra_1 <- replicate_4_extra_1[-double_4_extra_1]
}
double_4_extra_2 <- which(replicate_4_extra_2 %in% replicate_4_extra_2[which(replicate_4_extra_2>=150|replicate_4_extra_2<75)][which(replicate_4_extra_2[which(replicate_4_extra_2>=150|replicate_4_extra_2<75)]%in% replicate_4_extra_1|replicate_4_extra_2[which(replicate_4_extra_2>=150|replicate_4_extra_2<75)]%in% replicate_4_extra_3|replicate_4_extra_2[which(replicate_4_extra_2>=150|replicate_4_extra_2<75)]%in% replicate_4_extra_4)])
if(length(double_4_extra_2>0)){
  replicate_4_extra_2 <- replicate_4_extra_2[-double_4_extra_2]
}
double_4_extra_3 <- which(replicate_4_extra_3 %in% replicate_4_extra_3[which(replicate_4_extra_3>=225|replicate_4_extra_3<150)][which(replicate_4_extra_3[which(replicate_4_extra_3>=225|replicate_4_extra_3<150)]%in% replicate_4_extra_1|replicate_4_extra_3[which(replicate_4_extra_3>=225|replicate_4_extra_3<150)]%in% replicate_4_extra_2|replicate_4_extra_3[which(replicate_4_extra_3>=225|replicate_4_extra_3<150)]%in% replicate_4_extra_4)])
if(length(double_4_extra_3>0)){
  replicate_4_extra_3 <- replicate_4_extra_3[-double_4_extra_3]
}
double_4_extra_4 <- which(replicate_4_extra_4 %in% replicate_4_extra_4[which(replicate_4_extra_4<225)][which(replicate_4_extra_4[which(replicate_4_extra_4<225)]%in% replicate_4_extra_1|replicate_4_extra_4[which(replicate_4_extra_4<225)]%in% replicate_4_extra_2|replicate_4_extra_4[which(replicate_4_extra_4<225)]%in% replicate_4_extra_3)])
if(length(double_4_extra_4>0)){
  replicate_4_extra_4 <- replicate_4_extra_4[-double_4_extra_4]
}
replicate_4_extra_1 <- env_replicate_4$Site[order(env_replicate_4$Env)][replicate_4_extra_1]
replicate_4_extra_2 <- env_replicate_4$Site[order(env_replicate_4$Env)][replicate_4_extra_2]
replicate_4_extra_3 <- env_replicate_4$Site[order(env_replicate_4$Env)][replicate_4_extra_3]
replicate_4_extra_4 <- env_replicate_4$Site[order(env_replicate_4$Env)][replicate_4_extra_4]


####################################
#### -- TRAIT-BASED APPROACH -- ####
####################################
####### -- Traits Distribution Models for the Trait-based approach-- 
###### PREDICTION OF COMMUNITY WEIGHTED MEAN  ######
#### -- This computes path predicts the community weighted mean in the context of interpolation -- ####
gam_data <- left_join(CWM_simul, env_simul_obs, by="Site")
gam_data <- gam_data[order(gam_data$Site),]
### replicate 1:
## - For interpolation group 1 - ##
models_CWM_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_1),])
CWM_trait_based_inter1<-data.frame("Trait"=predict(models_CWM_trait_grow1, gam_data[which(gam_data$Site%in%replicate_1_inter_1),])$predicted)
rownames(CWM_trait_based_inter1) <- gam_data$Site[which(gam_data$Site%in%replicate_1_inter_1)]

## - For extrapolation group 1 - ##
models_CWM_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_1),])
CWM_trait_based_extra1<-data.frame("Trait"=predict(models_CWM_trait_grow1, gam_data[which(gam_data$Site%in%replicate_1_extra_1),])$predicted)
rownames(CWM_trait_based_extra1) <- gam_data$Site[which(gam_data$Site%in%replicate_1_extra_1)]

## - For interpolation group 2 - ##
models_CWM_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_2),])
CWM_trait_based_inter2<-data.frame("Trait"=predict(models_CWM_trait_grow2, gam_data[which(gam_data$Site%in%replicate_1_inter_2),])$predicted)
rownames(CWM_trait_based_inter2) <- gam_data$Site[which(gam_data$Site%in%replicate_1_inter_2)]

## - For extrapolation group 2 - ##
models_CWM_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_2),])
CWM_trait_based_extra2<-data.frame("Trait"=predict(models_CWM_trait_grow2, gam_data[which(gam_data$Site%in%replicate_1_extra_2),])$predicted)
rownames(CWM_trait_based_extra2) <- gam_data$Site[which(gam_data$Site%in%replicate_1_extra_2)]

## - For interpolation group 3 - ##
models_CWM_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_3),])
CWM_trait_based_inter3<-data.frame("Trait"=predict(models_CWM_trait_grow3, gam_data[which(gam_data$Site%in%replicate_1_inter_3),])$predicted)
rownames(CWM_trait_based_inter3) <- gam_data$Site[which(gam_data$Site%in%replicate_1_inter_3)]

## - For extrapolation group 3 - ##
models_CWM_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_3),])
CWM_trait_based_extra3<-data.frame("Trait"=predict(models_CWM_trait_grow3, gam_data[which(gam_data$Site%in%replicate_1_extra_3),])$predicted)
rownames(CWM_trait_based_extra3) <- gam_data$Site[which(gam_data$Site%in%replicate_1_extra_3)]

## - For interpolation group 4 - ##
models_CWM_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_4),])
CWM_trait_based_inter4<-data.frame("Trait"=predict(models_CWM_trait_grow4, gam_data[which(gam_data$Site%in%replicate_1_inter_4),])$predicted)
rownames(CWM_trait_based_inter4) <- gam_data$Site[which(gam_data$Site%in%replicate_1_inter_4)]

## - For extrapolation group 4 - ##
models_CWM_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_4),])
CWM_trait_based_extra4<-data.frame("Trait"=predict(models_CWM_trait_grow4, gam_data[which(gam_data$Site%in%replicate_1_extra_4),])$predicted)
rownames(CWM_trait_based_extra4) <- gam_data$Site[which(gam_data$Site%in%replicate_1_extra_4)]

## assembly in a single prediction table for the 4 groups :
CWM_AFPL_inter_replicate_1 <- rbind(CWM_trait_based_inter1, CWM_trait_based_inter2, CWM_trait_based_inter3, CWM_trait_based_inter4)
CWM_AFPL_inter_replicate_1 <- cbind("Site"=as.numeric(rownames(CWM_AFPL_inter_replicate_1)), CWM_AFPL_inter_replicate_1)
CWM_AFPL_inter_replicate_1 <- CWM_AFPL_inter_replicate_1[order(CWM_AFPL_inter_replicate_1$Site),]
rownames(CWM_AFPL_inter_replicate_1) <- CWM_AFPL_inter_replicate_1$Site

CWM_AFPL_extra_replicate_1 <- rbind(CWM_trait_based_extra1, CWM_trait_based_extra2, CWM_trait_based_extra3, CWM_trait_based_extra4)
CWM_AFPL_extra_replicate_1 <- cbind("Site"=as.numeric(rownames(CWM_AFPL_extra_replicate_1)), CWM_AFPL_extra_replicate_1)
CWM_AFPL_extra_replicate_1 <- CWM_AFPL_extra_replicate_1[order(CWM_AFPL_extra_replicate_1$Site),]
rownames(CWM_AFPL_extra_replicate_1) <- CWM_AFPL_extra_replicate_1$Site


### replicate 2
## - For interpolation group 1 - ##
models_CWM_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_1),])
CWM_trait_based_inter1<-data.frame("Trait"=predict(models_CWM_trait_grow1, gam_data[which(gam_data$Site%in%replicate_2_inter_1),])$predicted)
rownames(CWM_trait_based_inter1) <- gam_data$Site[which(gam_data$Site%in%replicate_2_inter_1)]

## - For extrapolation group 1 - ##
models_CWM_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_1),])
CWM_trait_based_extra1<-data.frame("Trait"=predict(models_CWM_trait_grow1, gam_data[which(gam_data$Site%in%replicate_2_extra_1),])$predicted)
rownames(CWM_trait_based_extra1) <- gam_data$Site[which(gam_data$Site%in%replicate_2_extra_1)]

## - For interpolation group 2 - ##
models_CWM_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_2),])
CWM_trait_based_inter2<-data.frame("Trait"=predict(models_CWM_trait_grow2, gam_data[which(gam_data$Site%in%replicate_2_inter_2),])$predicted)
rownames(CWM_trait_based_inter2) <- gam_data$Site[which(gam_data$Site%in%replicate_2_inter_2)]

## - For extrapolation group 2 - ##
models_CWM_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_2),])
CWM_trait_based_extra2<-data.frame("Trait"=predict(models_CWM_trait_grow2, gam_data[which(gam_data$Site%in%replicate_2_extra_2),])$predicted)
rownames(CWM_trait_based_extra2) <- gam_data$Site[which(gam_data$Site%in%replicate_2_extra_2)]

## - For interpolation group 3 - ##
models_CWM_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_3),])
CWM_trait_based_inter3<-data.frame("Trait"=predict(models_CWM_trait_grow3, gam_data[which(gam_data$Site%in%replicate_2_inter_3),])$predicted)
rownames(CWM_trait_based_inter3) <- gam_data$Site[which(gam_data$Site%in%replicate_2_inter_3)]

## - For extrapolation group 3 - ##
models_CWM_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_3),])
CWM_trait_based_extra3<-data.frame("Trait"=predict(models_CWM_trait_grow3, gam_data[which(gam_data$Site%in%replicate_2_extra_3),])$predicted)
rownames(CWM_trait_based_extra3) <- gam_data$Site[which(gam_data$Site%in%replicate_2_extra_3)]

## - For interpolation group 4 - ##
models_CWM_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_4),])
CWM_trait_based_inter4<-data.frame("Trait"=predict(models_CWM_trait_grow4, gam_data[which(gam_data$Site%in%replicate_2_inter_4),])$predicted)
rownames(CWM_trait_based_inter4) <- gam_data$Site[which(gam_data$Site%in%replicate_2_inter_4)]

## - For extrapolation group 4 - ##
models_CWM_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_4),])
CWM_trait_based_extra4<-data.frame("Trait"=predict(models_CWM_trait_grow4, gam_data[which(gam_data$Site%in%replicate_2_extra_4),])$predicted)
rownames(CWM_trait_based_extra4) <- gam_data$Site[which(gam_data$Site%in%replicate_2_extra_4)]

## assembly in a single prediction table for the 4 groups :
CWM_AFPL_inter_replicate_2 <- rbind(CWM_trait_based_inter1, CWM_trait_based_inter2, CWM_trait_based_inter3, CWM_trait_based_inter4)
CWM_AFPL_inter_replicate_2 <- cbind("Site"=as.numeric(rownames(CWM_AFPL_inter_replicate_2)), CWM_AFPL_inter_replicate_2)
CWM_AFPL_inter_replicate_2 <- CWM_AFPL_inter_replicate_2[order(CWM_AFPL_inter_replicate_2$Site),]
rownames(CWM_AFPL_inter_replicate_2) <- CWM_AFPL_inter_replicate_2$Site

CWM_AFPL_extra_replicate_2 <- rbind(CWM_trait_based_extra1, CWM_trait_based_extra2, CWM_trait_based_extra3, CWM_trait_based_extra4)
CWM_AFPL_extra_replicate_2 <- cbind("Site"=as.numeric(rownames(CWM_AFPL_extra_replicate_2)), CWM_AFPL_extra_replicate_2)
CWM_AFPL_extra_replicate_2 <- CWM_AFPL_extra_replicate_2[order(CWM_AFPL_extra_replicate_2$Site),]
rownames(CWM_AFPL_extra_replicate_2) <- CWM_AFPL_extra_replicate_2$Site

### replicate 3
## - For interpolation group 1 - ##
models_CWM_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_1),])
CWM_trait_based_inter1<-data.frame("Trait"=predict(models_CWM_trait_grow1, gam_data[which(gam_data$Site%in%replicate_3_inter_1),])$predicted)
rownames(CWM_trait_based_inter1) <- gam_data$Site[which(gam_data$Site%in%replicate_3_inter_1)]

## - For extrapolation group 1 - ##
models_CWM_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_1),])
CWM_trait_based_extra1<-data.frame("Trait"=predict(models_CWM_trait_grow1, gam_data[which(gam_data$Site%in%replicate_3_extra_1),])$predicted)
rownames(CWM_trait_based_extra1) <- gam_data$Site[which(gam_data$Site%in%replicate_3_extra_1)]

## - For interpolation group 2 - ##
models_CWM_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_2),])
CWM_trait_based_inter2<-data.frame("Trait"=predict(models_CWM_trait_grow2, gam_data[which(gam_data$Site%in%replicate_3_inter_2),])$predicted)
rownames(CWM_trait_based_inter2) <- gam_data$Site[which(gam_data$Site%in%replicate_3_inter_2)]

## - For extrapolation group 2 - ##
models_CWM_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_2),])
CWM_trait_based_extra2<-data.frame("Trait"=predict(models_CWM_trait_grow2, gam_data[which(gam_data$Site%in%replicate_3_extra_2),])$predicted)
rownames(CWM_trait_based_extra2) <- gam_data$Site[which(gam_data$Site%in%replicate_3_extra_2)]

## - For interpolation group 3 - ##
models_CWM_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_3),])
CWM_trait_based_inter3<-data.frame("Trait"=predict(models_CWM_trait_grow3, gam_data[which(gam_data$Site%in%replicate_3_inter_3),])$predicted)
rownames(CWM_trait_based_inter3) <- gam_data$Site[which(gam_data$Site%in%replicate_3_inter_3)]

## - For extrapolation group 3 - ##
models_CWM_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_3),])
CWM_trait_based_extra3<-data.frame("Trait"=predict(models_CWM_trait_grow3, gam_data[which(gam_data$Site%in%replicate_3_extra_3),])$predicted)
rownames(CWM_trait_based_extra3) <- gam_data$Site[which(gam_data$Site%in%replicate_3_extra_3)]

## - For interpolation group 4 - ##
models_CWM_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_4),])
CWM_trait_based_inter4<-data.frame("Trait"=predict(models_CWM_trait_grow4, gam_data[which(gam_data$Site%in%replicate_3_inter_4),])$predicted)
rownames(CWM_trait_based_inter4) <- gam_data$Site[which(gam_data$Site%in%replicate_3_inter_4)]

## - For extrapolation group 4 - ##
models_CWM_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_4),])
CWM_trait_based_extra4<-data.frame("Trait"=predict(models_CWM_trait_grow4, gam_data[which(gam_data$Site%in%replicate_3_extra_4),])$predicted)
rownames(CWM_trait_based_extra4) <- gam_data$Site[which(gam_data$Site%in%replicate_3_extra_4)]

## assembly in a single prediction table for the 4 groups :
CWM_AFPL_inter_replicate_3 <- rbind(CWM_trait_based_inter1, CWM_trait_based_inter2, CWM_trait_based_inter3, CWM_trait_based_inter4)
CWM_AFPL_inter_replicate_3 <- cbind("Site"=as.numeric(rownames(CWM_AFPL_inter_replicate_3)), CWM_AFPL_inter_replicate_3)
CWM_AFPL_inter_replicate_3 <- CWM_AFPL_inter_replicate_3[order(CWM_AFPL_inter_replicate_3$Site),]
rownames(CWM_AFPL_inter_replicate_3) <- CWM_AFPL_inter_replicate_3$Site

CWM_AFPL_extra_replicate_3 <- rbind(CWM_trait_based_extra1, CWM_trait_based_extra2, CWM_trait_based_extra3, CWM_trait_based_extra4)
CWM_AFPL_extra_replicate_3 <- cbind("Site"=as.numeric(rownames(CWM_AFPL_extra_replicate_3)), CWM_AFPL_extra_replicate_3)
CWM_AFPL_extra_replicate_3 <- CWM_AFPL_extra_replicate_3[order(CWM_AFPL_extra_replicate_3$Site),]
rownames(CWM_AFPL_extra_replicate_3) <- CWM_AFPL_extra_replicate_3$Site

### replicate 4
## - For interpolation group 1 - ##
models_CWM_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_1),])
CWM_trait_based_inter1<-data.frame("Trait"=predict(models_CWM_trait_grow1, gam_data[which(gam_data$Site%in%replicate_4_inter_1),])$predicted)
rownames(CWM_trait_based_inter1) <- gam_data$Site[which(gam_data$Site%in%replicate_4_inter_1)]

## - For extrapolation group 1 - ##
models_CWM_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_1),])
CWM_trait_based_extra1<-data.frame("Trait"=predict(models_CWM_trait_grow1, gam_data[which(gam_data$Site%in%replicate_4_extra_1),])$predicted)
rownames(CWM_trait_based_extra1) <- gam_data$Site[which(gam_data$Site%in%replicate_4_extra_1)]

## - For interpolation group 2 - ##
models_CWM_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_2),])
CWM_trait_based_inter2<-data.frame("Trait"=predict(models_CWM_trait_grow2, gam_data[which(gam_data$Site%in%replicate_4_inter_2),])$predicted)
rownames(CWM_trait_based_inter2) <- gam_data$Site[which(gam_data$Site%in%replicate_4_inter_2)]

## - For extrapolation group 2 - ##
models_CWM_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_2),])
CWM_trait_based_extra2<-data.frame("Trait"=predict(models_CWM_trait_grow2, gam_data[which(gam_data$Site%in%replicate_4_extra_2),])$predicted)
rownames(CWM_trait_based_extra2) <- gam_data$Site[which(gam_data$Site%in%replicate_4_extra_2)]

## - For interpolation group 3 - ##
models_CWM_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_3),])
CWM_trait_based_inter3<-data.frame("Trait"=predict(models_CWM_trait_grow3, gam_data[which(gam_data$Site%in%replicate_4_inter_3),])$predicted)
rownames(CWM_trait_based_inter3) <- gam_data$Site[which(gam_data$Site%in%replicate_4_inter_3)]

## - For extrapolation group 3 - ##
models_CWM_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_3),])
CWM_trait_based_extra3<-data.frame("Trait"=predict(models_CWM_trait_grow3, gam_data[which(gam_data$Site%in%replicate_4_extra_3),])$predicted)
rownames(CWM_trait_based_extra3) <- gam_data$Site[which(gam_data$Site%in%replicate_4_extra_3)]

## - For interpolation group 4 - ##
models_CWM_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_4),])
CWM_trait_based_inter4<-data.frame("Trait"=predict(models_CWM_trait_grow4, gam_data[which(gam_data$Site%in%replicate_4_inter_4),])$predicted)
rownames(CWM_trait_based_inter4) <- gam_data$Site[which(gam_data$Site%in%replicate_4_inter_4)]

## - For extrapolation group 4 - ##
models_CWM_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_4),])
CWM_trait_based_extra4<-data.frame("Trait"=predict(models_CWM_trait_grow4, gam_data[which(gam_data$Site%in%replicate_4_extra_4),])$predicted)
rownames(CWM_trait_based_extra4) <- gam_data$Site[which(gam_data$Site%in%replicate_4_extra_4)]

## assembly in a single prediction table for the 4 groups :
CWM_AFPL_inter_replicate_4 <- rbind(CWM_trait_based_inter1, CWM_trait_based_inter2, CWM_trait_based_inter3, CWM_trait_based_inter4)
CWM_AFPL_inter_replicate_4 <- cbind("Site"=as.numeric(rownames(CWM_AFPL_inter_replicate_4)), CWM_AFPL_inter_replicate_4)
CWM_AFPL_inter_replicate_4 <- CWM_AFPL_inter_replicate_4[order(CWM_AFPL_inter_replicate_4$Site),]
rownames(CWM_AFPL_inter_replicate_4) <- CWM_AFPL_inter_replicate_4$Site

CWM_AFPL_extra_replicate_4 <- rbind(CWM_trait_based_extra1, CWM_trait_based_extra2, CWM_trait_based_extra3, CWM_trait_based_extra4)
CWM_AFPL_extra_replicate_4 <- cbind("Site"=as.numeric(rownames(CWM_AFPL_extra_replicate_4)), CWM_AFPL_extra_replicate_4)
CWM_AFPL_extra_replicate_4 <- CWM_AFPL_extra_replicate_4[order(CWM_AFPL_extra_replicate_4$Site),]
rownames(CWM_AFPL_extra_replicate_4) <- CWM_AFPL_extra_replicate_4$Site

###### PREDICTION OF FUNCTIONAL DISPERSION  ######
#### -- This computes path predicts the community weighted mean in the context of interpolation -- ####
gam_data <- left_join(FDis_simul, env_simul_obs, by="Site")
gam_data <- gam_data[order(gam_data$Site),]
### replicate 1:
## - For interpolation group 1 - ##
models_FDis_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_1),])
FDis_trait_based_inter1<-data.frame("Trait"=predict(models_FDis_trait_grow1, gam_data[which(gam_data$Site%in%replicate_1_inter_1),])$predicted)
rownames(FDis_trait_based_inter1) <- gam_data$Site[which(gam_data$Site%in%replicate_1_inter_1)]

## - For extrapolation group 1 - ##
models_FDis_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_1),])
FDis_trait_based_extra1<-data.frame("Trait"=predict(models_FDis_trait_grow1, gam_data[which(gam_data$Site%in%replicate_1_extra_1),])$predicted)
rownames(FDis_trait_based_extra1) <- gam_data$Site[which(gam_data$Site%in%replicate_1_extra_1)]

## - For interpolation group 2 - ##
models_FDis_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_2),])
FDis_trait_based_inter2<-data.frame("Trait"=predict(models_FDis_trait_grow2, gam_data[which(gam_data$Site%in%replicate_1_inter_2),])$predicted)
rownames(FDis_trait_based_inter2) <- gam_data$Site[which(gam_data$Site%in%replicate_1_inter_2)]

## - For extrapolation group 2 - ##
models_FDis_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_2),])
FDis_trait_based_extra2<-data.frame("Trait"=predict(models_FDis_trait_grow2, gam_data[which(gam_data$Site%in%replicate_1_extra_2),])$predicted)
rownames(FDis_trait_based_extra2) <- gam_data$Site[which(gam_data$Site%in%replicate_1_extra_2)]

## - For interpolation group 3 - ##
models_FDis_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_3),])
FDis_trait_based_inter3<-data.frame("Trait"=predict(models_FDis_trait_grow3, gam_data[which(gam_data$Site%in%replicate_1_inter_3),])$predicted)
rownames(FDis_trait_based_inter3) <- gam_data$Site[which(gam_data$Site%in%replicate_1_inter_3)]

## - For extrapolation group 3 - ##
models_FDis_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_3),])
FDis_trait_based_extra3<-data.frame("Trait"=predict(models_FDis_trait_grow3, gam_data[which(gam_data$Site%in%replicate_1_extra_3),])$predicted)
rownames(FDis_trait_based_extra3) <- gam_data$Site[which(gam_data$Site%in%replicate_1_extra_3)]

## - For interpolation group 4 - ##
models_FDis_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_4),])
FDis_trait_based_inter4<-data.frame("Trait"=predict(models_FDis_trait_grow4, gam_data[which(gam_data$Site%in%replicate_1_inter_4),])$predicted)
rownames(FDis_trait_based_inter4) <- gam_data$Site[which(gam_data$Site%in%replicate_1_inter_4)]

## - For extrapolation group 4 - ##
models_FDis_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_4),])
FDis_trait_based_extra4<-data.frame("Trait"=predict(models_FDis_trait_grow4, gam_data[which(gam_data$Site%in%replicate_1_extra_4),])$predicted)
rownames(FDis_trait_based_extra4) <- gam_data$Site[which(gam_data$Site%in%replicate_1_extra_4)]

## assembly in a single prediction table for the 4 groups :
FDis_AFPL_inter_replicate_1 <- rbind(FDis_trait_based_inter1, FDis_trait_based_inter2, FDis_trait_based_inter3, FDis_trait_based_inter4)
FDis_AFPL_inter_replicate_1 <- cbind("Site"=as.numeric(rownames(FDis_AFPL_inter_replicate_1)), FDis_AFPL_inter_replicate_1)
FDis_AFPL_inter_replicate_1 <- FDis_AFPL_inter_replicate_1[order(FDis_AFPL_inter_replicate_1$Site),]
rownames(FDis_AFPL_inter_replicate_1) <- FDis_AFPL_inter_replicate_1$Site

FDis_AFPL_extra_replicate_1 <- rbind(FDis_trait_based_extra1, FDis_trait_based_extra2, FDis_trait_based_extra3, FDis_trait_based_extra4)
FDis_AFPL_extra_replicate_1 <- cbind("Site"=as.numeric(rownames(FDis_AFPL_extra_replicate_1)), FDis_AFPL_extra_replicate_1)
FDis_AFPL_extra_replicate_1 <- FDis_AFPL_extra_replicate_1[order(FDis_AFPL_extra_replicate_1$Site),]
rownames(FDis_AFPL_extra_replicate_1) <- FDis_AFPL_extra_replicate_1$Site

### replicate 2:
## - For interpolation group 1 - ##
models_FDis_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_1),])
FDis_trait_based_inter1<-data.frame("Trait"=predict(models_FDis_trait_grow1, gam_data[which(gam_data$Site%in%replicate_2_inter_1),])$predicted)
rownames(FDis_trait_based_inter1) <- gam_data$Site[which(gam_data$Site%in%replicate_2_inter_1)]

## - For extrapolation group 1 - ##
models_FDis_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_1),])
FDis_trait_based_extra1<-data.frame("Trait"=predict(models_FDis_trait_grow1, gam_data[which(gam_data$Site%in%replicate_2_extra_1),])$predicted)
rownames(FDis_trait_based_extra1) <- gam_data$Site[which(gam_data$Site%in%replicate_2_extra_1)]

## - For interpolation group 2 - ##
models_FDis_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_2),])
FDis_trait_based_inter2<-data.frame("Trait"=predict(models_FDis_trait_grow2, gam_data[which(gam_data$Site%in%replicate_2_inter_2),])$predicted)
rownames(FDis_trait_based_inter2) <- gam_data$Site[which(gam_data$Site%in%replicate_2_inter_2)]

## - For extrapolation group 2 - ##
models_FDis_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_2),])
FDis_trait_based_extra2<-data.frame("Trait"=predict(models_FDis_trait_grow2, gam_data[which(gam_data$Site%in%replicate_2_extra_2),])$predicted)
rownames(FDis_trait_based_extra2) <- gam_data$Site[which(gam_data$Site%in%replicate_2_extra_2)]

## - For interpolation group 3 - ##
models_FDis_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_3),])
FDis_trait_based_inter3<-data.frame("Trait"=predict(models_FDis_trait_grow3, gam_data[which(gam_data$Site%in%replicate_2_inter_3),])$predicted)
rownames(FDis_trait_based_inter3) <- gam_data$Site[which(gam_data$Site%in%replicate_2_inter_3)]

## - For extrapolation group 3 - ##
models_FDis_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_3),])
FDis_trait_based_extra3<-data.frame("Trait"=predict(models_FDis_trait_grow3, gam_data[which(gam_data$Site%in%replicate_2_extra_3),])$predicted)
rownames(FDis_trait_based_extra3) <- gam_data$Site[which(gam_data$Site%in%replicate_2_extra_3)]

## - For interpolation group 4 - ##
models_FDis_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_4),])
FDis_trait_based_inter4<-data.frame("Trait"=predict(models_FDis_trait_grow4, gam_data[which(gam_data$Site%in%replicate_2_inter_4),])$predicted)
rownames(FDis_trait_based_inter4) <- gam_data$Site[which(gam_data$Site%in%replicate_2_inter_4)]

## - For extrapolation group 4 - ##
models_FDis_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_4),])
FDis_trait_based_extra4<-data.frame("Trait"=predict(models_FDis_trait_grow4, gam_data[which(gam_data$Site%in%replicate_2_extra_4),])$predicted)
rownames(FDis_trait_based_extra4) <- gam_data$Site[which(gam_data$Site%in%replicate_2_extra_4)]

## assembly in a single prediction table for the 4 groups :
FDis_AFPL_inter_replicate_2 <- rbind(FDis_trait_based_inter1, FDis_trait_based_inter2, FDis_trait_based_inter3, FDis_trait_based_inter4)
FDis_AFPL_inter_replicate_2 <- cbind("Site"=as.numeric(rownames(FDis_AFPL_inter_replicate_2)), FDis_AFPL_inter_replicate_2)
FDis_AFPL_inter_replicate_2 <- FDis_AFPL_inter_replicate_2[order(FDis_AFPL_inter_replicate_2$Site),]
rownames(FDis_AFPL_inter_replicate_2) <- FDis_AFPL_inter_replicate_2$Site

FDis_AFPL_extra_replicate_2 <- rbind(FDis_trait_based_extra1, FDis_trait_based_extra2, FDis_trait_based_extra3, FDis_trait_based_extra4)
FDis_AFPL_extra_replicate_2 <- cbind("Site"=as.numeric(rownames(FDis_AFPL_extra_replicate_2)), FDis_AFPL_extra_replicate_2)
FDis_AFPL_extra_replicate_2 <- FDis_AFPL_extra_replicate_2[order(FDis_AFPL_extra_replicate_2$Site),]
rownames(FDis_AFPL_extra_replicate_2) <- FDis_AFPL_extra_replicate_2$Site

### replicate 3:
## - For interpolation group 1 - ##
models_FDis_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_1),])
FDis_trait_based_inter1<-data.frame("Trait"=predict(models_FDis_trait_grow1, gam_data[which(gam_data$Site%in%replicate_3_inter_1),])$predicted)
rownames(FDis_trait_based_inter1) <- gam_data$Site[which(gam_data$Site%in%replicate_3_inter_1)]

## - For extrapolation group 1 - ##
models_FDis_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_1),])
FDis_trait_based_extra1<-data.frame("Trait"=predict(models_FDis_trait_grow1, gam_data[which(gam_data$Site%in%replicate_3_extra_1),])$predicted)
rownames(FDis_trait_based_extra1) <- gam_data$Site[which(gam_data$Site%in%replicate_3_extra_1)]

## - For interpolation group 2 - ##
models_FDis_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_2),])
FDis_trait_based_inter2<-data.frame("Trait"=predict(models_FDis_trait_grow2, gam_data[which(gam_data$Site%in%replicate_3_inter_2),])$predicted)
rownames(FDis_trait_based_inter2) <- gam_data$Site[which(gam_data$Site%in%replicate_3_inter_2)]

## - For extrapolation group 2 - ##
models_FDis_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_2),])
FDis_trait_based_extra2<-data.frame("Trait"=predict(models_FDis_trait_grow2, gam_data[which(gam_data$Site%in%replicate_3_extra_2),])$predicted)
rownames(FDis_trait_based_extra2) <- gam_data$Site[which(gam_data$Site%in%replicate_3_extra_2)]

## - For interpolation group 3 - ##
models_FDis_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_3),])
FDis_trait_based_inter3<-data.frame("Trait"=predict(models_FDis_trait_grow3, gam_data[which(gam_data$Site%in%replicate_3_inter_3),])$predicted)
rownames(FDis_trait_based_inter3) <- gam_data$Site[which(gam_data$Site%in%replicate_3_inter_3)]

## - For extrapolation group 3 - ##
models_FDis_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_3),])
FDis_trait_based_extra3<-data.frame("Trait"=predict(models_FDis_trait_grow3, gam_data[which(gam_data$Site%in%replicate_3_extra_3),])$predicted)
rownames(FDis_trait_based_extra3) <- gam_data$Site[which(gam_data$Site%in%replicate_3_extra_3)]

## - For interpolation group 4 - ##
models_FDis_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_4),])
FDis_trait_based_inter4<-data.frame("Trait"=predict(models_FDis_trait_grow4, gam_data[which(gam_data$Site%in%replicate_3_inter_4),])$predicted)
rownames(FDis_trait_based_inter4) <- gam_data$Site[which(gam_data$Site%in%replicate_3_inter_4)]

## - For extrapolation group 4 - ##
models_FDis_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_4),])
FDis_trait_based_extra4<-data.frame("Trait"=predict(models_FDis_trait_grow4, gam_data[which(gam_data$Site%in%replicate_3_extra_4),])$predicted)
rownames(FDis_trait_based_extra4) <- gam_data$Site[which(gam_data$Site%in%replicate_3_extra_4)]

## assembly in a single prediction table for the 4 groups :
FDis_AFPL_inter_replicate_3 <- rbind(FDis_trait_based_inter1, FDis_trait_based_inter2, FDis_trait_based_inter3, FDis_trait_based_inter4)
FDis_AFPL_inter_replicate_3 <- cbind("Site"=as.numeric(rownames(FDis_AFPL_inter_replicate_3)), FDis_AFPL_inter_replicate_3)
FDis_AFPL_inter_replicate_3 <- FDis_AFPL_inter_replicate_3[order(FDis_AFPL_inter_replicate_3$Site),]
rownames(FDis_AFPL_inter_replicate_3) <- FDis_AFPL_inter_replicate_3$Site

FDis_AFPL_extra_replicate_3 <- rbind(FDis_trait_based_extra1, FDis_trait_based_extra2, FDis_trait_based_extra3, FDis_trait_based_extra4)
FDis_AFPL_extra_replicate_3 <- cbind("Site"=as.numeric(rownames(FDis_AFPL_extra_replicate_3)), FDis_AFPL_extra_replicate_3)
FDis_AFPL_extra_replicate_3 <- FDis_AFPL_extra_replicate_3[order(FDis_AFPL_extra_replicate_3$Site),]
rownames(FDis_AFPL_extra_replicate_3) <- FDis_AFPL_extra_replicate_3$Site

### replicate 4:
## - For interpolation group 1 - ##
models_FDis_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_1),])
FDis_trait_based_inter1<-data.frame("Trait"=predict(models_FDis_trait_grow1, gam_data[which(gam_data$Site%in%replicate_4_inter_1),])$predicted)
rownames(FDis_trait_based_inter1) <- gam_data$Site[which(gam_data$Site%in%replicate_4_inter_1)]

## - For extrapolation group 1 - ##
models_FDis_trait_grow1<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_1),])
FDis_trait_based_extra1<-data.frame("Trait"=predict(models_FDis_trait_grow1, gam_data[which(gam_data$Site%in%replicate_4_extra_1),])$predicted)
rownames(FDis_trait_based_extra1) <- gam_data$Site[which(gam_data$Site%in%replicate_4_extra_1)]

## - For interpolation group 2 - ##
models_FDis_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_2),])
FDis_trait_based_inter2<-data.frame("Trait"=predict(models_FDis_trait_grow2, gam_data[which(gam_data$Site%in%replicate_4_inter_2),])$predicted)
rownames(FDis_trait_based_inter2) <- gam_data$Site[which(gam_data$Site%in%replicate_4_inter_2)]

## - For extrapolation group 2 - ##
models_FDis_trait_grow2<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_2),])
FDis_trait_based_extra2<-data.frame("Trait"=predict(models_FDis_trait_grow2, gam_data[which(gam_data$Site%in%replicate_4_extra_2),])$predicted)
rownames(FDis_trait_based_extra2) <- gam_data$Site[which(gam_data$Site%in%replicate_4_extra_2)]

## - For interpolation group 3 - ##
models_FDis_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_3),])
FDis_trait_based_inter3<-data.frame("Trait"=predict(models_FDis_trait_grow3, gam_data[which(gam_data$Site%in%replicate_4_inter_3),])$predicted)
rownames(FDis_trait_based_inter3) <- gam_data$Site[which(gam_data$Site%in%replicate_4_inter_3)]

## - For extrapolation group 3 - ##
models_FDis_trait_grow3<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_3),])
FDis_trait_based_extra3<-data.frame("Trait"=predict(models_FDis_trait_grow3, gam_data[which(gam_data$Site%in%replicate_4_extra_3),])$predicted)
rownames(FDis_trait_based_extra3) <- gam_data$Site[which(gam_data$Site%in%replicate_4_extra_3)]

## - For interpolation group 4 - ##
models_FDis_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_4),])
FDis_trait_based_inter4<-data.frame("Trait"=predict(models_FDis_trait_grow4, gam_data[which(gam_data$Site%in%replicate_4_inter_4),])$predicted)
rownames(FDis_trait_based_inter4) <- gam_data$Site[which(gam_data$Site%in%replicate_4_inter_4)]

## - For extrapolation group 4 - ##
models_FDis_trait_grow4<- rfsrc(Trait ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_4),])
FDis_trait_based_extra4<-data.frame("Trait"=predict(models_FDis_trait_grow4, gam_data[which(gam_data$Site%in%replicate_4_extra_4),])$predicted)
rownames(FDis_trait_based_extra4) <- gam_data$Site[which(gam_data$Site%in%replicate_4_extra_4)]

## assembly in a single prediction table for the 4 groups :
FDis_AFPL_inter_replicate_4 <- rbind(FDis_trait_based_inter1, FDis_trait_based_inter2, FDis_trait_based_inter3, FDis_trait_based_inter4)
FDis_AFPL_inter_replicate_4 <- cbind("Site"=as.numeric(rownames(FDis_AFPL_inter_replicate_4)), FDis_AFPL_inter_replicate_4)
FDis_AFPL_inter_replicate_4 <- FDis_AFPL_inter_replicate_4[order(FDis_AFPL_inter_replicate_4$Site),]
rownames(FDis_AFPL_inter_replicate_4) <- FDis_AFPL_inter_replicate_4$Site

FDis_AFPL_extra_replicate_4 <- rbind(FDis_trait_based_extra1, FDis_trait_based_extra2, FDis_trait_based_extra3, FDis_trait_based_extra4)
FDis_AFPL_extra_replicate_4 <- cbind("Site"=as.numeric(rownames(FDis_AFPL_extra_replicate_4)), FDis_AFPL_extra_replicate_4)
FDis_AFPL_extra_replicate_4 <- FDis_AFPL_extra_replicate_4[order(FDis_AFPL_extra_replicate_4$Site),]
rownames(FDis_AFPL_extra_replicate_4) <- FDis_AFPL_extra_replicate_4$Site


######################################
#### -- SPECIES-BASED APPROACH -- ####
######################################
### replicate 1:
gam_data <- data.frame("Site"=as.numeric(rownames(com_replicate_1)), com_replicate_1)
gam_data <- left_join(gam_data,env_simul_obs, by="Site") # assembly of species recovery data and environmental variables.
####### -- Species Distribution Models for the Species-based approach-- ####
###### -- PREDICTION OF SPECIES ABUNDANCE ######
#### -- This computes path predicts the relative abundance of species in the context of interpolation -- ####
## - For interpolation group 1 - ##
species_models_recovery_grow1<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_1),])
species_models_recovery_inter1 <- predict(species_models_recovery_grow1, gam_data[which(gam_data$Site%in%replicate_1_inter_1),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_1),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter1 <- cbind(species_models_recovery_inter1, predict(model, gam_data[which(gam_data$Site%in%replicate_1_inter_1),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter1) <- colnames(com_replicate_1)
rownames(species_models_recovery_inter1) <- gam_data$Site[which(gam_data$Site%in%replicate_1_inter_1)]

## - For extrapolation group 1 - ##
species_models_recovery_grow1<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_1),])
species_models_recovery_extra1 <- predict(species_models_recovery_grow1, gam_data[which(gam_data$Site%in%replicate_1_extra_1),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_1),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra1 <- cbind(species_models_recovery_extra1, predict(model, gam_data[which(gam_data$Site%in%replicate_1_extra_1),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra1) <- colnames(com_replicate_1)
rownames(species_models_recovery_extra1) <- gam_data$Site[which(gam_data$Site%in%replicate_1_extra_1)]

## - For interpolation group 2 - ##
species_models_recovery_grow2<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_2),])
species_models_recovery_inter2 <- predict(species_models_recovery_grow2, gam_data[which(gam_data$Site%in%replicate_1_inter_2),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_2),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter2 <- cbind(species_models_recovery_inter2, predict(model, gam_data[which(gam_data$Site%in%replicate_1_inter_2),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter2) <- colnames(com_replicate_1)
rownames(species_models_recovery_inter2) <- gam_data$Site[which(gam_data$Site%in%replicate_1_inter_2)]

## - For extrapolation group 2 - ##
species_models_recovery_grow2<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_2),])
species_models_recovery_extra2 <- predict(species_models_recovery_grow2, gam_data[which(gam_data$Site%in%replicate_1_extra_2),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_2),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra2 <- cbind(species_models_recovery_extra2, predict(model, gam_data[which(gam_data$Site%in%replicate_1_extra_2),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra2) <- colnames(com_replicate_1)
rownames(species_models_recovery_extra2) <- gam_data$Site[which(gam_data$Site%in%replicate_1_extra_2)]

## - For interpolation group 3 - ##
species_models_recovery_grow3<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_3),])
species_models_recovery_inter3 <- predict(species_models_recovery_grow3, gam_data[which(gam_data$Site%in%replicate_1_inter_3),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_3),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter3 <- cbind(species_models_recovery_inter3, predict(model, gam_data[which(gam_data$Site%in%replicate_1_inter_3),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter3) <- colnames(com_replicate_1)
rownames(species_models_recovery_inter3) <- gam_data$Site[which(gam_data$Site%in%replicate_1_inter_3)]

## - For extrapolation group 3 - ##
species_models_recovery_grow3<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_3),])
species_models_recovery_extra3 <- predict(species_models_recovery_grow3, gam_data[which(gam_data$Site%in%replicate_1_extra_3),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_3),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra3 <- cbind(species_models_recovery_extra3, predict(model, gam_data[which(gam_data$Site%in%replicate_1_extra_3),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra3) <- colnames(com_replicate_1)
rownames(species_models_recovery_extra3) <- gam_data$Site[which(gam_data$Site%in%replicate_1_extra_3)]

## - For interpolation group 4 - ##
species_models_recovery_grow4<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_4),])
species_models_recovery_inter4 <- predict(species_models_recovery_grow4, gam_data[which(gam_data$Site%in%replicate_1_inter_4),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_inter_4),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter4 <- cbind(species_models_recovery_inter4, predict(model, gam_data[which(gam_data$Site%in%replicate_1_inter_4),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter4) <- colnames(com_replicate_1)
rownames(species_models_recovery_inter4) <- gam_data$Site[which(gam_data$Site%in%replicate_1_inter_4)]

## - For extrapolation group 4 - ##
species_models_recovery_grow4<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_4),])
species_models_recovery_extra4 <- predict(species_models_recovery_grow4, gam_data[which(gam_data$Site%in%replicate_1_extra_4),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_1_extra_4),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra4 <- cbind(species_models_recovery_extra4, predict(model, gam_data[which(gam_data$Site%in%replicate_1_extra_4),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra4) <- colnames(com_replicate_1)
rownames(species_models_recovery_extra4) <- gam_data$Site[which(gam_data$Site%in%replicate_1_extra_4)]

## assembly in a single prediction table for the 4 groups :
species_models_inter_replicate_1 <- data.frame(rbind(species_models_recovery_inter1, species_models_recovery_inter2, species_models_recovery_inter3, species_models_recovery_inter4))
colnames(species_models_inter_replicate_1)<- colnames(species_models_recovery_inter1)
species_models_inter_replicate_1 <- cbind("Site"=as.numeric(rownames(species_models_inter_replicate_1)), species_models_inter_replicate_1)
rownames(species_models_inter_replicate_1)=species_models_inter_replicate_1$Site
species_models_inter_replicate_1 <- species_models_inter_replicate_1[,-1]

CWM_PFAL_inter_replicate_1 <- CWM_simul(species_models_inter_replicate_1)
FDis_PFAL_inter_replicate_1 <- FDis_simul(species_models_inter_replicate_1)

species_models_extra_replicate_1 <- data.frame(rbind(species_models_recovery_extra1, species_models_recovery_extra2, species_models_recovery_extra3, species_models_recovery_extra4))
colnames(species_models_extra_replicate_1)<- colnames(species_models_recovery_extra1)
species_models_extra_replicate_1 <- cbind("Site"=as.numeric(rownames(species_models_extra_replicate_1)), species_models_extra_replicate_1)
rownames(species_models_extra_replicate_1)=species_models_extra_replicate_1$Site
species_models_extra_replicate_1 <- species_models_extra_replicate_1[,-1]

CWM_PFAL_extra_replicate_1 <- CWM_simul(species_models_extra_replicate_1)
FDis_PFAL_extra_replicate_1 <- FDis_simul(species_models_extra_replicate_1)

### Replicate 2:
gam_data <- data.frame("Site"=as.numeric(rownames(com_replicate_2)), com_replicate_2)
gam_data <- left_join(gam_data,env_simul_obs, by="Site") # assembly of species recovery data and environmental variables.
####### -- Species Distribution Models for the Species-based approach-- ####
###### -- PREDICTION OF SPECIES ABUNDANCE ######
#### -- This computes path predicts the relative abundance of species in the context of interpolation -- ####
## - For interpolation group 1 - ##
species_models_recovery_grow1<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_1),])
species_models_recovery_inter1 <- predict(species_models_recovery_grow1, gam_data[which(gam_data$Site%in%replicate_2_inter_1),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_1),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter1 <- cbind(species_models_recovery_inter1, predict(model, gam_data[which(gam_data$Site%in%replicate_2_inter_1),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter1) <- colnames(com_replicate_2)
rownames(species_models_recovery_inter1) <- gam_data$Site[which(gam_data$Site%in%replicate_2_inter_1)]

## - For extrapolation group 1 - ##
species_models_recovery_grow1<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_1),])
species_models_recovery_extra1 <- predict(species_models_recovery_grow1, gam_data[which(gam_data$Site%in%replicate_2_extra_1),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_1),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra1 <- cbind(species_models_recovery_extra1, predict(model, gam_data[which(gam_data$Site%in%replicate_2_extra_1),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra1) <- colnames(com_replicate_2)
rownames(species_models_recovery_extra1) <- gam_data$Site[which(gam_data$Site%in%replicate_2_extra_1)]

## - For interpolation group 2 - ##
species_models_recovery_grow2<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_2),])
species_models_recovery_inter2 <- predict(species_models_recovery_grow2, gam_data[which(gam_data$Site%in%replicate_2_inter_2),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_2),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter2 <- cbind(species_models_recovery_inter2, predict(model, gam_data[which(gam_data$Site%in%replicate_2_inter_2),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter2) <- colnames(com_replicate_2)
rownames(species_models_recovery_inter2) <- gam_data$Site[which(gam_data$Site%in%replicate_2_inter_2)]

## - For extrapolation group 2 - ##
species_models_recovery_grow2<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_2),])
species_models_recovery_extra2 <- predict(species_models_recovery_grow2, gam_data[which(gam_data$Site%in%replicate_2_extra_2),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_2),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra2 <- cbind(species_models_recovery_extra2, predict(model, gam_data[which(gam_data$Site%in%replicate_2_extra_2),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra2) <- colnames(com_replicate_2)
rownames(species_models_recovery_extra2) <- gam_data$Site[which(gam_data$Site%in%replicate_2_extra_2)]

## - For interpolation group 3 - ##
species_models_recovery_grow3<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_3),])
species_models_recovery_inter3 <- predict(species_models_recovery_grow3, gam_data[which(gam_data$Site%in%replicate_2_inter_3),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_3),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter3 <- cbind(species_models_recovery_inter3, predict(model, gam_data[which(gam_data$Site%in%replicate_2_inter_3),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter3) <- colnames(com_replicate_2)
rownames(species_models_recovery_inter3) <- gam_data$Site[which(gam_data$Site%in%replicate_2_inter_3)]

## - For extrapolation group 3 - ##
species_models_recovery_grow3<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_3),])
species_models_recovery_extra3 <- predict(species_models_recovery_grow3, gam_data[which(gam_data$Site%in%replicate_2_extra_3),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_3),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra3 <- cbind(species_models_recovery_extra3, predict(model, gam_data[which(gam_data$Site%in%replicate_2_extra_3),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra3) <- colnames(com_replicate_2)
rownames(species_models_recovery_extra3) <- gam_data$Site[which(gam_data$Site%in%replicate_2_extra_3)]

## - For interpolation group 4 - ##
species_models_recovery_grow4<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_4),])
species_models_recovery_inter4 <- predict(species_models_recovery_grow4, gam_data[which(gam_data$Site%in%replicate_2_inter_4),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_inter_4),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter4 <- cbind(species_models_recovery_inter4, predict(model, gam_data[which(gam_data$Site%in%replicate_2_inter_4),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter4) <- colnames(com_replicate_2)
rownames(species_models_recovery_inter4) <- gam_data$Site[which(gam_data$Site%in%replicate_2_inter_4)]

## - For extrapolation group 4 - ##
species_models_recovery_grow4<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_4),])
species_models_recovery_extra4 <- predict(species_models_recovery_grow4, gam_data[which(gam_data$Site%in%replicate_2_extra_4),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_2_extra_4),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra4 <- cbind(species_models_recovery_extra4, predict(model, gam_data[which(gam_data$Site%in%replicate_2_extra_4),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra4) <- colnames(com_replicate_2)
rownames(species_models_recovery_extra4) <- gam_data$Site[which(gam_data$Site%in%replicate_2_extra_4)]

## assembly in a single prediction table for the 4 groups :
species_models_inter_replicate_2 <- data.frame(rbind(species_models_recovery_inter1, species_models_recovery_inter2, species_models_recovery_inter3, species_models_recovery_inter4))
colnames(species_models_inter_replicate_2)<- colnames(species_models_recovery_inter1)
species_models_inter_replicate_2 <- cbind("Site"=as.numeric(rownames(species_models_inter_replicate_2)), species_models_inter_replicate_2)
rownames(species_models_inter_replicate_2)=species_models_inter_replicate_2$Site
species_models_inter_replicate_2 <- species_models_inter_replicate_2[,-1]

CWM_PFAL_inter_replicate_2 <- CWM_simul(species_models_inter_replicate_2)
FDis_PFAL_inter_replicate_2 <- FDis_simul(species_models_inter_replicate_2)

species_models_extra_replicate_2 <- data.frame(rbind(species_models_recovery_extra1, species_models_recovery_extra2, species_models_recovery_extra3, species_models_recovery_extra4))
colnames(species_models_extra_replicate_2)<- colnames(species_models_recovery_extra1)
species_models_extra_replicate_2 <- cbind("Site"=as.numeric(rownames(species_models_extra_replicate_2)), species_models_extra_replicate_2)
rownames(species_models_extra_replicate_2)=species_models_extra_replicate_2$Site
species_models_extra_replicate_2 <- species_models_extra_replicate_2[,-1]

CWM_PFAL_extra_replicate_2 <- CWM_simul(species_models_extra_replicate_2)
FDis_PFAL_extra_replicate_2 <- FDis_simul(species_models_extra_replicate_2)

### Replicate 3:
gam_data <- data.frame("Site"=as.numeric(rownames(com_replicate_3)), com_replicate_3)
gam_data <- left_join(gam_data,env_simul_obs, by="Site") # assembly of species recovery data and environmental variables.
####### -- Species Distribution Models for the Species-based approach-- ####
###### -- PREDICTION OF SPECIES ABUNDANCE ######
#### -- This computes path predicts the relative abundance of species in the context of interpolation -- ####
## - For interpolation group 1 - ##
species_models_recovery_grow1<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_1),])
species_models_recovery_inter1 <- predict(species_models_recovery_grow1, gam_data[which(gam_data$Site%in%replicate_3_inter_1),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_1),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter1 <- cbind(species_models_recovery_inter1, predict(model, gam_data[which(gam_data$Site%in%replicate_3_inter_1),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter1) <- colnames(com_replicate_3)
rownames(species_models_recovery_inter1) <- gam_data$Site[which(gam_data$Site%in%replicate_3_inter_1)]

## - For extrapolation group 1 - ##
species_models_recovery_grow1<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_1),])
species_models_recovery_extra1 <- predict(species_models_recovery_grow1, gam_data[which(gam_data$Site%in%replicate_3_extra_1),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_1),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra1 <- cbind(species_models_recovery_extra1, predict(model, gam_data[which(gam_data$Site%in%replicate_3_extra_1),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra1) <- colnames(com_replicate_3)
rownames(species_models_recovery_extra1) <- gam_data$Site[which(gam_data$Site%in%replicate_3_extra_1)]

## - For interpolation group 2 - ##
species_models_recovery_grow2<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_2),])
species_models_recovery_inter2 <- predict(species_models_recovery_grow2, gam_data[which(gam_data$Site%in%replicate_3_inter_2),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_2),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter2 <- cbind(species_models_recovery_inter2, predict(model, gam_data[which(gam_data$Site%in%replicate_3_inter_2),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter2) <- colnames(com_replicate_3)
rownames(species_models_recovery_inter2) <- gam_data$Site[which(gam_data$Site%in%replicate_3_inter_2)]

## - For extrapolation group 2 - ##
species_models_recovery_grow2<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_2),])
species_models_recovery_extra2 <- predict(species_models_recovery_grow2, gam_data[which(gam_data$Site%in%replicate_3_extra_2),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_2),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra2 <- cbind(species_models_recovery_extra2, predict(model, gam_data[which(gam_data$Site%in%replicate_3_extra_2),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra2) <- colnames(com_replicate_3)
rownames(species_models_recovery_extra2) <- gam_data$Site[which(gam_data$Site%in%replicate_3_extra_2)]

## - For interpolation group 3 - ##
species_models_recovery_grow3<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_3),])
species_models_recovery_inter3 <- predict(species_models_recovery_grow3, gam_data[which(gam_data$Site%in%replicate_3_inter_3),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_3),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter3 <- cbind(species_models_recovery_inter3, predict(model, gam_data[which(gam_data$Site%in%replicate_3_inter_3),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter3) <- colnames(com_replicate_3)
rownames(species_models_recovery_inter3) <- gam_data$Site[which(gam_data$Site%in%replicate_3_inter_3)]

## - For extrapolation group 3 - ##
species_models_recovery_grow3<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_3),])
species_models_recovery_extra3 <- predict(species_models_recovery_grow3, gam_data[which(gam_data$Site%in%replicate_3_extra_3),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_3),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra3 <- cbind(species_models_recovery_extra3, predict(model, gam_data[which(gam_data$Site%in%replicate_3_extra_3),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra3) <- colnames(com_replicate_3)
rownames(species_models_recovery_extra3) <- gam_data$Site[which(gam_data$Site%in%replicate_3_extra_3)]

## - For interpolation group 4 - ##
species_models_recovery_grow4<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_4),])
species_models_recovery_inter4 <- predict(species_models_recovery_grow4, gam_data[which(gam_data$Site%in%replicate_3_inter_4),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_inter_4),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter4 <- cbind(species_models_recovery_inter4, predict(model, gam_data[which(gam_data$Site%in%replicate_3_inter_4),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter4) <- colnames(com_replicate_3)
rownames(species_models_recovery_inter4) <- gam_data$Site[which(gam_data$Site%in%replicate_3_inter_4)]

## - For extrapolation group 4 - ##
species_models_recovery_grow4<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_4),])
species_models_recovery_extra4 <- predict(species_models_recovery_grow4, gam_data[which(gam_data$Site%in%replicate_3_extra_4),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_3_extra_4),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra4 <- cbind(species_models_recovery_extra4, predict(model, gam_data[which(gam_data$Site%in%replicate_3_extra_4),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra4) <- colnames(com_replicate_3)
rownames(species_models_recovery_extra4) <- gam_data$Site[which(gam_data$Site%in%replicate_3_extra_4)]

## assembly in a single prediction table for the 4 groups :
species_models_inter_replicate_3 <- data.frame(rbind(species_models_recovery_inter1, species_models_recovery_inter2, species_models_recovery_inter3, species_models_recovery_inter4))
colnames(species_models_inter_replicate_3)<- colnames(species_models_recovery_inter1)
species_models_inter_replicate_3 <- cbind("Site"=as.numeric(rownames(species_models_inter_replicate_3)), species_models_inter_replicate_3)
rownames(species_models_inter_replicate_3)=species_models_inter_replicate_3$Site
species_models_inter_replicate_3 <- species_models_inter_replicate_3[,-1]

CWM_PFAL_inter_replicate_3 <- CWM_simul(species_models_inter_replicate_3)
FDis_PFAL_inter_replicate_3 <- FDis_simul(species_models_inter_replicate_3)

species_models_extra_replicate_3 <- data.frame(rbind(species_models_recovery_extra1, species_models_recovery_extra2, species_models_recovery_extra3, species_models_recovery_extra4))
colnames(species_models_extra_replicate_3)<- colnames(species_models_recovery_extra1)
species_models_extra_replicate_3 <- cbind("Site"=as.numeric(rownames(species_models_extra_replicate_3)), species_models_extra_replicate_3)
rownames(species_models_extra_replicate_3)=species_models_extra_replicate_3$Site
species_models_extra_replicate_3 <- species_models_extra_replicate_3[,-1]

CWM_PFAL_extra_replicate_3 <- CWM_simul(species_models_extra_replicate_3)
FDis_PFAL_extra_replicate_3 <- FDis_simul(species_models_extra_replicate_3)

### Replicate 4:
gam_data <- data.frame("Site"=as.numeric(rownames(com_replicate_4)), com_replicate_4)
gam_data <- left_join(gam_data,env_simul_obs, by="Site") # assembly of species recovery data and environmental variables.
####### -- Species Distribution Models for the Species-based approach-- ####
###### -- PREDICTION OF SPECIES ABUNDANCE ######
#### -- This computes path predicts the relative abundance of species in the context of interpolation -- ####
## - For interpolation group 1 - ##
species_models_recovery_grow1<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_1),])
species_models_recovery_inter1 <- predict(species_models_recovery_grow1, gam_data[which(gam_data$Site%in%replicate_4_inter_1),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_1),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter1 <- cbind(species_models_recovery_inter1, predict(model, gam_data[which(gam_data$Site%in%replicate_4_inter_1),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter1) <- colnames(com_replicate_4)
rownames(species_models_recovery_inter1) <- gam_data$Site[which(gam_data$Site%in%replicate_4_inter_1)]

## - For extrapolation group 1 - ##
species_models_recovery_grow1<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_1),])
species_models_recovery_extra1 <- predict(species_models_recovery_grow1, gam_data[which(gam_data$Site%in%replicate_4_extra_1),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_1),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra1 <- cbind(species_models_recovery_extra1, predict(model, gam_data[which(gam_data$Site%in%replicate_4_extra_1),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra1) <- colnames(com_replicate_4)
rownames(species_models_recovery_extra1) <- gam_data$Site[which(gam_data$Site%in%replicate_4_extra_1)]

## - For interpolation group 2 - ##
species_models_recovery_grow2<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_2),])
species_models_recovery_inter2 <- predict(species_models_recovery_grow2, gam_data[which(gam_data$Site%in%replicate_4_inter_2),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_2),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter2 <- cbind(species_models_recovery_inter2, predict(model, gam_data[which(gam_data$Site%in%replicate_4_inter_2),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter2) <- colnames(com_replicate_4)
rownames(species_models_recovery_inter2) <- gam_data$Site[which(gam_data$Site%in%replicate_4_inter_2)]

## - For extrapolation group 2 - ##
species_models_recovery_grow2<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_2),])
species_models_recovery_extra2 <- predict(species_models_recovery_grow2, gam_data[which(gam_data$Site%in%replicate_4_extra_2),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_2),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra2 <- cbind(species_models_recovery_extra2, predict(model, gam_data[which(gam_data$Site%in%replicate_4_extra_2),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra2) <- colnames(com_replicate_4)
rownames(species_models_recovery_extra2) <- gam_data$Site[which(gam_data$Site%in%replicate_4_extra_2)]

## - For interpolation group 3 - ##
species_models_recovery_grow3<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_3),])
species_models_recovery_inter3 <- predict(species_models_recovery_grow3, gam_data[which(gam_data$Site%in%replicate_4_inter_3),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_3),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter3 <- cbind(species_models_recovery_inter3, predict(model, gam_data[which(gam_data$Site%in%replicate_4_inter_3),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter3) <- colnames(com_replicate_4)
rownames(species_models_recovery_inter3) <- gam_data$Site[which(gam_data$Site%in%replicate_4_inter_3)]

## - For extrapolation group 3 - ##
species_models_recovery_grow3<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_3),])
species_models_recovery_extra3 <- predict(species_models_recovery_grow3, gam_data[which(gam_data$Site%in%replicate_4_extra_3),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_3),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra3 <- cbind(species_models_recovery_extra3, predict(model, gam_data[which(gam_data$Site%in%replicate_4_extra_3),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra3) <- colnames(com_replicate_4)
rownames(species_models_recovery_extra3) <- gam_data$Site[which(gam_data$Site%in%replicate_4_extra_3)]

## - For interpolation group 4 - ##
species_models_recovery_grow4<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_4),])
species_models_recovery_inter4 <- predict(species_models_recovery_grow4, gam_data[which(gam_data$Site%in%replicate_4_inter_4),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_inter_4),])"))) #Abundance modeling for each species at sites in interpolation groups 2, 3, and 4. 
  species_models_recovery_inter4 <- cbind(species_models_recovery_inter4, predict(model, gam_data[which(gam_data$Site%in%replicate_4_inter_4),])$predicted) # prediction of abundance for each species at sites i interpolation group 1.
}
colnames(species_models_recovery_inter4) <- colnames(com_replicate_4)
rownames(species_models_recovery_inter4) <- gam_data$Site[which(gam_data$Site%in%replicate_4_inter_4)]

## - For extrapolation group 4 - ##
species_models_recovery_grow4<- rfsrc(colnames(gam_data)[2]~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_4),])
species_models_recovery_extra4 <- predict(species_models_recovery_grow4, gam_data[which(gam_data$Site%in%replicate_4_extra_4),])$predicted
for (i in 3 : (dim(gam_data)[2]-1)){ 
  eval(parse(text=paste0("model=rfsrc(",colnames(gam_data)[i]," ~ Env, data=gam_data[-which(gam_data$Site%in%replicate_4_extra_4),])"))) #Abundance modeling for each species at sites in extrapolation groups 2, 3, and 4. 
  species_models_recovery_extra4 <- cbind(species_models_recovery_extra4, predict(model, gam_data[which(gam_data$Site%in%replicate_4_extra_4),])$predicted) # prediction of abundance for each species at sites i extrapolation group 1.
}
colnames(species_models_recovery_extra4) <- colnames(com_replicate_4)
rownames(species_models_recovery_extra4) <- gam_data$Site[which(gam_data$Site%in%replicate_4_extra_4)]

## assembly in a single prediction table for the 4 groups :
species_models_inter_replicate_4 <- data.frame(rbind(species_models_recovery_inter1, species_models_recovery_inter2, species_models_recovery_inter3, species_models_recovery_inter4))
colnames(species_models_inter_replicate_4)<- colnames(species_models_recovery_inter1)
species_models_inter_replicate_4 <- cbind("Site"=as.numeric(rownames(species_models_inter_replicate_4)), species_models_inter_replicate_4)
rownames(species_models_inter_replicate_4)=species_models_inter_replicate_4$Site
species_models_inter_replicate_4 <- species_models_inter_replicate_4[,-1]

CWM_PFAL_inter_replicate_4 <- CWM_simul(species_models_inter_replicate_4)
FDis_PFAL_inter_replicate_4 <- FDis_simul(species_models_inter_replicate_4)

species_models_extra_replicate_4 <- data.frame(rbind(species_models_recovery_extra1, species_models_recovery_extra2, species_models_recovery_extra3, species_models_recovery_extra4))
colnames(species_models_extra_replicate_4)<- colnames(species_models_recovery_extra1)
species_models_extra_replicate_4 <- cbind("Site"=as.numeric(rownames(species_models_extra_replicate_4)), species_models_extra_replicate_4)
rownames(species_models_extra_replicate_4)=species_models_extra_replicate_4$Site
species_models_extra_replicate_4 <- species_models_extra_replicate_4[,-1]

CWM_PFAL_extra_replicate_4 <- CWM_simul(species_models_extra_replicate_4)
FDis_PFAL_extra_replicate_4 <- FDis_simul(species_models_extra_replicate_4)


### total predicted species communities : 
species_models_inter_tot <- full_join(species_models_inter_replicate_1, species_models_inter_replicate_2)
species_models_inter_tot <- full_join(species_models_inter_tot, species_models_inter_replicate_3)
species_models_inter_tot <- full_join(species_models_inter_tot, species_models_inter_replicate_4)
rownames(species_models_inter_tot)<- c(rownames(species_models_inter_replicate_1), rownames(species_models_inter_replicate_2),rownames(species_models_inter_replicate_3),rownames(species_models_inter_replicate_4))
species_models_inter_tot <- cbind("Site"=as.numeric(rownames(species_models_inter_tot)), species_models_inter_tot)
species_models_inter_tot <- species_models_inter_tot[order(species_models_inter_tot$Site), -1]
for (i in 1: 1200){
  for (j in 1 :dim(species_models_inter_tot)[2]){
    if (is.na(species_models_inter_tot[i,j])){
      species_models_inter_tot[i,j]=0
    }
  }
}

species_models_inter_tot <- cbind("Site"=as.numeric(rownames(species_models_inter_tot)), species_models_inter_tot)

species_models_extra_tot <- full_join(species_models_extra_replicate_1, species_models_extra_replicate_2)
species_models_extra_tot <- full_join(species_models_extra_tot, species_models_extra_replicate_3)
species_models_extra_tot <- full_join(species_models_extra_tot, species_models_extra_replicate_4)
rownames(species_models_extra_tot)<- c(rownames(species_models_extra_replicate_1), rownames(species_models_extra_replicate_2),rownames(species_models_extra_replicate_3),rownames(species_models_extra_replicate_4))
species_models_extra_tot <- cbind("Site"=as.numeric(rownames(species_models_extra_tot)), species_models_extra_tot)
species_models_extra_tot <- species_models_extra_tot[order(species_models_extra_tot$Site), -1]
for (i in 1: dim(species_models_extra_tot)[1]){
  for (j in 1 :dim(species_models_extra_tot)[2]){
    if (is.na(species_models_extra_tot[i,j])){
      species_models_extra_tot[i,j]=0
    }
  }
}

species_models_extra_tot <- cbind("Site"=as.numeric(rownames(species_models_extra_tot)), species_models_extra_tot)


#### Observed and predicted functional indices
FDis_extra_replicate_1 <- left_join(FDis_PFAL_extra_replicate_1, FDis_simul_env, by="Site")
FDis_extra_replicate_1 <- left_join(FDis_extra_replicate_1, FDis_AFPL_extra_replicate_1, by="Site")
colnames(FDis_extra_replicate_1)<-c("Site","PFAL", "obs","AFPL")

FDis_inter_replicate_1 <- left_join(FDis_PFAL_inter_replicate_1, FDis_simul_env, by="Site")
FDis_inter_replicate_1 <- left_join(FDis_inter_replicate_1, FDis_AFPL_inter_replicate_1, by="Site")
colnames(FDis_inter_replicate_1)<-c("Site","PFAL", "obs","AFPL")

CWM_extra_replicate_1 <- left_join(CWM_PFAL_extra_replicate_1, CWM_simul_env, by="Site")
CWM_extra_replicate_1 <- left_join(CWM_extra_replicate_1, CWM_AFPL_extra_replicate_1, by="Site")
colnames(CWM_extra_replicate_1)<-c("Site", "PFAL", "obs","AFPL")

CWM_inter_replicate_1 <- left_join(CWM_PFAL_inter_replicate_1, CWM_simul_env, by="Site")
CWM_inter_replicate_1 <- left_join(CWM_inter_replicate_1, CWM_AFPL_inter_replicate_1, by="Site")
colnames(CWM_inter_replicate_1)<-c("Site", "PFAL", "obs","AFPL")

FDis_extra_replicate_2 <- left_join(FDis_PFAL_extra_replicate_2, FDis_simul_env, by="Site")
FDis_extra_replicate_2 <- left_join(FDis_extra_replicate_2, FDis_AFPL_extra_replicate_2, by="Site")
colnames(FDis_extra_replicate_2)<-c("Site","PFAL", "obs","AFPL")

FDis_inter_replicate_2 <- left_join(FDis_PFAL_inter_replicate_2, FDis_simul_env, by="Site")
FDis_inter_replicate_2 <- left_join(FDis_inter_replicate_2, FDis_AFPL_inter_replicate_2, by="Site")
colnames(FDis_inter_replicate_2)<-c("Site","PFAL", "obs","AFPL")

CWM_extra_replicate_2 <- left_join(CWM_PFAL_extra_replicate_2, CWM_simul_env, by="Site")
CWM_extra_replicate_2 <- left_join(CWM_extra_replicate_2, CWM_AFPL_extra_replicate_2, by="Site")
colnames(CWM_extra_replicate_2)<-c("Site", "PFAL", "obs","AFPL")

CWM_inter_replicate_2 <- left_join(CWM_PFAL_inter_replicate_2, CWM_simul_env, by="Site")
CWM_inter_replicate_2 <- left_join(CWM_inter_replicate_2, CWM_AFPL_inter_replicate_2, by="Site")
colnames(CWM_inter_replicate_2)<-c("Site", "PFAL", "obs","AFPL")

FDis_extra_replicate_3 <- left_join(FDis_PFAL_extra_replicate_3, FDis_simul_env, by="Site")
FDis_extra_replicate_3 <- left_join(FDis_extra_replicate_3, FDis_AFPL_extra_replicate_3, by="Site")
colnames(FDis_extra_replicate_3)<-c("Site","PFAL", "obs","AFPL")

FDis_inter_replicate_3 <- left_join(FDis_PFAL_inter_replicate_3, FDis_simul_env, by="Site")
FDis_inter_replicate_3 <- left_join(FDis_inter_replicate_3, FDis_AFPL_inter_replicate_3, by="Site")
colnames(FDis_inter_replicate_3)<-c("Site","PFAL", "obs","AFPL")

CWM_extra_replicate_3 <- left_join(CWM_PFAL_extra_replicate_3, CWM_simul_env, by="Site")
CWM_extra_replicate_3 <- left_join(CWM_extra_replicate_3, CWM_AFPL_extra_replicate_3, by="Site")
colnames(CWM_extra_replicate_3)<-c("Site", "PFAL", "obs","AFPL")

CWM_inter_replicate_3 <- left_join(CWM_PFAL_inter_replicate_3, CWM_simul_env, by="Site")
CWM_inter_replicate_3 <- left_join(CWM_inter_replicate_3, CWM_AFPL_inter_replicate_3, by="Site")
colnames(CWM_inter_replicate_3)<-c("Site", "PFAL", "obs","AFPL")

FDis_extra_replicate_4 <- left_join(FDis_PFAL_extra_replicate_4, FDis_simul_env, by="Site")
FDis_extra_replicate_4 <- left_join(FDis_extra_replicate_4, FDis_AFPL_extra_replicate_4, by="Site")
colnames(FDis_extra_replicate_4)<-c("Site","PFAL", "obs","AFPL")

FDis_inter_replicate_4 <- left_join(FDis_PFAL_inter_replicate_4, FDis_simul_env, by="Site")
FDis_inter_replicate_4 <- left_join(FDis_inter_replicate_4, FDis_AFPL_inter_replicate_4, by="Site")
colnames(FDis_inter_replicate_4)<-c("Site","PFAL", "obs","AFPL")

CWM_extra_replicate_4 <- left_join(CWM_PFAL_extra_replicate_4, CWM_simul_env, by="Site")
CWM_extra_replicate_4 <- left_join(CWM_extra_replicate_4, CWM_AFPL_extra_replicate_4, by="Site")
colnames(CWM_extra_replicate_4)<-c("Site", "PFAL", "obs","AFPL")

CWM_inter_replicate_4 <- left_join(CWM_PFAL_inter_replicate_4, CWM_simul_env, by="Site")
CWM_inter_replicate_4 <- left_join(CWM_inter_replicate_4, CWM_AFPL_inter_replicate_4, by="Site")
colnames(CWM_inter_replicate_4)<-c("Site", "PFAL", "obs","AFPL")



####### FIGURE S9 

Figure_S9 <- data.frame("Assemble-First Approach"=c(CWM_inter_replicate_1$AFPL,
                                                    1500,521,
                                                    FDis_inter_replicate_1$AFPL,
                                                    98,5,
                                                    CWM_extra_replicate_1$AFPL,
                                                    1424, 607,
                                                    FDis_extra_replicate_1$AFPL,
                                                    121,40,
                                                    CWM_inter_replicate_2$AFPL,
                                                    1500,521,
                                                    FDis_inter_replicate_2$AFPL,
                                                    98,5,
                                                    CWM_extra_replicate_2$AFPL,
                                                    1424, 607,
                                                    FDis_extra_replicate_2$AFPL,
                                                    121,40,
                                                    CWM_inter_replicate_3$AFPL,
                                                    1500,521,
                                                    FDis_inter_replicate_3$AFPL,
                                                    98,5,
                                                    CWM_extra_replicate_3$AFPL,
                                                    1424, 607,
                                                    FDis_extra_replicate_3$AFPL,
                                                    121,40,
                                                    CWM_inter_replicate_4$AFPL,
                                                    1500,521,
                                                    FDis_inter_replicate_4$AFPL,
                                                    98,5,
                                                    CWM_extra_replicate_4$AFPL,
                                                    1424, 607,
                                                    FDis_extra_replicate_4$AFPL,
                                                    121,40),
                        "Predict-First Approach"=c(CWM_inter_replicate_1$PFAL,
                                                   1500,521,
                                                   FDis_inter_replicate_1$PFAL,
                                                   98,5,
                                                   CWM_extra_replicate_1$PFAL,
                                                   1424,607,
                                                   FDis_extra_replicate_1$PFAL,
                                                   121,40,
                                                   CWM_inter_replicate_2$AFPL,
                                                   1500,521,
                                                   FDis_inter_replicate_2$AFPL,
                                                   98,5,
                                                   CWM_extra_replicate_2$AFPL,
                                                   1424, 607,
                                                   FDis_extra_replicate_2$AFPL,
                                                   121,40,
                                                   CWM_inter_replicate_3$AFPL,
                                                   1500,521,
                                                   FDis_inter_replicate_3$AFPL,
                                                   98,5,
                                                   CWM_extra_replicate_3$AFPL,
                                                   1424, 607,
                                                   FDis_extra_replicate_3$AFPL,
                                                   121,40,
                                                   CWM_inter_replicate_4$AFPL,
                                                   1500,521,
                                                   FDis_inter_replicate_4$AFPL,
                                                   98,5,
                                                   CWM_extra_replicate_4$AFPL,
                                                   1424, 607,
                                                   FDis_extra_replicate_4$AFPL,
                                                   121,40),
                        "Indice"=c(rep("CWM", 302), rep("FDis", 302), rep("CWM",dim(CWM_extra_replicate_1)[2]+2), rep("FDis", dim(FDis_extra_replicate_1)[2]+2),
                                   rep("CWM", 302), rep("FDis", 302), rep("CWM",dim(CWM_extra_replicate_2)[2]+2), rep("FDis", dim(FDis_extra_replicate_2)[2]+2),
                                   rep("CWM", 302), rep("FDis", 302), rep("CWM",dim(CWM_extra_replicate_3)[2]+2), rep("FDis", dim(FDis_extra_replicate_3)[2]+2),
                                   rep("CWM", 302), rep("FDis", 302), rep("CWM",dim(CWM_extra_replicate_4)[2]+2), rep("FDis", dim(FDis_extra_replicate_4)[2]+2)),
                        "Crossvalidation"=c(rep("Interpolation",604), rep("Extrapolation",(dim(CWM_extra_replicate_1)[2]+2)*2),
                                            rep("Interpolation",604), rep("Extrapolation",(dim(CWM_extra_replicate_2)[2]+2)*2),
                                            rep("Interpolation",604), rep("Extrapolation",(dim(CWM_extra_replicate_3)[2]+2)*2),
                                            rep("Interpolation",604), rep("Extrapolation",(dim(CWM_extra_replicate_4)[2]+2)*2)),
                        "alpha"=c(rep(c(rep(1,300),c(0,0)),2), rep(c(rep(1,dim(CWM_extra_replicate_1)[2]),c(0,0)),2),
                                  rep(c(rep(1,300),c(0,0)),2), rep(c(rep(1,dim(CWM_extra_replicate_2)[2]),c(0,0)),2),
                                  rep(c(rep(1,300),c(0,0)),2), rep(c(rep(1,dim(CWM_extra_replicate_3)[2]),c(0,0)),2),
                                  rep(c(rep(1,300),c(0,0)),2), rep(c(rep(1,dim(CWM_extra_replicate_4)[2]),c(0,0)),2)),
                        "set"=c(rep(1,604+(dim(CWM_extra_replicate_1)[2]+2)*2), rep(2,604+(dim(CWM_extra_replicate_2)[2]+2)*2), rep(3,604+(dim(CWM_extra_replicate_3)[2]+2)*2), rep(4, 604+(dim(CWM_extra_replicate_4)[2]+2)*2))
)

Figure_S9$Crossvalidation <- factor(Figure_S9$Crossvalidation, levels = c("Interpolation",  "Extrapolation"))


CWM_inter_replicate<-rbind(CWM_inter_replicate_1, CWM_inter_replicate_2, CWM_inter_replicate_3, CWM_inter_replicate_4)
CWM_extra_replicate<-rbind(CWM_extra_replicate_1, CWM_extra_replicate_2, CWM_extra_replicate_3, CWM_extra_replicate_4)

FDis_inter_replicate<-rbind(FDis_inter_replicate_1, FDis_inter_replicate_2, FDis_inter_replicate_3, FDis_inter_replicate_4)
FDis_extra_replicate<-rbind(FDis_extra_replicate_1, FDis_extra_replicate_2, FDis_extra_replicate_3, FDis_extra_replicate_4)


cor(CWM_inter_replicate$PFAL, CWM_inter_replicate$AFPL)^2
cor(CWM_extra_replicate$PFAL, CWM_extra_replicate$AFPL)^2
cor(FDis_inter_replicate$PFAL, FDis_inter_replicate$AFPL)^2
cor(FDis_extra_replicate$PFAL, FDis_extra_replicate$AFPL)^2
data_lab <- data.frame(x=c(720, 750, 45,70), 
                       y=c(1500, 1500, 170, 160), 
                       lab=c("R2 > 0.99",
                             "R2 =0.98",
                             "R2 =0.92",
                             "R2 = 0.74"),
                       Indice=c("CWM","CWM","FDis","FDis"),
                       Crossvalidation = factor(rep(c("Interpolation","Extrapolation"),2), levels = c("Interpolation","Extrapolation")))
Figure_S9 %>%
  filter(alpha == 1) %>%
  ggplot() +
  aes(x = Predict.First.Approach, y = Assemble.First.Approach) +
  geom_hex()+
  scale_fill_gradient(low = "grey", high = "black")+
  xlab("Predict-First Approach")+
  ylab("Assemble-First Approach")+
  geom_point(data = Figure_S9[which(Figure_S9$alpha==0),], aes(x = Predict.First.Approach, y = Assemble.First.Approach), alpha=0)+
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    strip.text.x = element_text(size = 18),
    strip.text.y = element_text(size = 18)) + 
  geom_abline(slope=1, intercept = 0)+
  facet_nested(~ Indice + Crossvalidation, scales = "free", independent = "y")+
  geom_label(data=data_lab, mapping=aes(x=x, y=y, label=lab), size=3)


####### FIGURE S15
figure_S15 <- data.frame("RMSE"=c(rmse(CWM_inter_replicate_1$obs, CWM_inter_replicate_1$AFPL), # CWM trait based approach interpolation
                                 rmse(CWM_inter_replicate_2$obs, CWM_inter_replicate_2$AFPL),
                                 rmse(CWM_inter_replicate_3$obs, CWM_inter_replicate_3$AFPL),
                                 rmse(CWM_inter_replicate_4$obs, CWM_inter_replicate_4$AFPL),
                                 rmse(CWM_extra_replicate_1$obs, CWM_extra_replicate_1$AFPL),
                                 rmse(CWM_extra_replicate_2$obs, CWM_extra_replicate_2$AFPL),
                                 rmse(CWM_extra_replicate_3$obs, CWM_extra_replicate_3$AFPL),
                                 rmse(CWM_extra_replicate_4$obs, CWM_extra_replicate_4$AFPL),
                                 rmse(FDis_inter_replicate_1$obs, FDis_inter_replicate_1$AFPL),
                                 rmse(FDis_inter_replicate_2$obs, FDis_inter_replicate_2$AFPL),
                                 rmse(FDis_inter_replicate_3$obs, FDis_inter_replicate_3$AFPL),
                                 rmse(FDis_inter_replicate_4$obs, FDis_inter_replicate_4$AFPL),
                                 rmse(FDis_extra_replicate_1$obs, FDis_extra_replicate_1$AFPL),
                                 rmse(FDis_extra_replicate_2$obs, FDis_extra_replicate_2$AFPL),
                                 rmse(FDis_extra_replicate_3$obs, FDis_extra_replicate_3$AFPL),
                                 rmse(FDis_extra_replicate_4$obs, FDis_extra_replicate_4$AFPL),
                                 rmse(CWM_inter_replicate_1$obs, CWM_inter_replicate_1$PFAL),
                                 rmse(CWM_inter_replicate_2$obs, CWM_inter_replicate_2$PFAL),
                                 rmse(CWM_inter_replicate_3$obs, CWM_inter_replicate_3$PFAL),
                                 rmse(CWM_inter_replicate_4$obs, CWM_inter_replicate_4$PFAL),
                                 rmse(CWM_extra_replicate_1$obs, CWM_extra_replicate_1$PFAL),
                                 rmse(CWM_extra_replicate_2$obs, CWM_extra_replicate_2$PFAL),
                                 rmse(CWM_extra_replicate_3$obs, CWM_extra_replicate_3$PFAL),
                                 rmse(CWM_extra_replicate_4$obs, CWM_extra_replicate_4$PFAL),
                                 rmse(FDis_inter_replicate_1$obs, FDis_inter_replicate_1$PFAL),
                                 rmse(FDis_inter_replicate_2$obs, FDis_inter_replicate_2$PFAL),
                                 rmse(FDis_inter_replicate_3$obs, FDis_inter_replicate_3$PFAL),
                                 rmse(FDis_inter_replicate_4$obs, FDis_inter_replicate_4$PFAL),
                                 rmse(FDis_extra_replicate_1$obs, FDis_extra_replicate_1$PFAL),
                                 rmse(FDis_extra_replicate_2$obs, FDis_extra_replicate_2$PFAL),
                                 rmse(FDis_extra_replicate_3$obs, FDis_extra_replicate_3$PFAL),
                                 rmse(FDis_extra_replicate_4$obs, FDis_extra_replicate_4$PFAL)),
                        "R2"=c(cor(CWM_inter_replicate_1$obs, CWM_inter_replicate_1$AFPL)^2, # CWM trait based approach interpolation
                               cor(CWM_inter_replicate_2$obs, CWM_inter_replicate_2$AFPL)^2,
                               cor(CWM_inter_replicate_3$obs, CWM_inter_replicate_3$AFPL)^2,
                               cor(CWM_inter_replicate_4$obs, CWM_inter_replicate_4$AFPL)^2,
                               cor(CWM_extra_replicate_1$obs, CWM_extra_replicate_1$AFPL)^2,
                               cor(CWM_extra_replicate_2$obs, CWM_extra_replicate_2$AFPL)^2,
                               cor(CWM_extra_replicate_3$obs, CWM_extra_replicate_3$AFPL)^2,
                               cor(CWM_extra_replicate_4$obs, CWM_extra_replicate_4$AFPL)^2,
                               cor(FDis_inter_replicate_1$obs, FDis_inter_replicate_1$AFPL)^2,
                               cor(FDis_inter_replicate_2$obs, FDis_inter_replicate_2$AFPL)^2,
                               cor(FDis_inter_replicate_3$obs, FDis_inter_replicate_3$AFPL)^2,
                               cor(FDis_inter_replicate_4$obs, FDis_inter_replicate_4$AFPL)^2,
                               cor(FDis_extra_replicate_1$obs, FDis_extra_replicate_1$AFPL)^2,
                               cor(FDis_extra_replicate_2$obs, FDis_extra_replicate_2$AFPL)^2,
                               cor(FDis_extra_replicate_3$obs, FDis_extra_replicate_3$AFPL)^2,
                               cor(FDis_extra_replicate_4$obs, FDis_extra_replicate_4$AFPL)^2,
                               cor(CWM_inter_replicate_1$obs, CWM_inter_replicate_1$PFAL)^2,
                               cor(CWM_inter_replicate_2$obs, CWM_inter_replicate_2$PFAL)^2,
                               cor(CWM_inter_replicate_3$obs, CWM_inter_replicate_3$PFAL)^2,
                               cor(CWM_inter_replicate_4$obs, CWM_inter_replicate_4$PFAL)^2,
                               cor(CWM_extra_replicate_1$obs, CWM_extra_replicate_1$PFAL)^2,
                               cor(CWM_extra_replicate_2$obs, CWM_extra_replicate_2$PFAL)^2,
                               cor(CWM_extra_replicate_3$obs, CWM_extra_replicate_3$PFAL)^2,
                               cor(CWM_extra_replicate_4$obs, CWM_extra_replicate_4$PFAL)^2,
                               cor(FDis_inter_replicate_1$obs, FDis_inter_replicate_1$PFAL)^2,
                               cor(FDis_inter_replicate_2$obs, FDis_inter_replicate_2$PFAL)^2,
                               cor(FDis_inter_replicate_3$obs, FDis_inter_replicate_3$PFAL)^2,
                               cor(FDis_inter_replicate_4$obs, FDis_inter_replicate_4$PFAL)^2,
                               cor(FDis_extra_replicate_1$obs, FDis_extra_replicate_1$PFAL)^2,
                               cor(FDis_extra_replicate_2$obs, FDis_extra_replicate_2$PFAL)^2,
                               cor(FDis_extra_replicate_3$obs, FDis_extra_replicate_3$PFAL)^2,
                               cor(FDis_extra_replicate_4$obs, FDis_extra_replicate_4$PFAL)^2),
                        "Indices"=c(rep(c(rep("CWM", 8), rep("FDis", 8)),2)),
                        "Approaches"=c(rep("Assemble-First Approach", 4*4), rep("Predict-First Approach", 4*4)),
                        "Crossvalidation"=c(rep(c(rep("Interpolation", 4), rep("Extrapolation", 4)),4)))

figure_S15 <- data.frame("value"=c(figure_S15$R2, figure_S15$RMSE), 
                        "Indices"=rep(figure_S15$Indices,2),
                        "Approaches"=rep(figure_S15$Approaches,2),
                        "Crossvalidation"=rep(figure_S15$Crossvalidation,2),
                        "Metrics"=c(rep("R2", length(figure_S15$R2)), rep("RMSE", length(figure_S15$RMSE))))

figure_S15$Crossvalidation <- fct_relevel(figure_S15$Crossvalidation, c("Interpolation", "Extrapolation"))
figure_S15$Approaches <- fct_relevel(figure_S15$Approaches, c("Predict-First Approach","Assemble-First Approach"))


figure_S15 %>%
  ggplot() +
  aes(x=Crossvalidation, y = value, fill = Approaches) +
  scale_fill_manual(values=c("grey", "white"))+
  geom_boxplot() +
  theme_bw() +
  facet_nested( Metrics ~ Indices, scales = "free")+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 10))

###### FIGURE S18
order_species_inter_simul <- species_models_inter_tot
for (i in 1 : 1200){ 
  order_species_inter_simul[i,] <- order(-species_models_inter_tot[i,])
  print(i)
}

order_species_extra_simul<- species_models_extra_tot
for (i in 1 : dim(species_models_extra_tot)[1]){ 
  order_species_extra_simul[i,] <- order(-species_models_extra_tot[i,])
}

# calculation of functional indices and their correlation to those observed for each species addition
plant_inter_figS18 <- data.frame(matrix(0, ncol=dim(species_models_inter_tot)[2], nrow=dim(species_models_inter_tot)[1]), row.names = rownames(species_models_inter_tot))
colnames(plant_inter_figS18)=colnames(species_models_inter_tot)
R2_CWM_inter_figS18 <- data.frame(matrix(0, nrow=dim(species_models_inter_tot)[2], ncol=1))
rownames(R2_CWM_inter_figS18) <- 1:dim(species_models_inter_tot)[2]
colnames(R2_CWM_inter_figS18)<- c("Trait")
R2_FDis_inter_figS18 <- data.frame(matrix(0, nrow=dim(species_models_inter_tot)[2], ncol=1))
rownames(R2_FDis_inter_figS18)<- 1:dim(species_models_inter_tot)[2]
colnames(R2_FDis_inter_figS18)<- c("Trait")
for (m in 2:dim(species_models_inter_tot)[2]){ ## One new species is added each time
  for (i in 1 : dim(species_models_inter_tot)[1]){ ## we add the new species in the botanical survey for each site
    plant_inter_figS18[i,as.numeric(order_species_inter_simul[i,m])] = species_models_inter_tot[i,as.numeric(order_species_inter_simul[i,m])]
  }
  CWM_inter_figS18 <- CWM_simul(plant_inter_figS18)
  CWM_inter_figS18 <- left_join(CWM_inter_figS18, CWM_simul_env, by="Site")
  R2_CWM_inter_figS18[m,] <- c(cor(CWM_inter_figS18$Trait.x, CWM_inter_figS18$Trait.y)^2)
  
  FDis_inter_figS18 <- FDis(plant_inter_figS18)
  FDis_inter_figS18 <- left_join(FDis_inter_figS18, FDis_simul_env, by="Site")
  R2_FDis_inter_figS18[m,] <- c(cor(FDis_inter_figS18$Trait.x, FDis_inter_figS18$Trait.y)^2)
  print(m)
}
print("inter recovery OK")

# calculation of functional indices and their correlation to those observed for each species addition
plant_extra_figS18 <- data.frame(matrix(0, ncol=dim(species_models_extra_tot)[2], nrow=dim(species_models_extra_tot)[1]), row.names = rownames(species_models_extra_tot))
colnames(plant_extra_figS18)=colnames(species_models_extra_tot)
R2_CWM_extra_figS18 <- data.frame(matrix(0, nrow=dim(species_models_extra_tot)[2], ncol=1))
rownames(R2_CWM_extra_figS18)=1:dim(species_models_extra_tot)[2]
colnames(R2_CWM_extra_figS18)<- c("Trait")
R2_FDis_extra_figS18 <- data.frame(matrix(0, nrow=dim(species_models_extra_tot)[2], ncol=1))
rownames(R2_FDis_extra_figS18)=1:dim(species_models_extra_tot)[2]
colnames(R2_FDis_extra_figS18)<- c("Trait")
for (m in 2:dim(species_models_extra_tot)[2]){ ## One new species is added each time
  for (i in 1 : dim(species_models_extra_tot)[1]){ ## we add the new species in the botanical survey for each site
    plant_extra_figS18[i,as.numeric(order_species_extra_simul[i,m])] = species_models_extra_tot[i,as.numeric(order_species_extra_simul[i,m])]
  }
  CWM_extra_figS18 <- CWM(plant_extra_figS18)
  CWM_extra_figS18 <- left_join(CWM_extra_figS18, CWM_simul_env, by="Site")
  R2_CWM_extra_figS18[m,] <- c(cor(CWM_extra_figS18$Trait.x, CWM_extra_figS18$Trait.y)^2)
  
  FDis_extra_figS18 <- FDis(plant_extra_figS18)
  FDis_extra_figS18 <- left_join(FDis_extra_figS18, FDis_simul_env, by="Site")
  R2_FDis_extra_figS18[m,] <- c(cor(FDis_extra_figS18$Trait.x, FDis_extra_figS18$Trait.y)^2)
  print(m)
}
print("extra recovery OK")


figure_S18 <- data.frame("Number of species" =rep(1:dim(species_models_inter_tot)[2],4),
                         "R2"=c(R2_CWM_inter_figS18$Trait,
                                R2_CWM_extra_figS18$Trait,
                                R2_FDis_inter_figS18$Trait,
                                R2_FDis_extra_figS18$Trait),
                         "Indices"=c(rep("CWM", dim(species_models_inter_tot)[2]*2), rep("FDis", dim(species_models_inter_tot)[2]*2)), 
                         "Crossvalidation"=rep(c(rep("Interpolation", dim(species_models_inter_tot)[2]), rep("Extrapolation", dim(species_models_inter_tot)[2])), 2), 
                         "Trait_based models"=c(rep(cor(CWM_inter_replicate$obs, CWM_inter_replicate$AFPL)^2, dim(species_models_inter_tot)[2]),
                                                rep(cor(CWM_extra_replicate$obs, CWM_extra_replicate$AFPL)^2, dim(species_models_inter_tot)[2]),
                                                rep(cor(FDis_inter_replicate$obs, FDis_inter_replicate$AFPL)^2, dim(species_models_inter_tot)[2]),
                                                rep(cor(FDis_extra_replicate$obs, FDis_extra_replicate$AFPL)^2, dim(species_models_inter_tot)[2]))
)

figure_S18$Crossvalidation <- fct_relevel(figure_S18$Crossvalidation, c("Interpolation", "Extrapolation"))


figure_S18 %>%
  ggplot() +
  aes(x = Number.of.species, y = R2)+
  scale_x_continuous(trans='log2')+
  geom_point(shape = "circle", size = 1, col="black")+
  geom_hline(aes(yintercept = Trait_based.models),col="darkgrey", linetype = "dashed", size=2)+
  theme_bw() +  
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 10))+
  xlab("Number of species")+
  facet_nested(~ Indices + Crossvalidation)
## END OF THE SCRIPT
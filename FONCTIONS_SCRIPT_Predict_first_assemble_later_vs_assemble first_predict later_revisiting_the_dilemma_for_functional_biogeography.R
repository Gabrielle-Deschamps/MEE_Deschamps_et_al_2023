#This script hold some useful functions to compute all path analyses from :
#Deschamps G.; Giovanni P.; Brun P.; Thuiller, W. (2022) 
#Predict first - assemble later vs assemble first - predict later: revisiting the dilemma for functional biogeography.

#### Packages
library(ggplot2)
library(tidyverse)
library(readxl)
library(Metrics)
library(ggh4x)
library(randomForestSRC)
library(BIOMOD)


#### -- Data preparation -- ####
# This function returns an array of data that will be used to calculate CWM and FDis.
data_table <- function(data){
  data_table <- data.frame("Site"=rownames(data), data)
  colnames(data_table)<-c("Site", colnames(data))
  data_table <- data_table %>% pivot_longer(-Site)
  colnames(data_table)<- c("Site", "Taxa", "Abundance")
  data_table$Abundance <- as.numeric(data_table$Abundance)
  data_table <- left_join(data_table, traits, by="Taxa")
  data_table <- cbind(data_table, "A_LNC"=data_table$Abundance*data_table$LNC,
                   "A_SLA"=data_table$Abundance*data_table$SLA,
                   "A_PLH"=data_table$Abundance*data_table$PLH,
                   "pa_LNC"=data_table$Abundance*(data_table$LNC/data_table$LNC),
                   "pa_SLA"=data_table$Abundance*(data_table$SLA/data_table$SLA),
                   "pa_PLH"=data_table$Abundance*(data_table$PLH/data_table$PLH))
  data_table <- data_table %>% group_by(Site)
  return(data_table)
}
#### -- Calcul of the CWM -- ####
# This function returns an array of data containing the CWM for LNC, SLA and Plant Height.

CWM <- function(data){
  data <- data_table(data)
  LNC = data %>%
    summarise("LNC"=c(sum(A_LNC, na.rm = TRUE)/sum(pa_LNC, na.rm = TRUE)))
  SLA = data %>% 
    summarise("SLA"=c(sum(A_SLA, na.rm = TRUE)/sum(pa_SLA, na.rm = TRUE)))
  PLH = data %>%
    summarise("PLH"=c(sum(A_PLH, na.rm = TRUE)/sum(pa_PLH, na.rm = TRUE)))
  CWM = data.frame("Site"=LNC$Site, "LNC"=LNC$LNC, 'SLA'=SLA$SLA, "PLH"=PLH$PLH)
  return(CWM)
}

#### -- Calcul of the FDis -- ####
# This function returns an array of data containing the FDis for LNC, SLA and Plant Height.

FDis <-  function(data){
  data <- data_table(data)
  LNC = data %>%
    summarise("LNC"=c(sum(pa_LNC*abs(LNC-(sum(A_LNC,na.rm = TRUE)/sum(pa_LNC,na.rm = TRUE))),na.rm = TRUE)/sum(pa_LNC,na.rm = TRUE)))
  SLA = data %>% 
    summarise("SLA"=c(sum(pa_SLA*abs(SLA-(sum(A_SLA,na.rm = TRUE)/sum(pa_SLA,na.rm = TRUE))),na.rm = TRUE)/sum(pa_SLA,na.rm = TRUE)))
  PLH = data %>%
    summarise("PLH"=c(sum(pa_PLH*abs(PLH-(sum(A_PLH,na.rm = TRUE)/sum(pa_PLH,na.rm = TRUE))),na.rm = TRUE)/sum(pa_PLH,na.rm = TRUE)))
  FDis = data.frame("Site"=LNC$Site, "LNC"=LNC$LNC, 'SLA'=SLA$SLA, "PLH"=PLH$PLH)
  return(FDis)
}


#### -- Shuffle traits matrix -- ####
# This function returns an array of data that will be used to calculate CWM and FDis using the shuffle traits matrix.
data_table_shuffle <- function(data){
  data_table <- data.frame("Site"=rownames(data), data)
  data_table <- data_table %>% pivot_longer(-Site)
  colnames(data_table)<- c("Site", "Taxa", "Abundance")
  data_table$Taxa <- str_replace(data_table$Taxa, "\\.", " ")
  data_table <- left_join(data_table, shuffle_traits, by="Taxa")
  data_table <- cbind(data_table, 
                      "A_LNC"=data_table$Abundance*data_table$LNC,
                      "A_SLA"=data_table$Abundance*data_table$SLA,
                      "A_PLH"=data_table$Abundance*data_table$PLH,
                      "pa_LNC"=data_table$Abundance*(data_table$LNC/data_table$LNC),
                      "pa_SLA"=data_table$Abundance*(data_table$SLA/data_table$SLA),
                      "pa_PLH"=data_table$Abundance*(data_table$PLH/data_table$PLH))
  data_table <- data_table %>% group_by(Site)
  return(data_table)
}

# This function returns an array of data containing the CWM for LNC, SLA and Plant Height using the shuffle trait matrix.
CWM_shuffle <- function(data){
  data <- data_table_shuffle(data)
  LNC = data %>%
    summarise("LNC"=c(sum(A_LNC, na.rm = TRUE)/sum(pa_LNC, na.rm = TRUE)))
  SLA = data %>% 
    summarise("SLA"=c(sum(A_SLA, na.rm = TRUE)/sum(pa_SLA, na.rm = TRUE)))
  PLH = data %>%
    summarise("PLH"=c(sum(A_PLH, na.rm = TRUE)/sum(pa_PLH, na.rm = TRUE)))
  CWM = data.frame("Site"=LNC$Site, "LNC"=LNC$LNC, 'SLA'=SLA$SLA, "PLH"=PLH$PLH)
  return(CWM)
}

# This function returns an array of data containing the FDis for LNC, SLA and Plant Height using the shuffle traits matrix.
FDis_shuffle <- function(data){
  data <- data_table_shuffle(data)
  LNC = data %>%
    summarise("LNC"=c(sum(pa_LNC*abs(LNC-(sum(A_LNC,na.rm = TRUE)/sum(pa_LNC,na.rm = TRUE))),na.rm = TRUE)/sum(pa_LNC,na.rm = TRUE)))
  SLA = data %>% 
    summarise("SLA"=c(sum(pa_SLA*abs(SLA-(sum(A_SLA,na.rm = TRUE)/sum(pa_SLA,na.rm = TRUE))),na.rm = TRUE)/sum(pa_SLA,na.rm = TRUE)))
  PLH = data %>%
    summarise("PLH"=c(sum(pa_PLH*abs(PLH-(sum(A_PLH,na.rm = TRUE)/sum(pa_PLH,na.rm = TRUE))),na.rm = TRUE)/sum(pa_PLH,na.rm = TRUE)))
  FDis = data.frame("Site"=LNC$Site, "LNC"=LNC$LNC, 'SLA'=SLA$SLA, "PLH"=PLH$PLH)
  return(FDis)
}


#### -- Simulated trait functional indices -- ####
# This function returns an array of data that will be used to calculate CWM and FDis using the simulated trait.
data_table_simul <- function(data){
  data_table <- data.frame("Site"=rownames(data), data)
  data_table <- data_table %>% pivot_longer(-Site)
  colnames(data_table)<- c("Site", "Taxa", "Abundance")
  data_table <- left_join(data_table, traits_simul, by="Taxa")
  data_table <- cbind(data_table, "A_trait"=data_table$Abundance*data_table$Trait,
                      "pa_trait"=data_table$Abundance*(data_table$Trait/data_table$Trait))
  data_table <- data_table %>% group_by(Site)
  return(data_table)
}

# This function returns an array of data containing the CWM for the simulated trait.
CWM_simul <- function(data){
  data <- data_table_simul(data)
  Trait = data %>%
    summarise("Trait"=c(sum(A_trait, na.rm = TRUE)/sum(pa_trait, na.rm = TRUE)))
  CWM = data.frame("Site"=Trait$Site, "Trait"=Trait$Trait)
  return(CWM)
}

# This function returns an array of data containing the FDis for the simulated trait.
FDis_simul <- function(data){
  data <- data_table_simul(data)
  Trait = data %>%
    summarise("Trait"=c(sum(pa_trait*abs(Trait-(sum(A_trait,na.rm = TRUE)/sum(pa_trait,na.rm = TRUE))),na.rm = TRUE)/sum(pa_trait,na.rm = TRUE)))
  FDis = data.frame("Site"=Trait$Site, "Trait"=Trait$trait)
  return(FDis)
}
## END OF THE SCRIPT


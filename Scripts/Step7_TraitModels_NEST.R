##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Step 7: Phylogenetic Trait Models - NEST TRAITS
# Authors: Sarah L. Jennings, Emma M. Garrison
##########################################################
# The objective of this script is to run phylogenetic models using nest traits
# 4 life history traits: nest strategy, nest site, and nest safety


# load required packages 
library(tidyverse)
library(here)
library(nlme)
library(phytools)
library(ape)
library(geiger)
library(easystats)
library(phylolm)
library(logistf)

###################### Prep

# load in "Coastal_Species_Nest.rds"
C_Nest_dat <- readRDS(here("Outputs", "Coastal_Species_Nest.rds"))
str(C_Nest_dat)
nrow(C_Nest_dat)

C_Nest_dat2 <- C_Nest_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_Nest_dat2)


C_Nest_dat2$Urban <- ifelse(C_Nest_dat2$Urban == "U", 1, 0)
#View(C_Nest_dat2)
colnames(C_Nest_dat2)


######################## UAI and Nest Strategy ##########################
# 0 = enclosed
# 1 = open

# create a new data frame that contains only species with both UAI and nest strategy
UAI_NestStr <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(NestStr)) %>% as.data.frame()
length(UAI_NestStr$NestStr)
# 727 species with UAI and NestStr

###### add and pair tree

# add rownames to data
row.names(UAI_NestStr) <- UAI_NestStr$Species_Jetz

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_NestStr_phydat <- geiger::treedata(tree_out, UAI_NestStr, sort=T)

UAI_NestStr_phy <- UAI_NestStr_phydat$phy
UAI_NestStr_dat <- as.data.frame(UAI_NestStr_phydat$data)

str(UAI_NestStr_dat)
length(UAI_NestStr_dat$NestStr)

# convert traits of interest to numeric
UAI_NestStr_dat$aveUAI <- as.numeric(UAI_NestStr_dat$aveUAI)
UAI_NestStr_dat$Mass_log <- as.numeric(UAI_NestStr_dat$Mass_log)
UAI_NestStr_dat$NestStr <- as.numeric(UAI_NestStr_dat$NestStr)

# run phylogenetic linear model
UAI_GLS_neststr <- gls(aveUAI~ NestStr + Mass_log, data = UAI_NestStr_dat, 
                  correlation = corPagel(0.5, phy= UAI_NestStr_phy, fixed=F, form = ~Species_Jetz), 
                  method = "ML") 

# model summary and results
summary(UAI_GLS_neststr) 
confint(UAI_GLS_neststr)

# model diagnostics
check_model(UAI_GLS_neststr) 
qqnorm(resid(UAI_GLS_neststr)) 
qqline(resid(UAI_GLS_neststr)) 
hist(resid(UAI_GLS_neststr)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_neststr, here("Models/UAI", "UAI_GLS_neststr.rds"))


######################## MUTI and Nest Strategy ##########################
# 0 = enclosed
# 1 = open

# create a new data frame that contains only species with both MUTI and nest strategy
MUTI_NestStr <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(NestStr)) %>% as.data.frame()
length(MUTI_NestStr$NestStr)
# 115 species with MUTIscore and NestStr

###### add and pair tree

# add rownames to data
row.names(MUTI_NestStr) <- MUTI_NestStr$Species_Jetz

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_NestStr_phydat <- geiger::treedata(tree_out, MUTI_NestStr, sort=T)

MUTI_NestStr_phy <- MUTI_NestStr_phydat$phy
MUTI_NestStr_dat <- as.data.frame(MUTI_NestStr_phydat$data)

str(MUTI_NestStr_dat)
length(MUTI_NestStr_dat$NestStr)

# convert traits of interest to numeric
MUTI_NestStr_dat$MUTIscore <- as.numeric(MUTI_NestStr_dat$MUTIscore)
MUTI_NestStr_dat$Mass_log <- as.numeric(MUTI_NestStr_dat$Mass_log)
MUTI_NestStr_dat$NestStr <- as.numeric(MUTI_NestStr_dat$NestStr)


# run phylogenetic linear model
MUTI_GLS_neststr <- gls(MUTIscore~ NestStr + Mass_log, data = MUTI_NestStr_dat, 
                       correlation = corPagel(0.5, phy= MUTI_NestStr_phy, fixed = F, form = ~Species_Jetz), 
                       method = "ML") 

# model summary and results
summary(MUTI_GLS_neststr) 
confint(MUTI_GLS_neststr)

# model diagnostic
check_model(MUTI_GLS_neststr) 
qqnorm(resid(MUTI_GLS_neststr)) 
qqline(resid(MUTI_GLS_neststr)) 
hist(resid(MUTI_GLS_neststr)) 


# save model for easy retrieval 
saveRDS(MUTI_GLS_neststr, here("Models/MUTI", "MUTI_GLS_neststr.rds"))


######################## UN and Nest Strategy ##########################
# 0 = enclosed
# 1 = open

# create a new data frame that contains only species with both UN and nest strategy
UN_NestStr <- C_Nest_dat2 %>% filter(!is.na(Urban)) %>%
  filter(!is.na(NestStr)) %>% column_to_rownames(., var="Species_Jetz")
length(UN_NestStr$NestStr)
# 122 species with UN and NestStr

UN_NestStr %>% filter(!is.na(NestStr)) %>% 
  filter(!is.na(Urban)) %>% 
  group_by(NestStr) %>% count()

# do not run --> only 8 species for UN with "enclosed" = 0 Nest Strategy 
# too small of an N for this category 

##########################################################################
##########################################################################

######################## UAI and Nest Site Low vs High + Flexible ##########################
# 0 = low
# 1 = high and flexible species

UAI_NestHighplusFlex <- C_Nest_dat2 %>% 
  filter(!is.na(aveUAI)) %>% 
  filter(!is.na(NestSite_HighplusFlex)) %>% as.data.frame()
length(UAI_NestHighplusFlex$NestSite_HighplusFlex)
# 786 species 

###### add and pair tree

# add rownames to data
row.names(UAI_NestHighplusFlex) <- UAI_NestHighplusFlex$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_NestHighplusFlex_phydat <- geiger::treedata(tree_out, UAI_NestHighplusFlex, sort=T)

UAI_NestHighplusFlex_phy <- UAI_NestHighplusFlex_phydat$phy
UAI_NestHighplusFlex_dat <- as.data.frame(UAI_NestHighplusFlex_phydat$data)

str(UAI_NestHighplusFlex_dat)
length(UAI_NestHighplusFlex_dat$NestSite_HighplusFlex)

# convert traits of interest to numeric
UAI_NestHighplusFlex_dat$aveUAI <- as.numeric(UAI_NestHighplusFlex_dat$aveUAI)
UAI_NestHighplusFlex_dat$Mass_log <- as.numeric(UAI_NestHighplusFlex_dat$Mass_log)
UAI_NestHighplusFlex_dat$NestSite_HighplusFlex <- as.numeric(UAI_NestHighplusFlex_dat$NestSite_HighplusFlex)

# run a phylogenetic linear model
UAI_GLS_nest_highplusflex <- gls(aveUAI~ NestSite_HighplusFlex + Mass_log, data = UAI_NestHighplusFlex_dat, 
                       correlation = corPagel(0.5, phy = UAI_NestHighplusFlex_phy, fixed = F, form = ~Species_Jetz), 
                       method = "ML") 

# model summary and results
summary(UAI_GLS_nest_highplusflex) 
confint(UAI_GLS_nest_highplusflex)

# model diagnostics
check_model(UAI_GLS_nest_highplusflex) 
qqnorm(resid(UAI_GLS_nest_highplusflex)) 
qqline(resid(UAI_GLS_nest_highplusflex))
hist(resid(UAI_GLS_nest_highplusflex))

# save model
saveRDS(UAI_GLS_nest_highplusflex, here("Models/UAI", "UAI_GLS_nest_highplusflex.rds"))

######################## UAI and Nest Site Low + Flexible vs High ##########################
# 0 = low and flexible species
# 1 = high species

UAI_NestLowplusFlex <- C_Nest_dat2 %>% 
  filter(!is.na(aveUAI)) %>% 
  filter(!is.na(NestSite_LowplusFlex)) %>% as.data.frame()
length(UAI_NestLowplusFlex$NestSite_LowplusFlex)
# 786 species

###### add and pair tree

# add rownames to data
row.names(UAI_NestLowplusFlex) <- UAI_NestLowplusFlex$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_NestLowplusFlex_phydat <- geiger::treedata(tree_out, UAI_NestLowplusFlex, sort=T)

UAI_NestLowplusFlex_phy <- UAI_NestLowplusFlex_phydat$phy
UAI_NestLowplusFlex_dat <- as.data.frame(UAI_NestLowplusFlex_phydat$data)

str(UAI_NestLowplusFlex_dat)
length(UAI_NestLowplusFlex_dat$NestSite_LowplusFlex)

# convert traits of interest to numeric
UAI_NestLowplusFlex_dat$aveUAI <- as.numeric(UAI_NestLowplusFlex_dat$aveUAI)
UAI_NestLowplusFlex_dat$Mass_log <- as.numeric(UAI_NestLowplusFlex_dat$Mass_log)
UAI_NestLowplusFlex_dat$NestSite_LowplusFlex <- as.numeric(UAI_NestLowplusFlex_dat$NestSite_LowplusFlex)

# run a phylogenetic linear model
UAI_GLS_nest_lowplusflex <- gls(aveUAI~ NestSite_LowplusFlex + Mass_log, data = UAI_NestLowplusFlex_dat, 
                                 correlation = corPagel(0.5, phy = UAI_NestLowplusFlex_phy, fixed = F, form = ~Species_Jetz), 
                                 method = "ML") 

# model summary and results
summary(UAI_GLS_nest_lowplusflex) 
confint(UAI_GLS_nest_lowplusflex, level = 0.95)
confint(UAI_GLS_nest_lowplusflex, level = 0.85)

# model diagnostics
check_model(UAI_GLS_nest_lowplusflex) 
qqnorm(resid(UAI_GLS_nest_lowplusflex)) 
qqline(resid(UAI_GLS_nest_lowplusflex))
hist(resid(UAI_GLS_nest_lowplusflex))

# save model
saveRDS(UAI_GLS_nest_lowplusflex, here("Models/UAI", "UAI_GLS_nest_lowplusflex.rds"))

############################### UAI and High vs Low (flexible species excluded) ########################
# Run model with all flexible species excluded that focuses only on species that use either HIGH or LOW nest sites 

UAI_NestLowHigh <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(NestSite_LowHigh)) 
nrow(UAI_NestLowHigh)
# 565 species 

###### add and pair tree

# add rownames to data
row.names(UAI_NestLowHigh) <- UAI_NestLowHigh$Species_Jetz

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_NestLowHigh_phydat <- geiger::treedata(tree_out, UAI_NestLowHigh, sort=T)

UAI_NestLowHigh_phy <- UAI_NestLowHigh_phydat$phy
UAI_NestLowHigh_dat <- as.data.frame(UAI_NestLowHigh_phydat$data)

str(UAI_NestLowHigh_dat)
length(UAI_NestLowHigh_dat$NestSite_LowHigh) # 565

# convert traits of interest to numeric
UAI_NestLowHigh_dat$aveUAI <- as.numeric(UAI_NestLowHigh_dat$aveUAI)
UAI_NestLowHigh_dat$Mass_log <- as.numeric(UAI_NestLowHigh_dat$Mass_log)
UAI_NestLowHigh_dat$NestSite_LowHigh <- as.numeric(UAI_NestLowHigh_dat$NestSite_LowHigh)

# run a phylogenetic linear model
UAI_GLS_nest_lowhigh <- gls(aveUAI~ NestSite_LowHigh + Mass_log, data = UAI_NestLowHigh_dat, 
                          correlation = corPagel(0.5, phy = UAI_NestLowHigh_phy, fixed = F, form = ~Species_Jetz), 
                          method = "ML") 


# model summary and results
summary(UAI_GLS_nest_lowhigh) 
confint(UAI_GLS_nest_lowhigh, level = 0.95)

# model diagnostics
check_model(UAI_GLS_nest_lowhigh) 
qqnorm(resid(UAI_GLS_nest_lowhigh)) 
qqline(resid(UAI_GLS_nest_lowhigh))
hist(resid(UAI_GLS_nest_lowhigh))

# save model
saveRDS(UAI_GLS_nest_lowhigh, here("Models/UAI", "UAI_GLS_nest_lowhigh.rds"))

######################## MUTI and Nest Site Low vs High + Flexible ##########################
# 0 = low
# 1 = high and flexible species

MUTI_NestHighplusFlex <- C_Nest_dat2 %>% 
  filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(NestSite_HighplusFlex)) %>% as.data.frame()
length(MUTI_NestHighplusFlex$NestSite_HighplusFlex)
# 128 species 

###### add and pair tree

# add rownames to data
row.names(MUTI_NestHighplusFlex) <- MUTI_NestHighplusFlex$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_NestHighplusFlex_phydat <- geiger::treedata(tree_out, MUTI_NestHighplusFlex, sort=T)

MUTI_NestHighplusFlex_phy <- MUTI_NestHighplusFlex_phydat$phy
MUTI_NestHighplusFlex_dat <- as.data.frame(MUTI_NestHighplusFlex_phydat$data)

str(MUTI_NestHighplusFlex_dat)
length(MUTI_NestHighplusFlex_dat$NestSite_HighplusFlex)

# convert traits of interest to numeric
MUTI_NestHighplusFlex_dat$MUTIscore <- as.numeric(MUTI_NestHighplusFlex_dat$MUTIscore)
MUTI_NestHighplusFlex_dat$Mass_log <- as.numeric(MUTI_NestHighplusFlex_dat$Mass_log)
MUTI_NestHighplusFlex_dat$NestSite_HighplusFlex <- as.numeric(MUTI_NestHighplusFlex_dat$NestSite_HighplusFlex)

# run a phylogenetic linear model
MUTI_GLS_nest_highplusflex <- gls(MUTIscore ~ NestSite_HighplusFlex + Mass_log, data = MUTI_NestHighplusFlex_dat, 
                                 correlation = corPagel(0.5, phy = MUTI_NestHighplusFlex_phy, fixed = F, form = ~Species_Jetz), 
                                 method = "ML") 

# model summary and results
summary(MUTI_GLS_nest_highplusflex) 
confint(MUTI_GLS_nest_highplusflex)


# model diagnostics
check_model(MUTI_GLS_nest_highplusflex) 
qqnorm(resid(MUTI_GLS_nest_highplusflex)) 
qqline(resid(MUTI_GLS_nest_highplusflex))
hist(resid(MUTI_GLS_nest_highplusflex))

# save model
saveRDS(MUTI_GLS_nest_highplusflex, here("Models/MUTI", "MUTI_GLS_nest_highplusflex.rds"))

######################## MUTI and Nest Site Low + Flexible vs High ##########################
# 0 = low and flexible species
# 1 = high species

MUTI_NestLowplusFlex <- C_Nest_dat2 %>% 
  filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(NestSite_LowplusFlex)) %>% as.data.frame()
length(MUTI_NestLowplusFlex$NestSite_LowplusFlex)
# 128 species

###### add and pair tree

# add rownames to data
row.names(MUTI_NestLowplusFlex) <- MUTI_NestLowplusFlex$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_NestLowplusFlex_phydat <- geiger::treedata(tree_out, MUTI_NestLowplusFlex, sort=T)

MUTI_NestLowplusFlex_phy <- MUTI_NestLowplusFlex_phydat$phy
MUTI_NestLowplusFlex_dat <- as.data.frame(MUTI_NestLowplusFlex_phydat$data)

str(MUTI_NestLowplusFlex_dat)
length(MUTI_NestLowplusFlex_dat$NestSite_LowplusFlex)

# convert traits of interest to numeric
MUTI_NestLowplusFlex_dat$MUTIscore <- as.numeric(MUTI_NestLowplusFlex_dat$MUTIscore)
MUTI_NestLowplusFlex_dat$Mass_log <- as.numeric(MUTI_NestLowplusFlex_dat$Mass_log)
MUTI_NestLowplusFlex_dat$NestSite_LowplusFlex <- as.numeric(MUTI_NestLowplusFlex_dat$NestSite_LowplusFlex)

# run a phylogenetic linear model
MUTI_GLS_nest_lowplusflex <- gls(MUTIscore~ NestSite_LowplusFlex + Mass_log, data = MUTI_NestLowplusFlex_dat, 
                                correlation = corPagel(0.5, phy = MUTI_NestLowplusFlex_phy, fixed = F, form = ~Species_Jetz), 
                                method = "ML") 

# model summary and results
summary(MUTI_GLS_nest_lowplusflex) 
confint(MUTI_GLS_nest_lowplusflex, level = 0.95)
confint(MUTI_GLS_nest_lowplusflex, level = 0.85)

# model diagnostics
check_model(MUTI_GLS_nest_lowplusflex) 
qqnorm(resid(MUTI_GLS_nest_lowplusflex)) 
qqline(resid(MUTI_GLS_nest_lowplusflex))
hist(resid(MUTI_GLS_nest_lowplusflex))

# save model
saveRDS(MUTI_GLS_nest_lowplusflex, here("Models/MUTI", "MUTI_GLS_nest_lowplusflex.rds"))

############################### MUTI and High vs Low (flexible species excluded) ########################
# Run model with all flexible species excluded that focuses only on species that use either HIGH or LOW nest sites 

MUTI_NestLowHigh <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(NestSite_LowHigh)) 
nrow(MUTI_NestLowHigh)
# 89 species 

###### add and pair tree

# add rownames to data
row.names(MUTI_NestLowHigh) <- MUTI_NestLowHigh$Species_Jetz

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_NestLowHigh_phydat <- geiger::treedata(tree_out, MUTI_NestLowHigh, sort=T)

MUTI_NestLowHigh_phy <- MUTI_NestLowHigh_phydat$phy
MUTI_NestLowHigh_dat <- as.data.frame(MUTI_NestLowHigh_phydat$data)

str(MUTI_NestLowHigh_dat)
length(MUTI_NestLowHigh_dat$NestSite_LowHigh)

# convert traits of interest to numeric
MUTI_NestLowHigh_dat$MUTIscore <- as.numeric(MUTI_NestLowHigh_dat$MUTIscore)
MUTI_NestLowHigh_dat$Mass_log <- as.numeric(MUTI_NestLowHigh_dat$Mass_log)
MUTI_NestLowHigh_dat$NestSite_LowHigh <- as.numeric(MUTI_NestLowHigh_dat$NestSite_LowHigh)

# Model is not working with corPagel starting point = 0.5 and fixed = F
# we need to find a value of lambda value to fix in the model 

# Create an empty vector to store AIC values
AIC_values <- numeric()

# Loop through different values of the parameter for corPagel
for (i in seq(0, 1, by = 0.1)) {
  # Fit the gls model with the current value of i
  model <- gls(MUTIscore ~ NestSite_LowHigh + Mass_log, 
               data = MUTI_NestLowHigh_dat, 
               correlation = corPagel(i, phy = MUTI_NestLowHigh_phy, fixed = TRUE, form = ~Species_Jetz), 
               method = "ML")
  # Extract AIC value and store it in the vector
  AIC_values <- c(AIC_values, AIC(model))
}

# Print AIC values
print(AIC_values)
# 0.3 = best AIC score 

# run a phylogenetic linear model with corPagel fixed at 0.3
MUTI_GLS_nest_lowhigh <- gls(MUTIscore ~ NestSite_LowHigh + Mass_log, data = MUTI_NestLowHigh_dat, 
                            correlation = corPagel(0.3, phy = MUTI_NestLowHigh_phy, fixed = T, form = ~Species_Jetz), 
                            method = "ML") 

# model summary and results
summary(MUTI_GLS_nest_lowhigh) 
confint(MUTI_GLS_nest_lowhigh, level = 0.95)
confint(MUTI_GLS_nest_lowhigh, level = 0.85)

# model diagnostics
check_model(MUTI_GLS_nest_lowhigh) 
qqnorm(resid(MUTI_GLS_nest_lowhigh)) 
qqline(resid(MUTI_GLS_nest_lowhigh))
hist(resid(MUTI_GLS_nest_lowhigh))

# save model
saveRDS(MUTI_GLS_nest_lowhigh, here("Models/MUTI", "MUTI_GLS_nest_lowhigh.rds"))


######################## UN and Nest Site Low vs High plus Flexible ##########################
# 0 = low 
# 1 = high plus flexible

UN_NestHighplusFlex <- C_Nest_dat2 %>% 
  filter(!is.na(Urban)) %>%
  filter(!is.na(NestSite_HighplusFlex)) %>% column_to_rownames(., var = "Species_Jetz")
length(UN_NestHighplusFlex$NestSite_HighplusFlex)
# 128 species

###### add and pair tree

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_NestHighplusFlex_phydat <- geiger::treedata(tree_out, UN_NestHighplusFlex, sort=T)

UN_NestHighplusFlex_phy <- UN_NestHighplusFlex_phydat$phy
UN_NestHighplusFlex_dat <- as.data.frame(UN_NestHighplusFlex_phydat$data)

str(UN_NestHighplusFlex_dat)
length(UN_NestHighplusFlex_dat$NestSite_HighplusFlex)

# convert traits of interest to numeric
UN_NestHighplusFlex_dat$Urban <- as.numeric(UN_NestHighplusFlex_dat$Urban)
UN_NestHighplusFlex_dat$Mass_log <- as.numeric(UN_NestHighplusFlex_dat$Mass_log)
UN_NestHighplusFlex_dat$NestSite_Low <- as.numeric(UN_NestHighplusFlex_dat$NestSite_HighplusFlex)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
phyglm_UN_nest_highplusflex_scale <- phyloglm(Urban ~ NestSite_HighplusFlex + scale(Mass_log), 
                                      data = UN_NestHighplusFlex_dat, 
                                      phy = UN_NestHighplusFlex_phy,
                                      boot = 1000) 

# fails to converge

# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ NestSite_HighplusFlex + scale(Mass_log), 
                 data = UN_NestHighplusFlex_dat, 
                 phy = UN_NestHighplusFlex_phy,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylogenetic signal)


# try fixing alpha near upper bounds
set.seed(499)
phyglm_UN_nest_highplusflex_fix <- phyloglm(Urban ~ NestSite_HighplusFlex + scale(Mass_log), 
                                        data = UN_NestHighplusFlex_dat, 
                                        phy = UN_NestHighplusFlex_phy,
                                start.alpha = 0.55,
                                boot = 1000)
summary(phyglm_UN_nest_highplusflex_fix) # this model converges
confint(phyglm_UN_nest_highplusflex_fix)

# save model
saveRDS(phyglm_UN_nest_highplusflex_fix, here("Models/UN", "phyglm_UN_nest_highplusflex_fix.rds"))
# load model
phyglm_UN_nest_highplusflex_fix <- readRDS(here("Models/UN", "phyglm_UN_nest_highplusflex_fix.rds"))


# compare with non-phylogenetic model
glm_UN_nest_highplusflex <- logistf(Urban ~ NestSite_HighplusFlex + scale(Mass_log), 
                           data = UN_NestHighplusFlex)
summary(glm_UN_nest_highplusflex)
# some slight differences for coefficients but same conclusions reached


# get alpha, t, and half life for the model
(phyglm_UN_nest_highplusflex_scale$mean.tip.height) # t
(alpha_Nlow <- phyglm_UN_nest_highplusflex_scale$alpha) # alpha
(hl_NLow <- log(2)/alpha_Nlow) # half life
# compared to t, this is a small half life


######################## UN and Nest Site Low plus Flexible vs High ##########################
# 0 = low plus flexible
# 1 = high 

UN_NestLowplusFlex <- C_Nest_dat2 %>% 
  filter(!is.na(Urban)) %>%
  filter(!is.na(NestSite_LowplusFlex)) %>% column_to_rownames(., var = "Species_Jetz")
length(UN_NestLowplusFlex$NestSite_LowplusFlex)
# 128 species

###### add and pair tree

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_NestLowplusFlex_phydat <- geiger::treedata(tree_out, UN_NestLowplusFlex, sort=T)

UN_NestLowplusFlex_phy <- UN_NestLowplusFlex_phydat$phy
UN_NestLowplusFlex_dat <- as.data.frame(UN_NestLowplusFlex_phydat$data)

str(UN_NestLowplusFlex_dat)
length(UN_NestLowplusFlex_dat$NestSite_LowplusFlex)

# convert traits of interest to numeric
UN_NestLowplusFlex_dat$Urban <- as.numeric(UN_NestLowplusFlex_dat$Urban)
UN_NestLowplusFlex_dat$Mass_log <- as.numeric(UN_NestLowplusFlex_dat$Mass_log)
UN_NestLowplusFlex_dat$NestSite_LowplusFlex <- as.numeric(UN_NestLowplusFlex_dat$NestSite_LowplusFlex)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
phyglm_UN_nest_lowplusflex_scale <- phyloglm(Urban ~ NestSite_LowplusFlex + scale(Mass_log), 
                                              data = UN_NestLowplusFlex_dat, 
                                              phy = UN_NestLowplusFlex_phy,
                                              boot = 1000) 

# fails to converge

# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ NestSite_LowplusFlex + scale(Mass_log), 
                 data = UN_NestLowplusFlex_dat, 
                 phy = UN_NestLowplusFlex_phy,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylogenetic signal)


# try fixing alpha near upper bounds
set.seed(273)
phyglm_UN_nest_lowplusflex_fix <- phyloglm(Urban ~ NestSite_LowplusFlex + scale(Mass_log), 
                                       data = UN_NestLowplusFlex_dat, 
                                       phy = UN_NestLowplusFlex_phy,
                                      log.alpha.bound = 4.05,
                                       start.alpha = 0.55,
                                      boot = 1000)
summary(phyglm_UN_nest_lowplusflex_fix) # this model converges
confint(phyglm_UN_nest_lowplusflex_fix)

# save model
saveRDS(phyglm_UN_nest_lowplusflex_fix, here("Models/UN", "phyglm_UN_nest_lowplusflex_fix.rds"))
# load model
phyglm_UN_nest_lowplusflex_fix <- readRDS(here("Models/UN", "phyglm_UN_nest_lowplusflex_fix.rds"))


# compare with non-phylogenetic model
glm_UN_nest_lowplusflex <- logistf(Urban ~ NestSite_LowplusFlex + scale(Mass_log), 
                                    data = UN_NestLowplusFlex)
summary(glm_UN_nest_lowplusflex)
# some slight differences for coefficients but same conclusions reached

# get alpha, t, and half life for the model
(phyglm_UN_nest_lowplusflex_scale$mean.tip.height) # t
(alpha_Nlow <- phyglm_UN_nest_lowplusflex_scale$alpha) # alpha
(hl_NLow <- log(2)/alpha_Nlow) # half life
# compared to t, this is a small half life

############################### UN and High vs Low (flexible species excluded) ########################

# Run model with all flexible species excluded that focuses only on species that use either HIGH or LOW nest sites 
UN_NestLowHigh<- C_Nest_dat2 %>% 
  filter(!is.na(Urban)) %>%
  filter(!is.na(NestSite_LowHigh)) %>% 
  column_to_rownames(., var = "Species_Jetz")
length(UN_NestLowHigh$NestSite_LowHigh)
# 103 species

###### add and pair tree
tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_NestLowHigh_phydat <- geiger::treedata(tree_out, UN_NestLowHigh, sort=T)

UN_NestLowHigh_phy <- UN_NestLowHigh_phydat$phy
UN_NestLowHigh_dat <- as.data.frame(UN_NestLowHigh_phydat$data)

str(UN_NestLowHigh_dat)
length(UN_NestLowHigh_dat$NestSite_LowHigh)

# convert traits of interest to numeric
UN_NestLowHigh_dat$Urban <- as.numeric(UN_NestLowHigh_dat$Urban)
UN_NestLowHigh_dat$Mass_log <- as.numeric(UN_NestLowHigh_dat$Mass_log)
UN_NestLowHigh_dat$NestSite_LowHigh <- as.numeric(UN_NestLowHigh_dat$NestSite_LowHigh)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(712)
phyglm_UN_nest_lowhigh_scale <- phyloglm(Urban ~ NestSite_LowHigh + scale(Mass_log), 
                                             data = UN_NestLowHigh_dat, 
                                             phy = UN_NestLowHigh_phy,
                                             boot = 1000) 

summary(phyglm_UN_nest_lowhigh_scale)
confint(phyglm_UN_nest_lowhigh_scale)

# save model
saveRDS(phyglm_UN_nest_lowhigh_scale, here("Models/UN", "phyglm_UN_nest_lowhigh_scale.rds"))

# compare with non-phylogenetic model
glm_UN_nest_lowhigh <- logistf(Urban ~ NestSite_LowHigh + scale(Mass_log), 
                                   data = UN_NestLowHigh)
summary(glm_UN_nest_lowhigh)
# some slight differences for coefficients but same conclusions reached

# get alpha, t, and half life for the model
(phyglm_UN_nest_lowhigh_scale$mean.tip.height) # t
(alpha_Nlow <- phyglm_UN_nest_lowhigh_scale$alpha) # alpha
(hl_NLow <- log(2)/alpha_Nlow) # half life
# compared to t, this is a small half life

##########################################################################
##########################################################################

######################## UAI and Nest Safety ##########################

# create a new data frame that contains only species with both UAI and nest safety
UAI_NestSafety <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(nest.safety)) %>% as.data.frame()
length(UAI_NestSafety$nest.safety)
# 760 species with UAI and nest.safety

###### add and pair tree

# add rownames to data
row.names(UAI_NestSafety) <- UAI_NestSafety$Species_Jetz

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_NestSafety_phydat <- geiger::treedata(tree_out, UAI_NestSafety, sort=T)

UAI_NestSafety_phy <- UAI_NestSafety_phydat$phy
UAI_NestSafety_dat <- as.data.frame(UAI_NestSafety_phydat$data)

str(UAI_NestSafety_dat)
length(UAI_NestSafety_dat$nest.safety)

### convert traits of interest to numeric
UAI_NestSafety_dat$aveUAI <- as.numeric(UAI_NestSafety_dat$aveUAI)
UAI_NestSafety_dat$Mass_log <- as.numeric(UAI_NestSafety_dat$Mass_log)
UAI_NestSafety_dat$nest.safety <- as.numeric(UAI_NestSafety_dat$nest.safety)


# run phylogenetic linear model
UAI_GLS_nest_safety <- gls(aveUAI~ nest.safety + Mass_log, data = UAI_NestSafety_dat, 
                            correlation = corPagel(0.5, phy = UAI_NestSafety_phy, fixed=F, form = ~Species_Jetz), 
                            method = "ML") 

# model summary and results
summary(UAI_GLS_nest_safety) 
confint(UAI_GLS_nest_safety)

# model diagnostics
check_model(UAI_GLS_nest_safety) 
qqnorm(resid(UAI_GLS_nest_safety)) 
qqline(resid(UAI_GLS_nest_safety))
hist(resid(UAI_GLS_nest_safety))

# save model
saveRDS(UAI_GLS_nest_safety, here("Models/UAI", "UAI_GLS_nest_safety.rds"))


######################## MUTI and Nest Safety ##########################

# create a new data frame that contains only species with both MUTI and nest safety
MUTI_NestSafety <- C_Nest_dat2 %>% filter(!is.na(MUTIscore))  %>%
  filter(!is.na(nest.safety)) %>% as.data.frame()
length(MUTI_NestSafety$nest.safety)
# 125 species with MUTIscore and nest.safety

###### add and pair tree

# add rownames to data
row.names(MUTI_NestSafety) <- MUTI_NestSafety$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_NestSafety_phydat <- geiger::treedata(tree_out, MUTI_NestSafety, sort=T)

MUTI_NestSafety_phy <- MUTI_NestSafety_phydat$phy
MUTI_NestSafety_dat <- as.data.frame(MUTI_NestSafety_phydat$data)

str(MUTI_NestSafety_dat)
length(MUTI_NestSafety_dat$nest.safety)

# convert traits of interest to numeric
MUTI_NestSafety_dat$MUTIscore <- as.numeric(MUTI_NestSafety_dat$MUTIscore)
MUTI_NestSafety_dat$Mass_log <- as.numeric(MUTI_NestSafety_dat$Mass_log)
MUTI_NestSafety_dat$nest.safety <- as.numeric(MUTI_NestSafety_dat$nest.safety)


# run a phylogenetic linear model
MUTI_GLS_nest_safety <- gls(MUTIscore~ nest.safety + Mass_log, data = MUTI_NestSafety_dat, 
                           correlation = corPagel(0.5, phy = MUTI_NestSafety_phy, fixed = F, form = ~Species_Jetz), 
                           method = "ML") 

# model summary and results
summary(MUTI_GLS_nest_safety) 
confint(MUTI_GLS_nest_safety)

# model diagnostics
check_model(MUTI_GLS_nest_safety) 
qqnorm(resid(MUTI_GLS_nest_safety)) 
qqline(resid(MUTI_GLS_nest_safety))
hist(resid(MUTI_GLS_nest_safety))

# save model
saveRDS(MUTI_GLS_nest_safety, here("Models/MUTI", "MUTI_GLS_nest_safety.rds"))


######################## UN and Nest Safety ##########################


# create a new data frame that contains only species with both UN and nest safety
UN_NestSafety <- C_Nest_dat2 %>% filter(!is.na(Urban))  %>% 
  filter(!is.na(nest.safety)) %>% column_to_rownames(., var="Species_Jetz")
length(UN_NestSafety$nest.safety)
# 128 species with Urban and nest.safety

###### add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_NestSafety_phydat <- geiger::treedata(tree_out, UN_NestSafety, sort=T)

UN_NestSafety_phy <- UN_NestSafety_phydat$phy
UN_NestSafety_dat <- as.data.frame(UN_NestSafety_phydat$data)

str(UN_NestSafety_dat)
length(UN_NestSafety_dat$nest.safety)

### convert traits of interest to numeric
UN_NestSafety_dat$Urban <- as.numeric(UN_NestSafety_dat$Urban)
UN_NestSafety_dat$Mass_log <- as.numeric(UN_NestSafety_dat$Mass_log)
UN_NestSafety_dat$nest.safety <- as.numeric(UN_NestSafety_dat$nest.safety)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(753)
phyglm_UN_nest_safety_scale <- phyloglm( Urban ~ scale(nest.safety) + scale(Mass_log), 
                                         data = UN_NestSafety_dat, 
                                         phy = UN_NestSafety_phy, 
                                         boot = 1000) 
summary(phyglm_UN_nest_safety_scale) # this model converges
confint(phyglm_UN_nest_safety_scale)


# save model
saveRDS(phyglm_UN_nest_safety_scale, here("Models/UN", "phyglm_UN_nest_safety_scale.rds"))
# load model
phyglm_UN_nest_safety_scale <- readRDS(here("Models/UN", "phyglm_UN_nest_safety_scale.rds"))


# get alpha, t, and half life for the model
(phyglm_UN_nest_safety_scale$mean.tip.height) # t
(alpha_Nsafe <- phyglm_UN_nest_safety_scale$alpha) # alpha
(hl_Nsafe <- log(2)/alpha_Nsafe) # half life
# compared to t, this is small half life

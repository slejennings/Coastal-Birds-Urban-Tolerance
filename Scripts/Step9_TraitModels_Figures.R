##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Step 9: Trait model figures
# Authors: Sarah L. Jennings, Emma M. Garrison
##########################################################

# the purpose of this script is to create results figures for significant (or notable) trait models. 
# prioritizing models where results for a predictor trait are supported across multiple indexes, or differ across multiple indexes


# load packages
library(here)
library(tidyverse)
library(ape)
library(geiger)
library(nlme)
library(effects)
library(colorspace)
library(phylolm)
library(patchwork)
library(ggdist)
library(marginaleffects)
library(logistf)
library(scales)


##################################################################
################## LIFE HISTORY FIGURES ########################## 
##################################################################

################# UAI and Brood Value #######################
# load the model
UAI_GLS_bv <- readRDS(here("Models/UAI", "UAI_GLS_bv.rds"))

# load the data used for the model (this requires several steps. Target df is UAI_Brood_dat)
C_LifeHist_dat <- readRDS(here("Outputs", "Coastal_Species_LifeHistory.rds"))
C_LifeHist_dat2 <- C_LifeHist_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
C_LifeHist_dat2$Urban <- ifelse(C_LifeHist_dat2$Urban == "U", 1, 0)

# create a new data frame that contains only species with both UAI and brood values
UAI_Brood <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(brood_value)) %>% as.data.frame()

# add rownames to data
row.names(UAI_Brood) <- UAI_Brood$Species_Jetz

# add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))
UAI_Brood_phydat <- treedata(tree_out, UAI_Brood, sort=T)
UAI_Brood_phy <- UAI_Brood_phydat$phy
UAI_Brood_dat <- as.data.frame(UAI_Brood_phydat$data)

# convert traits of interest to numeric
UAI_Brood_dat$aveUAI <- as.numeric(UAI_Brood_dat$aveUAI)
UAI_Brood_dat$Mass_log <- as.numeric(UAI_Brood_dat$Mass_log)
UAI_Brood_dat$brood_value <- as.numeric(UAI_Brood_dat$brood_value)

#### make the plot for UAI and Brood Value ####

bv.UAIdb <- predictorEffect("brood_value" , UAI_GLS_bv)
BV.UAI.DF <-data.frame(bv.UAIdb)

BV_UAI_plot <- ggplot(data = BV.UAI.DF, aes(x = brood_value, y = fit)) +
  geom_line(color="#FFCDA1",lwd=1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#FFCDA1",alpha = .4, lwd = .1) +
  xlim(-1000,1000) +
  geom_point(data = UAI_Brood_dat, aes(x=jitter(brood_value, 1), y = aveUAI), color = "#ED6F00", bg = "#ED6F00", alpha = .4, size = 2, pch = 21) +
  coord_cartesian(ylim = c(-0, 3.5), xlim =c(-6,-1.7)) + theme_classic() +
  theme(axis.text.x =  element_text(color="black", size = 11),
        axis.text.y =  element_text(color="black", size = 11), 
        axis.title.x = element_text(size = 12, margin = margin(t = 8)), 
        axis.title.y = element_text(size =12), 
        legend.position = "none") +
  xlab("Brood Value") + 
  ylab("UAI") 

BV_UAI_plot

#######################################################################
######################## UAI and Clutch size ##########################

# load the model
UAI_GLS_clutch <- readRDS(here("Models/UAI", "UAI_GLS_clutch.rds"))

# load the data used for the model (this requires several steps. Target df is UAI_Clutch_dat)
# Note: need to have created C_LifeHist_dat2 in previous section
UAI_Clutch <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(clutch_size)) %>% as.data.frame()

# add rownames to data
row.names(UAI_Clutch) <- UAI_Clutch$Species_Jetz

# add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))
UAI_Clutch_phydat <- treedata(tree_out, UAI_Clutch, sort=T)
UAI_Clutch_phy <- UAI_Clutch_phydat$phy
UAI_Clutch_dat <- as.data.frame(UAI_Clutch_phydat$data)

# convert traits of interest to numeric
UAI_Clutch_dat$aveUAI <- as.numeric(UAI_Clutch_dat$aveUAI)
UAI_Clutch_dat$Mass_log <- as.numeric(UAI_Clutch_dat$Mass_log)
UAI_Clutch_dat$clutch_size <- as.numeric(UAI_Clutch_dat$clutch_size)

#### make the plot for UAI and clutch size ####

clutch.UAIdb <- predictorEffect("clutch_size", UAI_GLS_clutch)
clutch.UAI.DF <- data.frame(clutch.UAIdb)

clutch_UAI_plot<-ggplot(data = clutch.UAI.DF, aes(x = clutch_size, y = fit)) +
  geom_line(color = "#FFCDA1", lwd = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax=upper), fill="#FFCDA1", alpha=.4, lwd=.1)+
  xlim(-1000, 1000) +
  geom_point(data = UAI_Clutch_dat, aes(x = jitter(clutch_size, 1), y = aveUAI), color = "#ED6F00", bg = "#ED6F00", alpha=.4, size = 2, pch = 21) +
  coord_cartesian(ylim = c(-0, 4.1), xlim = c(0, 15.5)) + 
  theme_classic() +
  theme(axis.text.x =  element_text(color="black", size = 11),
        axis.text.y =  element_text(color="black", size = 11), 
        axis.title.x = element_text(size = 12, margin = margin(t = 8)), 
        axis.title.y = element_text(size = 12), 
        legend.position = "none") +
  xlab("Clutch Size") + 
  ylab("UAI") 

clutch_UAI_plot

#######################################################################
################### MUTI and Developmental Mode ##########################

# load the data used for the model (this requires several steps. Target df is MUTI_Develop_dat)
# Note: need to have created C_LifeHist_dat2 in previous section
MUTI_Develop <- C_LifeHist_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(developmental_mode)) %>% as.data.frame()

# add rownames to data
row.names(MUTI_Develop) <- MUTI_Develop$Species_Jetz

# add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))
MUTI_Develop_phydat <- treedata(tree_out, MUTI_Develop, sort=T)
MUTI_Develop_phy <- MUTI_Develop_phydat$phy
MUTI_Develop_dat <- as.data.frame(MUTI_Develop_phydat$data)

# convert traits of interest to numeric
MUTI_Develop_dat$MUTIscore <- as.numeric(MUTI_Develop_dat$MUTIscore)
MUTI_Develop_dat$Mass_log <- as.numeric(MUTI_Develop_dat$Mass_log)
#MUTI_Develop_dat$developmental_mode <- as.numeric(MUTI_Develop_dat$developmental_mode)
MUTI_Develop_dat$developmental_mode <- as.factor(MUTI_Develop_dat$developmental_mode)

# Run phylogenetic linear model with developmental mode as a factor
MUTI_GLS_develop <- gls(MUTIscore~ developmental_mode + Mass_log, data = MUTI_Develop_dat, 
                        correlation = corPagel(0.5, phy=MUTI_Develop_phy, fixed=F, form = ~Species_Jetz), 
                        method = "ML") 
summary(MUTI_GLS_develop)

# get predicted means and 95% CI
predMUTIDevelop <- avg_predictions(MUTI_GLS_develop, by="developmental_mode") %>%
  mutate(Developmental_Mode = as.factor(ifelse(developmental_mode == 0, "Precocial", "Altricial")))
colnames(predMUTIDevelop)

####### make the plot for MUTI and developmental mode ############

MUTI_Develop_dat_plot <- MUTI_Develop_dat %>%
  mutate(Developmental_Mode = as.factor(ifelse(developmental_mode == 0, "Precocial", "Altricial"))) # create factor version of developmental mode with altricial/precocial labels


# set the scale parameters for stat_slab and stat_dot to something around 0.5 (default is 0.9) so that the two slabs do not overlap
# scale the halfeye slab thickness by n (the number of observations in each group) so that the area of each slab represents sample size (and looks similar to the total area of its corresponding dotplot)
# this is achieved with the thickness argument
develop_MUTI_plot <- MUTI_Develop_dat_plot %>%
  ggplot(aes(y = MUTIscore, x = Developmental_Mode, fill = Developmental_Mode)) +
  stat_slab(aes(thickness = after_stat(pdf*n)), scale = 0.6, side = "right", alpha = 0.7, , show.legend = F) +
  stat_dots(side = "left", scale = 0.6, slab_linewidth = NA, show_interval=F, alpha = 0.7, show.legend = F) +
  scale_fill_manual(values=c("#009380", "#B5DD60")) +
  geom_pointinterval(data = predMUTIDevelop, aes(y = estimate, x = Developmental_Mode, ymin = conf.low, ymax = conf.high), 
                     size=4, lwd=3, show.legend = F) +
  theme_classic() +
  labs(x="Developmental Mode", y="MUTI") +
  theme(axis.text.x =  element_text(color="black", size = 11, margin = margin(t = 5)),
        axis.text.y =  element_text(color="black", size = 11), 
        axis.title.x = element_text(size =12, margin = margin(t = 8)), 
        axis.title.y = element_text(size =12)) 

develop_MUTI_plot


######################## UN and Developmental Mode ##########################
# load the data used for the model (this requires several steps. Target df is UN_Develop_dat)
# Note: need to have created C_LifeHist_dat2 in previous section
UN_Develop <- C_LifeHist_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(developmental_mode)) %>% column_to_rownames(., var = "Species_Jetz")

# add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))
UN_Develop_phydat <- treedata(tree_out, UN_Develop, sort=T)
UN_Develop_phy <- UN_Develop_phydat$phy
UN_Develop_dat <- as.data.frame(UN_Develop_phydat$data)

# convert traits of interest to numeric
UN_Develop_dat$Urban <- as.numeric(UN_Develop_dat$Urban)
UN_Develop_dat$Mass_log <- as.numeric(UN_Develop_dat$Mass_log)
UN_Develop_dat$developmental_mode <- as.factor(UN_Develop_dat$developmental_mode)

# at the time of this analysis, there was no way to extract predictions from a phyloglm model
# however, there was little/no support to suggest a phylogenetic signal and we reached very similar findings with a non-phylogenetic model
# we can extract predictions and plot marginal effects using this model
glm_UN_develop <- logistf(Urban ~ developmental_mode + scale(Mass_log), 
                          data = UN_Develop)

summary(glm_UN_develop)
phyglm_UN_develop_fix <- readRDS(here("Models/UN", "phyglm_UN_develop_fix.rds")) # import phyloglm version for comparison
summary(phyglm_UN_develop_fix)

# get predictions from the logistf version of the model
avg_predictions(glm_UN_develop, by="developmental_mode")
# this gives the probability that a species is Urban (Urban = 1) when their developmental mode is either Precocial (0) or Altricial (1)
# the probability that a species is classified as Urban is 29.2% when they are precocial, and therefore 70.8% for Non-urban and precocial
# the probability that a species is classified as Urban is 41% when they are altricial, and therefore 59% for Non-urban and altricial

# put these values into a data frame for plotting
develop_UN_summary <- data.frame(UN = as.factor(c("Urban", "Non-Urban", "Urban", "Non-Urban")),
                         Developmental_Mode = as.factor(c("Precocial", "Precocial", "Altricial", "Altricial")),
                         Percentage = as.numeric(c("29.2", "70.8", "41", "59"))) %>%
  mutate(Group = paste(UN, Developmental_Mode, sep="_"))

####### make the plot for UN and Developmental Mode #########

# make a stacked bar chart
develop_UN_plot <- ggplot(develop_UN_summary, aes(x = Developmental_Mode, y = Percentage, fill = Group)) + 
  geom_bar(stat = "identity", position = "fill", color = "black") + 
  scale_y_continuous(labels = scales::percent) +  # show percentages on the y-axis
  # Custom colors for each Group & Urban combination
  scale_fill_manual(
    values = c(
      "Non-Urban_Precocial" = "#C1DBEC",   # Color for "Precocial & Nonurban"
      "Urban_Precocial" = "#5C90C6",        # Color for "Precocial & Urban"
      "Non-Urban_Altricial" = "#C1DBEC",  # Color for "Altricial & Nonurban"
      "Urban_Altricial" = "#5C90C6"        # Color for "Altricial & Urban"
    )) + 
  labs(x = "Developmental Mode", y = "UN") + 
  theme_classic() + 
  theme(axis.text.x = element_text(color="black", size = 11, margin = margin(t = 5)),
        axis.text.y = element_text(color="black", size = 11),
        axis.title.x = element_text(size = 12, margin = margin(t = 8)), 
        axis.title.y = element_text(size = 12), 
        legend.position = "none") + 
  annotate(
    "text", 
    x = c(1, 1, 2, 2), 
    y = c(0.2, 0.7, 0.2, 0.7), 
    label = c("Urban", "Non\n Urban", "Urban", "Non\n Urban"),
    size = 11/.pt)

print(develop_UN_plot)

#######################################################################
############### Combine life history trait plots ######################
# plot all the life history results together

all_life_history_plots <- BV_UAI_plot + clutch_UAI_plot + develop_MUTI_plot + develop_UN_plot +
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)"))) & theme(plot.tag = element_text(size = 12))

all_life_history_plots

ggsave(here("Results", "LifeHistoryTraits_Plot.png"), width = 17, height = 11, units="cm", dpi=300)

##################################################################
################## NEST TRAIT FIGURES ############################ 
##################################################################

# Using nest site high models to plot the results
# For these models: 0 = Low, 1 = High plus flexible
#  0 = all species that exclusively use low nesting sites
# 1 = all species that use high nesting sites either exclusively or because they exhibit flexible in nest site height

################# UAI and Nest Site High #######################

# load the data used for the model (this requires several steps. Target df is UAI_NestHigh_dat)
# load in "Coastal_Species_Nest.rds"
C_Nest_dat <- readRDS(here("Outputs", "Coastal_Species_Nest.rds"))
C_Nest_dat2 <- C_Nest_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
C_Nest_dat2$Urban <- ifelse(C_Nest_dat2$Urban == "U", 1, 0)

# create a new data frame that contains only species with both UAI and nest site high
UAI_NestHigh <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(NestSite_High)) %>% as.data.frame()

# add rownames to data
row.names(UAI_NestHigh) <- UAI_NestHigh$Species_Jetz

# add and pair tree
tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))
UAI_NestHigh_phydat <- treedata(tree_out, UAI_NestHigh, sort=T)
UAI_NestHigh_phy <- UAI_NestHigh_phydat$phy
UAI_NestHigh_dat <- as.data.frame(UAI_NestHigh_phydat$data)

# convert traits of interest to numeric or factor
UAI_NestHigh_dat$aveUAI <- as.numeric(UAI_NestHigh_dat$aveUAI)
UAI_NestHigh_dat$Mass_log <- as.numeric(UAI_NestHigh_dat$Mass_log)
UAI_NestHigh_dat$NestSite_High <- as.factor(UAI_NestHigh_dat$NestSite_High)

# run phylogenetic linear model with nest site height as a factor
UAI_GLS_nest_high <- gls(aveUAI~ NestSite_High + Mass_log, data = UAI_NestHigh_dat, 
                         correlation = corPagel(0.5, phy = UAI_NestHigh_phy, fixed = F, form = ~Species_Jetz), 
                         method = "ML")


# get predicted means and 95% CI
predUAINest <- avg_predictions(UAI_GLS_nest_high, by="NestSite_High") %>%
  mutate(NestSiteHeight = as.factor(ifelse(NestSite_High == 0, "Low", "High & Flexible")),
         NestSiteHeight = fct_relevel(NestSiteHeight, "Low")) 
colnames(predUAINest)

################ make the plot for UAI and nest site height ####################

# set the scale parameters for stat_slab and stat_dot to something around 0.5 (default is 0.9) so that the two slabs do not overlap
# scale the halfeye slab thickness by n (the number of observations in each group) so that the area of each slab represents sample size (and looks similar to the total area of its corresponding dotplot)
# this is achieved with the thickness argument

nestsite_UAI_plot <- UAI_NestHigh_dat %>%
  mutate(NestSiteHeight = as.factor(ifelse(NestSite_High == 0, "Low", "High & Flexible")),
         NestSiteHeight = fct_relevel(NestSiteHeight, "Low")) %>%
  ggplot(aes(y = aveUAI, x = NestSiteHeight, fill = NestSiteHeight)) +
  stat_slab(aes(thickness = after_stat(pdf*n)), scale = 0.6, side = "right", alpha = 0.7, , show.legend = F) +
  stat_dots(side = "left", scale = 0.6, slab_linewidth = NA, show_interval=F, alpha = 0.7, show.legend = F) + # do not plot interval using raw data
  scale_fill_manual(values=c("#FDB16E", "#ED6F00")) + 
  geom_pointinterval(data = predUAINest, aes(y = estimate, x = NestSiteHeight, ymin = conf.low, ymax = conf.high), # use predicted means and 95% CI to plot intervals
                     size=4, lwd=3, show.legend = F) +
  theme_classic() +
  labs(y="UAI") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 13, margin = margin(t = 5)),
        axis.text.y = element_text(size = 13))


nestsite_UAI_plot

#################################################################
################# MUTI and Nest Site High #######################

# create a new data frame that contains only species with both MUTI and nest site high
MUTI_NestHigh <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(NestSite_High)) %>% as.data.frame()

# add rownames to data
row.names(MUTI_NestHigh) <- MUTI_NestHigh$Species_Jetz

# add and pair tree
tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))
MUTI_NestHigh_phydat <- treedata(tree_out, MUTI_NestHigh, sort=T)
MUTI_NestHigh_phy <- MUTI_NestHigh_phydat$phy
MUTI_NestHigh_dat <- as.data.frame(MUTI_NestHigh_phydat$data)

# convert traits of interest to numeric or factor
MUTI_NestHigh_dat$MUTIscore <- as.numeric(MUTI_NestHigh_dat$MUTIscore)
MUTI_NestHigh_dat$Mass_log <- as.numeric(MUTI_NestHigh_dat$Mass_log)
MUTI_NestHigh_dat$NestSite_High <- as.factor(MUTI_NestHigh_dat$NestSite_High)

# run phylogenetic linear model
MUTI_GLS_nest_high <- gls(MUTIscore~ NestSite_High + Mass_log, data = MUTI_NestHigh_dat, 
                          correlation = corPagel(0.5, phy = MUTI_NestHigh_phy, fixed = F, form = ~Species_Jetz), 
                          method = "ML") 

# get predicted means and 95% CI
predMUTINest <- avg_predictions(MUTI_GLS_nest_high, by="NestSite_High") %>%
  mutate(NestSiteHeight = as.factor(ifelse(NestSite_High == 0, "Low", "High & Flexible")),
         NestSiteHeight = fct_relevel(NestSiteHeight, "Low")) 
colnames(predMUTINest)

################ make the plot for MUTI and nest site height ####################

# set the scale parameters for stat_slab and stat_dot to something around 0.5 (default is 0.9) so that the two slabs do not overlap
# scale the halfeye slab thickness by n (the number of observations in each group) so that the area of each slab represents sample size (and looks similar to the total area of its corresponding dotplot)
# this is achieved with the thickness argument

nestsite_MUTI_plot <- MUTI_NestHigh_dat %>%
  mutate(NestSiteHeight = as.factor(ifelse(NestSite_High == 0, "Low", "High & Flexible")),
         NestSiteHeight = fct_relevel(NestSiteHeight, "Low")) %>%
  ggplot(aes(y = MUTIscore, x = NestSiteHeight, fill = NestSiteHeight)) +
  stat_slab(aes(thickness = after_stat(pdf*n)), scale = 0.5, side = "right", alpha = 0.7, , show.legend = F) +
  stat_dots(side = "left", scale = 0.5, slab_linewidth = NA, show_interval=F, alpha = 0.7, show.legend = F) + # do not plot interval using raw data
  scale_fill_manual(values=c("#009380", "#B5DD60")) + 
  geom_pointinterval(data = predMUTINest, aes(y = estimate, x = NestSiteHeight, ymin = conf.low, ymax = conf.high), # use predicted means and 95% CI to plot intervals
                     size=4, lwd=3, show.legend = F) +
  theme_classic() +
  labs(y="MUTI") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 13, margin = margin(t = 5)),
        axis.text.y = element_text(size = 13))


nestsite_MUTI_plot

###############################################################
################# UN and Nest Site High #######################
# load the data used for the model (this requires several steps. Target df is UN_Develop_dat)
# Note: need to have created C_LifeHist_dat2 in previous section
UN_NestHigh <- C_Nest_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(NestSite_High)) %>% as.data.frame()

# add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))
UN_NestHigh_phydat <- treedata(tree_out, UN_NestHigh, sort=T)
UN_NestHigh_phy <- UN_NestHigh_phydat$phy
UN_NestHigh_dat <- as.data.frame(UN_NestHigh_phydat$data)

# convert traits of interest to numeric
UN_NestHigh_dat$Urban <- as.numeric(UN_NestHigh_dat$Urban)
UN_NestHigh_dat$Mass_log <- as.numeric(UN_NestHigh_dat$Mass_log)
UN_NestHigh_dat$NestSite_High <- as.factor(UN_NestHigh_dat$NestSite_High)

# at the time of this analysis, there was no way to extract predictions from a phyloglm model
# however, there was little/no support to suggest a phylogenetic signal and we reached very similar findings with a non-phylogenetic model
# we can extract predictions and plot marginal effects using this model
glm_UN_nest_high <- logistf(Urban ~ NestSite_High + scale(Mass_log), 
                            data = UN_NestHigh)

summary(glm_UN_nest_high) 
phyglm_UN_nest_high_fix <- readRDS(here("Models/UN", "phyglm_UN_nest_high_fix.rds")) # import phyloglm version for comparison
summary(phyglm_UN_nest_high_fix)

# get predictions from the logistf version of the model
avg_predictions(glm_UN_nest_high, by="NestSite_High")
# this gives the probability that a species is Urban (Urban = 1) when their nest site is either Low (0) or High and Flexible (1)
# the probability that a species is classified as Urban is 25.8% when they use Low nest sites, and therefore 74.2% for Non-urban and Low nest sites
# the probability that a species is classified as Urban is 42.2% when they use High nest sites or are flexible, and therefore 57.8% for Non-urban

# put these values into a data frame for plotting
NestSite_UN_summary <- data.frame(UN = as.factor(c("Urban", "Non-Urban", "Urban", "Non-Urban")),
                                 NestSiteHeight = as.factor(c("Low", "Low", "High & Flexible", "High & Flexible")),
                                 Percentage = as.numeric(c("25.8", "74.2", "42.2", "57.8"))) %>%
  mutate(Group = paste(UN, NestSiteHeight, sep="_"),
         NestSiteHeight = fct_relevel(NestSiteHeight, "Low"))

####### make the plot for UN and Nest Site Height #########

# make a stacked bar chart
nestsite_UN_plot <- ggplot(NestSite_UN_summary, aes(x = NestSiteHeight, y = Percentage, fill = Group)) + 
  geom_bar(stat = "identity", position = "fill", color = "black") + 
  scale_y_continuous(labels = scales::percent) +  # show percentages on the y-axis
  # Custom colors for each Group & Urban combination
  scale_fill_manual(
    values = c(
      "Non-Urban_Low" = "#C1DBEC",   
      "Urban_Low" = "#5C90C6",        
      "Non-Urban_High & Flexible" = "#C1DBEC",  
      "Urban_High & Flexible" = "#5C90C6"        
    )) + 
  labs(y = "UN") + 
  theme_classic() + 
  theme(axis.text.x = element_text(color="black", size = 13, margin = margin(t = 5)),
        axis.text.y = element_text(color="black", size = 13),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14), 
        legend.position = "none") + 
  annotate(
    "text", 
    x = c(1, 1, 2, 2), 
    y = c(0.15, 0.65, 0.15, 0.65), 
    label = c("Urban", "Non\n Urban", "Urban", "Non\n Urban"),
    size = 11/.pt)

print(nestsite_UN_plot)


################################################################
############ Combine nest site height plots ####################
all_nestsite_plots <- nestsite_UAI_plot + nestsite_MUTI_plot + nestsite_UN_plot + 
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)")), caption = "Nest Site", ) & theme(plot.tag = element_text(size = 14), plot.caption = element_text(size = 14, hjust=0.5, margin = margin(t = 5)))

print(all_nestsite_plots)

ggsave(here("Results", "NestSite_Plot.png"), width = 25, height = 9, units="cm", dpi=300)

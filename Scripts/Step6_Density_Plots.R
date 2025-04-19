##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Step 6: Make density plots for the 3 UT indexes
# Authors: Sarah L. Jennings, Emma M. Garrison
##########################################################

# The objective of this script is to make a density plot to show the distribution 
# of Urban Tolerance scores across each of the three indexes (MUTI, UAI, and UN)


# load packages
library(tidyverse)
library(here)
library(ggdist)
library(patchwork)


# import data set of coastal species and urban tolerance index scores
AllIndexesCoastal <- readRDS(here("Outputs", "Coastal_Birds_List.rds"))
colnames(AllIndexesCoastal)
str(AllIndexesCoastal)
AllIndexesCoastal$Urban <- ifelse(AllIndexesCoastal$Urban == "U", 1, 0)

# how many species have ALL three index scores?
species_with_all_indexes <- AllIndexesCoastal %>% 
  filter(!is.na(aveUAI)) %>% 
  filter(!is.na(Urban)) %>% 
  filter(!is.na(MUTIscore)) 

nrow(species_with_all_indexes)
# there are only 49 species that are represented by all three indexes
unique(species_with_all_indexes$Species_Jetz)


# Also import the list of all bird species with urban tolerance scores to include as a comparison with coastal species
AllIndexes_AllBirds <- readRDS(here("Outputs", "UAI_MUTI_UN_final.rds"))
colnames(AllIndexes_AllBirds)
str(AllIndexes_AllBirds)
AllIndexes_AllBirds$Urban <- ifelse(AllIndexes_AllBirds$Urban == "U", 1, 0)


#######################################################################
##### Make density plot for UAI 

# simplify all birds list to only those with a UAI score
# first, add a group column
AllIndexes_AllBirds$Group <- "All"
colnames(AllIndexes_AllBirds)

# next, find all species with UAI scores
AllBirds_UAI <- AllIndexes_AllBirds %>% 
  filter(!is.na(aveUAI) & is.finite(aveUAI))
nrow(AllBirds_UAI) #4347

# find the average UAI score for all bird species
mean(AllBirds_UAI$aveUAI) # 1.182

AllBirds_UAI_2 <- AllBirds_UAI %>% 
  select(aveUAI, Group)

# simplify the coastal birds list to only those with a UAI score
# first, add a group column
AllIndexesCoastal$Group <- "Coastal"
colnames(AllIndexesCoastal)

# next, find all coastal species with UAI scores
CoastalBirds_UAI <- AllIndexesCoastal %>% 
  filter(!is.na(aveUAI) & is.finite(aveUAI))

# find the average UAI score for coastal bird species
mean(CoastalBirds_UAI$aveUAI) # 1.433

CoastalBirds_UAI_2 <- CoastalBirds_UAI %>% 
  select(aveUAI, Group)

### combine data for UAI to make plot
combined_UAI_for_density <- bind_rows(AllBirds_UAI_2, CoastalBirds_UAI_2)
nrow(combined_UAI_for_density) 

# Create the plot
UAI_density <- ggplot() + 
  # Density plots according to groups 
  geom_density(
    data = combined_UAI_for_density, 
    aes(x = aveUAI, fill = Group), 
    alpha = 0.5, 
    adjust = 0.5, 
  ) + 
  geom_boxplot(
    data = combined_UAI_for_density,
    aes(x = aveUAI, fill = Group), 
    width = 0.12, 
    color = 'black',
    alpha = 0.5, 
    position = position_nudge(y = -0.1), 
    show.legend = F
  ) +  
  scale_fill_manual(
    values = c("All" = "#FFCDA1", "Coastal" = "#ED6F00")
  ) +
  labs(
    x = "UAI", 
    y = "Density"
  ) + 
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 8)),  # Adjust axis title size
    axis.title.y = element_text(margin = margin(r = 10), size = 14), 
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.position = c(.95, .85),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),            
    legend.title = element_blank(),         # Remove the legend title
    legend.text = element_text(size = 10)
  )

print(UAI_density)


#######################################################################
##### Make density plot for MUTI 

# simplify all bird species list to only those with a MUTI score
AllBirds_MUTI <- AllIndexes_AllBirds %>% 
  filter(!is.na(MUTIscore) & is.finite(MUTIscore))
nrow(AllBirds_MUTI) # 431

AllBirds_MUTI_2 <- AllBirds_MUTI %>% 
  select(MUTIscore, Group)
nrow(AllBirds_MUTI_2) # 431 looks good! 

# find the average MUTI score for all bird species
mean(AllBirds_MUTI_2$MUTIscore) #-0.0051

# simplify the coastal species list to only those with a MUTI score
CoastalBirds_MUTI <- AllIndexesCoastal %>% 
  filter(!is.na(MUTIscore) & is.finite(MUTIscore))
nrow(CoastalBirds_MUTI) # 128

CoastalBirds_MUTI_2 <- CoastalBirds_MUTI %>% 
  select(MUTIscore, Group)
nrow(CoastalBirds_MUTI_2) # 128

# find the average MUTI score for coastal species
mean(CoastalBirds_MUTI_2$MUTIscore) # 0.0176

### combine data for MUTI to make plot
combined_MUTI_for_density <- bind_rows(AllBirds_MUTI_2, CoastalBirds_MUTI_2)
nrow(combined_MUTI_for_density) 

# Create the plot
MUTI_density <- ggplot() + 
  # Density plots according to groups 
  geom_density(
    data = combined_MUTI_for_density, 
    aes(x = MUTIscore, fill = Group), 
    alpha = 0.5, 
    adjust = 0.5, 
  ) + 
  # Boxplot for both groups i believe? 
  geom_boxplot(
    data = combined_MUTI_for_density,
    aes(x = MUTIscore, fill = Group), 
    width = 0.054, 
    color = 'black',
    alpha = 0.5, 
    position = position_nudge(y = -0.045), 
    show.legend = F
  ) +  
  scale_fill_manual(
    values = c("All" = "#B5DD60", "Coastal" = "#17A77E")
  ) +
  labs(
    x = "MUTI", 
    y = "Density"
  ) + 
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 8)),  # Adjust axis title size
    axis.title.y = element_text(margin = margin(r = 13), size = 14), 
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.position = c(.95, .85),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),            
    legend.title = element_blank(), # Remove the legend title
    legend.text = element_text(size = 10)
  )

print(MUTI_density)

###################################### UN DENSITY PLOT ###################################### 

### As UN is a binomial variable, use percentages in a stacked barplot rather than a density plot

# simplify all bird species list to only those with a UN score
AllBirds_UN <- AllIndexes_AllBirds %>% 
  filter(!is.na(Urban) & is.finite(Urban))
nrow(AllBirds_UN) # 533

AllBirds_UN_2 <- AllBirds_UN %>% 
  select(Urban, Group)
nrow(AllBirds_UN_2)

# How many species are Urban? 
sum(AllBirds_UN_2$Urban == 1, na.rm = TRUE) # 206 out of 533 total species
sum(AllBirds_UN_2$Urban == 1, na.rm = TRUE)/nrow(AllBirds_UN_2) # 38.6% are urban

# simplify the coastal species list to only those with a UN score
CoastalBirds_UN <- AllIndexesCoastal %>% 
  filter(!is.na(Urban) & is.finite(Urban))
nrow(CoastalBirds_UN) # 128

CoastalBirds_UN_2 <- CoastalBirds_UN %>% 
  select(Urban, Group)
nrow(CoastalBirds_UN_2) # 128

# How many coastal species are Urban? 
sum(CoastalBirds_UN_2$Urban == 1, na.rm = TRUE) # 41
sum(CoastalBirds_UN_2$Urban == 1, na.rm = TRUE)/nrow(CoastalBirds_UN_2) # 32% are urban for coastal species
1-(sum(CoastalBirds_UN_2$Urban == 1, na.rm = TRUE)/nrow(CoastalBirds_UN_2)) # 68% are non-urban for coastal species

### combine data for UN to make plot
combined_UN_for_density <- bind_rows(AllBirds_UN_2, CoastalBirds_UN_2)
nrow(combined_UN_for_density) # 661

combined_UN_for_density_summary <- combined_UN_for_density %>%
  group_by(Group, Urban) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(Group) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  mutate(Group_Urban = paste(Group, Urban, sep = " & "))  # Combine Group and Urban for legend


UN_stacked_barplot <- ggplot(combined_UN_for_density_summary, aes(x = Group, y = percentage, fill = Group_Urban)) + 
  geom_bar(stat = "identity", position = "fill", color = "black") + 
  scale_y_continuous(labels = scales::percent) +  # Show percentages on the y-axis
  # Custom colors for each Group & Urban combination
  scale_fill_manual(
    values = c(
      "All & 0" = "#DEEEF7",   # Color for "All & Nonurban"
      "All & 1" = "#A1C4E0",        # Color for "All & Urban"
      "Coastal & 0" = "#7FABD3",  # Color for "Coastal & Nonurban"
      "Coastal & 1" = "#5C90C6"        # Color for "Coastal & Urban" 
    )
  ) + 
  labs(
    x = "UN", 
    y = "Percentage", 
    fill = ""  # Legend title
  ) + 
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 8)), 
    axis.title.y = element_text(size = 14), 
    axis.text.x = element_text(size=13),
    axis.text.y = element_text(size=13),
    legend.position = "none"
  ) + 
  annotate(
    "text", 
    x = c(1, 1, 2, 2), 
    y = c(0.2, 0.7, 0.2, 0.7), 
    label = c("Urban", "Non\n Urban", "Urban", "Non\n Urban"),
    size = 13/.pt)

print(UN_stacked_barplot)

######################################################################
######################################################################
######################################################################

# arrange all three plots side by side 

all_density_3 <- UAI_density + MUTI_density + UN_stacked_barplot + 
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)"))) & theme(plot.tag = element_text(size = 14))

print(all_density_3)

ggsave(here("Results", "UT_Density_Plots.png"), width=25, height = 9, units="cm", dpi=300)

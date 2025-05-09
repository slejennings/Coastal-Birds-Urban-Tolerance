##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Step 3: Coastal Birds Selection
# Authors: Sarah L. Jennings, Emma M. Garrison
##########################################################

### Objective: refine list of all bird species (with UAI, MUTI, and/or UN scores) to COASTAL BIRDS ONLY 

# load packages
library(here)
library(tidyverse)
library(magrittr)

### Starting point - load in "UAI_MUTI_UN_final.rds" file 
# it is all species with at least one of the three urban tolerance index scores 

AllBirds <- readRDS(here("Outputs", "UAI_MUTI_UN_final.rds"))

print(AllBirds)
nrow(AllBirds) # 4433

# extract a data frame that contains a list of family names
unique_family_names <- AllBirds %>%
  distinct(Family_eBird) %>%
  tidyr::separate(., col = Family_eBird, # this splits the column Family_eBird into two columns
                  into = c("Family_Sci", "Family_English"), # the first column is called Family_Sci and the second Family_English
                  sep =  "^\\S*\\K\\s+") # split at the first space encountered

nrow(unique_family_names) # 205 families
print(unique_family_names)

# download this data frame as a csv, and assess all Families using Birds of the World

write.csv(unique_family_names, here("Notes", "families_all.csv"))

########################## Round 1 ###############################################

### Using Birds of the World (Cornell Lab or Ornithology), we classified each family represented by AllBirds dataframe as either: 
# "Coastal" or "Not-Coastal", based on whether the family habitat page mentioned coastal areas. 
# The column "Coastal" categorizes these families. 
# "Yes" = Coastal family 
# "No" = Non-coastal family 

# If the BOW family page mentioned use of coastal habitats (such as shorelines, beaches, estuaries, mangroves, etc.),
# then the family (and all representative species) was marked as "Yes" (coastal)
# If the BOW family page did not mention use of coastal habitats or resources, then the all representative species 
# within the family were marked as "No" (not-coastal).

####### For all species from a "Coastal" (Yes) family, mark them as a coastal species. 

#read in the csv that contains all families and their coastal status, as determined via investigation on Birds of the World (2022). 
List_of_Families <- read.csv(here("Notes", "families_all_coastal.csv")) %>% select(-X)
View(List_of_Families)
colnames(List_of_Families) 

# how many families were marked as "Coastal"? 
Families_Yes <- List_of_Families %>% 
  filter(Coastal == "Yes")
nrow(Families_Yes) 
# 37 families marked as Coastal 


# how many families were marked as "Non Coastal"? 
Families_No <- List_of_Families %>% 
  filter(Coastal == "No")
nrow(Families_No) 
# 168 families marked as non-coastal 

168 + 37
# 205 total


# reformat AllBirds Family_eBird column so that it is compatible to join with List_of_Families 
AllBirds$Family_Sci <- word(AllBirds$Family_eBird, 1)
head(AllBirds)

# join them together 
Coastal_Round_1 <- left_join(AllBirds, List_of_Families, by = "Family_Sci")
print(Coastal_Round_1)
nrow(Coastal_Round_1) # 4433 - correct number of rows in AllBirds, which was the left part of left_join 

#################################### Round 2 ############################################
# In this round, we will sort through the common names of species that were marked as "Yes" for Coastal 

Round_1_yes <- Coastal_Round_1 %>% 
  filter(Coastal == "Yes")

print(Round_1_yes)
nrow(Round_1_yes)
# 823 

# we will generate a list of descriptive terms from the common names of the 4433 species in AllBirds
# we will manually filter this list to identify words in the common names that relate to habitat and/or diet

twowordnames <- AllBirds %>%
  select(CommonName_eBird) %>%
  mutate(space = str_count(CommonName_eBird, " ")) %>% # identify how many words each species has in its names using number of spaces
  filter(space == 1) %>% # keep only species with names that contain two words (a single space)
  separate_wider_delim(CommonName_eBird, delim = " ", names = c("descriptor1", "species")) 
# some of these still have hyphenated terms in the species column that could be descriptive

threewordnames <- AllBirds %>%
  select(CommonName_eBird) %>%
  mutate(space = str_count(CommonName_eBird, " ")) %>% # identify how many words each species has in its names using number of spaces
  filter(space == 2) %>% # keep only species with names that contain two words (two spaces)
  separate_wider_delim(CommonName_eBird, delim = " ", names = c("descriptor1", "descriptor2", "species") )
# some of these still have hyphenated terms in the species column that could be descriptive

hyphennames <- data.frame(species = c(twowordnames$species, threewordnames$species)) %>%
  mutate(hyphen = str_count(species, "-")) %>% # identify species with hyphenated names
  filter(hyphen == 1) %>% # keep only species with names that contain a hyphen
  separate_wider_delim(species, delim = "-", names = c("descriptor1", "species") )

# put all the descriptive words from common names together into a single data frame and remove duplicates
alldescriptors <- data.frame( # make a data frame
  descriptors = sort( # arrange in alphabetical order
    c(twowordnames$descriptor1, threewordnames$descriptor1, threewordnames$descriptor2, hyphennames$descriptor1))) %>% 
  distinct() # remove duplicated descriptors

# export as csv 
# two researchers (EMG and SLJ) individually reviewed list and manually marked all terms that are associated with:
# a habitat (but not a specific place)
# a diet
# an Ocean or Sea (e.g., Pacific, Caribbean)
# or indicated a species association or interaction between the bird and another species (e.g., plant, insect, mammal)
write.csv(alldescriptors, here("Notes", "descriptors.csv"))

# import the marked up list
# use column coastal_Y.N.M to sort
descriptors_marked <- read.csv(here("Notes", "descriptors_marked.csv"), header=T)
head(descriptors_marked)

# first focus on non-coastal terms
noncoastal_terms <- descriptors_marked %>% 
  filter(coastal_Y.N.M =="N") %>% # keep all terms that are clearly non-coastal
  select(descriptors) # select the column of descriptive terms

# convert to a character vector 
NCterms_vec <- as.vector(noncoastal_terms[,1])
class(NCterms_vec)

# convert complete list of common names from Round_1_yes to a character vector 
head(Round_1_yes)
commname_vec <- as.vector(Round_1_yes[,5]) # common names are in 5th column
class(commname_vec)

# now, look through commname_vec which contains the species in Round_1_yes using the terms stored in NCterms_vec
# we are searching for key words that indicate NON-coastal habitats, diets or associations with other species that are non-coastal
# e.g., "freshwater", "alpine", "upland", "lake", "river", 
#"mountain", "prairie", "highland", "forest", "desert" 

R1_noncoastal <- map(NCterms_vec, str_detect, string = commname_vec) %>%
  reduce(`|`) %>% 
  magrittr::extract(commname_vec, .) %>%
  tibble()

colnames(R1_noncoastal) <- "CommonName_eBird"

print(R1_noncoastal, n=Inf)
nrow(R1_noncoastal)
# export data frame, look up species and mark ones to keep and remove
write.csv(R1_noncoastal, here("Notes", "Round_1_noncoastal.csv"))

# manually look up species and determine whether they should be marked as coastal or not
# read in edited csv
# use it to update coastal species list

# read in edited csv 
R1_noncoastal_edited <- read.csv(here("Notes", "Round_1_noncoastal_edited.csv"), header=T) %>%
  select(-X)

head(R1_noncoastal_edited)
colnames(R1_noncoastal_edited)
nrow(R1_noncoastal_edited) # 39
View(R1_noncoastal_edited)
# combine with Coastal_Round_1 to make the updates
# this takes a few steps


R1_noncoastal_edited_no <- R1_noncoastal_edited %>% 
  filter(Coastal == "No")
nrow(R1_noncoastal_edited_no) #18 species identified as NO from round 1 


# remind ourselves how many species are in the data
nrow(Coastal_Round_1) # 4433
colnames(Coastal_Round_1)


# get all the species that were not modified in any way in the past step (any species not in R1_noncoastal_edited)
Coastal_Round_1_edit1 <- anti_join(Coastal_Round_1, R1_noncoastal_edited, by="CommonName_eBird")
 

# does the number of rows in edit1 equal the number of rows in Round_1 minus the number of rows in noncoastal_edited?
# this should be "TRUE"
nrow(Coastal_Round_1_edit1) == nrow(Coastal_Round_1) - nrow(R1_noncoastal_edited)

# update the coastal classification for all species in R1_noncoastal_edited
Coastal_Round_1_edit2 <- Coastal_Round_1 %>% 
  select(-Coastal) %>% # removes this columns because it will be updated with new version in the next step
  left_join(R1_noncoastal_edited, ., by="CommonName_eBird")
nrow(Coastal_Round_1_edit2) # should be the same number as in R1_noncoastal_edited (39)

# bind everything back together
Coastal_Round_2 <- bind_rows(Coastal_Round_1_edit1, Coastal_Round_1_edit2)  

# make sure all the species are still present
nrow(Coastal_Round_2) == nrow(Coastal_Round_1)

# double check how many species are marked as "Coastal" now! 
Coastal_Round_2 %>% filter(Coastal == "Yes") %>% nrow()
# 823 - 18 = 805


#################################### Round 3 ###################################
# In this round, we will sort through the common names of species that were marked as "No" for Coastal
Round_1_no <- Coastal_Round_1 %>% 
  filter(Coastal == "No")

print(Round_1_no)
nrow(Round_1_no) # 3610 


# develop list of coastal descriptor terms to search the Common Names of species marked as not coastal
coastal_terms <- descriptors_marked %>% 
  filter(coastal_Y.N.M %in% c("M", "Y")) %>% # keep all terms that are clearly coastal (Y) and possibly coastal (M)
  select(descriptors) # select the column of descriptive terms

# convert to a character vector 
C.terms_vec <- as.vector(coastal_terms[,1])
class(C.terms_vec)

# convert complete list of common names from Round_1_yes to a character vector 
head(Round_1_no)
commname_vec_2 <- as.vector(Round_1_no[,5]) # common names are in 5th column

# now, edit this Round_1_no.csv file (containing all species from families marked as "No" in Round 1) -> 
#search through common names for coastal-identifier words stored in CTerms_vec (i.e. "Coastal", "Sea", "Tide", "Beach", "Mangrove", etc.) 
#for all flagged species, look on Birds of the World (2022) species page for mentions of coastal habitat/resource use. 

R1_coastal <- map(C.terms_vec, str_detect, string = commname_vec_2) %>%
  reduce(`|`) %>% 
  magrittr::extract(commname_vec_2, .) %>%
  tibble()

R1_coastal
nrow(R1_coastal) # 221

# export data frame, look up species and mark ones to keep and remove
write.csv(R1_coastal, here("Notes", "Round_1_coastal.csv"))

# # # # # # # # # # # # # # # # # # 

# upload edited csv 
R1_coastal_edited <- read.csv(here("Notes", "Round_1_coastal_edited.csv")) %>% select(-X)

# View to inspect 
View(R1_coastal_edited)
colnames(R1_coastal_edited)
nrow(R1_coastal_edited) # 221


# how many Coastal = Yes species were added? 
R1_coastal_edited_yes <- R1_coastal_edited %>% 
  filter(Coastal == "Yes")
nrow(R1_coastal_edited_yes) # 99 species identified as coastal! 


# combine with Coastal_Round_2 from the previous round to make the updates
# this takes a few steps

# remind ourselves how many species are in the data
nrow(Coastal_Round_2) # 4433

# get all the species that were not modified in any way in the past step (any species not in R1_coastal_edited)
Coastal_Round_2_edit1 <- anti_join(Coastal_Round_2, R1_coastal_edited, by="CommonName_eBird")

# does the number of rows in edit1 equal the number of rows in Round_2 minus the number of rows in coastal_edited?
# this should be "TRUE"
nrow(Coastal_Round_2_edit1) == nrow(Coastal_Round_2) - nrow(R1_coastal_edited)
#TRUE
nrow(R1_coastal_edited) # 221

# update the coastal classification for all species in R1_coastal_edited
Coastal_Round_2_edit2 <- Coastal_Round_2 %>% 
  select(-Coastal, -Notes_Species) %>% # remove these columns because they will be updated in the next step
  left_join(R1_coastal_edited, ., by="CommonName_eBird")
nrow(Coastal_Round_2_edit2) # should be the same number as in R1_coastal_edited (221) - yep 

# bind everything back together
Coastal_Round_3 <- bind_rows(Coastal_Round_2_edit1, Coastal_Round_2_edit2)  

# make sure all the species are still present
nrow(Coastal_Round_3) == nrow(Coastal_Round_2)

# how many species are on our coastal birds list now?
Coastal_Round_3 %>% filter(Coastal == "Yes") %>% nrow() 
# 904 species identified as Coastal! 

########################### Final Steps (rounds 4 and 5) ###########################################################
# make another pass to identify any additional coastal species or species that are not actually coastal 
# to do this, we will use the diet info columns that came from Wilman et al. 2014 (Elton traits) 

head(AllBirds)
# we need to search the original species list of 4433 species

# import Elton diet traits (from Wilman et al. 2014)
elton <- read.csv(here("Data", "elton.csv"), header=T)
head(elton)
colnames(elton) # look at column names to identify columns that could be useful

# these trait columns seem useful:
unique(elton$Diet.Vfish) # percentage of diet that is fish. Filter to retain any species with > 0
unique(elton$ForStrat.watbelowsurf) # percentage of time spent foraging below surf. Filter to retain any species with >0
unique(elton$ForStrat.wataroundsurf) # percentage of time spent foraging around surf. Filter to retain any species with >0
unique(elton$PelagicSpecialist) # Pelagic seabirds. Has values 0 or 1. Filter to retain species listed as 1


# join elton traits and AllBirds
AllBirds_elton <- elton %>%
  rename(Species_Jetz = Scientific) %>%
  select(Species_Jetz, Diet.Vfish, Diet.5Cat, ForStrat.watbelowsurf, # retain columns identified as useful above
         ForStrat.wataroundsurf, PelagicSpecialist) %>%
  left_join(AllBirds, ., by = "Species_Jetz")


# note that the filter requirements are OR statements 
# this means a bird should be retained if they have Diet.Vfish > 0 OR 
# if they do any of their foraging below surf OR
# if they forage around water OR
# if they are a pelagic specialist
# finally, we use a third filter to find all birds that meet the above requirements BUT are NOT currently classified as coastal
# we will want to investigate these species to see if they warrant inclusion 

coastaldiet <- AllBirds_elton %>% 
  filter(Diet.Vfish > 0 |   # Now apply a of filtering requirements. First, keep any species with some fish in diet OR
           ForStrat.watbelowsurf > 0 | # keep any species that do any of their foraging below surf OR
           ForStrat.wataroundsurf >0 | # keep any species that do any of their foraging around surf OR
           PelagicSpecialist == 1)  # keep any species that are classified as Pelagic Specialists

nrow(coastaldiet) # 688 species

# Are any birds in the coastaldiet list currently marked as non-coastal?
head(Coastal_Round_3)
R3_No <- Coastal_Round_3 %>% 
  filter(Coastal == "No")
R3_noncoastal <- inner_join(R3_No, coastaldiet)

nrow(R3_noncoastal)
# 58 species --> these are species that have been marked as "No" over the past 3 rounds, but their diet traits suggest a possible coastal association. 

# Export this list, look up all species and check whether they should be coastal
write.csv(R3_noncoastal, here("Notes", "Round_3_noncoastal.csv"))

#######################################

# Are any birds that are marked as coastal that don't have a diet with fish or forage in/around water?
# NOTE: species can still be coastal and not eat fish or foraging in/around water
# so this step may not result in the removal of any species from the coastal list
# we are using this as a final check to flag any species that may need a second look
R3_Yes <- Coastal_Round_3 %>% 
  filter(Coastal == "Yes") 

R3_coastal <- anti_join(R3_Yes, coastaldiet) %>%
  arrange(Family_eBird)

nrow(R3_coastal)
# 274
# many of these appear to be coastal
# we will export the list, double check their classification as coastal and re-import the list with any needed changes

# export to look up species in BOW
write.csv(R3_coastal, here("Notes", "Round_3_coastal.csv"))

# import edited coastal file
R3_coastal_edited <- read.csv(here("Notes", "Round_3_coastal_edited.csv"), header = T) %>% select(-X)
head(R3_coastal_edited)

# How many species that were originally classified as Coastal have diet traits that suggest they need to be re-classified? 
R3_coastal_edited_no <- R3_coastal_edited %>% 
  filter(Coastal == "No")
nrow(R3_coastal_edited_no) # 102 species re-classified as non-coastal


# import edited noncoastal file 
R3_noncoastal_edited <- read.csv(here("Notes", "Round_3_noncoastal_edited.csv"), header = T) %>% select(-X)

# How many species that were originally classified as non-coastal need to be reclassified after considering diet traits?
R3_noncoastal_edited_yes <- R3_noncoastal_edited %>% 
  filter(Coastal == "Yes")
nrow(R3_noncoastal_edited_yes) # 17 species identified as coastal 

# combine them into one object that contains all possible changes
R3_edits <- bind_rows(R3_coastal_edited, R3_noncoastal_edited)

# get all the species that were not modified in any way in the past step 
Coastal_Round_4_edit1 <- anti_join(Coastal_Round_3, R3_edits, by="CommonName_eBird")
# note: you need to specify to use CommonName_eBird for this join for it to work correctly

# does the number of rows in edit1 equal the number of rows in Round_3 minus the number of rows in R3_edits ?
# this should be "TRUE"
nrow(Coastal_Round_4_edit1) == nrow(Coastal_Round_3) - nrow(R3_edits)

colnames(R3_edits)
colnames(Coastal_Round_3)
# update the coastal classification for all species in R3_edits
Coastal_Round_4_edit2 <- Coastal_Round_3 %>% 
  select(-Coastal, -Notes_Species) %>% # removes these columns because they will be updated with new versions in the next step
  left_join(R3_edits, ., by="CommonName_eBird")

nrow(Coastal_Round_4_edit2) # 332 - should be the same number as in R3_edits
nrow(R3_edits) # 332

# bind everything back together
Coastal_All <- bind_rows(Coastal_Round_4_edit1, Coastal_Round_4_edit2)  

# make sure all the species are still present
nrow(Coastal_Round_3) == nrow(Coastal_All)
# TRUE 

# View Coastal_All for double-checking 
View(Coastal_All)

# how many species do we have? 
nrow(Coastal_All) # should equal 4433

Coastal_Only <- Coastal_All %>% 
  filter(Coastal == "Yes")

nrow(Coastal_Only)
# 819 birds are classified as Coastal

###########################################################################
# One last final check: we need all the species to have unique BirdTree names as we need to join them to the phylogeny that uses this naming scheme for the models
# Are the names unique for the BirdTree taxonomy? No
length(unique(Coastal_Only$Species_Jetz))

# look at the duplicates for BirdTree names
# these need to be resolved because we can only have one entry per name in the models
duplicates_Jetz <- Coastal_Only %>% count(Species_Jetz) %>%
  filter(n>1) %>%
  left_join(., Coastal_Only)

nrow(duplicates_Jetz) # 36

# how many of these have UAI?
duplicates_Jetz %>% filter(!is.na(aveUAI)) %>% nrow() # all of them

# for duplicated species that only have UAI and no other urban index
# we will just keep one entry

# pull out the UN species
UN_dups <- duplicates_Jetz %>% 
  select(Species_Jetz, Urban) %>% 
  filter(!is.na(Urban))  %>%
  select(Species_Jetz) %>%
  distinct() %>%
  mutate(UN = "Yes")

# pull out MUTI species
MUTI_dups <- duplicates_Jetz %>% 
  select(Species_Jetz, MUTIscore) %>% 
  filter(!is.na(MUTIscore))  %>%
  select(Species_Jetz) %>%
  distinct() %>%
  mutate(MUTI = "Yes")

# find duplicates that only have UAI and not UN or MUTI
UAIonly <- duplicates_Jetz %>% select(Species_Jetz) %>% distinct() %>%
  left_join(.,MUTI_dups) %>% left_join(., UN_dups) %>%
  filter(is.na(MUTI)) %>% filter(is.na(UN))

# for these species, keep the entry where BirdTree, BirdLife and eBird names are the same
duplicates_Jetz_UAIonly <- UAIonly %>%
  select(-UN, -MUTI) %>%
  left_join(., duplicates_Jetz)

fix1 <- duplicates_Jetz_UAIonly %>%
  mutate(match = 
    if_else(Species_BirdLife == Species_Jetz & Species_BirdLife == Species_eBird, "1", "0")) %>%
  filter(match == 1) %>%
  select(-n, -match)
# fix 1 is the first set of "resolved" duplicates. 

# now look at species that have both UAI and UN
UAI_UN <- duplicates_Jetz %>% select(Species_Jetz) %>% distinct() %>%
  left_join(.,MUTI_dups) %>% left_join(., UN_dups) %>%
  filter(!is.na(UN)) %>% filter(is.na(MUTI))
nrow(UAI_UN)

# for these species, keep the entry that has UN and UAI
duplicates_Jetz_UAI_UN <- UAI_UN %>%
  select(-UN, -MUTI) %>%
  left_join(., duplicates_Jetz)

fix2 <- duplicates_Jetz_UAI_UN %>%
  filter(!is.na(Urban)) %>%
  select(-n)
# this still has a single duplicated species, so we will keep the entry with shared species name across 3 naming schemes
fix2 <- fix2 %>%
  filter(Species_eBird != "Curruca cantillans")

# now look at species that have both UAI and MUTI
UAI_MUTI <- duplicates_Jetz %>% select(Species_Jetz) %>% distinct() %>%
  left_join(.,MUTI_dups) %>% left_join(., UN_dups) %>%
  filter(!is.na(MUTI)) %>% filter(is.na(UN))
nrow(UAI_MUTI)

# for these species, keep the entry that has UN and MUTI
duplicates_Jetz_UAI_MUTI <- UAI_MUTI %>%
  select(-UN, -MUTI) %>%
  left_join(., duplicates_Jetz)

fix3 <- duplicates_Jetz_UAI_MUTI %>%
  filter(!is.na(MUTIscore)) %>%
  select(-n)
# there is still one species that is duplicated because it has two MUTI scores
# look both species up in Birds of the World
# keeping Pacific Wren as its distribution is more coastal whereas Winter Wren occurs inland too
fix3 <- fix3 %>%
  filter(Species_eBird != "Troglodytes hiemalis")

# put the already fixed species together
# this will allow identification of all remaining duplicated Jetz names
fix123 <- bind_rows(fix1, fix2, fix3)

# which duplicates are still left? Use anti_join to find out
remaining_dups <- fix123 %>% select(Species_Jetz) %>%
  anti_join(duplicates_Jetz, .)

print(remaining_dups)

# making some decisions based on looking up the species on Birds of the World to resolve these final few
# the criteria used to decide which species to keep is described here:
# keep Clapper Rail and drop Ridgeway's Rail because Clapper rail overlaps more urban areas and has a larger range
# keep Eurasian moorhen and drop Common gaillnule because Eurasian moorhen has a wider distribution and occurs outside North America/Europe in areas less well represented in ornithological research
# keep Snowy plover and drop Kentish plover as it appears to use coastal habitats more based on range map
# keep Gray-cowled Wood-rail and drop Russet-napped Wood-rail. These were one species until recently. Gray-cowled has a large range with potential to overlap more urban areas

fix4 <- remaining_dups %>% 
  filter(CommonName_eBird %in% c("Clapper Rail", "Eurasian Moorhen", "Snowy Plover", "Gray-cowled Wood-Rail")) %>%
  select(-n)

allfixes <- bind_rows(fix123, fix4) # put all resolved duplicates together
nrow(allfixes) # 17

# get all the coastal species that did not have duplicate issues
Coastal_Only_nondups <- allfixes %>% select(Species_Jetz) %>%
  anti_join(Coastal_Only, .)
nrow(Coastal_Only_nondups)

# put together allfixes and Coastal_Only_nodups to get final and clean Coastal birds list
Coastal_Birds <- bind_rows(Coastal_Only_nondups, allfixes)

nrow(Coastal_Birds)# should be equal to the number of unique BirdTree species names (800)
length(unique(Coastal_Birds$Species_Jetz)) # 800

# View to confirm  
View(Coastal_Birds)

# save as rds 
saveRDS(Coastal_Birds, here("Outputs", "Coastal_Birds_List.rds"))

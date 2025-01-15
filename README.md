# GitHub Repository Read Me

 ## **Overview**

This repository contains data and code for the analyses in the manuscript Life history and nesting traits reflect urban tolerance in coastal birds. The analyses are organized into an R Project, which will reproduce the results, tables, and figures presented in the manuscript. This Read Me file describes the required software and the organization of the R Project and the associated files.

## **Correspondence**

Please direct questions about the data, analysis, and results to:

Sarah L. Jennings: sjenni02[at]calpoly.edu

Clinton D. Francis: cdfranci[at]calpoly.edu

## **Software**

R programming language version 4.3.2 (2023-10-31)

RStudio IDE version 2024.04.2+764

R packages with version number: ape (5.8), colorspace (2.1-0), broom (1.0.6), broom.mixed (0.2.9.5), effects (4.2-2), easystats (0.7.3), fuzzyjoin (0.1.6), geiger (2.0.11), ggdist (3.3.2), ggnewscale (0.4.10), ggtree (3.13.0), ggtreeExtra (1.15.0), here (1.0.1), logistf (1.26.0), nlme (3.1.164), openxlsx (4.2.5.2), patchwork (1.2.0), phylolm (2.6.4), phyr (1.1.0), phytools (2.1-1), scales (1.3.0), stringi (1.8.4), tidyverse (2.0.0), treeio (1.26.0)


## **Description of Folders and Files**
This repository contains an R project with various folders that are organized and named for their contents. Below we have provided a description of each folder, and the files contained within each.
Scripts Folder

***Description:*** contains R scripts to analyze the data, to generate the findings in the manuscript, and to create the tables and results figures. The scripts are sequential and are labeled accordingly (Step1 through Step 9).

***Contents:*** 15 files

1) Step1_UAINamestoSpecies.R

   *File description:* aligns common names from the UAI urban tolerance index with scientific names from the 3 avian naming schemes/taxonomies (Bird Tree, Bird Life, and eBird)


2) Step2_Join_MUTI_UAI_UN.R
   
   *File description:* join urban tolerance scores from MUTI and UN indexes to the list created in Step 1 


3) Step3_Coastal_Birds_Selection.R

    *File description:* filtering process to reduce complete species list to a refined list that contains on coastal bird species that will be used in subsequent analyses


4) Step4_CombineTraitswithSpecies.R

   *File description:* combine trait values from publicly available trait databases to the coastal species list


5) Step5_PhyloDiversityMetrics.R

    *File description:* obtain various phylogenetic diversity metrics for the 3 urban tolerance indexes and create a figure with the phylogeny and urban tolerance scores for coastal species


6) Step6_Density_Plots.R

    *File description:* create density plots to show the distribution of each urban tolerance metric within all bird species and the subset of coastal species


7) Step7_TraitModels_BODYMASS.R, Step7_TraitModels_DIET.R, Step7_TraitModels_LIFEHISTORY.R, Step7_TraitModels_NEST.R, Step7_TraitModels_SENSORY.R, Step7_TraitModels_SOCIAL.R, Step7_TraitModels_SSELECTION.R

    *File description:* fit phylogenetically-informed models to determine whether the urban tolerance of coastal birds is explained by ecological traits. There is one script for each category of traits as reflected in the script names.


8) Step8_TraitModels_ResultsTables.R

    *File description:* put all the results from the trait based models in to organized tables and export as an Excel spreadsheet


9) Step9_TraitModels_Figures.R

    *File description:* create figures that display the findings of the life history and nest traits models. 


### **Data Folder**

***Description:*** contains the datafiles that were used in the analysis, mostly stored as .csv. These files are avian taxonomic classification, species-specific urban tolerance classifications/scores, and/or species’ ecological traits from publicly available datasets. References for each source are provided.

***Contents:*** 17 files

1) amniote.csv

   *File description:* Life history traits for birds, mammals, and reptiles
  
    *Reference:* Myhrvold et al. 2015


2) avonet.csv

    *File description:* Morphological data for bird species

    *Reference:* Tobias et al. 2022


3) BirdNameConversion.csv

    *File description:* Naming schemes across 3 major avian taxonomies: BirdTree, BirdLife and eBird. From AVONET

    *Reference:* Tobias et al. 2022 


4) BirdTreeNames.csv

    *File description:* Scientific Family and Order information for species names in BirdTree taxonomy

    *Reference:* birdtree.org; Jetz et al. 2012

5) CoastalSppEyeGeometriesMain.csv

    *File description:* eye geometries for coastal bird species

    *Reference:* values come from multiple sources including this study (measured by members of our research group); Hall et al. 2009, Lindsey et al. 2013, and Ritland 1982. The values obtained
   by our research group represent averages taken from measurements of multiple museum specimens. See SEELabEyeMeasurements.csv for additional information.


7) SEELabEyeMeasurements.csv

    *File description:* eye geometries collected by our research group

    *Reference:* provides information about each measured specimen, including the morphometric values obtained, the museum/source where the specimen is located, the collection number, and      species. These values were averaged to produce the values provided in CoastalSppEyeGeometriesMain.csv. 


7) Delhey_2023_DS7.csv

    *File description:* Sexual selection and social traits for birds.Cooperative breeding information originally from Cockburn 2006. Territoriality originally from Tobias et al. 2016.        Developmental mode originally from Wang and Kimball 2016.

    *Reference:* Delhey et al. 2023


8) Dunn_2015_Dichromatism.csv

    *File description:* Plumage hue and brightness sexual dimorphism scores for bird species

    *Reference:* Dunn et al. 2015


9) ebird_taxonomy_v2022.csv

    *File description:* eBird/Clements taxonomy with scientific names for each species, their order, family, and common names

    *Reference:* Clements et al. 2022


10) elton.csv

    *File description:* Elton diet traits for bird species

    *Reference:* Wilman et al. 2014


11) Fanelli_Urban_Tolerance.csv

    *File description:* multivariate urban tolerance index (MUTI index) scores for bird species

    *Reference:* Fanelli et al. 2022


12) HuAndCardosoData.csv

    *File description:* urban/non-urban classification for bird species (UN index)

    *Reference:* Hu and Cardoso 2009


13) Jetz_ConsensusPhy.tre

    *File description:* avian phylogeny 

    *Reference:* birdtree.org; Jetz et al. 2012 


14) longevity.csv

    *File description:* generation lengths and maximum lifespan for bird species

    *Reference:* Bird et al. 2020


15) mikula_peakfreq.csv

    *File description:* peak vocal frequencies for bird species

    *Reference:* Mikula et al. 2021


16) nests.csv
 
    *File description:* nest traits for bird species

    *Reference:* Chia et al. 2023


17) UAI.csv

    *File description:* Urban association index (UAI) scores for bird species

    *Reference:* Neate-Clegg et al. 2023


### **Notes Folder**

***Description:*** contains various files stored as .csv files. Most contain lists of species that were exported from the R Project to be manually edited, often via research using Birds of the World (Billerman et al. 2022), to confirm their use of coastal/estuarine/marine habitats. One version of each file is the .cvs exported from R and the other version is the edited copy after our research process was complete. Annotations are provided in the associated scripts where these files are generated and used to explain their contents and the process that was undertaken to produce the edited/updated version.

***Contents:*** 18 files


### **Outputs Folder**

***Description:*** contains various files stored as .rds that contain organized data that is being moved between R scripts and used in subsequent analysis steps.

***Contents:*** 10 files


### **Models Folder**

***Description:*** contains the models produced by the R scripts stored as .rds files. There are 3 subfolders that refer to the different urban tolerance indexes.

***Contents for UAI subfolder:*** 20 files

***Contents for MUTI subfolder:*** 20 files

***Contents for UN subfolder:*** 20 files

### **Results Folder**

***Description:*** contains tables and figures that are presented in the manuscript that were generated by the scripts in this R Project

***Contents:*** 5 files


### ***References***

Bird, JP, R Martin, HR Akçakaya, J Gilroy, IJ Burfield, ST Garnett, A Symes, J Taylor, ÇH Şekercioğlu, and SHM Butchart. 2020. Generation Lengths of the World’s Birds and Their Implications for Extinction Risk. Conservation Biology 34 (5): 1252–61. https://doi.org/10.1111/COBI.13486.

Billerman, S, B Keeney, P Rodewalk, and T Schulenberg. 2022. Birds of the World. Ithaca, NY, USA: Cornell Laboratory of Ornithology. https://birdsoftheworld.org/bow/home

Clements, JF, TS Schulenberg, MJ IIiff, TA Fredericks, JA Gerbracht, D Lepage, SM Billerman, BL Sullivan, and CL Wood. 2022. The eBird/Clements checklist of Birds of the World: v2022. Downloaded from https://www.birds.cornell.edu/clementschecklist/introduction/updateindex/october-2022/2022-citation-checklist-download/ 

Chia, SY, YT Fang, YT Su, PY Tsai, C Hsieh, SH Tsao, JYJuang, CM Hung, and MN Tuanmu. 2023. A Global Database of Bird Nest Traits. Scientific Data 10 (1): 923. https://doi.org/10.1038/s41597-023-02837-1.

Cockburn, A. 2006. Prevalence of Different Modes of Parental Care in Birds. Proceedings of the Royal Society B: Biological Sciences 273 (1592): 1375–83. https://doi.org/10.1098/rspb.2005.3458.

Delhey, K, Mi Valcu, C Muck, J Dale, and B Kempenaers. 2023. Evolutionary Predictors of the Specific Colors of Birds. Proceedings of the National Academy of Sciences 120 (34): e2217692120. https://doi.org/10.1073/pnas.2217692120.

Dunn, PO, JK Armenta, and LA Whittingham. 2015. Natural and Sexual Selection Act on Different Axes of Variation in Avian Plumage Color. Science Advances 1 (2): e1400155. https://doi.org/10.1126/sciadv.1400155.

Fanelli, RE, PR Martin, OJ Robinson, and F Bonier. 2022. Estimates of Species-Level Tolerance of Urban Habitat in North American Birds. Ecology 103 (12): e3821. https://doi.org/10.1002/ECY.3821.
Hall, MI, AN Iwaniuk, and C Gutiérrez-Ibáñez. 2009. Optic Foramen Morphology and Activity Pattern in Birds. The Anatomical Record: Advances in Integrative Anatomy and Evolutionary Biology 292 (11): 1827–45. https://doi.org/10.1002/ar.21007.

Hu, Y, and GC Cardoso. 2009. Are Bird Species That Vocalize at Higher Frequencies Preadapted to Inhabit Noisy Urban Areas? Behavioral Ecology 20 (6): 1268–73. https://doi.org/10.1093/beheco/arp131.

Jetz, W, GH Thomas, JB Joy, K Hartmann, and AO Mooers. 2012. The Global Diversity of Birds in Space and Time. Nature 491 (7424): 444–48. https://doi.org/10.1038/nature11631.

Lisney, TJ, K Stecyk, J Kolominsky, BK Schmidt, JR Corfield, AN Iwaniuk, and DR Wylie. 2013. Ecomorphology of Eye Shape and Retinal Topography in Waterfowl (Aves: Anseriformes: Anatidae) with Different Foraging Modes. Journal of Comparative Physiology A 199 (5): 385–402. https://doi.org/10.1007/s00359-013-0802-1.

Mikula, P, M Valcu, H Brumm, M Bulla, W Forstmeier, T Petrusková, B Kempenaers, and T Albrecht. 2021. A Global Analysis of Song Frequency in Passerines Provides No Support for the Acoustic Adaptation Hypothesis but Suggests a Role for Sexual Selection. Ecology Letters 24 (3): 477–86. https://doi.org/10.1111/ele.13662.

Myhrvold, NP, E Baldridge, B Chan, D Sivam, DL Freeman, and SK Morgan Ernest. 2015. An Amniote Life-History Database to Perform Comparative Analyses with Birds, Mammals, and Reptiles. Ecology 96 (11): 3109–3109. https://doi.org/10.1890/15-0846R.1.

Neate-Clegg, MHC, BA Tonelli, C Youngflesh, JX Wu, GA Montgomery, ÇH Şekercioğlu, and MW Tingley. 2023. Traits Shaping Urban Tolerance in Birds Differ around the World. Current Biology 33 (9): 1677-1688.e6. https://doi.org/10.1016/J.CUB.2023.03.024.

Ritland, SM. 1982. The Allometry of the Vertebrate Eye. The University of Chicago. https://www.proquest.com/docview/303088065?pq-origsite=gscholar&fromopenview=true&sourcetype=Dissertations%20&%20Theses.

Tobias, JA, C Sheard, AL Pigot, AJM Devenish, J Yang, F Sayol, MHC Neate-Clegg, et al. 2022. AVONET: Morphological, Ecological and Geographical Data for All Birds. Ecology Letters 25 (3): 581–97. https://doi.org/10.1111/ele.13898.

Tobias, JA, C Sheard, N Seddon, A Meade, AJ Cotton, and S Nakagawa. 2016. Territoriality, Social Bonds, and the Evolution of Communal Signaling in Birds. Frontiers in Ecology and Evolution 4 (June). https://doi.org/10.3389/fevo.2016.00074.

Wang, N, and RT Kimball. 2016. Re-Evaluating the Distribution of Cooperative Breeding in Birds: Is It Tightly Linked with Altriciality? Journal of Avian Biology 47 (5): 724–30. https://doi.org/10.1111/jav.00869.

Wilman, H, J Belmaker, J Simpson, C de la Rosa, MM Rivadeneira, and W Jetz. 2014. EltonTraits 1.0: Species-Level Foraging Attributes of the World’s Birds and Mammals. Ecology 95 (7): 2027–2027. https://doi.org/10.1890/13-1917.1.






























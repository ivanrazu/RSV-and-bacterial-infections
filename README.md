# Sequential RSV and bacterial infections in the nasopharynx of Zambian infants and mothers

This repository contains the code and data analysis for the project investigating the interactions between Respiratory Syncytial Virus (RSV) and bacterial infections. The project specifically focuses on the effects sequential RSV-Bacterial (B) interactions with the bacterial species:
- **Streptococcus pneumoniae (SP)**
- **Moraxella catarrhalis (MC)**
- **Staphylococcus aureus (SA)**
- **Haemophilus influenzae (HI)**

## Table of Contents
1. [Project Overview](#project-overview)
2. [Directory Structure](#directory-structure)
3. [Data Description](#data-description)
6. [Scripts and Analysis](#scripts-and-analysis)
9. [License](#license)
10. [Contact](#contact)

## Project Overview
This project explores the longitudinal interactions between RSV and four bacterial pathogens (Streptococcus pneumoniae (SP), Moraxella catarrhalis (MC), Staphylococcus aureus (SA), and Haemophilus influenzae (HI)) using a birth cohort dataset from mother-infant pairs in Lusaka, Zambia. Individuals were categorized as either coinfected with RSV and one or more bacterial species, solely infected with one or more bacterial species, solely infected by RSV, or not infected by either bacteria or virus.

To simplify the analysis of potential interactions between RSV and each bacterial species, we examined sequential pairs of infections for each individual, focusing on one bacterial species at a time. Individuals were classified based on the first occurrence of RSV relative to a particular bacterial species (SP, MC, SA, or HI). For each potential viral-bacterial interaction, we categorized individuals into six groups:

1. **Bacteria First, RSV Later (B → RSV)**: Individuals initially infected with the specific bacterium and subsequently found infected with RSV at a later visit.
2. **Simultaneous Detection (B & RSV)**: Individuals in whom both RSV and the bacterium were detected for the first time at the same visit.
3. **RSV First, Bacteria Later (RSV → B)**: Individuals in whom RSV was detected first, followed by the specific bacterium at a later visit.
4. **Bacteria Only (B Only)**: Individuals who acquired the specific bacterium without any evidence of RSV infection (note that they may have been infected with other bacteria).
5. **RSV Only (RSV Only)**: Individuals infected with RSV but not the specific bacterium under consideration at the time of analysis (note that they may have been infected with other bacteria).
6. **Neither (Neither)**: Individuals who remained free of both RSV and the specific bacterium under consideration during the study duration (note that they may have been infected with other bacteria).

This classification scheme allows us to systematically investigate the interactions between RSV and each bacterial species within this population.

We use generalized additive models (GAM) to analyze and visualize the temporal relationship between RSV infection and subsequent bacterial carriage for each pathogen separately under these various scenarios.

## Directory Structure



## Data Description
- **Data Source**: The data analyzed in this project was collected during our previous project Southern Africa Mother Infant Pertusis Study (SAMIPS),  [Gill et al. 2016](https://academic.oup.com/cid/article/63/suppl_4/S154/2526406?login=true#google_vignette)
.
- **Data Files**: We provide a folder for each RSV-bacterial interaction labeled by their corresponding acronym (SP,MC,SA,HI). Within each folder you will find a folder named "GAM" and another folder named "Matlab files".  
- Within the "Matlab files" you will find the necessary data files containing the information of mothers and infants who acquired RSV and/or the corresponding bacteria and the Matlab scripts to generate the figures in the manuscript.
- Within the "Gam" folder, you will find the necessary data files and R scripts to generate the generalized additive model plots shown in Figures 5 and 6 in the main manuscript. 

## Installation
### Prerequisites
- Matlab (version 2021a or later)
- R (version 4.3.2)
- RStudio (optional but recommended)
- Required R packages:
  - `mgcv`: For GAM modeling.
  - `ggeffects`: For visualizing predictions.
  - `ggplot2`: For plotting.
  - `dplyr`: For data manipulation.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions, suggestions, please contact:

- Name: Ivan Ramirez-Zuniga
- Email: iramirezzuniga@uttyler.edu
- Affiliation: University of Texas at Tyler, Department of Mathematics


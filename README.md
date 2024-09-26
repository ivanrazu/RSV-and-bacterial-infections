# RSV and Bacterial Interactions

This repository contains the code and data analysis for the project investigating the interactions between Respiratory Syncytial Virus (RSV) and bacterial infections. The project specifically focuses on the effects sequential RSV-Bacterial (B) interactions with the bacterial species:
- **Streptococcus pneumoniae (SP)**
- **Moraxella catarrhalis (MC)**
- **Staphylococcus aureus (SA)**
- **Haemophilus influenzae (HI)**

## Table of Contents
1. [Project Overview](#project-overview)
2. [Directory Structure](#directory-structure)
3. [Data Description](#data-description)
5. [Usage](#usage)
6. [Scripts and Analysis](#scripts-and-analysis)
9. [License](#license)
10. [Contact](#contact)

## Project Overview
This project explores the longitudinal interactions between RSV and four bacterial pathogens (SP, MC, SA, and HI) using a birth cohort data set from mother-infant pairs from Lusaka, Zambia. We examined various scenarios, including:
- **RSV First, Bacteria Later (RSV->B)**: Infants who acquired RSV first and later the bacterial infection.
- **Bacteria First, RSV Later (B->RSV)**: Infants who acquired the bacterial infection first and later RSV.
- **Simultaneous Detection (RSV & B)**: Infants in whom both RSV and the bacterial pathogen were detected at the same time.
- **Virus-Only Cases(RSV)**: Infants who only acquired RSV without any bacterial infection.
- **Bacteria-Only Cases(B)**: Infants who only acquired one or more bacterial infections without RSV.

We use generalized additive models (GAM) to analyze and visualize the temporal relationship between RSV infection and subsequent bacterial carriage for each pathogen separately under these various scenarios.

## Directory Structure



## Data Description
- **Data Source**: The data analyzed in this project was collected during our previous project Southern Africa Mother Infant Pertusis Study (SAMIPS).
- **Data Files**: We provide a folder for each RSV-bacterial interaction labeled by their corresponding acronym (SP,MC,SA,HI). Within each folder you will find a folder named "GAM" and another folder named "Matlab files".  
- Within the "Matlab files" you will find the necessary data files containing the information of mothers and infants who acquired RSV and/or the corresponding bacteria and the Matlab scripts to generate the figures in the manuscript. 

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

Name: [Ivan Ramirez-Zuniga]
Email: [iramirezzuniga@uttyler.edu]
Affiliation: [University of Texas at Tyler, Department of Mathematics]


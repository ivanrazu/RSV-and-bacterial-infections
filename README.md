# RSV and Bacterial Interactions

This repository contains the code and data analysis for the project investigating the interactions between Respiratory Syncytial Virus (RSV) and bacterial infections. The project specifically focuses on the effects of prior RSV infection on the probability of subsequent carriage of four bacterial species:
- **Streptococcus pneumoniae (SP)**
- **Moraxella catarrhalis (MC)**
- **Staphylococcus aureus (SA)**
- **Haemophilus influenzae (HI)**

## Table of Contents
1. [Project Overview](#project-overview)
2. [Directory Structure](#directory-structure)
3. [Data Description](#data-description)
4. [Installation](#installation)
5. [Usage](#usage)
6. [Scripts and Analysis](#scripts-and-analysis)
7. [Results](#results)
8. [Contributing](#contributing)
9. [License](#license)
10. [Contact](#contact)

## Project Overview
This project explores the longitudinal interactions between RSV and four bacterial pathogens (SP, MC, SA, and HI) using datasets from different cohorts:
- **RSV-SP Infants**: Infants who acquired RSV first and later SP.
- **RSV-MC Infants**: Infants who acquired RSV first and later MC.
- **RSV-SA Infants**: Infants who acquired RSV first and later SA.
- **RSV-HI Infants**: Infants who acquired RSV first and later HI.
- **Control Infants**: Infants who acquired bacterial infections without prior RSV.

We use generalized additive models (GAM) to analyze and visualize the relationship between RSV infection and subsequent bacterial carriage for each pathogen separately.

## Directory Structure

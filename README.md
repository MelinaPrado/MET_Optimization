# MET_Optimization

This repository contains the project related to the paper "Optimizing multi-environment trials in the Southern US Rice belt via smart-climate-soil prediction based-models and economic importance" by Melina Prado, Adam Famoso, Kurt Guidry, and Roberto Fritsche-Neto (2024).

## Project Structure

- **RiceBelt_MET_OPT**: The main R project folder.
- **Scripts**: Contains all the R scripts used in the analysis.
- **Output**: Directory where all the generated output files will be saved.
- **Data**: Contains all the input data files required for the analysis.

## Scripts Description

### LSU_MET.R
This script performs a single trial analysis using linear mixed models to calculate grain yield BLUEs (Best Linear Unbiased Estimators) for each individual trial (location x year). It also identifies the environmental covariates that most significantly impact rice yield in the Southern U.S. Rice Belt and utilizes these covariates to cluster the LSU experimental sites into five major groups.

#### Key Sections:
- **First step**: Fitting a model for each trial
- **Second step**: EnvRtype
- **Third step**: SoilType
- **Fourth step**: Quality Control and feature selection
- **Fifth step**: Clustering and Plotting
- **Sixth step**: Map plotting
- **Seventh step**: PCA Biplot

### USA_TPE.R
This script delineates and characterizes the target population of environments (TPE) for the USA Rice Belt using predicted environmental covariates and economic data.

#### Key Sections:
- **First step**: EnvRtype
- **Second step**: SoilType
- **Third step**: Creating South and LA dataset
- **Fourth step**: Clustering and Plotting
- **Fifth step**: Map plotting
- **Sixth step**: PCA Biplot
- **Seventh step**: Comparing to LSU trials

### OPT_LSU_MET.R
This script optimizes the allocation of multi-environment trials (MET) by reducing the number of trials while maintaining accuracy using supervised learning. Additionally, this script compares the effect of the Enviromic-based kernel and Yield-based GxE matrix in the characterization and delineation of the TPEs.

#### Key Sections:
- **First step**: Cluster effect
- **Second step**: MET Model
- **Third step**: MET_EC Model
- **Fourth step**: WC_MET Model
- **Fifth step**: OPT_MET Model
- **Sixth step**: Plotting scenarios
- **Seventh step**: Yield-based environmental matrix
- **Eighth step**: Clustering and Plotting
- **Ninth step**: Map plotting
- **Tenth step**: Yield by cluster
- **Eleventh step**: Comparing environment matrices

## Usage
1. Clone this repository to your local machine.
2. Open the `RiceBelt_MET_OPT` project in RStudio.
3. Run the scripts in the following order:
   - `LSU_MET.R`
   - `OPT_LSU_MET.R`
   - `USA_TPE.R`
4. Check the `Output` folder for the generated results.

## Contact
For any questions or further information, please contact me at melinaprado@usp.br, or the corresponding author at rfneto@agcenter.lsu.edu,

Melina Prado.



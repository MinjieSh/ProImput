# ProImput
This is the repository of our article: Comparative assessment and outlook on methods for imputing proteomics data.

## Installation
1. Clone or download this GitHub repository
2. Set up your working directory in `Dependencies.R`
```
setwd("YOUR_PATH/ProImput")
```
3. Install all the packages required in `Dependencies.R` and then source `Dependencies.R`
4. Then play with the toy examples in the two R markdown files! 
```
Assessment Workflow for CAM (Setting 2).Rmd
Assessment Workflow for FRMF (Setting 2).Rmd
```
## Overview
In our paper, We first report an assessment of eight representative methods collectively targeting three typical missing mechanisms. The selected methods are compared on both realistic simulation and real proteomics datasets, and the performance is evaluated using three quantitative measures. 
![ProImput Workflow](images/ProImput%20Workflow.jpg)

To explore a more integrative strategy for improving imputation performance, we discuss a low-rank matrix factorization framework with fused regularization on sparsity and similarity – Fused Regularization Matrix factorization (FRMF), which can naturally integrate other-omics data such as gene expression or clinical variables. 

We also introduce a biologically-inspired latent variable modeling strategy - Convex analysis of Mixtures (CAM), which performs data imputation using the original intensity data (before log-transformation). 
![CAM Schematic](images/CAM%20Schematic.jpg)

The preliminary results on real proteomics data are provided together with an outlook into future development directions. 

############################################ R version 3.6.3

# if (!requireNamespace("BiocManager", quietly=TRUE))
# install.packages("BiocManager")
# BiocManager::install("pcaMethods")
# BiocManager::install("impute") 
# 
# install.packages("devtools")
# library(devtools)
# install_github("jeffwong/imputation")
# install_github("linxihui/NNLM")

############################################

# Assessment
library(directlabels)
library(ggplot2) 
library(devEMF)
source('Assessment/visualization.R')
source('Assessment/visualization_test.R')
source('Assessment/produce_simulated_MV.R')
source('Assessment/produce_simulated_MV_test.R')
source('Assessment/metrics.R')
source('Assessment/metrics_test.R')
source('Assessment/imputation_functions_wrappers.R')
source('Assessment/imputation_functions_wrappers_test.R')

# FRMF wrappers
library(lsa)
library(tilting)
source('FRMF/FRMF.R')
source('FRMF/FRMF_test.R')

# CAM (imputation) wrappers
library(nnls)
library(pcaMethods) 
library(imputation)
source('CAM_wrappers/CAM_cmplt.R')
source('CAM_wrappers/CAM_NIPALS.R')
source('CAM_wrappers/CAM_SVT.R')
source('CAM_wrappers/CAM_smpClus.R')
source('CAM_wrappers/CAM_helper.R')
source('CAM_wrappers/CAM_test.R')

# CAM source code
library(debCAM)
library(limSolve)
source("CAM/qpSum.R")
source("CAM/sbfsnnls.R")
source("CAM/sffsnnls.R")
source("CAM/sbfsqpSum.R")
source("CAM/sffsqpSum.R")
source("CAM/CAMPrep3.R")
source('CAM/CAMPrep3_impu.R')
source("CAM/CAMASestT.R")
source('CAM/CAMASest_impu.R')
source("CAM/reClustering.R")
source("CAM/CAMMGClusterG.R")

set.seed(1)

############################################ END

# save.image() 

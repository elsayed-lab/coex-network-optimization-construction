#######################################################################################
# list of R packages used in this analysis; this file exists to help renv detect all 
# needed packages without having to include them explicitly in coex_network_param_opt.R
#
# Packages with github sources that can't be picked up renv::hydrate():
#
# - slowkow/tftargets
# - tephenturner/annotables
# - elsayed-lab/hpgltools
# - elsayed-lab/Leishmania.major.Friedlin
# - elsayed-lab/Trypanosoma.cruzi.CLBrener.Esmeraldo
#
#######################################################################################

library(coop)
library(knitr)
library(matrixStats)
library(plyr)
#library(parallelDist)
library(doParallel)
library(printr)
library(reshape2)
library(RCurl)
library(flashClust)
library(gridExtra)
library(foreach)


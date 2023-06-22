if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scRNAseq")
library(devtools)
devtools::install_local("F:/bioinformation/SingleR-master(1).zip")
library(SingleR)

library(Seurat)
library(dplyr)
new_counts <- read.table(file="C:/Users/zhengchenlin/Desktop/GSM2884060_GW9_PFC1.UMI_TPM_no_ERCC.txt/GSM2884060_GW9_PFC1.UMI_TPM_no_ERCC.txt")
head(new_counts)
mydata <- CreateSeuratObject(counts = new_counts, min.cells = 3, project = "mydata_scRNAseq")
saveRDS(mydata, file = "C:/Users/zhengchenlin/Desktop/GSM2884060_GW9_PFC1.UMI_TPM_no_ERCC.txt/OPC20211202.rds")


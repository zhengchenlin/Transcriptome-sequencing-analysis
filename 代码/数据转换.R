#加载R包
library(Seurat)
library(dplyr)
library(patchwork)
library(mindr)
library(Matrix)
#读取文件

rawdata = "F:/bioinformation/data/皮层/GSE104276_RAW"

file_list = list.files(rawdata)
for (i in file_list) {
  assign(strsplit(i, "_", fixed = TRUE)[[1]][1], read.table(paste(rawdata,'/',i,sep = ''),comment.char = "#",header=T))
}






################

GSM4546857<-read.table("C:/Users/zhengchenlin/Downloads/GSM2884060_GW9_PFC1.UMI_TPM_no_ERCC.txt.gz",comment.char = "#",header=T)
dim(GSM4546857)
#GSM4546857[1:6,1:4]
rownames(GSM4546857)<-GSM4546857$Gene
#删除第一列
GSM4546857<-GSM4546857[,-1]
#转置
#GSM4546857<-t(GSM4546857)
object.size(GSM4546857)#2331231144 bytes
GSM4546857_sparse<-as(as.matrix(GSM4546857),"dgCMatrix")

GSM4546857_sparse[1:4,1:4]
object.size(GSM4546857_sparse)#166367952 bytes
save(GSM4546857_sparse,file = "GSM4546857_sparse.rds")
dim(GSM4546857_sparse)



#下面这段代码中，最重要的就是创建Seurat对象以及去除线粒体基因，其他都是对Seurat对象的可视化，其目的在于提高初学者对该对象的了解
## =============== 创建Seurat对象
OPC <- CreateSeuratObject(counts = GSM4546857_sparse, min.cells = 3)
OPC

## =============== 去除线粒体基因
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#此次数据检测到大量线粒体基因
grep(pattern = "^MT//.",rownames(OPC),value = T)
OPC[["percent.mt"]] <- PercentageFeatureSet(OPC, pattern = "^MT//.")
head(OPC@meta.data,5)
summary(OPC@meta.data$nCount_RNA)



seurat_object = CreateSeuratObject(counts = GSM4546857_sparse$`Gene Expression`)
seurat_object[['Protein']] = CreateAssayObject(counts = GSM4546857_sparse$`Antibody Capture`)

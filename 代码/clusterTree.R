#ѡȡPCΪ20
OPC <- FindNeighbors(OPC, dims = 1:16)
#???ò?ͬ??resolution
OPC <- FindClusters(OPC, resolution = c(0.01, 0.02, 0.03, 0.04,0.05, 0.06,0.07, 0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
OPC <- FindClusters(OPC, resolution = c(0.01, 0.02, 0.03, 0.04,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
#??ȡ????????
a=OPC@meta.data
b<-a[,c(1:17,39:43)]
##?鿴??????ͬ??resolution??ʱ?򣬾???????????
OPC@meta.data %>% View()
#install.packages('clustree')
library(ggraph)
library(clustree)
clustree(OPC, prefix = "RNA_snn_res.")
graph2ppt(file="10.clustree.pptx", width=9, aspectr=1.5)
DimPlot(OPC, reduction = "umap")

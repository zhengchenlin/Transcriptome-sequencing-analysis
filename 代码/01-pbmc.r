rm(list=ls())

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(MySeuratWrappers)

## =============1.Load the OPC dataset
data_dir <- "data/OPC3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19"
OPC.data <- Read10X(data.dir = data_dir)
OPC.data

# Initialize the Seurat object with the raw (non-normalized data)
# min.cell每个feature至少在多少个细胞中表达
# min.features每个细胞中至少有多少个feature被检测到
OPC <- CreateSeuratObject(counts = OPC.data, 
                          project = "OPC3k", 
                          min.cells = 3,
                          min.features = 200)
OPC



## =============3.QC
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
grep(pattern = "^MT-", rownames(OPC), value = T)
OPC[["percent.mt"]] <- PercentageFeatureSet(OPC, pattern = "^MT-")

# QC指标
head(OPC@meta.data, 5)
summary(OPC@meta.data$nCount_RNA)

# QC指标使用小提琴图可视化
VlnPlot(OPC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 密度曲线图
data <- OPC@meta.data
p <- ggplot(data = data, aes(x=nFeature_RNA)) + geom_density()
p

# 两两指标之间的相关性
plot1 <- FeatureScatter(OPC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OPC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# 过滤
OPC <- subset(OPC, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 25)




## =============4.标准化
# LogNormalize
OPC <- NormalizeData(OPC, normalization.method = "LogNormalize", scale.factor = 10000)

# 标准化后哦的值保存在：[["RNA"]]@data
normalized.data <- OPC[["RNA"]]@data
normalized.data[1:20,1:4]
dim(normalized.data)

## =============5.坚定高变基因
# 高变基因：在一些细胞中高表达，在另一些细胞中低表达的基因
# 变异指标： mean-variance relationship
# 默认返回两千个高变基因，用于下游PCA降维分析。
OPC <- FindVariableFeatures(OPC, selection.method = "vst", nfeatures = 2000)

# 提取前10的高变基因
top100 <- head(VariableFeatures(OPC), 100)
top100

# 展示高变基因
plot1 <- VariableFeaturePlot(OPC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2




## =============6.Scaling the data
#归一化处理：每一个基因在所有细胞中的均值变为0，方差为1，对于降维来说是必须的步骤
# 归一化后的值保存在：bmc[["RNA"]]@scale.data
OPC <- ScaleData(OPC)
scale.data <- OPC[["RNA"]]@scale.data
dim(scale.data)
scale.data[1:10,1:4]


# 可以选择全部基因归一化
all.genes <- rownames(OPC)
OPC <- ScaleData(OPC, features = all.genes)




## =============7.降维
# PCA降维，默认使用前面2000个高变基因，可使用features改变用于降维的基因集
OPC <- RunPCA(OPC, features = VariableFeatures(object = OPC))

# 可视化
VizDimLoadings(OPC, dims = 1:2, reduction = "pca")
DimPlot(OPC, reduction = "pca")
DimHeatmap(OPC, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(OPC, dims = 1:15, cells = 500, balanced = TRUE)





## =============8.确定使用PC个数
# each PC essentially representing a 'metafeature'
OPC <- JackStraw(OPC, num.replicate = 100)
OPC <- ScoreJackStraw(OPC, dims = 1:20)
JackStrawPlot(OPC, dims = 1:20)
ElbowPlot(OPC)




## =============9.对细胞聚类
# 首先基于PCA空间构建一个基于欧氏距离的KNN图
OPC <- FindNeighbors(OPC, dims = 1:20)

# 聚类并最优化
# resolution参数：值越大，细胞分群数越多，
# 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells
# Optimal resolution often increases for larger datasets. 
OPC <- FindClusters(OPC, resolution = 1.2)

# 查看聚类数ID
head(Idents(OPC), 5)

# 查看每个类别多少个细胞
head(OPC@meta.data)
table(OPC@meta.data$seurat_clusters)




## =============10.将细胞在低维空间可视化UMAP/tSNE
OPC <- RunUMAP(OPC, dims = 1:20)
OPC <- RunTSNE(OPC, dims = 1:10)

# 可视化
DimPlot(OPC, reduction = "umap", label = T, label.size = 5,pt.size = 1.2)
DimPlot(OPC, reduction = "tsne", label = T, label.size = 5,pt.size = 1.2)
saveRDS(OPC, file = "data/OPC_tutorial.rds")



## =============11.差异表达分析
# 在cluster2 vs else中差异表达
cluster1.markers <- FindMarkers(OPC, ident.1 = 0, min.pct = 0.25)
head(cluster1.markers, n = 100)


# 制定两个类cluster5 from cluster 0 and 3
cluster5.markers <- FindMarkers(OPC, ident.1 = 5, ident.2 = c(0, 1,2,3,4), min.pct = 0.25)
head(cluster5.markers, n = 10)

# 所有类的差异表达基因
# only.pos：只保留上调差异表达的基因
OPC.markers <- FindAllMarkers(OPC, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

getwd()
write.csv(cluster1.markers,file = "F:/bioinformation/data/整合数据/OPC.markers.csv")
head(OPC.markers)


#  筛选
OPC.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


# 选择一些基因进行可视化，作者这里根据自己的知识背景选择的相关基因
#小提琴图
VlnPlot(OPC, features = c("LGALS1", "PCLAF","MPZ"))
VlnPlot(OPC, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
#聚类图
FeaturePlot(OPC, features = c("RBFOX1"))
#气泡图
marker <- c("ZNF488","MYRF","E2F3","MBP","MAG","PLP1","TOP2A","MKI67","MYT1L","LHFPL3","PCNA","SLC1A2","HES1","CCND1","HES4","SNTG1","NXPH1","PDGFRA","PCDH15","ZFP36L2","HMGA1")
DotPlot(OPC, features = marker)+coord_flip()
#热图
top10 <- OPC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(OPC, features = top10$gene) + NoLegend()
#多基因小提琴图堆叠
plot_markers <- c("LINC01003","TIMP3","IRX1","ATP1B2","PLEKHH2","CAMK2D","GLRA3","TMCC3","SYNDIG1","LINC01170","RASGEF1B")
cols<- c('#FFB6C1','#FFF0F5','#FF00FF','#FF69B4','#FF1493')
VlnPlot(OPC, features = plot_markers,stacked=T,pt.size=0,cols =cols)
#格子热图
gene=c("CSPG4","PDGFRA","OLIG2","FABP7","MOBP","SOX10","ENPP6","BCAS1","MYRF","MBP","PLP1","MAG","MOG","CLDN11")
gene <- intersect(gene,rownames(as.matrix(OPC@assays$RNA@scale.data)))
project_avg<-Seurat::AverageExpression(nmGABA_new,return.seurat = T)
cm<-GetAssayData(project_avg,slot = 'data')[gene,]
cm<-as.data.frame(cm)
pheatmap::pheatmap(cm,scale = 'row',cluster_rows = F,cluster_cols = F,border_color = '#3f3f44',
                   color =colorRampPalette(c("#00cbff","#fffffb","#F27F0C"))(100)
)



## =============12.细胞类型鉴定
new.cluster.ids <- c("OPC1",   # IL7R, CCR7
                     "OL",    # CD14, LYZ
                     "OPC2",  # IL7R, S100A4
                     "OPC3",             # MS4A1
                     "OPC",         # CD8A
                     "5",  # FCGR3A, MS4A7
                     "NK",            # GNLY, NKG7
                     "DC",            # FCER1A, CST3
                     "Platelet")      # PPBP

names(new.cluster.ids) <- levels(OPC)
new.cluster.ids

# 原来的细胞聚类名称
Idents(OPC)

# 更改细胞聚类的名字
OPC <- RenameIdents(OPC, new.cluster.ids)

# 可视化
cols <- c('#ff8a80','#FFEFD5','#FF8C00','#FF4500')
DimPlot(OPC, reduction = "umap", label = TRUE, pt.size = 1.2,cols=cols) + NoLegend()
DimPlot(OPC, reduction = "tsne", label = TRUE, pt.size = 1.2) + NoLegend()

#提取特定群
cd4_OPC2 = OPC[, Idents(OPC) %in% c(4,13,18 )]
saveRDS(cd4_OPC2, file = "F:/bioinformation/data/基底神经节/OPC20211202.rds")

#亚群相关性热图
av <-AverageExpression(OPC,
                       group.by = "ident",
                       assays = "RNA")
av_jdsjj=av_jdsjj[[1]]
head(av_jdsjj)
cg=names(tail(sort(apply(FFF, 1, sd)),1000))
#pheatmap绘制热图
pheatmap::pheatmap(cor(FFF[cg,],method = 'spearman')) #默认是Pearson

# 保存结果
OPC@meta.data$cell_anno <- Idents(OPC)
write.csv(OPC@meta.data,file = "data/metadata.csv")
saveRDS(OPC, file = "F:/bioinformation/data/整合数据/OPC20220318.rds")

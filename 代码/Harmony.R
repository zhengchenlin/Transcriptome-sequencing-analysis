library(Seurat)
library(cowplot)
library(harmony)
#加载数据
load('OPC_Week.RData') #加载矩阵数据
#整合数据并降维
OPC <- CreateSeuratObject(counts = cbind( GW8_sparse ,
                                          GW9_sparse ,
                                          GW10_1_1_sparse ,
                                          GW10_1_2_sparse ,
                                          GW10_2_1_sparse ,
                                          GW10_2_2_sparse ,
                                          GW10_3_sparse ,
                                          GW12_sparse ,
                                          GW13_sparse ,
                                          GW16_1_1_sparse ,
                                          GW16_1_2_sparse ,
                                          GW16_1_3_sparse ,
                                          GW16_1_4_sparse ,
                                          GW16_1_5_sparse ,
                                          GW16_1_6_sparse ,
                                          GW16_1_7_sparse ,
                                          GW16_1_8_sparse ,
                                          GW16_1_9_sparse ,
                                          GW19_1_sparse ,
                                          GW19_2_sparse ,
                                          GW19_3_sparse ,
                                          GW23_1_1_sparse ,
                                          GW23_1_2_sparse ,
                                          GW23_1_3_sparse ,
                                          GW23_2_1_sparse ,
                                          GW23_2_2_sparse ,
                                          GW23_2_3_sparse ,
                                          GW23_2_4_sparse ,
                                          GW23_2_5_sparse ,
                                          GW26_1_1_sparse ,
                                          GW26_1_2_sparse ,
                                          GW26_1_3_sparse ,
                                          GW26_1_4_sparse ,
                                          GW26_1_5_sparse ,
                                          GW26_1_6_sparse ,
                                          GW26_1_7_sparse ,
                                          GW26_1_8_sparse ,
                                          GW26_1_9_sparse ,
                                          GW26_1_10_sparse ), 
                          project = "OPC", min.cells = 5) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = OPC@var.genes, npcs = 20, verbose = FALSE)

#赋值条件变量
OPC@meta.data$Week <- c(rep("GW10", ncol(GW10_1_1_sparse)),
                           rep("GW10", ncol(GW10_1_2_sparse)),
                           rep("GW10", ncol(GW10_2_1_sparse)),
                           rep("GW10", ncol(GW10_2_2_sparse)),
                           rep("GW10", ncol(GW10_3_sparse)),
                           rep("GW12", ncol(GW12_sparse)),
                           rep("GW13", ncol(GW13_sparse)),
                           rep("GW16", ncol(GW16_1_1_sparse)),
                           rep("GW16", ncol(GW16_1_2_sparse)),
                           rep("GW16", ncol(GW16_1_3_sparse)),
                           rep("GW16", ncol(GW16_1_4_sparse)),
                           rep("GW16", ncol(GW16_1_5_sparse)),
                           rep("GW16", ncol(GW16_1_6_sparse)),
                           rep("GW16", ncol(GW16_1_7_sparse)),
                           rep("GW16", ncol(GW16_1_7_sparse)),
                           rep("GW16", ncol(GW16_1_9_sparse)),
                           rep("GW19", ncol(GW19_1_sparse)),
                           rep("GW19", ncol(GW19_2_sparse)),
                           rep("GW19", ncol(GW19_3_sparse)),
                           rep("GW23", ncol(GW23_1_1_sparse)),
                           rep("GW23", ncol(GW23_1_2_sparse)),
                           rep("GW23", ncol(GW23_1_3_sparse)),
                           rep("GW23", ncol(GW23_2_1_sparse)),
                           rep("GW23", ncol(GW23_2_2_sparse)),
                           rep("GW23", ncol(GW23_2_3_sparse)),
                           rep("GW23", ncol(GW23_2_4_sparse)),
                           rep("GW23", ncol(GW23_2_5_sparse)),
                           rep("GW26", ncol(GW26_1_1_sparse)),
                           rep("GW26", ncol(GW26_1_10_sparse)),
                           rep("GW26", ncol(GW26_1_2_sparse)),
                           rep("GW26", ncol(GW26_1_3_sparse)),
                           rep("GW26", ncol(GW26_1_4_sparse)),
                           rep("GW26", ncol(GW26_1_5_sparse)),
                           rep("GW26", ncol(GW26_1_6_sparse)),
                           rep("GW26", ncol(GW26_1_7_sparse)),
                           rep("GW26", ncol(GW26_1_8_sparse)),
                           rep("GW26", ncol(GW26_1_9_sparse)),
                           rep("GW8", ncol(GW8_sparse)),
                           rep("GW9", ncol(GW9_sparse)))

#未经校正的PC中的数据集之间存在明显差异
options(repr.plot.height = 5, repr.plot.width = 12)
DimPlot(object = OPC, reduction = "pca", pt.size = .1, group.by = "Week")
VlnPlot(object = OPC, features = "PC_1", group.by = "Week",  pt.size = .1)
plot_grid(p1,p2)

#Run Harmony
options(repr.plot.height = 2.5, repr.plot.width = 6)
OPC <- OPC %>%
RunHarmony("Week", plot_convergence = TRUE)

#查看确认数据集在Harmony运行之后的前两个维度中得到很好的整合。
options(repr.plot.height = 5, repr.plot.width = 12)
DimPlot(object = OPC, reduction = "harmony", pt.size = .1, group.by = "Week")
VlnPlot(object = OPC, features = "harmony_1", group.by = "Week",pt.size = .1)
plot_grid(p1,p2)

#下游分析
OPC <- OPC
OPC <-RunTSNE(OPC,reduction = "harmony", dims = 1:20)
OPC <-RunUMAP(OPC,reduction = "harmony", dims = 1:20)
OPC <-FindNeighbors(OPC,reduction = "harmony", dims = 1:20)
OPC <-FindClusters(OPC,resolution = 0.07) %>%
  identity()

#鉴定混合情况
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(OPC, reduction = "tsne", group.by = "Week", pt.size = 1.1)
DimPlot(OPC, reduction = "tsne", pt.size = 1.1)
DimPlot(OPC, reduction = "tsne", group.by = "orig.ident", pt.size = 1.1, split.by = 'orig.ident')
DimPlot(OPC, reduction = "umap", group.by = "Week", pt.size = 1.1)
DimPlot(OPC, reduction = "umap",  pt.size = 1.1)
DimPlot(OPC, reduction = "umap", group.by = "Week", pt.size = 1.1, split.by = 'Week')
DimPlot(OPC, reduction = "umap", group.by = "orig.ident", pt.size = 1.1)
DimPlot(OPC, reduction = "umap", group.by = "orig.ident", pt.size = 1.1, split.by = 'orig.ident')

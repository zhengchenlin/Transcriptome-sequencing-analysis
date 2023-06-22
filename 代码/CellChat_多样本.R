library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

#提取细胞群
#GW9
GW9<-subset(OPC,cells = row.names(subset(OPC@meta.data,
                                         OPC@meta.data$Week=="BG_GW09"|
                                         OPC@meta.data$Week=="BS_GW09"|
                                         OPC@meta.data$Week=="SC_GW09"
)))
#GW9
GW12<-subset(OPC,cells = row.names(subset(OPC@meta.data,
                                         OPC@meta.data$Week=="BG_GW12"|
                                         OPC@meta.data$Week=="BS_GW12"|
                                         OPC@meta.data$Week=="SC_GW12"
)))

#########################分别创建CellChat对象#########################
#########################GW9
#提取细胞的基因表达数据
data.input  <- GW9@assays$RNA@data
# create a dataframe consisting of the cell labels
identity = data.frame(group =GW9@active.ident,row.names = names(GW9@active.ident)) 
unique(identity$group) # check the cell labels

#创建cellchat对象
cellchat <- createCellChat(data.input)
cellchat

#将metadata信息加到CellChat对象中
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

#分析样本cco.pbmc的细胞通讯网络
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
cellchat_GW9 <- cellchat
cellchat_GW9
saveRDS(cellchat_GW9, "cellchat_GW9.rds")



#########################GW12
#提取细胞的基因表达数据
data.input  <- GW12@assays$RNA@data
# create a dataframe consisting of the cell labels
identity = data.frame(group =GW12@active.ident,row.names = names(GW12@active.ident)) 
unique(identity$group) # check the cell labels

#创建cellchat对象
cellchat <- createCellChat(data.input)
cellchat

#将metadata信息加到CellChat对象中
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

#分析样本cco.pbmc的细胞通讯网络
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
cellchat_GW12 <- cellchat
cellchat_GW12
saveRDS(cellchat_GW12, "cellchat_GW12.rds")


#########################合并#########################
object.list <- list(GW9 = cellchat_GW9, GW12 = cellchat_GW12)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

#比较交互总数和交互强度
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#不同细胞群之间的相互作用数量或强度的差异
#圆图
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
#热图
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
#细胞之间交互数量
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#根据信号组的功能相似性识别信号组
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
#计算和可视化通路距离
rankSimilarity(cellchat, type = "functional")
#比较每个信号通路的整体信息流
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

#比较与每个细胞群相关的传出（或传入）信号
library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#识别上调和下调的信号配体对
netVisual_bubble(cellchat, sources.use = 4,  comparison = c(1, 2), angle.x = 45)

#使用圆图可视比较细胞-细胞通信
pathways.show <- c("MK") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

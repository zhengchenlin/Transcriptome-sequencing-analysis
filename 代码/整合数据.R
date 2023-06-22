library(Seurat)
library(harmony)
library(ggplot2)
#删除多余数据
JDSJJ@meta.data<-JDSJJ@meta.data[,c(-11:-28)]
NG@meta.data<-NG@meta.data[,c(-11:-32)]
JS@meta.data<-JS@meta.data[,c(-8:-26)]
#更改列名
colnames(JS@meta.data)[4] <- "Week"
#添加部位信息
JDSJJ@meta.data[,1 ] = "subpallium"
NG@meta.data[,1 ] = "brainstem"
JS@meta.data[,1 ] = "spinal_cord"
#更改时期名
JDSJJ@meta.data$Week = c('GW09'="SP_GW09",'GW10'="SP_GW10",'GW11'="SP_GW11",'GW12'="SP_GW12")[ as.character(JDSJJ@meta.data$Week)]
NG@meta.data$Week = c('GW09'="BS_GW09",'GW10'="BS_GW10",'GW11'="BS_GW11",'GW12'="BS_GW12")[ as.character(NG@meta.data$Week)]
JS@meta.data$Week = c('SC8'="SC_GW08",'SC9'="SC_GW09",'SC10'="SC_GW10",'SC12'="SC_GW12")[ as.character(JS@meta.data$Week)]
#整合降维
OPC <- merge(JDSJJ, y = c(NG,JS), add.cell.ids = c("basal ganglion", "brainstem","spinal cord"), project = "OPC",merge.data = TRUE)
red.genes <- c("HBA1","HBA2","HBB",'HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ') 
sex.genes<-c('DDX3Y','EIF2S3Y','UTY','KDM5D','XIST','TSIX')
mt.genes <- grep(pattern = "^MT-", rownames(OPC), value = T)
OPC<-subset(OPC,features = setdiff(row.names(OPC),c(red.genes,mt.genes,sex.genes)))
OPC <- NormalizeData(OPC, normalization.method = "LogNormalize", scale.factor = 10000)
OPC <- FindVariableFeatures(OPC, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(OPC)
OPC <- ScaleData(OPC, features = all.genes)
OPC <- RunPCA(OPC, features = VariableFeatures(object = OPC))
OPC <- OPC %>%
  RunHarmony("Week", plot_convergence = TRUE)
#下游分析
OPC <-RunUMAP(OPC,reduction = "harmony", dims = 1:20)
OPC <-FindNeighbors(OPC,reduction = "harmony", dims = 1:20)
OPC <-FindClusters(OPC,resolution = 0.1) %>%
  identity()
marker <- c("PDGFRA","ZFP36L2","FABP7","PLP1","CLDN11","MBP","COL20A1","MKI67","TOP2A","HMGA1","HEY1","MYT1L")
DotPlot(OPC, features = marker)+coord_flip()
new.cluster.ids <- c("OPC1",
                     "OL",
                     "OPC2",
                     "OPC3")
names(new.cluster.ids) <- levels(OPC)
OPC <- RenameIdents(OPC, new.cluster.ids)
cols <- c('#ff8a80','#FFEFD5','#FF8C00','#FF4500')
DimPlot(OPC, reduction = "umap", label = F, pt.size = 1.2,cols=cols)
DimPlot(OPC, reduction = "umap", group.by = "Week", pt.size = 1.1)
DimPlot(OPC, reduction = "umap", group.by = "orig.ident", pt.size = 1.1)
JDSJJ  = OPC[, OPC@meta.data$orig.ident=="subpallium"]
NG  = OPC[, OPC@meta.data$orig.ident=="brainstem"]
JS  = OPC[, OPC@meta.data$orig.ident=="spinal_cord"]
table(JDSJJ@meta.data$seurat_clusters)
table(NG@meta.data$seurat_clusters)
table(JS@meta.data$seurat_clusters)



#################monocle#################
library(monocle)
expr_matrix <- as(as.matrix(OPC@assays$RNA@counts), 'sparseMatrix')
p_data <- OPC@meta.data 
p_data$celltype <- OPC@active.ident 
f_data <- data.frame(gene_short_name = row.names(OPC),row.names = row.names(OPC))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
deg.cluster <- FindAllMarkers(OPC)
express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds, express_genes)
diff<- differentialGeneTest(cds[ express_genes, ], fullModelFormulaStr="~ celltype")
deg <- subset(diff, qval < 0.01)
deg <- deg[order(deg$qval,decreasing=F),]
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
ordering_genes <- row.names(subset(diff, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method='DDRTree')
cds <- orderCells(cds)
cds <- orderCells(cds,root_state = 4)
plot_complex_cell_trajectory(cds, color_by="celltype")+
  scale_color_manual( values=cols)





#按部位查找差异基因
nmGABA_new<-OPC
Idents(nmGABA_new)<-nmGABA_new$orig.ident
State_idents<-c('basal_ganglion','brainstem','spinal_cord')
names(State_idents)<-levels(nmGABA_new)
nmGABA_new<-RenameIdents(nmGABA_new,State_idents)
OPC.markers <- FindAllMarkers(nmGABA_new, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(OPC.markers,"place.csv")
jdsjj_marker <- subset(OPC.markers,OPC.markers$cluster=="basal_ganglion")[,7]
ng_marker <- subset(OPC.markers,OPC.markers$cluster=="brainstem")[,7]
js_marker <- subset(OPC.markers,OPC.markers$cluster=="spinal_cord")[,7]
jdsjj_ng <- intersect(ng_marker,jdsjj_marker)
length(intersect(ng_marker,js_marker)
length(intersect(js_marker,jdsjj_marker)
length(intersect(intersect(ng_marker,jdsjj_marker),js_marker))
table(OPC.markers$cluster)

library(venneuler)
vd <- venneuler(c(A=606, B=546, C=584, "A&B"=214, "B&C"=0 ,"A&C"=1 ,"A&B&C"=0))
plot(vd)


# 载入包dian
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)

jdsjj <- setdiff(setdiff(jdsjj_marker,ng_marker),js_marker)
ng <- setdiff(setdiff(ng_marker,jdsjj_marker),js_marker)
js <- setdiff(setdiff(js_marker,ng_marker),jdsjj_marker)
#####################基底神经节
eg = bitr(jdsjj, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
erich.go.BP = enrichGO(gene = eg$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.BP)
#####################脑干
eg = bitr(ng, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

erich.go.BP = enrichGO(gene = eg$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.BP)
#####################脊髓
eg = bitr(js, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
erich.go.BP = enrichGO(gene = eg$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.BP)
#####################基底神经节和脑干
eg = bitr(jdsjj_ng, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

erich.go.BP = enrichGO(gene = eg$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.BP)



library(dorothea)
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
TF <- intersect(jdsjj_ng,dorothea_regulon_human)
DotPlot(nmGABA_new, features = TF)+coord_flip()
project<-subset(OPC,features = ng)

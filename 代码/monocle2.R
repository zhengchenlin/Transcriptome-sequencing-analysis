##################################### 1 #####################################
library(monocle)
library(Seurat)


##################################### 2 #####################################
#创建CellDataSet
##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
expr_matrix <- as(as.matrix(OPC@assays$RNA@counts), 'sparseMatrix')
##提取表型信息到p_data(phenotype_data)里面 
p_data <- OPC@meta.data 
p_data$celltype <- OPC@active.ident 
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(OPC),row.names = row.names(OPC))
#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())


##################################### 3 #####################################
#估计size factor和离散度
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


##################################### 4 #####################################
#因为Seurat已经完成细胞过滤，此步可省略
cds <- detectGenes(cds, min_expr=0.1)
#print(head(fData(cds)))


##################################### 5 #####################################
#轨迹定义基因选择及可视化和构建轨迹
##使用seurat选择的高变基因
express_genes <- VariableFeatures(OPC)
cds <- setOrderingFilter(cds, express_genes)
#plot_ordering_genes(cds)
length(express_genes)
#print(head(pData(cds)))
##使用clusters差异表达基因
deg.cluster <- FindAllMarkers(OPC)
express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds, express_genes)
#plot_ordering_genes(cds)
##使用monocle选择的高变基因
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
#plot_ordering_genes(cds)

#这一步输入的expressed_genes来自于步骤4。
#后续分析使用的是该方法
#也可输入seurat筛选出的高变基因：expressed_genes <- VariableFeatures(pbmc) 
diff<- differentialGeneTest(cds[ express_genes, ], fullModelFormulaStr="~ celltype")




deg <- subset(diff, qval < 0.01) #选出2829个基因
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)
##差异基因的结果文件保存
#write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
## 轨迹构建基因可视化
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  




ordering_genes <- row.names(subset(diff, qval < 0.01))
length(ordering_genes)
cds <- setOrderingFilter(cds, ordering_genes)

#降维
cds <- reduceDimension(cds, max_components = 2, method='DDRTree')
# 拟时间轴轨迹构建和在拟时间内排列细胞
cds <- orderCells(cds)
##使用root_state参数可以设置拟时间轴的根，此处以state2为根
cds <- orderCells(cds,root_state = 2)


##################################### 6 #####################################
#可视化
#以细胞类型上色的拟时轨迹
plot_cell_trajectory(cds, color_by="celltype")+
  scale_color_manual( values=cols)
#树状图拟时轨迹
plot_complex_cell_trajectory(cds, color_by="celltype")+
  scale_color_manual( values=cols)
#以细胞状态上色（拆分）“分面"轨迹图
plot_cell_trajectory(cds, color_by = "celltype") + facet_wrap("~celltype", nrow = 1)+
  scale_color_manual( values=cols)
#细胞拟时轨迹
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)
#以细胞state上色
cell_colors <- c('lightblue','red','yellow','purple','green','blue')
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
#沿时间轴的细胞密度图
library(ggpubr)
df <- pData(cds) 
## pData(cds)取出的是cds对象中cds@phenoData@data的内容
View(df)
ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
#手动设置颜色注意
ClusterName_color_panel <- c("OPC1" = "#DC143C", 
                             "OL" = "#0000FF", 
                             "NPC" = "#20B2AA",
                             "OPC2" = "#FFA500",
                             "CD8 T" = "#9370DB"
                             )
ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+
  theme_classic2()+ 
  scale_fill_manual(name = "", values = ClusterName_color_panel)+
  scale_color_manual(name = "", values = ClusterName_color_panel)


##################################### 7 #####################################
#提取感兴趣的细胞（进行后续分析）
pdata <- Biobase::pData(cds)
ordergene <- rownames(ordering_genes)
s.cells <- subset(pdata, State=="1") %>% rownames()
save(s.cells, file = "Monocle_state7.rda")
#保存结果
write.csv(pData(cds), "pseudotime.csv")
save(cds, file = "cds20220325.rda")

#指定基因的可视化
##选择前4个top基因并将其对象取出
keygenes <- head(ordering_genes,4)
cds_subset <- cds[keygenes,]
##可视化：以state/celltype/pseudotime进行
plot_genes_in_pseudotime(cds_subset, color_by = "State",cell_size = 1.5)
plot_genes_in_pseudotime(cds_subset, color_by = "celltype")
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
#ggsave("Genes_pseudotimeplot.pdf", plot = plotc, width = 16, height = 8)

#指定基因
s.genes <- c("KIF20A","NPY","CDK1")
plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State",cell_size = 1.5)
plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
plot_genes_in_pseudotime(cds[s.genes,], color_by = "State",cell_size = 1.5)

#拟时序展示单个基因表达量
colnames(pData(cds))
pData(cds)$MBP = log2( exprs(cds)['MBP',]+1)
p1=plot_cell_trajectory(cds, color_by = "MBP")
pData(cds)$TF = log2(exprs(cds)['TF',]+1)
p2=plot_cell_trajectory(cds, color_by = "TF")
library(patchwork)
p1+p2


##################################### 8 #####################################
#寻找拟时相关的基因（拟时差异基因）
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(cds, cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
#write.csv(Time_diff, "1.csv", row.names = F)
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=5, show_rownames=T, return_heatmap=F)
#ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

#提取前面热图中每个cluster的基因
p$tree_row
clusters <- cutree(p$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.csv(clustering, "Time_clustering_all1.csv", row.names = T)

#拟时差异基因热图绘制（提取了前100个）
Time_genes <- top_n(Time_diff, n = 200, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=5, show_rownames=T, return_heatmap=T)
#ggsave("Time_heatmapTop100.pdf", p, width = 5, height = 10)

#显著差异基因按热图结果排序并保存
hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(Time_diff_sig, "Time_diff_sig.csv", row.names = F)

#也可手动选择基因来绘制热图，查看其表达模式
marker_genes <- row.names(subset(fData(cds),gene_short_name %in% c("PLP1","MOG","MBP","QKI","PCDH9","NRXN3","GRIN1","GRIN2A","CNTNAP2","SYT1","TENM2","RBFOX1","MEF2C","DPP10","NRG3","GRIK2","GRIN2B","GRIK1","ITPR2","SOX6","VCAN","CSPG4","PDGFRA","DSCAM","NEU4")))
diff_test_res <- differentialGeneTest(cds[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(cds[sig_gene_names,],
                        num_clusters = 6,
                        cores = 1,
                        show_rownames = T)



##################################### 9 #####################################
#分支点分析
BEAM_res <- BEAM(cds,
                 fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch",
                 reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
                 branch_states = NULL,
                 branch_point = 1,
                 relative_expr = TRUE,
                 branch_labels = NULL,
                 verbose = FALSE,
                 cores = 1,)
#全部差异基因
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
q <- plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T,
                            return_heatmap = T)
#热图展示TOP100基因
BEAM_genes <- top_n(BEAM_res, n = 20, desc(qval)) %>% pull(gene_short_name) %>% as.character()
q <- plot_genes_branched_heatmap(cds[BEAM_genes,],  branch_point = 1, 
                                 num_clusters = 4, show_rownames = T, return_heatmap = T)
ggsave("BEAM_heatmap.pdf", p$ph_res, width = 6.5, height = 10)
#提取前面热图中每个cluster的基因
q$ph_res$tree_row
clusters <- cutree(q$ph_res$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.csv(clustering, "Time_clustering_all1.csv", row.names = T)

#对感兴趣的基因可视化
genes <- row.names(subset(fData(cds),
                          gene_short_name %in% c("HES5","SOX6","SOX8","SOX10","TCF7L2","E2F3","SOX4")))

plot_genes_branched_pseudotime(cds[genes,],
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1,
                               cell_size = 1.5)

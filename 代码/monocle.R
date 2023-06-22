if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scRNAseq")
library(devtools)
devtools::install_local("F:/bioinformation/SingleR-master(1).zip")
library(SingleR)

library(monocle)
library(Seurat)

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
cds <- detectGenes(cds, min_expr=0.1)
print(head(fData(cds)))
##ʹ??seuratѡ???ĸ߱????????
express_genes <- VariableFeatures(OPC)
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)
length(express_genes)
print(head(pData(cds)))
##ʹ??clusters????????????
deg.cluster <- FindAllMarkers(OPC)
express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene

cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)
##ʹ??monocleѡ???ĸ߱????????
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)
#??????
diff<- differentialGeneTest(cds[ express_genes, ], fullModelFormulaStr="~ celltype")
ordering_genes <- row.names(subset(diff, qval < 0.01))
length(ordering_genes)
cds <- setOrderingFilter(cds, ordering_genes)
#??ά
cds <- reduceDimension(cds, max_components = 2, method='DDRTree')
cds <- orderCells(cds)
cds <- orderCells(cds,root_state = 2)
plot_cell_trajectory(cds, color_by="celltype")
plot_complex_cell_trajectory(cds, color_by="celltype")

plot_cell_trajectory(cds, color_by = "celltype") + facet_wrap("~celltype", nrow = 1)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)
#????ɫ
cell_colors <- c('lightblue','red','yellow','purple','green','blue')
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)

#??????State7??ϸ??????Ȥ
pdata <- Biobase::pData(cds)
ordergene <- rownames(ordering_genes)
s.cells <- subset(pdata, State=="1") %>% rownames()
save(s.cells, file = "Monocle_state7.rda")
##ѡ??ǰ4??top???򲢽???????ȡ??
keygenes <- head(ordering_genes,4)
cds_subset <- cds[keygenes,]
##???ӻ?????state/celltype/pseudotime????
plot_genes_in_pseudotime(cds_subset, color_by = "State",cell_size = 1.5)
plot_genes_in_pseudotime(cds_subset, color_by = "celltype")
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
ggsave("Genes_pseudotimeplot.pdf", plot = plotc, width = 16, height = 8)

#ָ??????
s.genes <- c("APLP1","PLP1","MBP","BCAS1")
plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State",cell_size = 1.5)
plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
plot_genes_in_pseudotime(cds[s.genes,], color_by = "State",cell_size = 1.5)

#??֧??????
BEAM_res <- BEAM(cds,
                 fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch",
                 reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
                 branch_states = NULL,
                 branch_point = 2,
                 relative_expr = TRUE,
                 branch_labels = NULL,
                 verbose = FALSE,
                 cores = 1,)
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2) #??????
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                 qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)#??632??gene??̫????


#ѡǰ100?????????ӻ?
BEAM_genes <- top_n(BEAM_res, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(cds[BEAM_genes,],  branch_point = 2, 
                                 num_clusters = 3, show_rownames = T, return_heatmap = T)
ggsave("BEAM_heatmap.pdf", p$ph_res, width = 6.5, height = 10)

genes <- row.names(subset(fData(cds),
                          gene_short_name %in% c( "FOS", "LRRC4C", "NOVA1","PNISR","RPS24")))

plot_genes_branched_pseudotime(cds[genes,],
                               branch_point = 2,
                               color_by = "State",
                               ncol = 1,
                               cell_size = 1.5)

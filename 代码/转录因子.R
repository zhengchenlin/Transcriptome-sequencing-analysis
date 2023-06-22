library(ggpubr)
#脑干
ClusterName_color_panel <- c("OPC1" = "#87cefa", 
                             "OPC2" = "#90ee90",
                             "OL1" = "#0ae1ff", 
                             "OL2" = "#1e90ff"
)
#脊髓
ClusterName_color_panel <- c("OPC1" = "#9acd32", 
                             "OPC2" = "#00FFFF", 
                             "OPC3" = "#3CB371",
                             "OL1" = "#00FF7F"
)
#基底神经节
ClusterName_color_panel <- c("OPC1" = "#FFB6C1", 
                             "OPC2" = "#D8BFD8", 
                             "OPC3" = "#FF00FF",
                             "OL1" = "#FF69B4",
                             "OL2" = "#FF1493"
)
genes <- row.names(subset(fData(cds),
                          gene_short_name %in% TF))

plot_genes_branched_pseudotime(cds[genes,],
                                   branch_point = 1,
                                   color_by = "celltype",
                                   ncol = 5,
                                   cell_size = 1.5)+
  theme_classic2()+ 
  scale_fill_manual(name = "", values = ClusterName_color_panel)+
  scale_color_manual(name = "", values = ClusterName_color_panel)




genes <- row.names(subset(fData(cds),
                          gene_short_name %in% c("TCF7L2","SOX8","ZNF704","SOX11","EEF1A1","SOX10")))

p2<-plot_genes_branched_pseudotime(cds[genes,],
                                   branch_point = 2,
                                   color_by = "State",
                                   ncol = 1,
                                   cell_size = 1.5)
genes <- row.names(subset(fData(cds),
                          gene_short_name %in% c("OLIG2","HOXB8","JUND","HES6","FOS","SOX9")))

p3<-plot_genes_branched_pseudotime(cds[genes,],
                                   branch_point = 2,
                                   color_by = "State",
                                   ncol = 1,
                                   cell_size = 1.5)
genes <- row.names(subset(fData(cds),
                          gene_short_name %in% c("TCF4","HIPK2","LMO4","JUN","NKX6-2","SOX6")))


p4<-plot_genes_branched_pseudotime(cds[genes,],
                                   branch_point = 2,
                                   color_by = "State",
                                   ncol = 1,
                                   cell_size = 1.5)
p1+p2
p3+p4
a+b
a

library(ggsci)
#拟时序展示单个基因表达量
colnames(pData(cds))
pData(cds)$TCF4 = log2(exprs(cds)['TCF4',]+1)
p <- plot_cell_trajectory(cds, color_by = "TCF4")+scale_color_gsea()
p+scale_color_gradient(low="#00e0e8",high="#ff6d00")

pData(cds)$ASCL1 = log2(exprs(cds)['ASCL1',]+1)
p2=plot_cell_trajectory(cds, color_by = "ASCL1")+scale_color_gsea()

pData(cds)$HES5 = log2( exprs(cds)['HES5',]+1)
p3=plot_cell_trajectory(cds, color_by = "HES5")+scale_color_gsea()

pData(cds)$SALL3 = log2(exprs(cds)['SALL3',]+1)
p4=plot_cell_trajectory(cds, color_by = "SALL3")+scale_color_gsea()

pData(cds)$E2F3 = log2(exprs(cds)['E2F3',]+1)
p5=plot_cell_trajectory(cds, color_by = "E2F3")+scale_color_gsea()

pData(cds)$OLIG1 = log2(exprs(cds)['OLIG1',]+1)
p6=plot_cell_trajectory(cds, color_by = "OLIG1")+scale_color_gsea()

pData(cds)$TCF7L2 = log2(exprs(cds)['TCF7L2',]+1)
p7=plot_cell_trajectory(cds, color_by = "TCF7L2")+scale_color_gsea()

pData(cds)$SOX8 = log2(exprs(cds)['SOX8',]+1)
p8=plot_cell_trajectory(cds, color_by = "SOX8")+scale_color_gsea()

pData(cds)$ZNF704 = log2(exprs(cds)['ZNF704',]+1)
p9=plot_cell_trajectory(cds, color_by = "ZNF704")+scale_color_gsea()

pData(cds)$SOX11 = log2(exprs(cds)['SOX11',]+1)
p10=plot_cell_trajectory(cds, color_by = "SOX11")+scale_color_gsea()

pData(cds)$EEF1A1 = log2(exprs(cds)['EEF1A1',]+1)
p11=plot_cell_trajectory(cds, color_by = "EEF1A1")+scale_color_gsea()

pData(cds)$SOX10 = log2(exprs(cds)['SOX10',]+1)
p12=plot_cell_trajectory(cds, color_by = "SOX10")+scale_color_gsea()

pData(cds)$OLIG2 = log2(exprs(cds)['OLIG2',]+1)
p13=plot_cell_trajectory(cds, color_by = "OLIG2")+scale_color_gsea()

pData(cds)$HOXB8 = log2(exprs(cds)['HOXB8',]+1)
p14=plot_cell_trajectory(cds, color_by = "HOXB8")+scale_color_gsea()

pData(cds)$JUND = log2(exprs(cds)['JUND',]+1)
p15=plot_cell_trajectory(cds, color_by = "JUND")+scale_color_gsea()

pData(cds)$HES6 = log2(exprs(cds)['HES6',]+1)
p16=plot_cell_trajectory(cds, color_by = "HES6")+scale_color_gsea()

pData(cds)$SOX9 = log2(exprs(cds)['SOX9',]+1)
p17=plot_cell_trajectory(cds, color_by = "SOX9")+scale_color_gsea()

pData(cds)$TCF4 = log2(exprs(cds)['TCF4',]+1)
p18=plot_cell_trajectory(cds, color_by = "TCF4")+scale_color_gsea()

pData(cds)$HIPK2 = log2(exprs(cds)['HIPK2',]+1)
p19=plot_cell_trajectory(cds, color_by = "HIPK2")+scale_color_gsea()

pData(cds)$LMO4 = log2(exprs(cds)['LMO4',]+1)
p20=plot_cell_trajectory(cds, color_by = "LMO4")+scale_color_gsea()

pData(cds)$JUN = log2(exprs(cds)['JUN',]+1)
p21=plot_cell_trajectory(cds, color_by = "JUN")+scale_color_gsea()

pData(cds)$"MBP" = log2(exprs(cds)['MBP',]+1)
p22=plot_cell_trajectory(cds, color_by = "MBP")+scale_color_gsea()

pData(cds)$SOX6 = log2(exprs(cds)['SOX6',]+1)
p23=plot_cell_trajectory(cds, color_by = "SOX6")+scale_color_gsea()

library(patchwork)
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12

p13+p14+p15+p16+p17+p18+p19+p20+p21+p23

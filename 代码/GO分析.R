
# 载入包dian
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)

#转换基因ID
eg = bitr(aa, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

erich.go.BP = enrichGO(gene = eg$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

##分析完成后，作图
#气泡图
dotplot(erich.go.BP)
#柱状图
barplot(erich.go.BP)
#热图展示特定GO
library(ComplexHeatmap)
col_fun = colorRamp2(c(1, 5), c("#ff7900", "#005676"))
Heatmap(AA, name=" ", border=F,col=col_fun,rect_gp=gpar(col="grey",lwd=2),cluster_rows =T)

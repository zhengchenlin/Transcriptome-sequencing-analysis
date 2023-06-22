library(org.Hs.eg.db)
library(stringr)
library(BiocGenerics)
library(clusterProfiler)
library(ggplot2)
require(DOSE)
library(Hmisc)
library(future)
library(future.apply)
library(clusterProfiler)
library(limma)
library(enrichplot)
library(psych)
library(RColorBrewer)
library(GSVA)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')

#
#利用相关性进行单基因GSEA####
#写循环分析KEGG通路的GSEA
## 需要进行分析的基因 "ZSCAN4", "LCOR",  "LCN2", "ADH6" 
matrix_uniq.filt.tumor <- OPC@assays$RNA@data
###
data01 <- c("SOX4")
## 提取基因的表达量
tar.exp <- matrix_uniq.filt.tumor[data01,]
## 转换数据格式
# class(tar.exp)
# "data.frame"
y <- as.numeric(tar.exp)
# class(y)
# [1] "numeric"
## 与其他基因进行相关性计算
data1 <- data.frame()
for (i in rownames(matrix_uniq.filt.tumor)){
  dd <- corr.test(as.numeric(matrix_uniq.filt.tumor[i,]), y, method = "pearson", adjust = "fdr")
  data1 <- rbind(data1, data.frame(gene =i, cor = dd$r, p.value=dd$p))  ## 筛选相关性高的基因
}
data1 <- data1[order(data1$cor,decreasing = T),]
gene1 <- data1$cor
names(gene1) <- mapIds(org.Hs.eg.db,keys = data1$gene,column = 'ENTREZID',
                       keytype = 'SYMBOL',multiVals='filter')
gene1 <- na.omit(gene1)
save(gene1,file = paste0('ZSCAN4.cor.Rdata'))
#分析
options(digits = 4)#输出结果的小数点后保留4位
kegggsea <-  gseKEGG(
  gene1,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = TRUE)
write.table(kegggsea,'ZSCAN4.KEGGgsea.txt',sep = '\t',quote = F,row.names = F)
#作图
pdf('ZSCAN4.GSEA.KEGG.pdf',width = 15,height = 10)
gseaplot2(kegggsea, geneSetID = 1:5,pvalue_table = FALSE,
          base_size = 13,  ## 字体大小
          rel_heights = c(1.5, 0.5, 0.5),color= brewer.pal(10,'Paired'))

dev.off()
#
pdf('ZSCAN4.GSEA.KEGG.dotplot.pdf', width = 8, height = 6)
dotplot(kegggsea, showCategory = 10)

dev.off()
## ---------------------------------
### GO富集分析
file <- 'ZSCAN4.cor.Rdata'
#load相关性结果
load(file)
#分析
options(digits = 4)#输出结果的小数点后保留4位
gogsea <-  gseGO(
  gene1,
  OrgDb = org.Hs.eg.db,
  ont = 'ALL',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  seed = TRUE)

write.table(gogsea,'ZSCAN4.GOgsea.txt',sep = '\t',quote = F,row.names = F)
gogsea[1:5,1:5]
#作图
pdf('ZSCAN4.GSEA-GO.pdf',width = 15,height = 10)
gseaplot2(gogsea, geneSetID = 1:10,pvalue_table = FALSE,## 不带P值
          base_size = 13,  ## 字体大小
          rel_heights = c(1.5, 0.5, 0.5))

dev.off()

pdf("ZSCAN4.GSEA.GO.dotplot.pdf", width = 10, height = 6)
dotplot(gogsea, showCategory = 10)

dev.off()


##-------------------------------------------------------
## 
### 多几基因一起分析，写一个循环体即可

for (ii in c("ZSCAN4", "LCOR",  "LCN2", "ADH6")){
  
  tar.exp <- matrix_uniq.filt.tumor[ii,]
  y <- as.numeric(tar.exp)
  data1 <- data.frame()
  for (i in rownames(matrix_uniq.filt.tumor)) {
    dd  <- corr.test(as.numeric(matrix_uniq.filt.tumor[i,]), y , method="pearson",adjust = "fdr")
    data1 = rbind(data1,data.frame(gene=i,cor=dd$r,p.value=dd$p))
  }
  data1 <- data1[order(data1$cor,decreasing = T),]
  gene1 <- data1$cor
  names(gene1) <- mapIds(org.Hs.eg.db,keys = data1$gene,column = 'ENTREZID',
                         keytype = 'SYMBOL',multiVals='filter')
  gene1 <- na.omit(gene1)
  save(gene1,file = paste0(ii,'.cor.Rdata'))
  #分析
  options(digits = 4)#输出结果的小数点后保留4位
  kegggsea <-  gseKEGG(
    gene1,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE,
    seed = TRUE)
  
  write.table(kegggsea,paste0(ii,'.KEGGgsea.txt'),sep = '\t',quote = F,row.names = F)
  #作图
  pdf(paste0(ii,'.KEGG.pdf'),width = 15,height = 10)
  p1 <- gseaplot2(kegggsea, geneSetID = 1:10,pvalue_table = FALSE,
                  base_size = 13,  ## 字体大小
                  rel_heights = c(1.5, 0.5, 0.5),color= brewer.pal(10,'Paired'))
  print(p1)
  dev.off()
  #
  pdf(paste0(ii, '.KEGG.dotplot.pdf'), width = 8, height = 6)
  p2 <- dotplot(kegggsea, showCategory = 10)
  print(p2)
  dev.off()
  #}
  
  #写循环分析GO的GSEA####
  file <- paste0(ii,'.cor.Rdata')
  #load相关性结果
  load(file)
  #分析
  options(digits = 4)#输出结果的小数点后保留4位
  gogsea <-  gseGO(
    gene1,
    OrgDb = org.Hs.eg.db,
    ont = 'ALL',
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    seed = TRUE)
  
  write.table(gogsea,paste0(ii,'.GOgsea.txt'),sep = '\t',quote = F,row.names = F)
  #作图
  pdf(paste0(ii,'.GO.pdf'),width = 15,height = 10)
  p3 <- gseaplot2(gogsea, geneSetID = 1:10,pvalue_table = FALSE,## 不带P值
                  base_size = 13,  ## 字体大小
                  rel_heights = c(1.5, 0.5, 0.5))
  print(p3)
  dev.off()
  
  pdf(paste0(ii,".GO.dotplot.pdf"), width = 10, height = 6)
  p4 <- dotplot(gogsea, showCategory = 10)
  print(p4)
  dev.off()
}

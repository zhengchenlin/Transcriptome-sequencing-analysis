library(Seurat)
library(GENIE3)
library(SeuratData)
library(tidyverse)

deg.cluster <-FindAllMarkers(OPC)
(deg.cluster %>% filter(cell == c("s.cell")) )$gene   -> OPC2
(deg.cluster  %>% top_n(5000,avg_log2FC))$gene   -> NG
(deg.cluster )$gene   -> OPC2
OPC2 <- OPC2[-c(grep("^RP",OPC2))]
OPC2<-unique(OPC2)
exprMatr <- as.matrix(OPC@assays$RNA@data[OPC2,])

# Genes that are used as candidate regulators
regulators <- c("SALL3","HES5","ZNF704","JUND","SOX4")
regulators<-unique(regulators)
weightMat <- GENIE3(exprMatr, regulators=regulators)


linkList <- getLinkList(weightMat, reportMax=1000)
write.csv(linkList, "F:/bioinformation/data/基底神经节/links1.csv", row.names = T)



JS<-GetAssayData(OPC_JS,slot = 'counts')
JS_ALL<-as.matrix(JS)
regulators<-c("SALL3","HES5","ZNF704","SOX6","JUND","SOX4")
weightMat <- GENIE3(as.matrix(JS), regulators=regulators,nCores = 4)
linkList <- getLinkList(weightMat)
##export table and visable by Cytoscope
write.csv(linkList, "F:/bioinformation/data/基底神经节/linkList.csv", row.names = T)


a <- ls()
rm(list=a[which(a!='OPC'&a!='gene_all_jdsjj')])

write.table(b,'CGE_linlist.txt',row.names = F,col.names = T,sep = ',')

##performe Lin1 and Lin2 differential analysis
nmGABA_new<-OPC
Idents(nmGABA_new)<-cds$State
State_idents<-c('s1','s2','s3')
names(State_idents)<-levels(nmGABA_new)
nmGABA_new<-RenameIdents(nmGABA_new,State_idents)
s3<-FindMarkers(nmGABA_new,ident.1 = 's3',ident.2 = c("s1","s2"))
s3$type<-ifelse(s3$p_val_adj<0.05&abs(s3$avg_log2FC)>0.25,
                   ifelse(s3$avg_logFC>0.25,'s3','s2','s1'),'NO')

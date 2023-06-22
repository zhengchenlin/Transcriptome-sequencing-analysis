
npcs=NPCS
npcs <- readRDS("~/mid-count/npcs.rds")
npcs<-subset(npcs,idents=c('ESCO2+/CENPM+','CENPF+/TOP2A+','HES1+/ASCL1+','PAX6+/NEUROD1+'))

#提取特定时期数据
npcs<-subset(npcs,cells = row.names(subset(npcs@meta.data,
                                           npcs@meta.data$Week=='GW09'
)))
DimPlot(npcs, reduction = "tsne",
        split.by = 'Week',cols=cols, label = F)
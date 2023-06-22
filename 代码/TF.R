## We compute Viper Scores 

# 0.安装R包 ----
devtools::install_github('caleblareau/BuenColors')
utils::install.packages(pkgs = "ggstatsplot")
# InstallData("pbmc3k") 

# 1.加载R包和测试数据 ----
rm(list = ls())
library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout = 10000)
data("pbmc3k")  
sce <- pbmc
library(Seurat)
# 一个seurat对象
library(Seurat)
library(dorothea)
library(tidyverse)

# 获取包自带数据库
## We read Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C","D"))

sce <- run_viper(OPC, regulon,
                 options = list(method = "scale", minsize = 4, 
                                eset.filter = FALSE, cores = 1,  
                                verbose = FALSE))

Assays(sce)



DefaultAssay(object = sce) <- "dorothea"
table(Idents(sce))

library(future)
# check the current active plan
plan()
plan("multiprocess", workers = 4)
plan()

sce.markers <- FindAllMarkers(object = sce, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
DT::datatable(sce.markers)
pro='dorothea-markers-for-pbmc3k'
write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
save(sce.markers,file = paste0(pro, '_sce.markers.Rdata'))

library(dplyr) 
sce.markers$fc = sce.markers$pct.1 - sce.markers$pct.2
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, fc)
sce@assays$dorothea@data[1:4,1:4]
top10 =top10[top10$fc > 0.5,] 
sce <- ScaleData(sce )

sce@assays$dorothea@scale.data[1:4,1:4]
DoHeatmap(sce,top10$gene,size=3,slot = 'scale.data')



## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df <- GetAssayData(sce, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()
viper_scores_df[1:4,1:4]


## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Idents(sce)), 
                            cell_type = as.character(Idents(sce)),
                            check.names = F)
head(CellsClusters)

## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))
# For visualization purposes, we select the 20 most variable TFs across clusters according to our scores.
head(summarized_viper_scores)

## We select the 20 most variable TFs. (20*9 populations = 180)
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(180, var) %>%
  distinct(tf)
highly_variable_tfs

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 
palette_length = 100
my_color = colorRampPalette(c("#11998e","#38ef7d"))(palette_length)
colnames(summarized_viper_scores_df)



summarized_viper_scores_df[1:4,1:4]




my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))
library(pheatmap)
viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       main = "DoRothEA (ABC)", angle_col = 45,
                       treeheight_col = 0,  border_color = NA) 



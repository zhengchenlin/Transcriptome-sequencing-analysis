########################去除列########################
#加载R包
library(Seurat)
#加载数据
data(OPC)
#显示mate.date中的列名
names(x = OPC[[]])
#去除特定列
toOPC[['CC.Difference']] <- NULL
names(x = OPC[[]])
#删除多列
OPC@meta.data<-OPC@meta.data[,c(-11:-17)]


########################改列名########################
#修改第四列名字为number
colnames(OPC@meta.data)[4] <- "Week"


########################替换列内容########################
#将矩阵第1列的数据都修改为basal_ganglion
OPC@meta.data[,1 ] = "basal_ganglion"
#替换某列中的特定数据
OPC@meta.data$orig.ident[which(OPC@meta.data$orig.ident =='basal_ganglion')] <- 'subpallium'

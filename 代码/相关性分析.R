av <-AverageExpression(OPC,
                       group.by = "ident",
                       assays = "RNA")
av=av[[1]]

merge_avg <- merge(av, AAA, by = 0, all = TRUE) # by = 0 表示按行名合并，all = TRUE 表示添加缺失值0

merge_avg[is.na(merge_avg)] <- 0
rownames(merge_avg) <- merge_avg$Row.names
merge_avg[,1] <- NULL
cg=names(tail(sort(apply(merge_avg, 1, sd)),1000))
avg_data <- cor(merge_avg[cg,],method = 'spearman')
avg_data <- avg_data[14:17,1:13]
pheatmap::pheatmap(t(avg_data))

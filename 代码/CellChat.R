devtools::install_github("sqjin/CellChat")
library(CellChat)
library(ggplot2)
library(ggalluvial)
options(stringsAsFactors = FALSE)

#提取细胞的基因表达数据
data.input  <- OPC@assays$RNA@data
# create a dataframe consisting of the cell labels
identity = data.frame(group =OPC@active.ident,row.names = names(OPC@active.ident)) 
unique(identity$group) # check the cell labels

#创建cellchat对象
cellchat <- createCellChat(data.input)
cellchat

#将metadata信息加到CellChat对象中
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
#计算每个群的细胞数量
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

#导入受配体数据ku
CellChatDB <- CellChatDB.human 

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object

#预处理
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel  这里似乎有一些bug，在Linux上居然不行。de了它。
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  

#相互作用推断
cellchat <- computeCommunProb(cellchat)  #注意这个函数如果你可以用就用，这个是作者的。
mycomputeCommunProb <-edit(computeCommunProb)  # computeCommunProb内部似乎有一些bug，同一套数据在window10上没事，到了Linux上有报错。发现是computeExpr_antagonist这个函数有问题，(matrix(1, nrow = 1, ncol = length((group))))，中应为(matrix(1, nrow = 1, ncol = length(unique(group))))？ 不然矩阵返回的不对。de了它。
environment(mycomputeCommunProb) <- environment(computeCommunProb)
cellchat <- mycomputeCommunProb(cellchat)  # 这儿是我de过的。

#通过计算链路的数量或汇总通信概率来计算细胞间的聚合通信网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

#查看结果
cellchat@netP$pathways


vertex.receiver = seq(1,4) # a numeric vector
# check the order of cell identity to set suitable vertex.receiver
cellchat@LR$LRsig$pathway_name
#cellchat@LR$LRsig$antagonist
#定义展示的pathway
pathways.show <- "PDGF"
#可视化
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)   # 原函数


netVisual_aggregate(cellchat, signaling = c("PTN"), layout = "circle", vertex.size = groupSize,pt.title=20)
#展示通路中受配体的贡献
netAnalysis_contribution(cellchat, signaling = pathways.show)

#识别细胞间通信网络中的主要发送者、接收者、调解者和影响者
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 15, height = 6, font.size = 10)

#识别分泌细胞外向交流模式
nPatterns = 3
# 同样在这里遇到了bug，难道说是我没有安装好吗，de了它。
cellchat <- myidentifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)  
myidentifyCommunicationPatterns <- edit(identifyCommunicationPatterns)
environment(myidentifyCommunicationPatterns) <- environment(identifyCommunicationPatterns)
#使用热图可视化
cellchat <- myidentifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# 使用桑吉图可视化
netAnalysis_river(cellchat, pattern = "outgoing")
# 使用气泡图可视化
netAnalysis_dot(cellchat, pattern = "outgoing")

#识别目标细胞的传入(incoming)通信模式 
cellchat <- myidentifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)  
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

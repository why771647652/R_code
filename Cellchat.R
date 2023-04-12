rm(list=ls())

library("Seurat")
library("tidyverse")
library("patchwork")
library("dplyr")
library("harmony")
library("ggrepel")
library("ggplot2")
library("reshape2")
library("ggalluvial")
library("usethis")
library("devtools")
library("remotes")
library("ComplexHeatmap")
library("ggsci")
library("GSEABase")
library("enrichplot")
library("msigdbr")
library("stringr")
library("RColorBrewer")
library("tibble")
library("ComplexHeatmap")
library("future")
library("GSVA")
library("GSEABase")
library("clusterProfiler")
library("AUCell")
library("BiocParallel")
library("org.Hs.eg.db")
library("fgsea")
library("circlize")
library("CellChat")

scrna1<-read_rds("Res.Mono.rds")

scrna1<-scrna1[,scrna1@meta.data$Res %in% c("Response.After","NonResponse.After")]

scrna2<-read_rds("Total.T.cells.rds")

scrna2<-scrna2[,scrna2@meta.data$Res %in% c("Response.After","NonResponse.After")]

scrna<-merge(scrna1,y=c(scrna2))

DefaultAssay(scrna) <- 'RNA'

#表达量数据标准化

#LogNormalize 的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000 ) 

scrna <- NormalizeData(scrna, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_object_copy <- scrna

seurat_object_copy <- SetIdent(seurat_object_copy,value=seurat_object_copy$Res)

levels <- levels(factor(seurat_object_copy$Res))

for(x in 1:length(levels)){
  value <- levels[x]
  seurat_object <- subset(seurat_object_copy,idents=value)
  cellchat <- createCellChat(object=seurat_object@assays$RNA@data,meta=seurat_object@meta.data,group.by = 'celltype')
  groupSize <- as.numeric(table(cellchat@idents))
  CellChatDB <- CellChatDB.human  #CellChatDB.human
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)#过滤出所需基因
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)#结果保存在cellchat@LR$LRsig
  cellchat <- projectData(cellchat,PPI.human)#慢，找到配体受体关系后，将配体手提对的表达之投射到PPI上，来对@data.signaling中的表达值进行矫正。结果保存在@data.project
  #根据表达值推测细胞互作的概率
  cellchat <- computeCommunProb(cellchat,raw.use = F,population.size = T)#慢,如果不想用上一步PPI矫正的结果，raw.use=T
  #cellchat <- filterCommunication(cellchat,min.cells=10)#过滤细胞少的组别
  df.net <- subsetCommunication(cellchat)
  file.name <- paste(value,'_net_lr.csv',sep = '')
  write.csv(df.net,file = file.name)
  #推断信号通路水平的细胞通讯网络，结果储存在@netP下面，有一个概率值和对应pvalue
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp <- subsetCommunication(cellchat,slot.name='netP')
  file.name <- paste(value,'_net_pathway.csv',sep = '')
  write.csv(df.netp,file = file.name)
  cellchat <- aggregateNet(cellchat)
  file.name <- paste(value,'_cellchat.Rdata',sep = '')
  save(cellchat,file = file.name)
}


#merge cellchat

load("~/Final.CRS/另外的Monocyte analysis/NonResponse.After_cellchat.Rdata")

cellchat1<-cellchat

load("~/Final.CRS/另外的Monocyte analysis/Response.After_cellchat.Rdata")

cellchat2<-cellchat

object.list <- list(NonRes = cellchat1, Res = cellchat2)

cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)


cellchat <- cellchat_comparisonAnalysis_Res_NonRes

###比较交互总数和交互强度
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


### heatmap
pathways.show <- c("MHC-I") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

pathways.show <- c("MHC-II") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

### bubble MHC-II pathway
df <- net
selected_pathway <- c('MHC-II')
selected_source <- c("CD16 Monocyte","DN Monocyte","VCAN Monocyte","LGALS2 Monocyte")
selected_target <- c('T cell')

df_selected <- df[which(df$source %in% selected_source),]
df_selected <- df_selected[which(df_selected$target %in% selected_target),]
df_selected <- df_selected[which(df_selected$pathway_name %in% selected_pathway),]
df_selected$interaction <- paste(df_selected$source,'|',df_selected$target,sep = '')
df_selected$log.prob <- log10(df_selected$prob)

write_csv2(df_selected,"MHC-II.Bubble.csv")

midpoint <- (max(df_selected$log.prob)+min(df_selected$log.prob))/2
bubble <- ggplot(df_selected,aes(x=interaction,y=interaction_name))+
  geom_point(aes(size = -log10(pval+0.001),colour = log.prob))+
  facet_wrap(vars(sample_location))+
  scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="firebrick3",mid = "white",low ="navy",midpoint = midpoint)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90,vjust = 0,hjust = 0),
        text = element_text(size = 10))
bubble

### bubble MHC-I pathway
df <- net
selected_pathway <- c('MHC-I')
selected_source <- c("CD16 Monocyte","DN Monocyte","VCAN Monocyte","LGALS2 Monocyte")
selected_target <- c('T cell')

df_selected <- df[which(df$source %in% selected_source),]
df_selected <- df_selected[which(df_selected$target %in% selected_target),]
df_selected <- df_selected[which(df_selected$pathway_name %in% selected_pathway),]
df_selected$interaction <- paste(df_selected$source,'|',df_selected$target,sep = '')
df_selected$log.prob <- log10(df_selected$prob)

write_csv2(df_selected,"MHC-I.Bubble.csv")
midpoint <- (max(df_selected$log.prob)+min(df_selected$log.prob))/2
bubble <- ggplot(df_selected,aes(x=interaction,y=interaction_name))+
  geom_point(aes(size = -log10(pval+0.001),colour = log.prob))+
  facet_wrap(vars(sample_location))+
  scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="firebrick3",mid = "white",low ="navy",midpoint = midpoint)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90,vjust = 0,hjust = 0),
        text = element_text(size = 10))
bubble
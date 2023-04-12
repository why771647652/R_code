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

################# 颜色设定================================================

mycol <- c(brewer.pal(5,"Set1"), brewer.pal(8,"Set2"),
           brewer.pal(11,"Set3"), brewer.pal(12,"Paired"),
           brewer.pal(8,"Accent"), brewer.pal(11,"Spectral"),
           brewer.pal(11,"BrBG"), brewer.pal(11,"PiYG"),
           brewer.pal(11,"PuOr"),brewer.pal(11,"RdBu"))

### Monocyte 提取

Monocyte<- read_rds("Final.Monocyte.rds")

### recluster new monocytes

Monocyte <-FindNeighbors(Monocyte, reduction = "harmony",dims = 1:20,k.param = 40) %>% FindClusters(resolution=0.5,n.iter = 50)

Monocyte <- RunTSNE(Monocyte,reduction = "harmony",dims = 1:20,n.neighbors = 20, min.dist = 0.2)


### Response anatotaion

Monocyte@meta.data$res.ano="other"
Monocyte$res.ano[Monocyte$orig.ident=="Patient1.Before"]="Res.Before"
Monocyte$res.ano[Monocyte$orig.ident=="Patient1.Normal"]="Response.Inter"
Monocyte$res.ano[Monocyte$orig.ident=="Patient2.Before"]="NR.Before"
Monocyte$res.ano[Monocyte$orig.ident=="Patient2.CRS"]="NR.Inter"
Monocyte$res.ano[Monocyte$orig.ident=="Patient3.Before"]="NR.Before"
Monocyte$res.ano[Monocyte$orig.ident=="Patient3.CRS"]="NR.Inter"
Monocyte$res.ano[Monocyte$orig.ident=="Patient4.Before"]="NR.Before"
Monocyte$res.ano[Monocyte$orig.ident=="Patient4.CRS"]="NR.Inter"

### cell type tsne 

plot1 <- DimPlot(Monocyte, pt.size = 0.2, cols = mycol, label=T, repel=T,
                 raster = FALSE, label.size = 5, reduction = "tsne",group.by = "celltype")
plot1

ggsave("Monocyte_cell_type.pdf",plot1,height = 5,width = 7)

### stage.ano tsne 

plot2 <- DimPlot(Monocyte, pt.size = 0.2, cols = mycol, label=T, repel=T,
                 raster = FALSE, label.size = 5, reduction = "tsne",group.by = "stage.ano")
plot2

ggsave("Monocyte_stage.ano.pdf",plot2,height = 5,width = 7)

### seurat_clusters tsne 

plot3 <- DimPlot(Monocyte, pt.size = 0.2, cols = mycol, label=T, repel=T,
                 raster = FALSE, label.size = 5, reduction = "tsne",group.by = "seurat_clusters")
plot3

ggsave("Monocyte_seurat_clusters.pdf",plot3,height = 5,width = 5.5)

### seurat_clusters tsne 

plot4 <- DimPlot(Monocyte, pt.size = 0.2, cols = mycol, label=T, repel=T,
                 raster = FALSE, label.size = 5, reduction = "tsne",group.by = "pat.ano")
plot4

ggsave("Monocyte_pat.ano.pdf",plot4,height = 5,width = 6.2)

### marker gene 表达图谱

mycolor <- c('lightgrey',"red")#设置颜色  

FeaturePlot(Monocyte,features = c("LGALS2"),min.cutoff = 0.3, max.cutoff = 1,cols = mycolor,pt.size = 0.5,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("LGALS2.pdf", width = 5.5, height = 5)

FeaturePlot(Monocyte,features = c("CD36"),min.cutoff = 2, max.cutoff =2.5,cols = mycolor,pt.size = 0.5,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
ggsave("CD36.pdf", width = 5.5, height = 5)

FeaturePlot(Monocyte,features = c("FCGR3A"),min.cutoff = 0.5, max.cutoff = 2,cols = mycolor,pt.size = 0.5,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("FCGR3A.pdf", width = 5.5, height = 5)

FeaturePlot(Monocyte,features = c("FCGR3B"),min.cutoff = 0, max.cutoff = 2,cols = mycolor,pt.size = 0.5,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("FCGR3B.pdf", width = 5.5, height = 5)

dim(Monocyte)

### TOP5 marker gene heatmap cell gene

seurat_object <- Monocyte
seurat_object <- SetIdent(seurat_object,value = seurat_object@meta.data$celltype)
markers <- FindAllMarkers(seurat_object,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25,return.thresh = 0.01)
top10 <- markers%>%group_by(cluster)%>%slice_max(n=10,order_by=avg_log2FC)

pal <- c(pal_npg('nrc')(10),pal_jco('default')(10),pal_aaas('default')(10),pal_npg('nrc',alpha = 0.5)(10))
new.cluster.ids <- data.frame(cluster=seurat_object@meta.data$seurat_clusters,
                              celltype=seurat_object@meta.data$celltype)
new.cluster.ids<- new.cluster.ids[!duplicated(new.cluster.ids$cluster),]
new.cluster.ids$cluster <- as.character(new.cluster.ids$cluster)
top10$cluster <- as.character(top10$cluster)
top10 <- merge(top10,new.cluster.ids,by.x='cluster',by.y='celltype',all.x=T)
top10 <- top10[order(top10$cluster),]
top10 <- top10$gene
top10 <- top10[!duplicated(top10)]

information <- data.frame(seurat_clusters=seurat_object@meta.data$seurat_clusters,
                          celltype=seurat_object@meta.data$celltype)
information <- information[!duplicated(information$seurat_clusters),]
rownames(information) <- information$seurat_clusters
information <- information[order(information$seurat_clusters),]
information <- information[order(information$celltype),]

set_color <- pal[1:length(levels(factor(information$celltype)))]
names(set_color) <- levels(factor(information$celltype))
set_color2 <- pal[1:length(levels(factor(information$seurat_clusters)))]
names(set_color2) <- levels(factor(information$seurat_clusters))	
column_annotation <- HeatmapAnnotation(df = information, 
                                       col = list(seurat_clusters=set_color2,
                                                  celltype=set_color))
seurat_object<- SetIdent(seurat_object,value=seurat_object@meta.data$seurat_clusters) 
mat <- AverageExpression(seurat_object,
                         group.by = 'ident')
mat <- mat$RNA
mat <- as.matrix(mat[top10,rownames(information)])
mat <- t(scale(t(mat)))

max(mat)
min(mat)
mat[mat>3] = 3
mat[mat<(-3)] = (-3)

colnames(mat)

mat<-mat[,c(4,5,1,6,0,3,2)]

complexheatmap <- Heatmap(mat, 
                          col=colorRamp2(c(-3,0,3),c('navy','white','firebrick3')), 
                          column_names_side = "top",  
                          border=T,
                          cluster_rows=F,
                          cluster_columns =F,
                          show_column_names = T, 
                          show_row_names = T,
                          column_split = information$celltype,
                          row_names_gp = gpar(fontsize = 8), 
                          column_names_gp = gpar(fontsize = 8), 
                          top_annotation = column_annotation,width = 5,height = 5)
complexheatmap

pdf('Top.10.marker_heatmap_complex.pdf',width = 6,height=5)
print(complexheatmap)
dev.off()

### CRS stack plot

Idents(Monocyte)=Monocyte@meta.data$celltype #y轴

tab<-table(Idents(Monocyte),Monocyte$stage.ano)

tab

tab<-prop.table(tab,2)*100

tab ###1是行，二是列

tab<-as.data.frame(tab)

write.csv(tab,"Monocyte.fraction1.csv")

ggplot(tab,aes(x=Freq,y=Var2,fill=Var1)) +
  geom_col(position = position_fill()) +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        # panel.background = element_rect(fill = 'grey90'),
        axis.line=element_line(size=1,colour="black"),
        aspect.ratio = 2) +
  scale_fill_brewer(palette = "Set3") +
  coord_cartesian(expand = 0) +
  xlab('Cell Proportion') + ylab("Patient ID")

ggsave("Monocyte_stage.ano_cellpro.pdf", width = 7, height = 5)

### 提取 CRS patient 资料 每个细胞的boxplot

Idents(Monocyte)=Monocyte@meta.data$celltype #y轴

tab1<-table(Idents(Monocyte),Monocyte$orig.ident)

tab1

tab1<-prop.table(tab1,2)*100

tab1 ###1是行，二是列

tab1<-as.data.frame(tab1)

mab<- tab1[tab1$Var2 %in% c("Patient2.Before","Patient2.CRS","Patient3.Before","Patient3.CRS","Patient4.Before","Patient4.CRS"),]

mab$Var3="other"
mab$Var3[mab$Var2=="Patient2.Before"]="Before.CRS"
mab$Var3[mab$Var2=="Patient2.CRS"]="CRS"
mab$Var3[mab$Var2=="Patient3.Before"]="Before.CRS"
mab$Var3[mab$Var2=="Patient3.CRS"]="CRS"
mab$Var3[mab$Var2=="Patient4.Before"]="Before.CRS"
mab$Var3[mab$Var2=="Patient4.CRS"]="CRS"

write.csv(mab,"Monocyte.fraction2.csv")

boxplot=ggboxplot(mab, x="Var3", y="Freq",
                  xlab="",
                  ylab="Cell Fraction",
                  legend.title="",
                  palette = "bioCo",
                  bxp.errorbar = T,
                  bxp.errorbar.width = 0.2,
                  fill = "Var3",
                  size = 1,short.panel.labs = FALSE)+ylim(0, 80)+facet_wrap(~Var1)+
  theme(
    plot.title    = element_text(color = "black", size   = 16, hjust = 0.5),
    axis.text.x   = element_text(color = "black", size = 12, angle = 0),
    axis.text.y   = element_text(color = "black", size = 16, angle = 0),
    axis.title.x  = element_text(color = "black", size = 20, angle = 0),
    axis.title.y  = element_text(color = "black", size = 16, angle = 90),
    legend.title  = element_text(color = "black", size  = 20),
    legend.text   = element_text(color = "black", size   = 16),
    axis.line.y = element_line(color = "black", linetype = "solid"), # y轴线特征
    axis.line.x = element_line (color = "black",linetype = "solid"), # x轴线特征
    panel.border = element_rect(linetype = "solid", size = 1.2,fill = NA)) # 图四周框起来

boxplot

ggsave("MonoBoxplot.pdf",boxplot, width = 8, height = 5)

### 定义CRS.Score

genes = str_to_upper(c("IL1B","IL6","TNF","CD93","CXCL9","CXCL10","CCL3L1","CXCL2","JUN","IL18","CCL5","PF4","CXCL8"))
genes <- as.data.frame(genes)
genes <- as.list(genes)

Monocyte<- AddModuleScore(Monocyte,features = genes,name = 'CRS')

### 定义TLR.Score

gene2 = str_to_upper(c("CXCL8","CCL3L1","CCL3","PPBP","CCL4L2","CXCL2","CCL4","GNG11","CXCL3","NFKBIA","CXCR4","CCL5","CDC42","PF4","GNG5","GNAI3","PF4V1"))

gene2 <- as.data.frame(gene2)

gene2 <- as.list(gene2)

Monocyte<- AddModuleScore(Monocyte,features = gene2,name = 'TLR')

### CRS TLR correlation
Monocyte1 <- Monocyte [,Monocyte @meta.data$celltype %in% c("CD36+ Monocyte")]

Idents(Monocyte1) <-"celltype"
allCells=names(Idents(Monocyte1))
allType = levels(Idents(Monocyte1))

choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(Monocyte1)== x ]
  cg=sample(cgCells,300)
  cg
}))

scrna1 = Monocyte1[, allCells %in% choose_Cells]

mat_cor <- scrna1@meta.data
mat_cor <- mat_cor[,c("CRS1","TLR1")]
color1 <- colorRampPalette(rev(brewer.pal(1,"YlGn")))#设置连续颜色
p1=ggplot(mat_cor, aes(x=CRS1, y=TLR1))+geom_point(size = 1,color = '#EC0101',alpha = 0.5)+#点大小和颜色设置
  theme_bw() +
  # 主题细节调整
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
  ) +
  # 添加回归线
  geom_smooth(method = 'lm',se = T,size = 2,fill = '#FEA82F',formula = y ~ x,color="#F9B208") +
  # 添加相关性系数及p值
  stat_cor(method = "pearson",digits = 3,size=8)
p1 
ggsave("TLR.CRS.plot.pdf",p1,height = 4, width = 4.3)


### Dotplot展示CRS发生过程中一些趋化因子细胞因子的变化

PBMC1 <- Monocyte[,Monocyte@meta.data$seurat_clusters %in% c(0)]

markers <- c("IL1B","IL1A","IL18","IL33","CCL3","CCL5","IL10","TNF","CCL4","CCL23","CCL7","CCL12","CXCL1","CXCL5","CXCL6","CXCL3","CXCL4","CXCL8","CXCL2")

scn1 <- PBMC1[,PBMC1@meta.data$stage.ano %in% c("CRS.Before","CRS.After")]

Idents(scn1) <-"stage.ano"

DotPlot(scn1, features = markers,cols = c("grey", "red"),scale = TRUE,scale.by = "size",cluster.idents = T)+RotatedAxis()+ theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size = 15,color = "black"),axis.text.y = element_text(size = 15,color = "black"))+#文字90度呈现 
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置 
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) +scale_size(range = c(8,15))

ggsave("CytokineChange_Dotplot.pdf",width = 10,height =4)


### CRS.Before以及CRS.After GSEA.pathway分析

Idents(PBMC1) <-"stage.ano"

DEGCRS <- FindMarkers(PBMC1,ident.1 = "CRS.Before",ident.2 = "CRS.After",min.pct = 0.01,logfc.threshold = 0.25)

DefaultAssay(Monocyte) <-"RNA"

genelist <- DEGCRS$avg_log2FC

names(genelist) = toupper(rownames(DEGCRS))

genelist <-sort(genelist,decreasing = T)

geneset<-msigdbr(species = "Homo sapiens",category = "C2",subcategory = "KEGG")

geneset<-geneset[,c("gs_name","gene_symbol")]

length(unique(geneset$gs_name))

egmt <- GSEA(genelist,TERM2GENE = geneset,minGSSize = 1,pvalueCutoff = 0.5)

mat <- data.frame(egmt)

mat <- mat[which(mat$pvalue < 0.05),]

mat <- mat[order(mat$NES,decreasing = T),]

write.csv(mat,"GSEA.CD36.CRS.CRS.Before.csv")

gseaplot2(egmt,geneSetID = c(9),pvalue_table = T,color = "red",base_size = 10,ES_geom = "line")

ggsave("Cytokine=Cytokine.pdf",width = 4.5,height = 4)

gseaplot2(egmt,geneSetID = c(14),pvalue_table = T,color = "red",base_size = 10,ES_geom = "line")

ggsave("TLR.pdf",width = 4.5,height = 4)

### Pathway 画图

plot_need<- mat[c(1:8,16:23),]

plot_need$pathway <- as.factor(plot_need$ID)
plot_need$threshold <- as.factor(ifelse(plot_need$NES>=0 ,'Up','Down'))
up <- nrow(plot_need[which(plot_need$threshold=='Up'),])
down <- nrow(plot_need[which(plot_need$threshold=='Down'),])

### GSEA 好看的图片展示
gsea_barplot <- ggplot(plot_need,aes(x = reorder(pathway,NES),y = NES))+
  geom_bar(aes(fill=threshold),stat='identity') + 
  scale_fill_manual(values = c('Up'='firebrick3','Down'='navy'))+
  coord_flip()+
  guides(fill='none')+
  xlab('')+
  ylab('')+
  labs(title = "")+
  geom_text(data = plot_need[1:up,],aes(x=pathway,y=-4,label=pathway),hjust=0,color='black',size = 3)+
  geom_text(data = plot_need[(up+1):(up+down),],aes(x=pathway,y=4,label=pathway),hjust=1,color='black',size = 3)+
  theme_bw()+ #鑳屾櫙鍙樹负鐧借壊
  theme(panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"), #鍘婚櫎榛樿濉厖鐨勭伆鑹诧紝骞跺皢x=0杞村拰y=0杞村姞绮楁樉绀?(size=1)
        panel.grid.major = element_blank(),   #涓嶆樉绀虹綉鏍肩嚎
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
gsea_barplot

ggsave("gsea_Monocyte_KEGG_CRS_Before.pdf",gsea_barplot,width = 8,height = 12)


### APC和TLR的correlation

mat_cor<-read.csv2("Correlatkion.APC.TLR.csv")

mat_cor <- mat_cor[,c("APC1","TLR1")]

color1 <- colorRampPalette(rev(brewer.pal(1,"YlGn")))#设置连续颜色

p1=ggplot(mat_cor, aes(x=APC1, y=TLR1))+geom_point(size = 1,color = '#EC0101',alpha = 0.5)+#点大小和颜色设置
  theme_bw() +
  # 主题细节调整
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
  ) +
  # 添加回归线
  geom_smooth(method = 'lm',se = T,size = 2,fill = '#FEA82F',formula = y ~ x,color="#F9B208") +
  # 添加相关性系数及p值
  stat_cor(method = "pearson",digits = 3,size=8)
p1 
ggsave("TLR.APC.plot.pdf",p1,height = 4, width = 4)

### CRS和TLR的correlation
Monocyte1 <- Monocyte [,Monocyte @meta.data$celltype %in% c("MGST1 Monocyte")]

Idents(Monocyte1) <-"celltype"
allCells=names(Idents(Monocyte1))
allType = levels(Idents(Monocyte1))

choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(Monocyte1)== x ]
  cg=sample(cgCells,500)
  cg
}))

scrna1 = Monocyte1[, allCells %in% choose_Cells]

saveRDS(scrna1,"Correlatkion.CRS.TLR.rds")

mat_cor <- scrna1@meta.data
mat_cor <- mat_cor[,c("CRS1","TLR1")]

write.csv2(mat_cor,"Correlatkion.CRS.TLR.csv")

color1 <- colorRampPalette(rev(brewer.pal(1,"YlGn")))#设置连续颜色
p1=ggplot(mat_cor, aes(x=CRS1, y=TLR1))+geom_point(size = 1,color = '#EC0101',alpha = 0.5)+#点大小和颜色设置
  theme_bw() +
  # 主题细节调整
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
  ) +
  # 添加回归线
  geom_smooth(method = 'lm',se = T,size = 2,fill = '#FEA82F',formula = y ~ x,color="#F9B208") +
  # 添加相关性系数及p值
  stat_cor(method = "pearson",digits = 3,size=8)
p1 

### APC res gene

scn2 <- Monocyte [,Monocyte @meta.data$res.ano %in% c("NR.Inter","Response.Inter")]

seurat_object  <- scn2

expr_matrix <-as.data.frame(t(as.matrix(seurat_object@assays$RNA@data)))
expr_matrix$rename_clusters <- seurat_object@meta.data$celltype
expr_matrix$stage.ano <- seurat_object@meta.data$res.ano

genes <- c("HLA-DPB1","CD74","HLA-A","HLA-F")

violin_list <- list()

for(x in 1:length(genes)){
  gene <- genes[x]
  expr_df <- data.frame(gene=expr_matrix[,gene],
                        stage.ano=expr_matrix$stage.ano,
                        rename_clusters=expr_matrix$rename_clusters)
  violin <- ggplot(expr_df,aes(x=stage.ano,y=gene))+
    geom_violin(trim = FALSE, aes(fill = factor(stage.ano)))+
    #geom_jitter(aes(fill = factor(res.day),color = factor(res.day)))+
    geom_boxplot(width=0.1,fill='white')+
    #stat_compare_means(method = "anova")+
    scale_fill_manual(values = c("Normal.Before" = "steelblue1","CRS.Before"="firebrick1"))+
    scale_color_manual(values = c("Normal.Before" = "steelblue1", "CRS.Before"="firebrick1"))+
    theme_bw()+ #背景变为白色
    theme(axis.text.x=element_text(colour="black"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.x=element_text(colour="black"), #设置x轴标题的字体属性
          axis.title.y=element_text(colour="black"), #设置y轴标题的字体属性
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.text=element_text(colour="black"),  #设置图例的子标题的字体属性
          legend.title=element_text(colour="black"), #设置图例的总标题的字体属性
          panel.grid.major = element_blank(),   #不显示网格线
          panel.grid.minor = element_blank())+
    #facet_wrap(~rename_clusters)+
    stat_summary_bin(fun = mean)+
    ylab(gene)
  
  violin_list[[x]] <- violin
}
pdf('Monocyte.Res.APC.pdf',width = 5,height = 3)

for(x in 1:length(violin_list)){
  print(violin_list[[x]])
}
dev.off()

### Dotplot展示CRS发生过程中一些趋化因子细胞因子的变化

PBMC1 <- Monocyte[,Monocyte@meta.data$seurat_clusters %in% c(0)]

markers <- c("IL1B","IL1A","IL18","CCL3","CCL5","IL10","TNF","CCL4","CCL23","CCL7","CXCL1","CXCL5","CXCL6","CXCL3","CXCL8","CXCL2")

scn1 <- PBMC1[,PBMC1@meta.data$stage.ano %in% c("CRS.Before","CRS.After")]

Idents(scn1) <-"stage.ano"

DotPlot(scn1, features = markers,cols = c("grey", "red"),scale = TRUE,scale.by = "size",cluster.idents = T)+RotatedAxis()+ theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size = 15,color = "black"),axis.text.y = element_text(size = 15,color = "black"))+#文字90度呈现 
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置 
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) +scale_size(range = c(8,15))

ggsave("CytokineChange_Dotplot.pdf",width = 10,height =4)








### Response stack plot

scn2 <- Monocyte [,Monocyte @meta.data$res.ano %in% c("NR.Inter","Response.Inter")]

Idents(scn2)=scn2@meta.data$celltype #y轴

tab<-table(Idents(scn2),scn2$res.ano)

tab

tab<-prop.table(tab,2)*100

tab ###1是行，二是列

tab<-as.data.frame(tab)

ggplot(tab,aes(x=Freq,y=Var2,fill=Var1)) +
  geom_col(position = position_fill()) +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        # panel.background = element_rect(fill = 'grey90'),
        axis.line=element_line(size=1,colour="black"),
        aspect.ratio = 2) +
  scale_fill_brewer(palette = "Set3") +
  coord_cartesian(expand = 0) +
  xlab('Cell Proportion') + ylab("Patient ID")

ggsave("Monocyte_res.ano_cellpro.pdf", width = 7, height = 3)

ggplot(tab, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = pal_nejm()(7))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.9,hjust = 0.9))+
  xlab("")+ylab("percentage")

ggsave("Monocyte_stage.ano_cellpro.pdf", width = 5, height = 3)



### 每群细胞的GSVA分析
av <- AverageExpression(Monocyte,group.by = "celltype",assays = "RNA")

all_gene_sets <-msigdbr(species = "Homo sapiens",category = "C2",subcategory = "KEGG")

length(unique(table(all_gene_sets$gs_name)))

tail(table(all_gene_sets$gs_name))

tail(sort(table(all_gene_sets$gs_name)))

gs <- split(all_gene_sets$gene_symbol,all_gene_sets$gs_name)

gs <- lapply(gs,unique)

gsc <- GSEABase::GeneSetCollection(mapply(function(geneIds,keggId){
  GeneSet(geneIds,geneIdType=EntrezIdentifier(),collectionType=KEGGCollection(keggId),setName=keggId)},gs,names(gs)))
gsc

geneset<-gsc

X=data.frame(av)

X<- as.matrix(X)

es.max <- gsva(expr = X,geneset,kcdf="Gaussian",verbose=T,parallel.sz=8)

class(es.max)

cg<- names(tail(sort(apply(es.max,1,sd)),10))

### 选取Top5genelist
df<-do.call(rbind,
            lapply(1:ncol(es.max),function(i){
              dat=data.frame(
                path =rownames(es.max),
                cluster =colnames(es.max)[i],
                sd.1=es.max[,i],
                sd.2=apply(es.max[,-i],1,sd)
              )
            }))

df$fc = df$sd.1-df$sd.2

Top5 <- df %>%group_by(cluster) %>%top_n(8,fc)

pheatmap(es.max[Top5$path,],show_rownames = T,scale="row",cluster_cols = F,cellwidth = 10,cellheight = 10,cluster_rows = F,legend = F)

ggsave("MonocyteGSVA.heatmap.pdf",width = 8,height = 15)



### 定义APC Score

gene3 = str_to_upper(c("HLA-DRB5","HLA-DQB1","CD74","HSPA8","PSME1","HLA-C","HLA-DQA1","HSPA5","GZMB"))
gene3 <- as.data.frame(gene3)
gene3 <- as.list(gene3)

Monocyte<- AddModuleScore(Monocyte,features = gene3,name = 'APC')

###计算APC.Score和TLR.Score的相关性(耿哥，需要帮助重新画图 )

### 随机抽取100个细胞
Idents(Monocyte) <-"celltype"
allCells=names(Idents(Monocyte))
allType = levels(Idents(Monocyte))

choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(Monocyte)== x ]
  cg=sample(cgCells,200)
  cg
}))

scrna1 = Monocyte[, allCells %in% choose_Cells]

mat_cor <- scrna1@meta.data
mat_cor <- mat_cor[,c("APC1","TLR1")]
color1 <- colorRampPalette(rev(brewer.pal(1,"YlGn")))#设置连续颜色
p1=ggplot(mat_cor, aes(x=APC1, y=TLR1))+geom_point(size = 1,color = '#EC0101',alpha = 0.5)+#点大小和颜色设置
  theme_bw() +
  # 主题细节调整
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
  ) +
  # 添加回归线
  geom_smooth(method = 'lm',se = T,size = 2,fill = '#FEA82F',formula = y ~ x,color="#F9B208") +
  # 添加相关性系数及p值
  stat_cor(method = "pearson",digits = 3,size=8)
p1 
ggsave("TLR.APC.plot.pdf",p1,height = 4, width = 4)

###计算CRS.Score和TLR.Score的相关性(耿哥，需要帮助重新画图 )
mat_cor <- Monocyte@meta.data

mat_cor <- mat_cor[,c("CRS1","TLR1")]

color1 <- colorRampPalette(rev(brewer.pal(1,"YlGn")))#设置连续颜色

p1=ggplot(mat_cor, aes(x=APC1, y=TLR1))+geom_point(size = 1,color = '#EC0101',alpha = 0.5)+#点大小和颜色设置
  theme_bw() +
  # 主题细节调整
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
  ) +
  # 添加回归线
  geom_smooth(method = 'lm',se = T,size = 2,fill = '#FEA82F',formula = y ~ x,color="#F9B208") +
  # 添加相关性系数及p值
  stat_cor(method = "pearson",digits = 3,size=8)

p1 

### correlation 第一种做图
p1 <-ggplot(data, aes(x = MGST1_Monocyte, y = TRDC_NKT))+geom_point(color = ifelse(pearson_rho >0, "#f15e4c","#009bc7"))+
  geom_smooth(formula = y ~ x,method = "lm", color = ifelse(pearson_rho >0, "#f15e4c","#009bc7"), fill = ifelse(pearson_rho >0, "#f15e4c","#009bc7")) +
  stat_cor(method = "pearson")  + theme_bw() 
p1
ggsave("NKT.Monocyte.plot.pdf",p1,height = 4, width = 4)

### correlation 第二种做图
p1=ggplot(mat_cor, aes(x=CRS1, y=TLR1))+geom_point(size = 1,color = '#EC0101',alpha = 0.5)+#点大小和颜色设置
  theme_bw() +
  # 主题细节调整
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
  ) +
  # 添加回归线
  geom_smooth(method = 'lm',se = T,size = 2,fill = '#FEA82F',formula = y ~ x,color="#F9B208") +
  # 添加相关性系数及p值
  stat_cor(method = "pearson",digits = 3,size=8)

p1 




###Normal.Before以及CRS.Before GSEA.pathway分析

DEGCRS <- FindMarkers(PBMC1,ident.1 = "CRS.Before",ident.2 = "Normal.Before",min.pct = 0.01,logfc.threshold = 0.25)

DefaultAssay(PBMC1) <-"RNA"

genelist <- DEGCRS$avg_log2FC

names(genelist) = toupper(rownames(DEGCRS))

genelist <-sort(genelist,decreasing = T)

geneset<-msigdbr(species = "Homo sapiens",category = "C2",subcategory = "PID")

geneset<-geneset[,c("gs_name","gene_symbol")]

length(unique(geneset$gs_name))

egmt <- GSEA(genelist,TERM2GENE = geneset,minGSSize = 1,pvalueCutoff = 0.5)

mat <- data.frame(egmt)

write.csv(mat,"GSEA.MGST1.Normal.Before.CRS.Before.csv")

gseaplot2(egmt,geneSetID = c(2),pvalue_table = T,color = "red",base_size = 10,ES_geom = "line")

ggsave("CD16CD36 CXCR4.pdf",width = 4.5,height = 4)

### CRS.before 以及 Normal.Before Heatmap

mat=log1p(AverageExpression(PBMC1,verbose = FALSE,group.by = "stage.ano")$RNA)

mab <- mat[,c("CRS.Before","Normal.Before")]

mab = as.data.frame(mab)

### antigen presentation pathway gene

MHC.Gene <- c("PSME2","HLA-DRB1","IFI30","HLA-DMA","TAP1","HLA-DQB1","HLA-DPB1","HLA-DPA1","CD74","TAPBP","PDIA3","HSP90AB1","TAP2","HLA-DRA","HLA-A","HSPA8","CALR","CIITA","HSP90AA1","RFXANK","HLA-F","PSME1","B2M","HLA-C","CTSL","HLA-DQA1","HLA-B","HSPA5")

mab_need <- mab[c("PSME2","HLA-DRB1","IFI30","HLA-DMA","TAP1","HLA-DQB1","HLA-DPB1","HLA-DPA1","CD74","TAPBP","PDIA3","HSP90AB1","TAP2","HLA-DRA","HLA-A","HSPA8","CALR","CIITA","HSP90AA1","RFXANK","HLA-F","PSME1","B2M","HLA-C","CTSL","HLA-DQA1","HLA-B","HSPA5"),]

p <- pheatmap(mab_need,cluster_rows = F,scale="row",cluster_cols = F,col=colorRamp2(c(-1.0,0,1.0),c('navy','white','firebrick3')),cellwidth = 15,cellheight = 15)

p

### CRS.Before 以及 CRS.After Heatmap

mat=log1p(AverageExpression(PBMC1,verbose = FALSE,group.by = "stage.ano")$RNA)

mab <- mat[,c("CRS.Before","CRS.After")]

mab = as.data.frame(mab)

### TLR 通路差异gene list

TLR<-c("CXCL8","IL1B","CCL3L1","FOS","CCL4","CCL5","NFKBIA","TLR4","PF4","PF4V1","CCL3","PPBP","CCL4L2","RAC1","TIRAP")

mab_need <- mab[c("CXCL8","IL1B","CCL3L1","FOS","CCL4","CCL5","NFKBIA","TLR4","PF4","PF4V1","CCL3","PPBP","CCL4L2","RAC1","TIRAP"),]

p <- pheatmap(mab_need,cluster_rows = F,scale="row",cluster_cols = F,col=colorRamp2(c(-1.0,0,1.0),c('navy','white','firebrick3')),cellwidth = 15,cellheight = 15)

p

### Cytokine-Cytokine Interaction 通路差异gene list

CCI<-c("CXCL8","IL1B","CCL3L1","CCL3","CCL4","PPBP","CXCR4","CCL4L2","CCL4","CXCL2","CXCL3","CCL5","IL18","PF4","TNFRSF1B")

mab_need <- mab[c("CXCL8","IL1B","CCL3L1","CCL3","CCL4","PPBP","CXCR4","CCL4L2","CCL4","CXCL2","CXCL3","CCL5","IL18","PF4","TNFRSF1B"),]

p <- pheatmap(mab_need,cluster_rows = F,scale="row",cluster_cols = F,col=colorRamp2(c(-1.0,0,1.0),c('navy','white','firebrick3')),cellwidth = 15,cellheight = 15)

p


### TLR pathway violin plot

seurat_object<-Monocyte [,Monocyte @meta.data$stage.ano %in% c("CRS.After","CRS.Before")]

seurat_object <- seurat_object[,seurat_object@meta.data$seurat_clusters %in% c(0)]

expr_matrix <-as.data.frame(t(as.matrix(seurat_object@assays$RNA@data)))
expr_matrix$rename_clusters <- seurat_object@meta.data$celltype
expr_matrix$stage.ano <- seurat_object@meta.data$stage.ano

genes <- c("CXCL8","IL1B","CCL3L1","FOS","CCL3","CCL4","CCL5","NFKBIA","RGS2","PF4","PF4V1","PPBP","CCL4L2","TNFRSF1B","PTGS2","IL18")

violin_list <- list()

for(x in 1:length(genes)){
  gene <- genes[x]
  expr_df <- data.frame(gene=expr_matrix[,gene],
                        stage.ano=expr_matrix$stage.ano,
                        rename_clusters=expr_matrix$rename_clusters)
  violin <- ggplot(expr_df,aes(x=stage.ano,y=gene))+
    geom_violin(trim = FALSE, aes(fill = factor(stage.ano)))+
    #geom_jitter(aes(fill = factor(res.day),color = factor(res.day)))+
    geom_boxplot(width=0.1,fill='white')+
    #stat_compare_means(method = "anova")+
    scale_fill_manual(values = c("CRS.Before" = "steelblue1","CRS.After"="firebrick1"))+
    scale_color_manual(values = c("CRS.Before" = "steelblue1", "CRS.After"="firebrick1"))+
    theme_bw()+ #背景变为白色
    theme(axis.text.x=element_text(colour="black"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.x=element_text(colour="black"), #设置x轴标题的字体属性
          axis.title.y=element_text(colour="black"), #设置y轴标题的字体属性
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.text=element_text(colour="black"),  #设置图例的子标题的字体属性
          legend.title=element_text(colour="black"), #设置图例的总标题的字体属性
          panel.grid.major = element_blank(),   #不显示网格线
          panel.grid.minor = element_blank())+
    #facet_wrap(~rename_clusters)+
    stat_summary_bin(fun = mean)+
    ylab(gene)
  
  violin_list[[x]] <- violin
}
pdf('Monocyte.TLR.pdf',width = 5,height = 3)

for(x in 1:length(violin_list)){
  print(violin_list[[x]])
}
dev.off()

## CRS Score different monocyte Vlnplot

Idents(Monocyte) <-"celltype"

genes = str_to_upper(c("IL1B","IL6","TNF","CD93","CXCL9","CXCL10","CCL3L1","CXCL2","JUN","IL18","CCL5","PF4","CXCL8"))
genes <- as.data.frame(genes)
genes <- as.list(genes)

Monocyte<- AddModuleScore(Monocyte,features = genes,name = 'CRS')

seurat_object1<-Monocyte [,Monocyte @meta.data$stage.ano %in% c("CRS.After")]

seurat_object2<-seurat_object1 [,seurat_object1 @meta.data$celltype %in% c("MGST1 Monocyte","LGALS2 Monocyte","CD14-CD16- Monocyte","CD16 Monocyte")]

CRS_matrix <- data.frame(seurat_object2@meta.data)

expr_df <- data.frame(gene=CRS_matrix[,"CRS1"],celltype=CRS_matrix$celltype)

violin <- ggplot(expr_df ,aes(x=celltype,y=gene))+ylim(-0.5, 1.5)+
  geom_violin(trim = FALSE, aes(fill = factor(celltype)))+
  #geom_jitter(aes(fill = factor(res.day),color = factor(res.day)))+
  geom_boxplot(width=0.1,fill='white')+
  #stat_compare_means(method = "anova")+
  scale_fill_manual(values = c("MGST1 Monocyte" = "steelblue1","LGALS2 Monocyte"="firebrick1","CD14-CD16- Monocyte"="firebrick4","CD16 Monocyte"="steelblue2"))+
  scale_color_manual(values = c("MGST1 Monocyte" = "steelblue1","LGALS2 Monocyte"="firebrick1","CD14-CD16- Monocyte"="firebrick4","CD16 Monocyte"="steelblue2"))+
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(colour="black",angle = 45), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(colour="black"), #设置x轴标题的字体属性
        axis.title.y=element_text(colour="black"), #设置y轴标题的字体属性
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(colour="black"),  #设置图例的子标题的字体属性
        legend.title=element_text(colour="black"), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+
  #facet_wrap(~rename_clusters)+
  stat_summary_bin(fun = mean)+
  ylab("gene")+
  scale_x_discrete(limits=c("CD14-CD16- Monocyte","CD16 Monocyte","LGALS2 Monocyte","MGST1 Monocyte"))
violin

ggsave("Monocyte.CRS.Vln.pdf",violin,height = 4,width = 5)

### APC pathway violin plot

seurat_object<-Monocyte [,Monocyte @meta.data$stage.ano %in% c("Normal.Before","CRS.Before")]

seurat_object <- seurat_object[,seurat_object@meta.data$seurat_clusters %in% c(0)]

expr_matrix <-as.data.frame(t(as.matrix(seurat_object@assays$RNA@data)))
expr_matrix$rename_clusters <- seurat_object@meta.data$celltype
expr_matrix$stage.ano <- seurat_object@meta.data$stage.ano

genes <- c("PSME2","HLA-DRB1","IFI30","HLA-DMA","TAP1","HLA-DQB1","HLA-DPB1","HLA-DPA1","CD74","TAPBP","PDIA3","HSP90AB1","TAP2","HLA-DRA","HLA-A","HSPA8","CALR","CIITA","HSP90AA1","RFXANK","HLA-F","PSME1","B2M","HLA-C","CTSL","HLA-DQA1","HLA-B","HSPA5")

violin_list <- list()

for(x in 1:length(genes)){
  gene <- genes[x]
  expr_df <- data.frame(gene=expr_matrix[,gene],
                        stage.ano=expr_matrix$stage.ano,
                        rename_clusters=expr_matrix$rename_clusters)
  violin <- ggplot(expr_df,aes(x=stage.ano,y=gene))+
    geom_violin(trim = FALSE, aes(fill = factor(stage.ano)))+
    #geom_jitter(aes(fill = factor(res.day),color = factor(res.day)))+
    geom_boxplot(width=0.1,fill='white')+
    #stat_compare_means(method = "anova")+
    scale_fill_manual(values = c("Normal.Before" = "steelblue1","CRS.Before"="firebrick1"))+
    scale_color_manual(values = c("Normal.Before" = "steelblue1", "CRS.Before"="firebrick1"))+
    theme_bw()+ #背景变为白色
    theme(axis.text.x=element_text(colour="black"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.x=element_text(colour="black"), #设置x轴标题的字体属性
          axis.title.y=element_text(colour="black"), #设置y轴标题的字体属性
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.text=element_text(colour="black"),  #设置图例的子标题的字体属性
          legend.title=element_text(colour="black"), #设置图例的总标题的字体属性
          panel.grid.major = element_blank(),   #不显示网格线
          panel.grid.minor = element_blank())+
    #facet_wrap(~rename_clusters)+
    stat_summary_bin(fun = mean)+
    ylab(gene)
  
  violin_list[[x]] <- violin
}
pdf('Monocyte.APC.pdf',width = 5,height = 3)

for(x in 1:length(violin_list)){
  print(violin_list[[x]])
}
dev.off()

### Response and non Response pathway

PBMC2 <- scn2

Idents(PBMC2) <-"res.ano"

DEGCRS <- FindMarkers(PBMC2,ident.1 = "Response.Inter",ident.2 = "NR.Inter",min.pct = 0.01,logfc.threshold = 0.25)

DefaultAssay(PBMC2) <-"RNA"

genelist <- DEGCRS$avg_log2FC

names(genelist) = toupper(rownames(DEGCRS))

genelist <-sort(genelist,decreasing = T)

geneset<-msigdbr(species = "Homo sapiens",category = "C2",subcategory = "KEGG")

geneset<-geneset[,c("gs_name","gene_symbol")]

length(unique(geneset$gs_name))

egmt <- GSEA(genelist,TERM2GENE = geneset,minGSSize = 1,pvalueCutoff = 0.5)

mat <- data.frame(egmt)

mat <- mat[which(mat$pvalue < 0.05),]

mat <- mat[order(mat$NES,decreasing = T),]

write.csv(mat,"GSEA.MGST1.Normal.Before.CRS.Before.csv")

gseaplot2(egmt,geneSetID = c(2),pvalue_table = T,color = "red",base_size = 10,ES_geom = "line")

ggsave("APC_RES_vs_NorRES.pdf",width = 4.5,height = 4)

gseaplot2(egmt,geneSetID = c(16),pvalue_table = T,color = "red",base_size = 10,ES_geom = "line")

ggsave("TLR_RES_vs_NorRES.pdf",width = 4.5,height = 4)
### 画图
### Pathway 画图

plot_need<- mat[c(1:8,24:31),]

plot_need$pathway <- as.factor(plot_need$ID)
plot_need$threshold <- as.factor(ifelse(plot_need$NES>=0 ,'Up','Down'))
up <- nrow(plot_need[which(plot_need$threshold=='Up'),])
down <- nrow(plot_need[which(plot_need$threshold=='Down'),])

### GSEA 好看的图片展示
gsea_barplot <- ggplot(plot_need,aes(x = reorder(pathway,NES),y = NES))+
  geom_bar(aes(fill=threshold),stat='identity') + 
  scale_fill_manual(values = c('Up'='firebrick3','Down'='navy'))+
  coord_flip()+
  guides(fill='none')+
  xlab('')+
  ylab('')+
  labs(title = "")+
  geom_text(data = plot_need[1:up,],aes(x=pathway,y=-4,label=pathway),hjust=0,color='black',size = 3)+
  geom_text(data = plot_need[(up+1):(up+down),],aes(x=pathway,y=4,label=pathway),hjust=1,color='black',size = 3)+
  theme_bw()+ #鑳屾櫙鍙樹负鐧借壊
  theme(panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"), #鍘婚櫎榛樿濉厖鐨勭伆鑹诧紝骞跺皢x=0杞村拰y=0杞村姞绮楁樉绀?(size=1)
        panel.grid.major = element_blank(),   #涓嶆樉绀虹綉鏍肩嚎
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
gsea_barplot

ggsave("gsea_Monocyte_KEGG_Res.pdf",gsea_barplot,width = 8,height = 12)

seurat_object <- PBMC2

expr_matrix <-as.data.frame(t(as.matrix(seurat_object@assays$RNA@data)))
expr_matrix$rename_clusters <- seurat_object@meta.data$celltype
expr_matrix$stage.ano <- seurat_object@meta.data$res.ano

genes <- c("HLA-DRB5","HLA-DRB1","IFI30","HLA-DMA","TAP1","HLA-DQB1","HLA-DPB1","HLA-DPA1","CD74","TAPBP","PDIA3","HSP90AB1","TAP2","HLA-DRA","HLA-A","HSPA8","CALR","CIITA","HSP90AA1","RFXANK","HLA-F","PSME1","B2M","HLA-C","CTSL","HLA-DQA1","HLA-B","HSPA5","GZMB")

violin_list <- list()

for(x in 1:length(genes)){
  gene <- genes[x]
  expr_df <- data.frame(gene=expr_matrix[,gene],
                        stage.ano=expr_matrix$stage.ano,
                        rename_clusters=expr_matrix$rename_clusters)
  violin <- ggplot(expr_df,aes(x=stage.ano,y=gene))+
    geom_violin(trim = FALSE, aes(fill = factor(stage.ano)))+
    #geom_jitter(aes(fill = factor(res.day),color = factor(res.day)))+
    geom_boxplot(width=0.1,fill='white')+
    #stat_compare_means(method = "anova")+
    scale_fill_manual(values = c("Normal.Before" = "steelblue1","CRS.Before"="firebrick1"))+
    scale_color_manual(values = c("Normal.Before" = "steelblue1", "CRS.Before"="firebrick1"))+
    theme_bw()+ #背景变为白色
    theme(axis.text.x=element_text(colour="black"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.x=element_text(colour="black"), #设置x轴标题的字体属性
          axis.title.y=element_text(colour="black"), #设置y轴标题的字体属性
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.text=element_text(colour="black"),  #设置图例的子标题的字体属性
          legend.title=element_text(colour="black"), #设置图例的总标题的字体属性
          panel.grid.major = element_blank(),   #不显示网格线
          panel.grid.minor = element_blank())+
    #facet_wrap(~rename_clusters)+
    stat_summary_bin(fun = mean)+
    ylab(gene)
  
  violin_list[[x]] <- violin
}
pdf('Monocyte.Res.APC.pdf',width = 5,height = 3)

for(x in 1:length(violin_list)){
  print(violin_list[[x]])
}
dev.off()

### 提取VCAN，DN，LGSL2 and CD16 for pseudotime

### cell anatotaion

scn3 <- Monocyte [,Monocyte @meta.data$celltype %in% c("MGST1 Monocyte","LGALS2 Monocyte","CD14-CD16- Monocyte","CD16 Monocyte")]

scn3  <-FindNeighbors(scn3 , reduction = "harmony",dims = 1:20,k.param = 40) %>% FindClusters(resolution=0.5,n.iter = 50)

scn3  <- RunTSNE(scn3 ,reduction = "harmony",dims = 1:20,n.neighbors = 20, min.dist = 0.2)

### cell type tsne 

plot1 <- DimPlot(scn3 , pt.size = 0.2, cols = mycol, label=T, repel=T,
                 raster = FALSE, label.size = 5, reduction = "tsne",group.by = "celltype")
plot1

ggsave("Pesudo_Monocyte_cell_type.pdf",plot1,height = 5,width = 7)

### 绘制APC及TOLL pathway的分布

all_gene_sets <-msigdbr(species = "Homo sapiens",category = "C2")

all_gene_sets <- all_gene_sets[all_gene_sets$gs_exact_source %in% c("hsa04612","hsa04060","R-HSA-168898","hsa04620"),]

av <- as.matrix(scn3@assays$RNA@data)

length(unique(table(all_gene_sets$gs_name)))

tail(table(all_gene_sets$gs_name))

tail(sort(table(all_gene_sets$gs_name)))

gs <- split(all_gene_sets$gene_symbol,all_gene_sets$gs_name)

gs <- lapply(gs,unique)

gsc <- GSEABase::GeneSetCollection(mapply(function(geneIds,keggId){
  GeneSet(geneIds,geneIdType=EntrezIdentifier(),collectionType=KEGGCollection(keggId),setName=keggId)},gs,names(gs)))
gsc

geneset<-gsc

X=data.frame(av)

X<- as.matrix(X)

es.max <- gsva(expr = X,geneset,kcdf="Gaussian",verbose=T,parallel.sz=8)

es <- data.frame(t(es.max),stringsAsFactors = F)

colnames(es) <- c("APC","Cytokine","TOLL","RTOLL")

### Metadata转入scn3

samples <- intersect(rownames(es),colnames(scn3))

gsub("\\.","-",rownames(es))

samples <- intersect(gsub("\\.","-",rownames(es)),colnames(scn3))

rownames(es) <- gsub("\\.","-",rownames(es))

scn3<-AddMetaData(scn3,es)

write_rds(scn3,"Monocyte.Monocle.scn3.rds")

### Metadata转入scn3

mycolor <- c('lightgrey',"red")#设置颜色  

FeaturePlot(scn3,features = c("APC"),min.cutoff = 0, max.cutoff = 0.1,cols = mycolor,pt.size = 0.2,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("APC_t_sne.pdf", width = 5.5, height = 5)

FeaturePlot(scn3,features = c("RTOLL"),min.cutoff = -0.2, max.cutoff = -0.1,cols = mycolor,pt.size = 0.2,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("TOLL_t_sne.pdf", width = 5.5, height = 5)


FeaturePlot(scn3,features = c("TOLL"),min.cutoff = -0.2, max.cutoff = -0.1,cols = mycolor,pt.size = 0.2,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("TOLL_t_sne.pdf", width = 5.5, height = 5)
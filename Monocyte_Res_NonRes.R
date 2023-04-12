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

### re-analysis monocyte


################# 颜色设定================================================

mycol <- c(brewer.pal(5,"Set1"), brewer.pal(8,"Set2"),
           brewer.pal(11,"Set3"), brewer.pal(12,"Paired"),
           brewer.pal(8,"Accent"), brewer.pal(11,"Spectral"),
           brewer.pal(11,"BrBG"), brewer.pal(11,"PiYG"),
           brewer.pal(11,"PuOr"),brewer.pal(11,"RdBu"))

##################

scrna_harmony<-read_rds("Final.scrna.harmony.rds")

scn<- read_rds("Res.Mono.rds")

scn <- scn[,scn@meta.data$seurat_clusters %in% c("0","1","2","3")]

scn <-FindNeighbors(scn, reduction = "harmony",dims = 1:30,k.param = 40) %>% FindClusters(resolution=0.5,n.iter = 50)

scn<- RunTSNE(scn,reduction = "harmony",dims = 1:30,n.neighbors = 20, min.dist = 0.2)

scn<- RunUMAP(scn,reduction = "harmony",dims = 1:20,n.neighbors = 20, min.dist = 0.5)

plot1 <- DimPlot(scn, pt.size = 0.2, cols = mycol, label=T, repel=T,
                 raster = FALSE, label.size = 5, reduction = "tsne",group.by = "seurat_clusters")
plot1

write_rds(scn,file = "Res.Mono.rds")


#### cell anatation

scn@meta.data$celltype="other"
scn@meta.data$celltype[scn@meta.data$seurat_clusters=="0"]="LGALS2 Monocyte"
scn@meta.data$celltype[scn@meta.data$seurat_clusters=="1"]="VCAN Monocyte"
scn@meta.data$celltype[scn@meta.data$seurat_clusters=="2"]="DN Monocyte"
scn@meta.data$celltype[scn@meta.data$seurat_clusters=="3"]="CD16 Monocyte"


### TOP5 marker gene heatmap cell gene

seurat_object <- scn
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

####

mycolor <- c('lightgrey',"red")#设置颜色  

FeaturePlot(scn,features = c("LGALS2"),min.cutoff = 0.3, max.cutoff = 1,cols = mycolor,pt.size = 0.2,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

FeaturePlot(scn,features = c("FCGR3A"),min.cutoff = 0.5, max.cutoff = 2,cols = mycolor,pt.size = 0.2,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

FeaturePlot(scn,features = c("FCGR3B"),min.cutoff = 0.5, max.cutoff = 2,cols = mycolor,pt.size = 0.2,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

FeaturePlot(scn,features = c("CD1C"),min.cutoff = 0.5, max.cutoff = 2,cols = mycolor,pt.size = 0.2,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

FeaturePlot(scn,features = c("VCAN"),min.cutoff = 3.0, max.cutoff = 3.3,cols = mycolor,pt.size = 0.2,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

FeaturePlot(scn,features = c("CD14"),min.cutoff = 1.5, max.cutoff = 2,cols = mycolor,pt.size = 0.2,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

### stack plot

scn@meta.data$Res="other"
scn@meta.data$Res[scn@meta.data$res.ano=="SD.Inter"]="Response.After"
scn@meta.data$Res[scn@meta.data$res.ano=="SD.D28"]="Response.After"
scn@meta.data$Res[scn@meta.data$res.ano=="SD.Before"]="Response.Before"
scn@meta.data$Res[scn@meta.data$res.ano=="PR.D28"]="Response.After"
scn@meta.data$Res[scn@meta.data$res.ano=="NR.Inter"]="NonResponse.After"
scn@meta.data$Res[scn@meta.data$res.ano=="NR.D28"]="NonResponse.After"
scn@meta.data$Res[scn@meta.data$res.ano=="NR.Before"]="NonResponse.Before"

Idents(scn)=scn@meta.data$celltype #y轴

tab<-table(Idents(scn),scn$Res)

tab

tab<-prop.table(tab,2)*100

tab ###1是行，二是列

tab<-as.data.frame(tab)

write.csv2(tab,"Monocyte_Res_stage.ano_cellpro.csv")

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

ggsave("Monocyte_stage.ano_cellpro.pdf", width = 7, height = 3)

####

PBMC2 <- scn[,scn@meta.data$celltype %in% c("CD16 Monocyte")]

Idents(PBMC2) <-"Res"

DEGCRS <- FindMarkers(PBMC2,ident.1 = "Response.After",ident.2 = "NonResponse.After",min.pct = 0.01,logfc.threshold = 0.25)

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

write.csv2(mat,"CD16.Res_NonRes.csv")

gseaplot2(egmt,geneSetID = c(2),pvalue_table = T,color = "red",base_size = 10,ES_geom = "line")

ggsave("TLR_CRSB_vs_NorB_-2.75.pdf",width = 4.5,height = 4)

gseaplot2(egmt,geneSetID = c(23),pvalue_table = T,color = "red",base_size = 10,ES_geom = "line")

ggsave("APC_CRSB_vs_NorB_1.88.pdf",width = 4.5,height = 4)

### APC and TOLL 

scn3<-scn

all_gene_sets <-msigdbr(species = "Homo sapiens",category = "C2")

all_gene_sets <- all_gene_sets[all_gene_sets$gs_exact_source %in% c("hsa04612","R-HSA-168898"),]

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

scn3<-readRDS("Monocyte.Monocle.scn3.rds")

### Metadata转入scn3

mycolor <- c('lightgrey',"red")#设置颜色  

FeaturePlot(scn3,features = c("APC"),min.cutoff = 0, max.cutoff = 0.1,cols = mycolor,pt.size = 0.2,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("APC_t_sne.pdf", width = 5.5, height = 5)

FeaturePlot(scn3,features = c("RTOLL"),min.cutoff = -0.2, max.cutoff = -0.1,cols = mycolor,pt.size = 0.2,reduction = "tsne")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("TOLL_t_sne.pdf", width = 5.5, height = 5)

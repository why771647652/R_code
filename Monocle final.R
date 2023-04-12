
library("Seurat")
library("dplyr")
library("tidyverse")
library("patchwork")
library("monocle")

HSMM<-read_rds("Mono.HSMM.rds")

###HSMM<-read_rds("T.HSMM.rds")

scrna<-read_rds("Monocyte.Monocle.scrna.rds")

scrna<- subset(scrna)

Idents(scrna) <-"celltype"

scrna1 <- scrna

as.data.frame(table(Idents(scrna1)))

scrna1<- subset(scrna1)

Idents(scrna1) <-"celltype"

marker_gene<- FindAllMarkers(scrna1,only.pos = FALSE,min.pct = 0.25,logfc.threshold = 0.25)

### 构建matrix，单细胞count矩阵
data <- as(as.matrix(scrna1@assays$RNA@counts),"sparseMatrix")

pd <- new("AnnotatedDataFrame",data=scrna1@meta.data)

fData <- data.frame(gene_short_name=row.names(data),row.names = row.names(data))

fd <- new("AnnotatedDataFrame",data=fData)

monocle_cds <-newCellDataSet(data,
                             phenoData = pd,
                             featureData = fd,
                             lowerDetectionLimit = 0.5,
                             expressionFamily = negbinomial.size())

###2.2 估计size factors 和 dispersions （需要）

monocle_cds <-estimateSizeFactors(monocle_cds)

monocle_cds <-estimateDispersions(monocle_cds)

###. 过滤低质量细胞
monocle_cds <- detectGenes(monocle_cds,min_expr = 3)

monocle_cds<-monocle_cds[fData(monocle_cds)$num_cells_expressed>10,]

### Gene 筛选
HSMM=monocle_cds

### 使用Seurat差异gene进行分析
marker_gene1 <- marker_gene[abs(marker_gene$avg_log2FC)>0.1,]

marker_gene1 <- subset(marker_gene,p_val<0.05)

marker_gene1 <-marker_gene1[order(marker_gene1$avg_log2FC),]

test_ordering_genes<-unique(marker_gene1$gene)

HSMM<- setOrderingFilter(HSMM,ordering_genes = test_ordering_genes)

plot_ordering_genes(HSMM)

HSMM <- reduceDimension(HSMM, max_components = 2,method="DDRTree")

HSMM <- orderCells(HSMM)

### root 更换图片起点的根

### 使用Seurat差异gene进行分析

marker_gene1 <- marker_gene[abs(marker_gene$avg_log2FC)>0.1,]

marker_gene1 <- subset(marker_gene,p_val<0.05)

marker_gene1 <-marker_gene1[order(marker_gene1$avg_log2FC),]

test_ordering_genes<-unique(marker_gene1$gene)


plot_cell_trajectory(HSMM,color_by = "celltype")+scale_color_nejm()
ggsave("celltype.pdf",width = 5,height = 5)

plot_cell_trajectory(HSMM,color_by = "State") +scale_color_nejm()
ggsave("state.pdf",width = 5,height = 5)

plot_cell_trajectory(HSMM,color_by = "Pseudotime")+scale_color_gsea()
ggsave("pesudotime.pdf",width = 5,height = 5)

plot_cell_trajectory(HSMM,color_by = "RTOLL")+scale_color_gsea()
ggsave("pesudotime.TOLL.pdf",width = 5,height = 5)

plot_cell_trajectory(HSMM,color_by = "TLR4")+scale_color_gsea()
ggsave("TLR4.pdf",width = 5,height = 5)

plot_cell_trajectory(HSMM,color_by = "APC")+scale_color_gsea()
ggsave("pesudotime.APC.pdf",width = 5,height = 5)
colnames(pData(HSMM))

### cluster差异gene

diff.genes <- marker_gene

sig_diff.genes <- subset(diff.genes,p_val<0.0001&abs(avg_log2FC)>0.8)

sig_diff.genes <- unique(as.character(sig_diff.genes$gene))

marker_genes <- row.names(subset(fData(HSMM),gene_short_name %in% sig_diff.genes))

diff_test_res <- differentialGeneTest(HSMM[marker_genes,],fullModelFormulaStr = "~sm.ns(Pseudotime)")

sig_gene_names <- rownames(subset(diff_test_res,qval<0.1))

sig_gene_names100 <-top_n(diff_test_res,n=50,desc(qval)) %>% pull(gene_short_name) %>% as.character()

P <- plot_pseudotime_heatmap(HSMM[sig_gene_names100,],
                             num_clusters = 4,
                             cores = 1,
                             show_rownames = T,return_heatmap = T)
P

ggsave("Monocyte_Time_HeatmapTop1001.pdf",P,height = 8,width = 5)


### 单基因变化

pData(HSMM)$KLF6 <- log2(exprs(HSMM)["KLF6",]+1)

p1=plot_cell_trajectory(HSMM,color_by = "KLF6")+scale_color_gsea()

p1

pData(HSMM)$LST1<- log2(exprs(HSMM)["LST1",]+1)

p1=plot_cell_trajectory(HSMM,color_by = "LST1")+scale_color_gsea()

p1

pData(HSMM)$CD74<- log2(exprs(HSMM)["CD74",]+1)

p1=plot_cell_trajectory(HSMM,color_by = "CD74")+scale_color_gsea()

p1

ggsave("CD74.pdf",p1,width = 5,height = 5)

pData(HSMM)$HLA_A<- log2(exprs(HSMM)["HLA-A",]+1)

p1=plot_cell_trajectory(HSMM,color_by = "HLA_A")+scale_color_gsea()

p1

ggsave("HLA_A.pdf",p1,width = 5,height = 5)

pData(HSMM)$HLA_DRA<- log2(exprs(HSMM)["HLA-DRA",]+1)

p1=plot_cell_trajectory(HSMM,color_by = "HLA_DRA")+scale_color_gsea()

p1

ggsave("HLA_DRA.pdf",p1,width = 5,height = 5)

pData(HSMM)$TLR2<- log2(exprs(HSMM)["TLR2",]+1)

p1=plot_cell_trajectory(HSMM,color_by = "SPP1")+scale_color_gsea()

p1

ggsave("HLA_A.pdf",p1,width = 5,height = 5)


library("Seurat")
library("tidyverse")
library("patchwork")
library("dplyr")
library("harmony")
library("ggrepel")
library("ggplot2")
library("ggpubr")
library("cowplot")
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


exp <- read.csv("FPKM_Rest.csv")

mrna<-read.table("mRNA_list.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",fill = T)

mrna_list<-c(mrna$gene_id)

exp1<-exp[exp$geneid %in% unlist(mrna_list),]
  
exp1=na.omit(exp1)

exp2=distinct(exp1,Symbol,.keep_all = T)

rownames(exp2)<-exp2$Symbol


exp3<-exp2[,-c(1:3)]


expMatrix<-exp3
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}

tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,]

colSums(tpms)

exp3<-tpms

a<-c(1,2,3)
a2<-c(4,5,6)

exp4<-exp3[,c(a,a2)]

group_list<-c(rep("shCtrl",3),rep("TKD",3))

group_list<-factor(group_list)

group_list<-relevel(group_list,ref = "shCtrl")

library(limma)

exp5<-normalizeBetweenArrays(exp4)

boxplot(exp4,col=group_list,las=2,outline=FALSE)

boxplot(exp5,col=group_list,las=2,outline=FALSE)

#### limma

data<-exp5  
data <- log2(data+1)

group=c(rep("shCtrl",3),rep("TKD",3))   

design <- model.matrix(~0+factor(group))
colnames(design)= levels(factor(group))
rownames(design)=colnames(data)

contrast.matrix <- makeContrasts(TKD-shCtrl,levels =design)

paste0(unique(group),collapse = "-")

fit<- lmFit(data,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2,coef = 1,n=Inf)
allDiff1 = na.omit(tempOutput)

write.table(allDiff1,file="Rest_DEG.txt",sep="\t",quote=F)

allDiff<-allDiff1[allDiff1$AveExpr>0.1,]

library(ggplot2)                         

library(ggrepel)                               

data<-allDiff
data$significant="stable"
log2(1.5)
data$significant[data$logFC>=0.58 & data$P.Value <0.01]="up"
data$significant[data$logFC<= -0.58 & data$P.Value <0.01]="down"

label.deg<-c("ID3","CCR7","BACH2","IL7R","TOX","BCL6","BCL2","BCL2L1","BCL2L11","PDCD1","LAG3","ITGA4","EOMES","ITGA4","ITGAL","SELPLG")

p<-ggplot(data,aes(logFC,-1*log10(P.Value)))+xlim(-5,5)+ylim(0,6)+
  geom_point(aes(color=significant),size=0.8)+theme_classic()+
  scale_color_manual(values = c("#2878B5","#E0E0E0","#BB2525"))+
  geom_hline(yintercept = 2,linetype=4,size=0.3)+
  geom_vline(xintercept = c(-0.58,0.58),linetype=4,size=0.3)+
  theme(title=element_text(size = 18),text = element_text(size=18))+
  labs(x="log2(foldchange)",y="-log10(p_value)")

data_selected <- data[label.deg,]

p1<-p + geom_label_repel(data=data_selected,
                     aes(label=rownames(data_selected)))
p1

ggsave("rest_vlocano_plot.pdf",p1,height = 5,width = 7)

annotation_col1 = data.frame(CellType =c(rep("CAR shCtrl",3),rep("CAR TKD",3)) )
rownames(annotation_col1)=colnames(exp5)
exp6=as.data.frame(exp5)

exprSet.map=exp6[c("ID3","CCR7","SLAMF6","LEF1","BACH2","TCF7","IL7R","LY6E","BCL6","BCL2","BCL2L1","BCL2L11","BAK","CASP3","ZAP70","LCK","CD69","FAS","TOX","PDCD1","LAG3","CTLA4","CD96","EOMES"),]
####聚类
pheatmap::pheatmap(exprSet.map, #热图的数据
                   cluster_rows =F,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   annotation_col =annotation_col1,
                   show_colnames=T,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("#035490", "white","#C25454"))(100))

####聚类
p<-pheatmap::pheatmap(exprSet.map, #热图的数据
                   cluster_rows =F,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   show_colnames=T,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("#035490", "white","#C25454"))(100))
p

ggsave("heatmap_resting.pdf",p,height = 5,width = 4)

exprSet.map=exp6[c("ITGA4","ITGAL","SELPLG","ANXA8","GSTM3","MT1H","KCNK1","TUBB6","BST2","MT1E","FABP5","S100A6"),]
####聚类
pheatmap::pheatmap(exprSet.map, #热图的数据
                   cluster_rows =F,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   annotation_col =annotation_col1,
                   show_colnames=T,
                   scale = "row", #以行来标准化，这个功能很不错
                   color =colorRampPalette(c("#2878B5", "white","#BB2525"))(100))

p<-pheatmap::pheatmap(exprSet.map, #热图的数据
                      cluster_rows =F,#行聚类
                      cluster_cols =T,#列聚类，可以看出样本之间的区分度
                      show_colnames=T,
                      scale = "row", #以行来标准化，这个功能很不错
                      color =colorRampPalette(c("#006cbe", "white","#c60000"))(100))


library(msigdbr)
all_gene_sets = msigdbr(species = "Homo sapiens",category='C2')
length(unique(table(all_gene_sets$gs_name)))
tail(table(all_gene_sets$gs_name))
head(all_gene_sets[,1:5])

gcSample = split(all_gene_sets$gene_symbol,
                 all_gene_sets$gs_name)
names(gcSample)
lapply(gcSample[1:3],function(x)head(x,5))

names(gcSample)

file="sink-examp.txt"
gs=gcSample
write.gmt <- function(gs,file){
  sink(file)
  lapply(names(gs), function(i){
    cat( paste(c(i,'tmp',gs[[i]]),collapse='\t') )
    cat('\n')
  })
  sink()
}
write.gmt(gs,file) 

### geneset
geneset<-all_gene_sets[,c("gs_name","gene_symbol")]

length(unique(geneset$gs_name))


geneList <- allDiff1$logFC
names(geneList) <- rownames(allDiff1)
geneList <- sort(geneList, decreasing = T)
geneList <- geneList[geneList != 0]
head(geneList)

library(clusterProfiler)
set.seed(123456)

egmt <- GSEA(geneList,TERM2GENE = geneset,minGSSize = 1,pvalueCutoff = 0.5)

mat <- data.frame(egmt)

mat <- mat[which(mat$pvalue < 0.05),]

mat <- mat[order(mat$NES,decreasing = T),]

write.csv(egmt@result,"GSEA.C2.TKDvsshctrl.Before.csv")

gseaplot2(egmt,geneSetID = c(693),pvalue_table = T,color = "red",base_size = 10,ES_geom = "line")

gseaplot2(egmt,geneSetID = c(1626),pvalue_table = T,color = "red",base_size = 10,ES_geom = "line")

gseaplot2(egmt,geneSetID = c(1250),pvalue_table = T,color = "red",base_size = 10,ES_geom = "line")

gseaplot2(egmt,geneSetID = c(24),pvalue_table = T,color = "red",base_size = 10,ES_geom = "line")

gseaplot2(egmt,geneSetID = c(3901),pvalue_table = T,color = "red",base_size = 10,ES_geom = "line")

gseaplot2(egmt,geneSetID = c(1023),pvalue_table = T,color = "red",base_size = 10,ES_geom = "line")





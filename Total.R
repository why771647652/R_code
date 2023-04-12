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

################# 颜色设定================================================

mycol <- c(brewer.pal(5,"Set1"), brewer.pal(8,"Set2"),
           brewer.pal(11,"Set3"), brewer.pal(12,"Paired"),
           brewer.pal(8,"Accent"), brewer.pal(11,"Spectral"),
           brewer.pal(11,"BrBG"), brewer.pal(11,"PiYG"),
           brewer.pal(11,"PuOr"),brewer.pal(11,"RdBu"))

#################=========================================================

scrna_harmony<-read_rds("Final.scrna.harmony.rds")

plot1 <- DimPlot(scrna_harmony, pt.size = 0.2, cols = mycol, label=T, repel=T,
                 raster = FALSE, label.size = 5, reduction = "umap",
                 group.by = "celltype")
plot1

ggsave("TotalCelltype.pdf", plot1, width = 6, height = 5)

plot2 <- DimPlot(scrna_harmony, pt.size = 0.2, cols = mycol, label=F, repel=T,
                 raster = FALSE, label.size = 5, reduction = "umap",
                 group.by = "pat.ano")
plot2

ggsave("Totalpatano.pdf", plot2, width = 6.5, height = 5)

plot2 <- DimPlot(scrna_harmony, pt.size = 0.2, cols = mycol, label=F, repel=T,
                 raster = FALSE, label.size = 5, reduction = "umap",
                 group.by = "orig.ident")
plot2

plot3 <- DimPlot(scrna_harmony, pt.size = 0.2, cols = mycol, label=F, repel=T,
                 raster = FALSE, label.size = 5, reduction = "umap",
                 group.by = "stage.ano")
plot3

ggsave("Totalstageano.pdf", plot3, width = 7, height = 5)

### orig.ano 注释stage
scrna_harmony@meta.data$stage.ano="other"
scrna_harmony$stage.ano[scrna_harmony$orig.ident=="Patient1.Before"]="NonCRS.Before"
scrna_harmony$stage.ano[scrna_harmony$orig.ident=="Patient1.Normal"]="NonCRS"
scrna_harmony$stage.ano[scrna_harmony$orig.ident=="Patient1.Remission"]="NonCRS.Remission"
scrna_harmony$stage.ano[scrna_harmony$orig.ident=="Patient2.Before"]="CRS.Before"
scrna_harmony$stage.ano[scrna_harmony$orig.ident=="Patient2.CRS"]="CRS"
scrna_harmony$stage.ano[scrna_harmony$orig.ident=="Patient2.Remission"]="CRS.Remission"
scrna_harmony$stage.ano[scrna_harmony$orig.ident=="Patient3.Before"]="CRS.Before"
scrna_harmony$stage.ano[scrna_harmony$orig.ident=="Patient3.CRS"]="CRS"
scrna_harmony$stage.ano[scrna_harmony$orig.ident=="Patient3.Remission"]="CRS.Remission"
scrna_harmony$stage.ano[scrna_harmony$orig.ident=="Patient4.Before"]="CRS.Before"
scrna_harmony$stage.ano[scrna_harmony$orig.ident=="Patient4.CRS"]="CRS"
scrna_harmony$stage.ano[scrna_harmony$orig.ident=="Patient5.Remission"]="CRS.Remission"

### patient 注释
scrna_harmony@meta.data$pat.ano="other"
scrna_harmony$pat.ano[scrna_harmony$orig.ident=="Patient1.Before"]="Patient1"
scrna_harmony$pat.ano[scrna_harmony$orig.ident=="Patient1.Normal"]="Patient1"
scrna_harmony$pat.ano[scrna_harmony$orig.ident=="Patient1.Remission"]="Patient1"
scrna_harmony$pat.ano[scrna_harmony$orig.ident=="Patient2.Before"]="Patient2"
scrna_harmony$pat.ano[scrna_harmony$orig.ident=="Patient2.CRS"]="Patient2"
scrna_harmony$pat.ano[scrna_harmony$orig.ident=="Patient2.Remission"]="Patient2"
scrna_harmony$pat.ano[scrna_harmony$orig.ident=="Patient3.Before"]="Patient3"
scrna_harmony$pat.ano[scrna_harmony$orig.ident=="Patient3.CRS"]="Patient3"
scrna_harmony$pat.ano[scrna_harmony$orig.ident=="Patient3.Remission"]="Patient3"
scrna_harmony$pat.ano[scrna_harmony$orig.ident=="Patient4.Before"]="Patient4"
scrna_harmony$pat.ano[scrna_harmony$orig.ident=="Patient4.CRS"]="Patient4"
scrna_harmony$pat.ano[scrna_harmony$orig.ident=="Patient5.Remission"]="Patient5"

### Response
scrna_harmony@meta.data$res.ano="other"
scrna_harmony$res.ano[scrna_harmony$orig.ident=="Patient1.Before"]="SD.Before"
scrna_harmony$res.ano[scrna_harmony$orig.ident=="Patient1.Normal"]="SD.Inter"
scrna_harmony$res.ano[scrna_harmony$orig.ident=="Patient1.Remission"]="SD.D28"
scrna_harmony$res.ano[scrna_harmony$orig.ident=="Patient2.Before"]="NR.Before"
scrna_harmony$res.ano[scrna_harmony$orig.ident=="Patient2.CRS"]="NR.Inter"
scrna_harmony$res.ano[scrna_harmony$orig.ident=="Patient2.Remission"]="NR.D28"
scrna_harmony$res.ano[scrna_harmony$orig.ident=="Patient3.Before"]="NR.Before"
scrna_harmony$res.ano[scrna_harmony$orig.ident=="Patient3.CRS"]="NR.Inter"
scrna_harmony$res.ano[scrna_harmony$orig.ident=="Patient3.Remission"]="NR.D28"
scrna_harmony$res.ano[scrna_harmony$orig.ident=="Patient4.Before"]="NR.Before"
scrna_harmony$res.ano[scrna_harmony$orig.ident=="Patient4.CRS"]="NR.Inter"
scrna_harmony$res.ano[scrna_harmony$orig.ident=="Patient5.Remission"]="PR.D28"

### 定义CRS related gene set

genes = str_to_upper(c("IL1B","IL6","TNF","CD93","CXCL9","CXCL10","CCL3L1","CXCL2","JUN","IL18","CCL5","PF4","CXCL8"))

genes <- as.data.frame(genes)

genes <- as.list(genes)

scrna_harmony<- AddModuleScore(scrna_harmony,features = genes,name = 'CRS')

colnames(scrna_harmony@meta.data)

scn1 <- scrna_harmony[,scrna_harmony@meta.data$stage.ano %in% c("CRS")]

FeaturePlot(scn1,features = c("CRS1"),min.cutoff = 0, max.cutoff = 0.5,cols = mycolor,pt.size = 0.2)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("CRSscorefea.pdf", width = 5.5, height = 5)

### cytokine 表达图谱

mycolor <- c('lightgrey',"red")#设置颜色  

FeaturePlot(scn1,features = c("IFNG"),min.cutoff = 0, max.cutoff = 1,cols = mycolor,pt.size = 0.2)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("IFNG.pdf", width = 5, height = 5)

FeaturePlot(scn1,features = c("IL1B"),min.cutoff = 0.2, max.cutoff = 1,cols = mycolor,pt.size = 0.2)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("IL1B.pdf", width = 5, height = 5)

FeaturePlot(scn1,features = c("CXCL8"),min.cutoff = 1, max.cutoff = 3,cols = mycolor,pt.size = 0.2)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("CXCL8.pdf", width = 5, height = 5)

FeaturePlot(scn1,features = c("CXCL2"),min.cutoff = 0, max.cutoff = 0.1,cols = mycolor,pt.size = 0.2)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("CXCL2.pdf", width = 5, height = 5)

FeaturePlot(scn1,features = c("IL18"),min.cutoff = 0, max.cutoff = 0.1,cols = mycolor,pt.size = 0.2)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("IL18.pdf", width = 5, height = 5)

FeaturePlot(scn1,features = c("IL6"),min.cutoff = 0, max.cutoff = 0.1,cols = mycolor,pt.size = 0.2)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("IL6.pdf", width = 5, height = 5)

FeaturePlot(scn1,features = c("CCL5"),min.cutoff = 0, max.cutoff = 0.1,cols = mycolor,pt.size = 0.2)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("CCL2.pdf", width = 5, height = 5)

FeaturePlot(scn1,features = c("TNF"),min.cutoff = 0, max.cutoff = 0.1,cols = mycolor,pt.size = 0.2)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsave("TNF.pdf", width = 5, height = 5)

## CRS Score different monocyte Vlnplot

seurat_object1<-scrna_harmony[,scrna_harmony @meta.data$stage.ano %in% c("CRS")]

CRS_matrix <- data.frame(seurat_object1@meta.data)

expr_df <- data.frame(gene=CRS_matrix[,"CRS1"],celltype=CRS_matrix$celltype)

violin <- ggplot(expr_df ,aes(x=celltype,y=gene))+ylim(-0.5, 1.5)+
  geom_violin(trim = FALSE, aes(fill = factor(celltype)))+
  #geom_jitter(aes(fill = factor(res.day),color = factor(res.day)))+
  geom_boxplot(width=0.1,fill='white')+
  #stat_compare_means(method = "anova")+
  scale_fill_manual(values = c("MGST1 Monocyte" = "steelblue1","LGALS2 Monocyte"="firebrick1","Neutrophils"="firebrick2","mDC"="firebrick3","CD14-CD16- Monocyte"="firebrick4"))+
  scale_color_manual(values = c("MGST1 Monocyte" = "steelblue1","LGALS2 Monocyte"="firebrick1","Neutrophils"="firebrick2","mDC"="firebrick3","CD14-CD16- Monocyte"="firebrick4"))+
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
  scale_x_discrete(limits=c("B cell","T cell","mDC","NK cell","NKT cell","Platelet","Neutrophil","Monocyte"))
violin

ggsave("Total.CRS.Vln.pdf",violin,height = 4,width = 6)


### CRS gene vln

library(MySeuratWrappers) 

scn <- scrna_harmony[,scrna_harmony@meta.data$stage.ano %in% c("CRS.Before","CRS")]

VlnPlot(scn, features =c("CRS1"), split.by="stage.ano",stacked=T,pt.size=0,cols = my36colors,direction = "vertical", x.lab = '', y.lab = '')+theme()#不显示坐标刻度 

ggsave("CRSvln.pdf", width = 6, height = 4)

####. scrna_harmony stack plot

scn <- scrna_harmony[,scrna_harmony@meta.data$stage.ano %in% c("CRS.Before","CRS","NonCRS","NonCRS.Before")]

Idents(scn)=scn@meta.data$celltype #y轴

tab<-table(Idents(scn),scn$stage.ano)

tab

tab<-prop.table(tab,2)*100

tab ###1是行，二是列

tab<-as.data.frame(tab)

write.csv2(tab,"Cell_propotion.csv")

ggplot(tab, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = pal_nejm()(8))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.9,hjust = 0.9))+
  xlab("")+ylab("percentage")

ggsave("Total_stage.ano_cellpro.pdf", width = 5, height = 3)


plot3 <- DimPlot(scn, pt.size = 0.2, cols = mycol, label=F, repel=T,
                 raster = FALSE, label.size = 5, reduction = "umap",
                 group.by = "stage.ano")
plot3

ggsave("UMAP.pat.ano.pdf", plot3, width = 7, height = 5)

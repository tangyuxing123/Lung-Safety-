rm(list = ls())
setwd("D:\\GG市\\单细胞\\GSE149614\\sc0306")

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
# install_github("immunogenomics/harmony")
library(ggplot2)

#读取----------------------------------------------------------------------------------------
scRNAlist <- fread('GSE149614_HCC.scRNAseq.S71915-NT.count.txt', data.table = FALSE)
rownames(scRNAlist) <- scRNAlist[,1]
scRNAlist <- scRNAlist[,-1]
scRNAlist <- CreateSeuratObject(counts = scRNAlist)


scRNA <- merge(x=scRNAlist[[1]],y=scRNAlist[-1])
scRNA <- JoinLayers(scRNA)#连接layers层的count数据

#-------------------------------------------------------------------------


scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
#计算红细胞比例-----------------------------------------------------------
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
col.num <- length(levels(scRNA@active.ident))
#质控前---------------------------------------------------------------------
violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0)
violin

#第一步 数据过滤------------------------------------------------------------
# 计算线粒体、血红蛋白基因比例
scRNA <- PercentageFeatureSet(scRNA,
                              pattern = "^MT-", 
                              col.name = "pMT")

scRNA <- PercentageFeatureSet(scRNA,
                              pattern = "^HBA|^HBB",
                              col.name = "pHB")

# 画图-----------------------------------------------------------------------
VlnPlot(scRNA, 
        features = c("nCount_RNA", "nFeature_RNA", "pMT", "pHB"), 
        log = T,
        pt.size = 0)

# 设置参数
nFeature_lower <- 200
nFeature_upper <- 8000
pMT_upper <- 20
pHB_upper <- 5


# 过滤，提取子集很重要----------------------------------------------------
scRNA <- subset(scRNA,
                subset =nFeature_RNA > nFeature_lower &
                  nFeature_RNA < nFeature_upper &
                  pMT < pMT_upper&
                  pHB < pHB_upper )

# 第二步 标准化（简略版）------------------------------------------------------------
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
#第三步  降维聚类 --------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)

# 高变基因
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) ###或者2000
#只对高变基因进行scale
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)



gc()
scRNA<- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
library(harmony)

scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident")

plot2 <- ElbowPlot(scRNA, ndims=50, reduction="pca") 
plot2



#PC选取---------------------------------------------------------------
pc.num=1:23
scRNA <- FindNeighbors(scRNA, dims = pc.num)
scRNA = FindClusters(scRNA ,resolution = c(seq(0.05,0.5,0.05)))#分辨率可以自己调

save(scRNA,file = "data.rds")


install.packages("ggplot2")
load("data_ann.rds")
library(clustree)
clustree(scRNA@meta.data, prefix = "RNA_snn_res.")
scRNA <- FindClusters(scRNA,resolution = 0.01)
#tSNE------------------------------------------------------
# scRNA = RunTSNE(scRNA, dims = pc.num)
scRNA= RunUMAP(scRNA, dims = pc.num)
plot1 = DimPlot(scRNA, reduction = "umap",label = T) 
plot1







#####数据归一化、筛选高变基因与PCA降维（详细版）---------------------------------------------
#harmony整合
#数据归一化
scRNA_harmony <- scRNA

scRNAlist <- NormalizeData(scRNAlist) %>% #数据归一化
  FindVariableFeatures(selection.method = "vst",nfeatures = 2000) %>% #筛选高变基因2000
  ScaleData() %>% #数据标准化
  RunPCA(npcs = 50, verbose = T)#npcs计算存储PC数（默认50）
a=DimPlot(scRNAlist,reduction = "pca",group.by = "orig.ident")
#PCA存在批次效应
a

#提取前15个高变基因-------------------------------------------------------------------
top15 <- head(VariableFeatures(scRNAlist), 15) 
plot1 <- VariableFeaturePlot(scRNAlist) 
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=3) 
# 合并图片
feat_15 <- CombinePlots(plots = list(plot1,plot2),legend = "bottom")
feat_15
#保存图片
ggsave(file = "feat_15.pdf",plot = feat_15,he = 10,wi = 15 )



scRNA_harmony <- scRNA

#####RunHarmony去批次#####-----------------------------------------------------
#整合指定Seurat对象和metadata中的变量名
scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "orig.ident")
scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
b=DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")

b
#合并图片
pca_harmony_integrated <- CombinePlots(list(a,b),ncol=1)
pca_harmony_integrated




#####聚类、umap/tsne降维（方法1）-----------------------------------------------------
ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony")
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:10) %>% FindClusters(resolution = 0.05)#设置分辨率
##umap/tsne降维
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:10)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:10)

#绘图
umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "orig.ident") 
tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)
#合并
umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,umap_integrated1,umap_integrated2),ncol=2)
# 输出
umap_tsne_integrated
#保存图片
ggsave("umap_tsne_integrated.pdf",umap_tsne_integrated,wi=25,he=15)
#保存数据
save(scRNA_harmony,scRNAlist,file = "scdata2.Rdata")
load("scdata_harmony.Rdata")




####方法2----------------------------------------------------------------------------------------------
library(gplots)
library(ggplot2)
library(Seurat) 

ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony")#根据拐点
scRNA_harmony=sce.all.filt
####调试dims####
# 默认resolution = 0.2
for (i in c(6, 11, 15, 20, 23, 25, 30, 35, 40)) {
  scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.2)
  scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(scRNA_harmony, reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5) + labs(title = paste0("dims: ", i)))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)#assign() 函数 plot_i赋值给plotname
  print(plot_i)
}
p = plot_6 + plot_11 + plot_15 + plot_20 + plot_23 + plot_25 + plot_30 + plot_35 + plot_40 +plot_layout( ncol = 3 )
ggsave( p, filename = "dims.06-40.pdf", width = 15, height = 8 )
print(plot_40)


####调试resolution####------------------------------------------------------

# BiocManager::install("rlang")
# install.packages("ggplot2")
library(clustree)
set.seed(123)#设置非随机
#设置dim=30
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:23)
scRNA_harmony <- FindClusters(
  object = scRNA_harmony,
  resolution = c( seq( 0.1, 1, 0.1) ) # 0.1,-1之间，间隔 0.1
)
#删除RNA_snn_res
# scRNA_harmony@meta.data <- scRNA_harmony@meta.data[, -which(colnames(scRNA_harmony@meta.data) == "RNA_snn_res.1.5")]

#可视化
clustree(scRNA_harmony@meta.data, prefix = "RNA_snn_res.")
ggsave( filename = "分群clustree.pdf", width = 12, height = 9 )

#查看分群效果
colnames(scRNA_harmony@meta.data)
table( scRNA_harmony$RNA_snn_res.0.3 )

p1 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.1", 
              label = T, repel = F, shuffle = T )
#
p2 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.2", 
              label = T, repel = F, shuffle = T )
#
p3 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.3", 
              label = T, repel = F, shuffle = T )
#
p4 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.4", 
              label = T, repel = F, shuffle = T )
#
p5 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.5", 
              label = T, repel = F, shuffle = T )
#
p6 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.6", 
              label = T, repel = F, shuffle = T )
#
p7 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.7", 
              label = T, repel = F, shuffle = T )
#
p8 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.8", 
              label = T, repel = F, shuffle = T )
#
p9 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.9", 
              label = T, repel = F, shuffle = T )
#
p10 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.1", 
               label = T, repel = F, shuffle = T )
#
p = p1 + p2 + p3 + p4 + p5 + p6 + p7+ p8+ p9+ p10 +plot_layout( ncol = 4 )
ggsave( p, filename = "RNA_snn_res.0.1-1.pdf", width = 25, height = 8 )
p
p3 = DimPlot( scRNA_harmony, reduction = "umap", group.by = "RNA_snn_res.0.1", 
              label = T, repel = F, shuffle = T ,pt.size = 1,label.size = 5)  
p3
# plot1 =DimPlot(scRNA_harmony, group.by = "RNA_snn_res.0.2",reduction = "umap",label = T,pt.size = 1,repel = F,label.size = 5) 

#####正式降维####------------------------------------------------------------------------------
#dim30，reso0.3
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:23) %>% FindClusters(resolution = 0.1)
table(scRNA_harmony@meta.data$seurat_clusters)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:23)#dims30

#UMAP----------------------------------------------------------------------------------------------
plot1 =DimPlot(scRNA_harmony, reduction = "umap",label = T, repel = F, shuffle = T ,pt.size = 1,label.size = 8) 
# theme(text = element_text(size = 20))  # 字体大小
plot1
ggsave( filename = "Umap_dim23-reso0.1.pdf", width = 12, height = 9 )
save(scRNA_harmony,file = "HCC_scRNA_harmony_dim_23_rseo0.1.Rdata")
load("HCC_scRNA_harmony_dim_23_rseo0.1.Rdata")
save(scRNA_harmony,file = "HCC_scRNA_harmony_dim_23_rseo0.1.rds")




#注释----------------------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scCustomize)
library(viridis)
library(RColorBrewer)
library(gridExtra)
plot1 = DimPlot(scRNA_harmony, reduction = "umap",label = T) 
plot1

Idents(scRNA_harmony)= "seurat_clusters"
cell.markers <- FindAllMarkers(object = scRNA_harmony, 
                               only.pos = FALSE, # 是否只保留表达相对上调的基因，设置FALSE则会保留下调的
                               test.use = "wilcox", # 默认使用 wilcox 非参数检验，其它选项可以查看说明
                               slot = "data", # 需要注意的是，默认使用 data，而不是 counts
                               min.pct = 0.25, # 设置表达比例的阈值，没有统一标准，差异基因很多的情况下可以把阈值调高，差异基因有1000个够用
                               logfc.threshold = 0.4) # 设置 log2FC 即差异倍率的阈值，没有统一标准，差异基因很多的情况下可以把阈值调高

write.table(cell.markers,"marker.txt",sep="\t")


cell_mark <- c( "TAGLN","ACTA2","SYNPO2",#SMCs\Smooth Muscle Cell
                "DCN","CFD","MMP2",#Fibro 
                "MLANA","CDH19","PCSK2",#Melano/Schwann
                "KRT5","KRT14",#BasalKeras
                "TPSAB1","CTSG",#	Mast_cells
                "CCL21","TFF3","PROX1",#LymphEndo
                "KLRB1","NCR1",#NK
                "APOC3","FABP1","APOA1",#hepa
                'CD19', 'CD79A', 'MS4A1',#B
                'CD8B', 'ZNF683','FCGR3A','FCGR3B','NCAM1',#NKT
                'KRT19','ITGAE','VCAM1','IL6R','ANPEP','CD24','CDH1',#Epithelial cell
                'CD3D','CD3E','CD8A',## Tcells
                'PECAM1','VWF','MCAM','CD34','ESM1', ## Endothelial
                'PTPRC','CD14','AIF1','TYROBP','CD163')#Myeloid 





DotPlot(scRNA_harmony, features =cell_mark) + 
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ 
  scale_color_gradientn(colours = c('#336699','#66CC66','#FFCC33'))



#################手动注释---------------------------------------
cell_type <- c("0"= "T",
               "1"= "Hepatocyte",
               "2"= "Myeloid", 
               "3"= "Myeloid", 
               "4"= "Endothelial", 
               "5"= "Hepatocyte",
               "6"= "Epithelial", 
               "7"= "SMCs", 
               "8"= "B",
               "9"= "Hepatocyte")


scRNA_harmony@meta.data$cell_type = unname(cell_type[scRNA_harmony@meta.data$seurat_clusters])
Idents(scRNA_harmony)=scRNA_harmony@meta.data$cell_type

DimPlot(scRNA_harmony, reduction = 'umap', 
        label = TRUE, pt.size=1)



DotPlot(scRNA_harmony, features =cell_mark) + 
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ 
  scale_color_gradientn(colours = c('#336699','#66CC66','#FFCC33'))

save(scRNA_harmony,file = "data_marker_cell.rds")


# 功能：进行单细胞转录组基本分析流程
# 主要问题：
    # 1.细胞、基因过滤。数据标准化。
    # 2.线性降维和非线性降维。
    # 3.高变化基因鉴定及可视化。
    # 4.细胞聚类及可视化。
# 输入：
    # barcodes.tsv：样本/细胞名称
    # genes.tsv：基因名称/探针号
    # matrix.mtx：表达矩阵
# 输出：
    # MarkerCluster.png
    # MarkerHeatMap.png
    # MarkerVln.png
    # PCA1.png
    # PCA2.png
    # PCA3.png
    # PCSelectElbow.png
    # PCSelectJackStraw.png
    # QCPlot.png
    # TypeDefine.png
    # Umap.png
    # highlyVarGenesWithTop10.png
    # pbmc3k_final.rds
    # pbmc_tutorial.rds

    # 
    # 


# 0工具包
rm(list = ls())
gc()
library(modifyTools)
library(writexl)
library(ggplot2)
# install.packages('Seurat')
# 系统命令行pip install umap-learn
library(Seurat)
library(patchwork)
library(dplyr)
# 0函数


# 1目录
# 获取当前脚本文件路径全称
wholeName = parent.frame(2)$filename
# 创建输出文件夹，获得输出目录
outputDir = envirPrepare(wholeName)


# 2数据
# 当前目录文件
fileNames = dir()
# 数据：读入数据；生成Seurat对象
# 目录：数据所在位置
pbmc.data = Read10X(data.dir = getwd())#一种矩阵
pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", 
                            min.cells = 3, min.features = 200)#参数基础筛选


# 3数据质量控制：去除低质量细胞
# 计算线粒体基因组比例
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因标记
# 质控的可视化
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = paste(outputDir,"QCPlot.png",sep = ""),
       width = 10,height = 8)
# 质控的实现
pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# 特征相关性可视化：额外提供的函数，非必要步骤
# FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")


# 4数据标准化（Normalizing）:数据位置pbmc[["RNA"]]@data
pbmc = NormalizeData(pbmc)#默认参数normalization.method = "LogNormalize", scale.factor = 10000


# 5高变化基因鉴定（Highly variable features）
# 前2000高变化
pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)#筛选方法及基因数
# 标记前10
top10 = head(VariableFeatures(pbmc), 10)
plot1 = VariableFeaturePlot(pbmc)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave(filename = paste(outputDir,"highlyVarGenesWithTop10.png",sep = ""),
       width = 10,height = 8)


# 6数据缩放(Scaling):数据存储在pbmc[["RNA"]]@scale.data
all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)


# 7PCA线性降维:pbmc[["pca"]]产看信息
# PCA分析：默认只使用高变化基因，可手动选择其它子集
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# 三种可视化方法:推荐热图
# VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
# ggsave(filename = paste(outputDir,"PCA1.png",sep = ""),
#        width = 10,height = 8)
# DimPlot(pbmc, reduction = "pca")
# ggsave(filename = paste(outputDir,"PCA2.png",sep = ""),
#        width = 10,height = 8)
DimHeatmap(pbmc, dims = 1:2, cells = 500, balanced = TRUE,fast = F)#fast禁用ggplot2功能以提高速度，导致不可定制化保存不方便
ggsave(filename = paste(outputDir,"PCA3.png",sep = ""),
       width = 10,height = 8)


# 8数据集维度确定（PC个数，对应基因降维）
# JackStraw算法
pbmc = JackStraw(pbmc, num.replicate = 100)
pbmc = ScoreJackStraw(pbmc, dims = 1:20)
# 可视化：确定PC个数
JackStrawPlot(pbmc, dims = 1:15)
ggsave(filename = paste(outputDir,"PCSelectJackStraw.png",sep = ""),
       width = 10,height = 8)
ElbowPlot(pbmc)
ggsave(filename = paste(outputDir,"PCSelectElbow.png",sep = ""),
       width = 10,height = 8)


# 9细胞聚类:Idents(pbmc)查看聚类情况
# 构建KNN图
pbmc = FindNeighbors(pbmc, dims = 1:10)#PC的选择
# 聚类
pbmc = FindClusters(pbmc, resolution = 0.5)#resulution0.4-1.2对3000细胞效果最好，一般细胞越多取值设置高一些


# 10非线性降维
# 调用python中的Umap
pbmc = RunUMAP(pbmc, dims = 1:10,umap.method = "umap-learn",metric = "correlation")
DimPlot(pbmc, reduction = "umap")
ggsave(filename = paste(outputDir,"Umap.png",sep = ""),
       width = 10,height = 8)
# 保存所有计算结果：后续可直接导入
saveRDS(pbmc, file = paste(outputDir,"pbmc_tutorial.rds",sep = ""))


# 11类别标志物
# 鉴定某类标记物:ident.1为目标类别，pct任意类别中的比例
cluster2.markers = FindMarkers(pbmc, 
                               ident.1 = 2,#目标细胞类型 
                               # ident.2 = c(0, 3),#对照类别，默认为其它所有
                               min.pct = 0.25,#基因在任意细胞类型中的检测比例
                               # 
                               )
# 鉴定所有标记物
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, #只有正向标志物
                               min.pct = 0.25, 
                               logfc.threshold = 0.25,#LFC阈值
                               # test.use = "roc",#ROC方式
                               )
# 查看所有标志物信息
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
# 小提琴图可视化
head(cluster2.markers, n = 5)
VlnPlot(pbmc, 
        features = c("MS4A1", "CD79A"),
        # slot = "counts",
        # log = T
        )
ggsave(filename = paste(outputDir,"MarkerVln.png",sep = ""),
       width = 10,height = 8)
# 聚类图可视化
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14"))
ggsave(filename = paste(outputDir,"MarkerCluster.png",sep = ""),
       width = 10,height = 8)
# 热图可视化
pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
ggsave(filename = paste(outputDir,"MarkerHeatMap.png",sep = ""),
       width = 10,height = 8)


# 12细胞类别标注
new.cluster.ids = c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) = levels(pbmc)
pbmc = RenameIdents(pbmc, new.cluster.ids)
saveRDS(pbmc, file = paste(outputDir,"pbmc3k_final.rds",sep = ""))
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(filename = paste(outputDir,"TypeDefine.png",sep = ""),
       width = 10,height = 8)

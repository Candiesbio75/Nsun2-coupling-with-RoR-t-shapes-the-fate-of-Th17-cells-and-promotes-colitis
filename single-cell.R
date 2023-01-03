source activate r4_py37_env
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(DoubletFinder)
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

################################step1:将mapping之后的数据进行merge
setwd('analysis')
file <- c('WT_H20','WT_DSS','KO_H2O','KO_DSS')
path <- './'
###################读取原始数据
for(i in 1:length(file)){
  input <- paste(path,file[i],'filtered_feature_bc_matrix',sep='/')
  tmp <- Read10X(input)
  seurat_obj <- CreateSeuratObject(counts = tmp, 
                                   min.features = 250, 
                                  project = file[i])
 assign(file[i], seurat_obj)
}

###################合并数据
merged_seurat <- merge(x = WT_H20, 
                       y = c(WT_DSS, KO_H2O, KO_DSS), 
                       add.cell.id = c('WT_H20','WT_DSS','KO_H2O','KO_DSS'))
#查看是数据head(WT_H20@meta.data)

#---------------------------以下是质量评估----------------------------------
#-------------------------https://zhuanlan.zhihu.com/p/134135859------------
# 计算每个细胞每个UMI的基因数目并添加到元数据中,用于评估文库复杂度
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# 计算线粒体比率,用于评估细胞活性
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# 创建元数据数据框
metadata <- merged_seurat@meta.data
# 为元数据添加细胞ID
metadata$cells <- rownames(metadata)

# 重命名列
#orig.ident：通常包含所知的样品名，默认为我们赋给project的值
#nCount_RNA：每个细胞的UMI数目
#nFeature_RNA：每个细胞所检测到的基因数目
metadata <- metadata %>%
dplyr::rename(seq_folder = orig.ident,
              nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# 创建样本列
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^WT_H20_"))] <- "WT_H20"
metadata$sample[which(str_detect(metadata$cells, "^WT_DSS_"))] <- "WT_DSS"
metadata$sample[which(str_detect(metadata$cells, "^KO_H2O_"))] <- "KO_H2O"
metadata$sample[which(str_detect(metadata$cells, "^KO_DSS_"))] <- "KO_DSS"


# 将元数据添加回Seurat对象中
merged_seurat@meta.data <- metadata
# 任何时候都要创建.RData对象保存进度
save(merged_seurat, file="merged_filtered_seurat.RData")


setwd('./analysis/UMAP_mito20_fitering4500')
load('./analysis/merged_filtered_seurat.RData')
metadata <- merged_seurat@meta.data
pdf('raw.pdf')
# 可视化每个样本的细胞计数
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
 theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


# 可视化每个细胞的测序深度
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
 geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000)

# 通过频数图可视化每个细胞检测出的基因数分布
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# 通过箱线图可视化每个细胞检测到的基因的分布
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# 可视化每个细胞检测到的线粒体基因表达分布
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)



# 通过可视化每一个UMI检测到的基因数来可视化基因表达的整体复杂性
metadata %>%
   ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# 可视化检测到的基因数和UMI数之间的关系，并且观察是否存在大量低数目的基因数/UMI数的细胞
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
dev.off()




################################step2:现在我们已经生成了要评估的各种指标，我们可以用可视化的方式来探索它们。我们将评估各种指标，然后决定哪些细胞质量较低，而应从分析中删除：
#细胞计数 (Cell counts)
#每个细胞的UMI计数 (UMI counts per cell)
#每个细胞检测到的基因 (Genes detected per cell)
#检测到的UMI数对比基因数 (UMIs vs. genes detected)
#线粒体计数比率 (Mitochondrial counts ratio)
#复杂度 (Novelty)
setwd('./analysis/UMAP_mito20_fitering4500')
load('./analysis/merged_filtered_seurat.RData')
###################过滤线粒体比例20%
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500 & nGene <= 4500) & 
                            (nGene >= 250) & 
                            (mitoRatio <= 0.2))		
save(filtered_seurat, file="filtered_seurat.RData")
metadata <- filtered_seurat@meta.data

pdf('filter.pdf')
# 可视化每个样本的细胞计数
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
 geom_bar() +
 theme_classic() +
 theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
 ggtitle("NCells")


# 可视化每个细胞的测序深度
metadata %>% 
 ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
   geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000)

# 通过频数图可视化每个细胞检测出的基因数分布
metadata %>% 
 ggplot(aes(color=sample, x=nGene, fill= sample)) + 
 geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# 通过箱线图可视化每个细胞检测到的基因的分布
metadata %>% 
 ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
 theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# 可视化每个细胞检测到的线粒体基因表达分布
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)



# 通过可视化每一个UMI检测到的基因数来可视化基因表达的整体复杂性
metadata %>%
 ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# 可视化检测到的基因数和UMI数之间的关系，并且观察是否存在大量低数目的基因数/UMI数的细胞
metadata %>% 
ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
 geom_point() + 
 scale_colour_gradient(low = "gray90", high = "black") +
 stat_smooth(method=lm) +
  scale_x_log10() + 
 scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
dev.off()


##############################去除double########################################
setwd('./analysis/UMAP_mito20_fitering4500')
load('filtered_seurat.RData')
split.list <- SplitObject(filtered_seurat, split.by = "seq_folder")
split.list <- lapply(X = split.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x)
  seurat_integrated.1 <- ScaleData(object = x, verbose = FALSE)
  seurat_integrated.2 <- RunPCA(object = seurat_integrated.1, npcs = 20, verbose = FALSE)
  seurat_integrated.3 <- FindNeighbors(object = seurat_integrated.2, reduction = "pca", dims = 1:20)
  seurat_integrated.4 <- FindClusters(seurat_integrated.3, resolution = 0.5)
  seurat_integrated.6 <- RunUMAP(object = seurat_integrated.4, reduction = "pca", dims = 1:20)
  seurat_integrated.7 <- paramSweep_v3(seurat_integrated.6, PCs = 1:20,sct = FALSE)
  
 seurat_integrated.8 <- summarizeSweep(seurat_integrated.7, GT = FALSE)
  seurat_integrated.9 <- find.pK(seurat_integrated.8)
  
  mpK<-as.numeric(as.vector(seurat_integrated.9$pK[which.max(seurat_integrated.9$BCmetric)]))
  annotations <- seurat_integrated.6@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           
  nExp_poi <- round(0.075*length(seurat_integrated.6@active.ident)) 
  
  mpK <-as.numeric(as.vector(seurat_integrated.9$pK[which.max(seurat_integrated.9$BCmetric)]))
  doublet  <- doubletFinder_v3(seurat_integrated.6, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
  
  save(doublet,file=paste(x@meta.data$seq_folder[1],'_doublet.Robj',sep='')) 
  #########################返回数值
})




library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(scDblFinder)
library(SeuratWrappers)
library(SingleCellExperiment)

#load all data

count_data <- read.table(file = "../C_P2-5-11-25-30/version2/C_P2-5-11-25-30_merged_counts.dge.txt.gz", header = TRUE, row.names = 1)
sobj_c <- CreateSeuratObject(counts = count_data, min.cells = 3, min.features = 200, project = "C")

count_data <- read.table(file = "../D_P7-9-12-26-27/version2/D_P7-9-12-26-27_merged_counts.dge.txt.gz", header = TRUE, row.names = 1)
sobj_d <- CreateSeuratObject(counts = count_data, min.cells = 3, min.features = 200, project = "D")

count_data <- read.table(file = "../E_P1-2-6-26-28/Version2/E_P1-2-6-26-28_merged_counts.dge.txt.gz", header = TRUE, row.names = 1)
sobj_e <- CreateSeuratObject(counts = count_data, min.cells = 3, min.features = 200, project = "E")

count_data <- read.table(file = "../F_P3-9-11-24-30/version2/F_P3-9-11-24-30_merged_counts.dge.txt.gz", header = TRUE, row.names = 1)
sobj_f <- CreateSeuratObject(counts = count_data, min.cells = 3, min.features = 200, project = "F")

#save backups
backup_c <- sobj_c
backup_d <- sobj_d
backup_e <- sobj_e
backup_f <- sobj_f

#load all demultiplexing information
data_c <- read.csv("../C_P2-5-11-25-30/version2/assigned/result.csv")
data_d <- read.csv("../D_P7-9-12-26-27/version2/assigned/result.csv")
data_e <- read.csv("../E_P1-2-6-26-28/Version2/assigned/result.csv")
data_f <- read.csv("../F_P3-9-11-24-30/version2/assigned/result.csv")


#add demultiplexing information
sobj_c$demultiplexed <- data_c$sample[match(colnames(sobj_c), data_c$X)]
sobj_c$type <- data_c$type[match(colnames(sobj_c), data_c$X)]
sobj_c$type[sobj_c$type == "unknown"] <- "singlet"

sobj_d$demultiplexed <- data_d$sample[match(colnames(sobj_d), data_d$X)]
sobj_d$type <- data_d$type[match(colnames(sobj_d), data_d$X)]
sobj_d$type[sobj_d$type == "unknown"] <- "singlet"

sobj_e$demultiplexed <- data_e$sample[match(colnames(sobj_e), data_e$X)]
sobj_e$type <- data_e$type[match(colnames(sobj_e), data_e$X)]
sobj_e$type[sobj_e$type == "unknown"] <- "singlet"

sobj_f$demultiplexed <- data_f$sample[match(colnames(sobj_f), data_f$X)]
sobj_f$type <- data_f$type[match(colnames(sobj_f), data_f$X)]
sobj_f$type[sobj_f$type == "unknown"] <- "singlet"

#start downstream
sobj_c[["percent.mt"]] <- PercentageFeatureSet(sobj_c, pattern = "^mt-")
sobj_c <- subset(sobj_c, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
sobj_c <- NormalizeData(sobj_c, normalization.method = "LogNormalize", scale.factor = 10000)
sobj_c <- FindVariableFeatures(sobj_c, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sobj_c)
sobj_c <- ScaleData(sobj_c, features = all.genes)

sobj_d[["percent.mt"]] <- PercentageFeatureSet(sobj_d, pattern = "^mt-")
sobj_d <- subset(sobj_d, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
sobj_d <- NormalizeData(sobj_d, normalization.method = "LogNormalize", scale.factor = 10000)
sobj_d <- FindVariableFeatures(sobj_d, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sobj_d)
sobj_d <- ScaleData(sobj_d, features = all.genes)

sobj_e[["percent.mt"]] <- PercentageFeatureSet(sobj_e, pattern = "^mt-")
sobj_e <- subset(sobj_e, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
sobj_e <- NormalizeData(sobj_e, normalization.method = "LogNormalize", scale.factor = 10000)
sobj_e <- FindVariableFeatures(sobj_e, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sobj_e)
sobj_e <- ScaleData(sobj_e, features = all.genes)

sobj_f[["percent.mt"]] <- PercentageFeatureSet(sobj_f, pattern = "^mt-")
sobj_f <- subset(sobj_f, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
sobj_f <- NormalizeData(sobj_f, normalization.method = "LogNormalize", scale.factor = 10000)
sobj_f <- FindVariableFeatures(sobj_f, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sobj_f)
sobj_f <- ScaleData(sobj_f, features = all.genes)


#######seurat doublet detection

cds <- as.cell_data_set(sobj_c)
colLabels(cds) <- sobj_c$demultiplexed
cds <- scDblFinder(cds)
sobj_c$doublet_score <- cds$scDblFinder.score
sobj_c$doublet_class <- cds$scDblFinder.class

cds <- as.cell_data_set(sobj_d)
colLabels(cds) <- sobj_d$demultiplexed
cds <- scDblFinder(cds)
sobj_d$doublet_score <- cds$scDblFinder.score
sobj_d$doublet_class <- cds$scDblFinder.class

cds <- as.cell_data_set(sobj_e)
colLabels(cds) <- sobj_e$demultiplexed
cds <- scDblFinder(cds)
sobj_e$doublet_score <- cds$scDblFinder.score
sobj_e$doublet_class <- cds$scDblFinder.class

cds <- as.cell_data_set(sobj_f)
colLabels(cds) <- sobj_f$demultiplexed
cds <- scDblFinder(cds)
sobj_f$doublet_score <- cds$scDblFinder.score
sobj_f$doublet_class <- cds$scDblFinder.class

######################integrate data###############

sobj_all <- merge(sobj_c, y = c(sobj_d, sobj_e, sobj_f), add.cell.ids = c("C", "D", "E", "F"), project = "vierbuchen")

saveRDS(sobj_all, "sobj_all.RDS")
sobj_all <- readRDS("sobj_all.RDS")


pdf("plots/VlnPlot.pdf", width = 18, height = 12)
VlnPlot(sobj_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


plot1 <- FeatureScatter(sobj_all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sobj_all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf("plots/QC_plot.pdf", width = 18, height = 12)
plot1 + plot2
dev.off()

#sobj_all <- subset(sobj_all, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
#sobj_all

sobj_all <- NormalizeData(sobj_all, normalization.method = "LogNormalize", scale.factor = 10000)

sobj_all <- FindVariableFeatures(sobj_all, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sobj_all)
sobj_all <- ScaleData(sobj_all, features = all.genes)

sobj_all <- RunPCA(sobj_all, features = VariableFeatures(object = sobj_all))

pdf("plots/Elbow_plot.pdf", width = 10, height = 8)
ElbowPlot(sobj_all)
dev.off()

#sobj_all <- FindNeighbors(sobj_all, dims = 1:20)
#sobj_all <- FindClusters(sobj_all, resolution = 0.5)
sobj_all <- RunUMAP(sobj_all, dims = 1:20)


pdf("plots/UMAP_plot.pdf", width = 10, height = 8)
DimPlot(sobj_all, reduction = "umap", label = T)
dev.off()


pdf("plots/UMAP_plot_nfeature.pdf", width = 10, height = 8)
FeaturePlot(sobj_all, features = "nFeature_RNA",reduction = "umap", label = T)
dev.off()

pdf("plots/UMAP_plot_demultiplexed.pdf", width = 10, height = 8)
DimPlot(sobj_all, reduction = "umap", group.by = "demultiplexed", label = T)
dev.off()

sobj_all$low_qual <- ifelse(sobj_all$nFeature_RNA > 2000, "false", "true")

pdf("plots/UMAP_plot_low_qual.pdf", width = 10, height = 8)
DimPlot(sobj_all, reduction = "umap", group.by = "low_qual", label = T)
dev.off()

sobj_all$is_doublet <- ifelse(sobj_all$type == "doublet" | sobj_all$doublet_class == "doublet", "true", "false")
table(sobj_all$is_doublet)

pdf("plots/UMAP_plot_doublet.pdf", width = 10, height = 8)
DimPlot(sobj_all, reduction = "umap", group.by = "is_doublet", label = T)
dev.off()

#filter low quality cells

sobj_all_filt <- subset(sobj_all, subset = nFeature_RNA > 2000 & is_doublet == "false")


sobj_all_filt <- NormalizeData(sobj_all_filt, normalization.method = "LogNormalize", scale.factor = 10000)

sobj_all_filt <- FindVariableFeatures(sobj_all_filt, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sobj_all_filt)
sobj_all_filt <- ScaleData(sobj_all_filt, features = all.genes)

sobj_all_filt <- RunPCA(sobj_all_filt, features = VariableFeatures(object = sobj_all_filt))

sobj_all_filt <- RunUMAP(sobj_all_filt, dims = 1:20)

pdf("plots/filt_UMAP_plot.pdf", width = 10, height = 8)
DimPlot(sobj_all_filt, reduction = "umap", label = T)
dev.off()


#integration
sobj_all_filt <- IntegrateLayers(object = sobj_all_filt, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                 verbose = FALSE)
sobj_all_filt <- RunUMAP(sobj_all_filt, dims = 1:20, reduction = "integrated.cca")

pdf("plots/filt_UMAP_plot_integrated.pdf", width = 10, height = 8)
DimPlot(sobj_all_filt, reduction = "umap", label = T)
dev.off()

pdf("plots/filt_UMAP_plot_integrated_demultiplexed.pdf", width = 10, height = 8)
DimPlot(sobj_all_filt, group.by = "demultiplexed", reduction = "umap", label = T)
dev.off()

#feature plots of genes of interest
goi <- c("Pou5f1", "Otx2", "Foxa2", "Sox17", "T", "Sox1")
goi %in% rownames(sobj_all_filt)

for(gene in goi){
  pdf(paste0("plots/UMAP_plot_",gene ,".pdf"), width = 10, height = 8)
  plot <- FeaturePlot(sobj_all_filt, features = gene, reduction = "umap", label = T)
  print(plot)
  dev.off()
  
}



saveRDS(sobj_all_filt, "sobj_all_filt.RDS")
#sobj_all_filt <- readRDS("sobj_all_filt.RDS")

# GSE104556
# Invadolysin (LmLn)
# mouse testis
# use pbmc_0.3

library(dplyr)
library(Seurat)


pbmc.data <- ReadMtx(mtx = "count_matrix.mtx.gz", 
                     features = "features.tsv.gz",
                     cells = "barcodes.tsv.gz")

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GSE104556", min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- subset(pbmc, subset = nFeature_RNA > 2000 & nFeature_RNA < 7500 & percent.mt < 2)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc_0.5 <- FindClusters(pbmc, resolution = 0.5)
pbmc_0.4 <- FindClusters(pbmc, resolution = 0.4)
pbmc_0.3 <- FindClusters(pbmc, resolution = 0.3)
pbmc_0.2 <- FindClusters(pbmc, resolution = 0.2)

pbmc_0.4 <- RunUMAP(pbmc_0.4, dims = 1:15)
pbmc_0.3 <- RunUMAP(pbmc_0.3, dims = 1:15)
pbmc_0.2 <- RunUMAP(pbmc_0.2, dims = 1:15)

DimPlot(pbmc_0.4, reduction = "umap")
DimPlot(pbmc_0.3, reduction = "umap", label = TRUE)
DimPlot(pbmc_0.2, reduction = "umap")

# saveRDS(pbmc_0.3, file = "../pbmc_03.rds")

cd_genes <- c("Spag17", # SP
              "Hspa1l", # CS
              "Dyrk4",  # ES
              "Esx1",   # RS
              "Prss54", # SC
              "Ccna1",  # SC
              "Gpat2",  # SC
              "Crabp1", # SG
              "Wt1",    # Ser
              "Cyp17a1" # Ley
              )

DotPlot(object = pbmc_0.3, features = cd_genes)

new.cluster.ids <- c("SC3", "ES", "CS", "SC2", "SC2", "RS", "SC4", "SC1", "SC1", "SP")
names(new.cluster.ids) <- levels(pbmc_0.3)
pbmc_0.3_anno <- RenameIdents(pbmc_0.3, new.cluster.ids)
DimPlot(pbmc_0.3_anno, reduction = "umap", label = TRUE) + NoLegend()

DotPlot(object = pbmc_0.3_anno, features = "Lmln")
VlnPlot(object = pbmc_0.3_anno, features = "Lmln")
FeaturePlot(pbmc_0.3_anno, feature = "Lmln")





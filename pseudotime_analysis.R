library(Seurat)
library(dplyr)
library(slingshot)
library(SingleCellExperiment) # for setting sce object for slingshot

#### Loading WT datasets ####
for (x in list("WT1","WT2","WT3")) {
  name = paste(x,"data",sep=".") # assign the name for Seurat object, it will look like "WT1.data"
  data =  Read10X(data.dir = x)
  # project name specified the "active.ident" of Seurat object
  assign(name, CreateSeuratObject(counts = data, project = "WT", min.cells = 3, min.features = 200))
}

WT = merge(WT1.data,WT2.data, add.cell.ids = c("WT1", "WT2"))
WT = merge(WT,WT3.data, add.cell.ids = c("WT1_2", "WT3"))

#### Performing normalization, clustering and dimensionality reduction ####
WT = NormalizeData(WT)
WT = FindVariableFeatures(WT, selection.method = "vst", nfeatures = 2000)
all.genes = rownames(WT)
WT = ScaleData(WT, features = all.genes)
WT = RunPCA(WT, features = VariableFeatures(object = WT))

ElbowPlot(WT)

WT = FindNeighbors(WT, dims = 1:10) # dimension is determined by ElbowPlot
WT = FindClusters(WT, resolution = 0.5)
WT = RunUMAP(WT, dims = 1:10, reduction = "pca")
DimPlot(WT, reduction = "umap", group.by = "ident") 

# Filtering away noise
WT = subset(WT,idents = c(6,7,8,9,11,12),invert = TRUE)

#### Setting the sce object ####
sce = as.SingleCellExperiment(WT, assay = "RNA")

#### Running slingshot ####
# slingshot integrates "getLineages" and "getCurves" functions into one
# Set UMAP for reduced dimensional matrix of coordinates and clustering information from Seurat 

# Set the palette for plotcol
library(grDevices)
library(RColorBrewer)
colors = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

# with UMAP embeddings
sce_UMAP = slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP')
summary(sce_UMAP$slingPseudotime_1)

plotcol = colors[cut(sce_UMAP$slingPseudotime_1, breaks=100)]

# Plot the dimensionality reduction results and the Slingshot output on a plot  
plot(reducedDims(sce_UMAP)$UMAP, col = plotcol, pch=16, asp = 1, xlab = "UMAP 1", ylab = "UMAP 2", main = "with UMAP embeddings")
lines(SlingshotDataSet(sce_UMAP), lwd=2, col='black') 

# with PCA embeddings
sce_PCA = slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'PCA')
summary(sce_PCA$slingPseudotime_1)

plotcol = colors[cut(sce_PCA$slingPseudotime_1, breaks=100)]

# Plot the dimensionality reduction results and the Slingshot output on a plot  
plot(reducedDims(sce_PCA)$PCA, col = plotcol, pch=16, asp = 1, main = "with PCA embeddings")
lines(SlingshotDataSet(sce_PCA), lwd=2, col='black')

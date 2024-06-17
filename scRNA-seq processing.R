#### Task 3 ####

library(Seurat)
library(patchwork) # for combing plots
library(metap) # for "FindConservedMarkers"
library(dplyr) # for pipe operator
library(harmony)

#### Loading Data ####
# Load the dataset
for (x in list("KO1","WT1","WT2","WT3")) {
  name = paste(x,"data",sep=".") # assign the name for Seurat object, it will look like "KO1.data"
  data =  Read10X(data.dir = x)
  # set the project name as "KO" and "WT" (first two characters)
  # project name specified the "active.ident" of Seurat object, which will be used in subsequently "group.by" argument
  assign(name, CreateSeuratObject(counts = data, project = substr(name,1,2), min.cells = 3, min.features = 200))
}

#### Quality Control ####
# If the droplet contains cell doublets or the droplet is empty, then the unique feature counts 
# (stored in @meta.data[["nFeature_RNA"]]) ) will be too many or too little

# Mitochondrial contamination presents in low quality or dying cells (add percent.mt for filtering)

# Using nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 as cutoffs

# Add the percent.mt matrics
KO1.data[["percent.mt"]] = PercentageFeatureSet(KO1.data, pattern = "^mt-")
WT1.data[["percent.mt"]] = PercentageFeatureSet(WT1.data, pattern = "^mt-")
WT2.data[["percent.mt"]] = PercentageFeatureSet(WT2.data, pattern = "^mt-")
WT3.data[["percent.mt"]] = PercentageFeatureSet(WT3.data, pattern = "^mt-")

# Using nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 for filtering
KO1.data = subset(KO1.data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
WT1.data = subset(WT1.data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
WT2.data = subset(WT2.data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
WT3.data = subset(WT3.data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#### Setup the Seurat objects ####
# Create a list for Seurat object
merge.object =  list(KO1.data,WT1.data,WT2.data,WT3.data)

#### Normalize and identify variable features for each dataset in merge.object independently ####
gamma_delta.list <- lapply(X = merge.object, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = gamma_delta.list)

#### Perform integration (Seurat default) ####
# Finding anchors across 4 datasets and using these anchors to integrate them into 1 Seurat object
# identify anchors (time consuming)
immune.anchors = FindIntegrationAnchors(object.list = gamma_delta.list, anchor.features = features)
# this command creates an 'integrated' data assay
immune.combined = IntegrateData(anchorset = immune.anchors)

#### Perform an integrated analysis ####
# specify that we will perform downstream analysis on the corrected data (store in @assays)
DefaultAssay(immune.combined) = "integrated"

#### Standard workflow for PCA, visualization and clustering ####
## Scaling the data
immune.combined = ScaleData(immune.combined, verbose = FALSE)

## PCA, linear dimensional reduction
# npcs specifies Total Number of PCs to compute and store (50 by default)
immune.combined = RunPCA(immune.combined, npcs = 30, verbose = FALSE)

# Visualize the result of PCA
# VizDimLoadings : Visualize top genes associated with reduction components
VizDimLoadings(immune.combined, dims = 1:2, reduction = "pca")

# DimPlot: each point is a cell 
DimPlot(immune.combined, reduction = "pca")
# DimHeatmap : draw cell and features(genes) 
# Both cells and genes are sorted by their principal component scores
DimHeatmap(immune.combined, dims = 1:2, cells = 500, balanced = TRUE)

## Determine the ‘dimensionality’ of the dataset that contain useful information
# From the result of ElbowPlot and JackStrawPlot, 14 components are choosen to include  
ElbowPlot(immune.combined)
# JackStraw is time-consuming
# immune.combined = JackStraw(immune.combined, num.replicate = 100)
# immune.combined = ScoreJackStraw(immune.combined, dims = 1:20)
# JackStrawPlot(immune.combined, dims = 1:20)

## Clustering, graph-based clustering approach
# dimension is determined by ElbowPlot or JackStrawPlot 
# resolution sets the ‘granularity’ of the downstream clustering, affect the number of clusters
immune.combined = FindNeighbors(immune.combined, reduction = "pca", dims = 1:14)
immune.combined = FindClusters(immune.combined, resolution = 0.5)

## UMAP, non-linear dimensional reduction
# reduction can specified PCA as input and therefore speed up the UMAP  
immune.combined = RunUMAP(immune.combined, reduction = "pca", dims = 1:30)

## Visualization of UMAP
# p1 using @meta.data[["orig.ident"]] to group, this attribute stores the "WT" and "KO" identifiers 
# from corresponding datasets
# p2 using default  "ident" to group
p1 = DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident") 
p2 = DimPlot(immune.combined, reduction = "umap", group.by = "ident", label = TRUE, repel = TRUE) 
p1 + p2

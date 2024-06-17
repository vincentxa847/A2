# Exploring scRNA-seq datasets with Seurat and Slingshot
## Method
[Standard Seurat pipeline](https://doi.org/10.1016/j.cell.2021.04.048) was performed for integration analysis of datasets with two conditions (WT & KO). 

First, **quality control** was performed to remove low-quality cells and doublets.
```
# Add the percent.mt matrics for all the four datasets
KO1.data[["percent.mt"]] = PercentageFeatureSet(KO1.data, pattern = "^mt-")
WT1.data[["percent.mt"]] = PercentageFeatureSet(WT1.data, pattern = "^mt-")
WT2.data[["percent.mt"]] = PercentageFeatureSet(WT2.data, pattern = "^mt-")
WT3.data[["percent.mt"]] = PercentageFeatureSet(WT3.data, pattern = "^mt-")

# Using nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 to filter droplet contains cell doublets, empty droplet and droplet with Mitochondrial contamination
KO1.data = subset(KO1.data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
WT1.data = subset(WT1.data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
WT2.data = subset(WT2.data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
WT3.data = subset(WT3.data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```
Then **integration** was performed using Seurat default method, which [identifies cell pairwise 
correspondences between single cells across datasets (anchors)](https://doi.org/10.1016/j.cell.2019.05.031) and uses these anchors to integrate multiple datasets. 
```
immune.anchors = FindIntegrationAnchors(object.list = gamma_delta.list, anchor.features = features)
# this command creates an 'integrated' data assay
immune.combined = IntegrateData(anchorset = immune.anchors)
# specify that we will perform downstream analysis on the corrected data (store in 'integrated' assays)
DefaultAssay(immune.combined) = "integrated"
```
**Scaling and PCA** were performed as a preparation for clustering.
```
# Scaling the data
immune.combined = ScaleData(immune.combined, verbose = FALSE)
# PCA, linear dimensional reduction
# npcs specifies Total Number of PCs to compute and store (50 by default)
immune.combined = RunPCA(immune.combined, npcs = 30, verbose = FALSE)
```
Nearest neighbours (NN) were built on dimensionally reduced form of the data, which is PCA data here.
From the result of Elbow Plot, first 14 principle components hold the majority of variation and therefore selected to build NN. 
Resolution parameter was set to 0.5, since clusters will be too many to granular clusters with similar features if parameter is larger than 0.5 was set.
```
# Clustering
# resolution sets the ‘granularity’ of the downstream clustering, affect the number of clusters
immune.combined = FindNeighbors(immune.combined, reduction = "pca", dims = 1:14)
immune.combined = FindClusters(immune.combined, resolution = 0.5)
```
![image](https://github.com/vincentxa847/A2/assets/118545004/fb2fc3cf-bcc4-490c-ab39-2ad195d447e5)\
*Elbow Plot to determine the number of PCs*

From the UMAP coordinates and the markers of clusters, 7 small clusters were removed before generating this annotated UMAP and the remaining cluster 0-6 were annotated using genes that have significant differential positive  expression in the assigned cluster.\
![image](https://github.com/vincentxa847/A2/assets/118545004/846042fd-5117-489b-b5f0-48db45ad4812)\
*Annotated UMAP*

## Result


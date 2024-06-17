# Exploring scRNA-seq datasets of immune cells under two biological conditions using Seurat and Slingshot
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
[Harmony](https://doi.org/10.1038/s41592-019-0619-0) was used as an alternative integration method to compare with the default Seurat integration approach
```
# Run Harmony
# the result of harmony will be stored in @reduction
# Note for UMAP and clustering, using "reduction = harmony"
four_samples_combined.harmony = RunHarmony(four_samples_combined, group.by.vars = "orig.ident", plot_convergence = TRUE)
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

### Datasets used
10X Chromium dataset of a γδ T cell dataset with three wild type (WT) runs and one knock out (KO). The KO is from a β2 integrin deficient mouse. 

## Result
Cluster in the upper left corner with 3 subclusters annotated with “Xcl1” ,“Smc4” and “Ikzf2” has featured markers “Cd27” ,“Cd28” ,“Sell” ,“Ccr7” ,“Plac8” and 
“Ly6c2”. These markers are associated with a subset of [CD27+ IFNγ-producing γδ T cells](https://doi.org/10.1101/2020.09.08.287854). 
Cluster in the lower right corner with 3 subclusters annotated with “Scart1”, “Cd40lg” and “Scart2” has featured markers “Il17a”,“Il17f”,“Fos”,“Zfp36” and “Gadd45b”. These markers are associated with [CD27- IL-17-producing γδ T cells](https://doi.org/10.3389/fimmu.2018.00796). Cluster annotated with “Klra7” has marker genes “Gzma”, “Klra7” and “Klre1” associated with NK cells. In conclusion, 3 main clusters were identified to be associated with CD27+ IFNγ-producing γδ T cell, CD27-IL-17-producing γδ T cells and NK cells.

It can be seen from the frequency of cell in the clusters in two conditions that cell numbers in cluster associated with the CD27+ IFNγ-producing γδ T cells decrease in β2 integrin KO condition. However, the number of cells in cluster associated with the IL-17-producing γδ T cells increase in β2 integrin KO condition. 
The result reveals that [β2 integrins are important regulators of γδ T cell homeostasis that promote the development of IFNγ-producing γδ T cells and inhibit the survival of IL-17-producing γδ T cells](https://doi.org/10.1073/pnas.1921930117).
![image](https://github.com/vincentxa847/Identify_Different_Cell_Types_by_Exploring_scRNA-seq_datasets_with_Seurat_and_Slingshot/assets/118545004/e624b210-72c7-4456-92b9-61f876744db0)
*Frequency of cell in the clusters in two conditions*

DE analysis of cluster Cd40lg between WT and KO was conducted using Seurat default Wilcoxon test. The number of genes that show up-regulation in the KO with the corrected p-values <0.01 and logFC >0.5 (target genes) is 2. The integration method was changed to Harmony, which groups by cell types to better compatible with cell subpopulations identification, the number of target genes increases to 67. 

WT and KO datasets integrate better using Harmony compared to the default Seurat integration method.
![image](https://github.com/vincentxa847/Identify_Different_Cell_Types_by_Exploring_scRNA-seq_datasets_with_Seurat_and_Slingshot/assets/118545004/845a726b-38c6-4f70-8588-3f57426d7efd)
*The performance of integration (mixing of datasets) of Harmony. Harmony (left) and  Default Seurat integration method (right)* 

[MAST](https://doi.org/10.1186/s13059-015-0844-5), was adopted to test DE genes, which is a two-part generalized linear model to model the rate of expression of various transcripts and the positive expression mean. It also takes the cellular detection rate into account to deal with biological factors such as cell volume of different cells and technical assay variability. 
The number of DE genes remains the same comparing with default Wilcoxon test, but the P-values and the order of DE genes change. Overall, using Harmony as integration method is better in this dataset. 


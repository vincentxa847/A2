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
# specify that we will perform downstream analysis on the corrected data (store in 'integrated' assays)
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

#### Finding biomarkers of cluster ####
# DE analysis is performed here
# find markers for every cluster that compares to all remaining cells, report only the positive ones
# set the pct cutoff to only test genes that are detected in a minimum fraction of cells 
# set the logfc cutoff to test genes with at least X-fold difference (log-scale) between groups 
cluster.markers = FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# take the genes show highest logFC different between groups
cluster.markers_top_8 = cluster.markers %>%
  group_by(cluster) %>%
  slice_max(n = 8, order_by = avg_log2FC)

## Visualize the markers
# select markers to visualize for each clusters
for (x in c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)) {
  name = paste("cluster.",".top8",sep=as.character(x)) # assign the variable name, it will look like "cluster.0.top8"
  assign(name, as.factor(subset(cluster.markers_top_8, cluster == as.character(x))[,'gene']))
}

# VlnPlot -- showing gene names to better investigate the function of genes
# cluster 0
VlnPlot(immune.combined, features = c("S100a4", "Il17a", "Tmem176a", "Tmem176b", "S100a6", "Cxcr6", "Ccr2", "Maf"))
# cluster 1
VlnPlot(immune.combined, features = c("5830411N06Rik", "Rgcc", "Ccr6", "Abi3bp", "Acsbg1", "Il7r", "S100a4", "Crip1") )
# cluster 2
VlnPlot(immune.combined, features = c("Ly6c2", "Ctla2a", "Klrc1", "Xcl1", "Klrc2", "Ms4a4b", "Ccl5", "Hopx") )
# cluster 3
VlnPlot(immune.combined, features = c("Gm8369", "Ccr7", "Dapl1", "Smc4", "Ms4a6b", "Slamf6", "Cd8b1", "Sell") )
# cluster 4
VlnPlot(immune.combined, features = c("Trdv4", "Cd163l1", "Ckb", "Cxcr6", "Serpinb1a", "Krt83", "Capg", "Tcrg-C1") )
# cluster 5
VlnPlot(immune.combined, features = c("Gzma", "Ccl5", "Klra7", "Klre1", "Nkg7", "Ccl4", "Klrd1", "Zeb2") )
# cluster 6
VlnPlot(immune.combined, features = c("Ikzf2", "Sox4", "Smc4", "Slamf6", "Ccr9", "Ifi27l2a", "Lef1", "Cd27")  )
# cluster 7 -- lysozyme 2 (Lyz2), apolipoprotein E (Apoe), cytochrome b-245 beta chain (Cybb), related to macrophage
VlnPlot(immune.combined, features = c("Lyz2", "Apoe", "Ifitm3", "Cst3", "Lst1", "Cybb", "Cd74", "Gngt2")  )
# cluster 8 -- protein tyrosine phosphatase family (Ptprb), endothelial PAS domain protein 1 (Epas1), related to endothelial cell
VlnPlot(immune.combined, features = c("Gsn", "Cd36", "Calcrl", "Ptprb", "Ramp2", "Epas1", "Egfl7", "Cd93")  )
# cluster 9 -- granzyme C (Gzmc), expressed by activated T cells
VlnPlot(immune.combined, features = c("Xcl1", "Gzmb", "Gzmc", "Tyrobp", "Fcer1g", "Cd7", "Rgs1", "Cd160")  )
# cluster 10 -- chromatin binding activity (Pclaf),inhibitor of apoptosis (Birc5), H2A clustered histone 23 (Hist1h2ap)
VlnPlot(immune.combined, features = c("Stmn1", "Hmgb2", "Gzma", "Hist1h2ap", "Birc5", "Pclaf", "Hmgn2", "Klra7") )
# cluster 11 -- S100 family of proteins containing 2 EF-hand calcium-binding motifs (S100a8,S100a9), family of interferon induced antiviral proteins (Ifitm1)
VlnPlot(immune.combined, features = c("S100a8", "S100a9", "Ifitm1", "Il1b", "Ccl6", "Ccl3", "Ccl4", "Tyrobp")  )
# cluster 12 -- B lymphocyte antigen receptor (Cd79a), immunoglobulin kappa constant, located in blood microparticle and extracellular exosome (Igkc)
VlnPlot(immune.combined, features = c("Igkc", "Cd74", "H2-Eb1", "Cd79a", "Ly6d", "Ighm", "Ebf1", "H2-Ab1")  )
# cluster 13 -- Enables adrenomedullin binding activity (Calcrl), integral membrane proteins (Cldn5), related to endothelial cells
VlnPlot(immune.combined, features = c("Ramp2", "Calcrl", "Cldn5", "Igfbp7", "Cdh5", "Tmem100", "Cd36", "Cyp4b1")  )

## Filtering away noise
# Based on the result of UMAP (small clusters) and cluster' markers, cluster 7,8,9,10,11,12,13 are filtered away
immune.combined.filtered = subset(immune.combined,idents = c(7,8,9,10,11,12,13),invert = TRUE)
# Visualize the result
p1 = DimPlot(immune.combined.filtered, reduction = "umap", group.by = "orig.ident") 
p2 = DimPlot(immune.combined.filtered, reduction = "umap", group.by = "ident", label = TRUE, repel = TRUE) 
p1 + p2

# From the result of Vlnplot and UMAP, it is spectulated that cluster 0,1,4 , cluster 2,3,6 are related cells
# and need to conduct in-cluster comparison to identify the sub-populations
# Using "FindMarkers" to do (finding markers between two groups)

#### Markers of cluster 0,1,4 ####

# Combining the logFC and visualization result (expression of markers in every cluster) to determine featured markers

# Finding markers of cluster 0.1.4
# Featured markers of cluster 0.1.4 were selected, which are "Il17a" "Il17f" "Fos" "Zfp36" "Gadd45b"
cluster.0.1.4.markers = FindMarkers(immune.combined.filtered, ident.1 = c(0,1,4), ident.2 = NULL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Finding markers of cluster0, cluster1 and cluster4 separately
for (x in c(0,1,4)){
  vector = c(0,1,4)
  name = paste("cluster.",".markers",sep = as.character(x))
  assign(name, FindMarkers(immune.combined.filtered, ident.1 = x, ident.2 = vector[!vector == x], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25))
}

# Featured markers of cluster 0 were selected, which are "Areg" "Stmn2" "Cd40lg"
# But Areg and Stmn2 are highly express in cluster 0.1.4, it seems that there no specific marker for cluster 0

# Featured markers of cluster 1 were selected, which are "5830411N06Rik"(Scart2) "Ccr6" "Ccr10" "Tsc22d1"
# "5830411N06Rik" "Ccr6" "Ccr10" are also highly expressed in cluster 0, and "Tsc22d1" is also highly expressed in cluster 6

# Featured markers of cluster 4 were selected, which are "Cd163l1" "Trdv4" "Ckb" "Krt83"


VlnPlot(immune.combined.filtered, features = c("Il17a", "Il17f", "Fos", "Zfp36", "Gadd45b"), n = 5)

#### Markers of cluster 2,3,6 ####

# Combining the logFC and visualization result (expression of markers in every cluster) to determine featured markers

# Finding markers of cluster 2,3,6 
# Featured markers of cluster 2,3,6 were selected, which are "Cd27" "Sell" "Cd28" "Ccr7" "Plac8" and "Ly6c2"
cluster.2.3.6.markers = FindMarkers(immune.combined.filtered, ident.1 = c(2,3,6), ident.2 = NULL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

for (x in c(2,3,6)){
  vector = c(2,3,6)
  name = paste("cluster.",".markers",sep = as.character(x))
  assign(name, FindMarkers(immune.combined.filtered, ident.1 = x, ident.2 = vector[!vector == x], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25))
}

# Featured markers of cluster 2 were selected, which are "Ctla2a" "Xcl1" "Klrc1" "Klrc2" and "Tyrobp"
# From the analysis of cluster2 marker, it was observed cluster 2 and 5 share similiar markers

# Featured markers of cluster 3 were selected, which are "Ccr7" "Dapl1" "Cd8b1" "Smc4" and "Themis"

# Featured markers of cluster 6 were selected, which are "Ikzf2" "Sox4" "Ccr9" "Trbc2" and "Tsc22d1"
# From the analysis of cluster2 marker, it was observed cluster 6 and 3 share similiar markers


VlnPlot(immune.combined.filtered, features = c("Cd27", "Sell", "Cd28", "Ccr7", "Plac8", "Ly6c2"))


#### cluster 5 ####
cluster.5.markers = FindMarkers(immune.combined, ident.1 = 5, ident.2 = NULL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Featured markers of cluster 5 were selected, which are "Gzma"  "Ccl5"  "Ccl5"  "Ccl4" "Klra7" "Klre1" 

# From the analysis of cluster5 markers, it was observed that markers in cluster 5 also have higher expression in cluster 2
VlnPlot(immune.combined.filtered, features = c("Gzma","Ccl5","Ccl4","Klra7","Klre1"))

#### Assigning markergenes to clusters for annotation ####
# marker selected for each cluster was determined by p value and LogFC of the DE analysis (FindMarkers)
new.cluster.ids = c("Cd40lg", "Scart2", "Xcl1", "Smc4", "Scart1", "Klra7","Ikzf2")
names(new.cluster.ids) = levels(immune.combined.filtered)
immune.combined.filtered.renamed = RenameIdents(immune.combined.filtered, new.cluster.ids)

# Visualize the result
p1 = DimPlot(immune.combined.filtered.renamed, reduction = "umap", group.by = "orig.ident") 
p2 = DimPlot(immune.combined.filtered.renamed, reduction = "umap", group.by = "ident", label = TRUE, repel = TRUE) 
p1 + p2

# From the plot below, it can be seen that cluster "Scart1", "Cd40lg" and "Scart2" contain cells from WT and KO datasets
# It is decided to use cluster "Cd40lg" for further DE analysis
DimPlot(immune.combined.filtered.renamed, reduction = "umap", split.by = "orig.ident")

#### DE analysis with cluster "Cd40lg" ####
# Using subset to select cells from identity class "Cd40lg"
Cd40lg = subset(immune.combined.filtered.renamed, idents = "Cd40lg")

# Change the identity class from "Cd40lg" to "WT" and "KO"
Idents(Cd40lg) = "orig.ident"

# Peform DE analysis using FindMarkers function (select only up-regulate genes)
# Note that ident.1 is KO to find up-regulate (only.pos = TRUE) genes in KO when comparing with WT 
# It will return a dataframe with gene name as row name and 5 columns, 
# using avg_log2FC (log fold-change of the average expression between the two groups) and
# p_val_adj (Adjusted p-value, based on Bonferroni correction) for requirement p-values <0.01 and logFC >0.5
Cd40lg.KO_affected = FindMarkers(Cd40lg, ident.1 = "KO", ident.2 = "WT",only.pos = TRUE)
# Two genes, Areg and Gzmb
Cd40lg.KO_affected_filtered = subset(Cd40lg.KO_affected,avg_log2FC > 0.5 & p_val_adj <0.01) 

#### Using Harmony as integration method ####
# Merge four Seurat objects into one (Note that some cell names are duplicated across objects, therefore 
# using add.cell.ids to distinguish cell names from different datasets)
# The identifiers of two conditions (KO & WT) are stored in meta.data "orig.ident"
four_samples_combined = merge(WT1.data,WT2.data, add.cell.ids = c("WT1", "WT2"))
four_samples_combined = merge(four_samples_combined,WT3.data, add.cell.ids = c("WT1_2", "WT3"))
four_samples_combined = merge(four_samples_combined,KO1.data, add.cell.ids = c("WT", "KO"))

# QC for combined Seurat
four_samples_combined[["percent.mt"]] = PercentageFeatureSet(four_samples_combined, pattern = "^mt-")
four_samples_combined = subset(four_samples_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize, FindVariableFeatures, ScaleData and RunPCA for combined Seurat
four_samples_combined = NormalizeData(four_samples_combined)
four_samples_combined = FindVariableFeatures(four_samples_combined,selection.method = "vst", nfeatures = 2000)
four_samples_combined = ScaleData(four_samples_combined, features = rownames(four_samples_combined))
four_samples_combined = RunPCA(four_samples_combined, npcs = 30, verbose = FALSE)
DimPlot(object = four_samples_combined, reduction = "pca",group.by = "orig.ident")

# Run Harmony
# the result of harmony will be stored in @reduction
# Note for UMAP and clustering, using "reduction = harmony"
four_samples_combined.harmony = RunHarmony(four_samples_combined, group.by.vars = "orig.ident", plot_convergence = TRUE)
DimPlot(object = four_samples_combined.harmony, reduction = "harmony",group.by = "orig.ident")
four_samples_combined.harmony = RunUMAP(four_samples_combined.harmony, reduction = "harmony",dims = 1:30)
# Downstream analysis
ElbowPlot(four_samples_combined.harmony)
four_samples_combined.harmony = FindNeighbors(four_samples_combined.harmony, reduction = "harmony", dims = 1:14)
four_samples_combined.harmony = FindClusters(four_samples_combined.harmony, resolution = 0.5)
DimPlot(four_samples_combined.harmony, reduction = "umap", label = TRUE)
DimPlot(four_samples_combined.harmony, reduction = "umap", group.by = "orig.ident", label = TRUE)

# Filter away noise, and DE analysis on cluster 5 (largest cluster that contains cells from the KO and the WT)
four_samples_combined.harmony.filtered = subset(four_samples_combined.harmony,idents = c(7,8,9,10,11),invert = TRUE)

cluster.5.harmony = subset(four_samples_combined.harmony.filtered, idents = "5")
Idents(cluster.5.harmony) = "orig.ident"
cluster.5.harmony.KO_affected = FindMarkers(cluster.5.harmony, ident.1 = "KO", ident.2 = "WT", verbose = FALSE,only.pos = TRUE)
# 67 genes, top 5 are "Itgb1","AC163354.1","Gzmb" ,"Klrk1","Tagln2" 
cluster.5.harmony.KO_affected_filtered = subset(cluster.5.harmony.KO_affected,avg_log2FC > 0.5 & p_val_adj <0.01)
row.names(cluster.5.harmony.KO_affected_filtered)

#### Using MAST as DE method ####
# For default integration method
Cd40lg.MAST.KO_affected = FindMarkers(Cd40lg, ident.1 = "KO", ident.2 = "WT", test.use = "MAST", only.pos = TRUE)
Cd40lg.MAST.KO_affected_filtered = subset(Cd40lg.MAST.KO_affected,avg_log2FC > 0.5 & p_val_adj <0.01)
row.names(Cd40lg.MAST.KO_affected_filtered) # 2 GENES , "Gzmb" "Areg"

# For Harmony intergation method
cluster.5.MAST.harmony.KO_affected = FindMarkers(cluster.5.harmony, ident.1 = "KO", ident.2 = "WT", test.use = "MAST",only.pos = TRUE)
cluster.5.MAST.harmony.KO_affected_filtered = subset(cluster.5.MAST.harmony.KO_affected,avg_log2FC > 0.5 & p_val_adj <0.01)
# 67 genes, top 5 are "Itgb1" "Gzmb" "Cd163l1" "AC163354.1" "Klrk1"
row.names(cluster.5.MAST.harmony.KO_affected)

#### Save the result of DE analysis with four different methods ####
write.table(Cd40lg.KO_affected_filtered, file='default_integration_and_DEmethod.tsv', quote=FALSE, sep='\t',col.names = NA)
write.table(cluster.5.harmony.KO_affected_filtered, file='Harmony_integration_and_Default_DEmethod.tsv', quote=FALSE, sep='\t')
write.table(Cd40lg.MAST.KO_affected_filtered, file='default_integration_and_MAST_DEmethod.tsv', quote=FALSE, sep='\t',col.names = NA)
write.table(cluster.5.MAST.harmony.KO_affected, file='Harmony_integration_and_MAST_DEmethod.tsv', quote=FALSE, sep='\t',col.names = NA)


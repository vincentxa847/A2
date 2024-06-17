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

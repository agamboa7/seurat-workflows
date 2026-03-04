####################################################################
#### Seurat - Guided Clustering Tutorial
#### Andrea AG
#### R.5.0

#
#####--- SETUP ----#####
#
setwd("C:/Users/andyg/Desktop/Thesis/Tutorials/Seurat/Clustering")

# ---- Install Seurat ----
# From Github:
# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
# From CRAN:
# install.packages('Seurat')

library(dplyr)
library(Seurat)
library(patchwork)

# ---- Load input data ----
input_data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19") # Download from the repo

# ---- Initialize the Seurat object with raw data ----
pbmc <- CreateSeuratObject(counts = input_data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Lets examine a few genes in the first thirty cells
input_data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# The . values in the matrix represent 0s (no molecules detected). 
# Since most values in an scRNA-seq matrix are 0, Seurat uses a sparse-matrix 
# representation whenever possible. This results in significant memory 
# and speed savings for Drop-seq/inDrop/10x data.

dense.size <- object.size(as.matrix(input_data))
dense.size

sparse.size <- object.size(input_data)
sparse.size

#
#####--- PREPROCESSING ----#####
#

# Selection and filtering of cells based on QC metrics, data normalization
# and scaling, and the detection of highly variable features

## QC and selecting cells for further analysis:

# Commonly used QC Metrics:
# - Number of unique genes detected in each cell (low quality cells or empty droplets often have very few genes and cell doublets or multiplets may exhibit an aberrantly high gene count)
# - Total number of molecules detected within a cell (correlates strongly with unique genes)
# - Percentage of reads that map to the mitochondrial genome (Low-quality / dying cells often exhibit extensive mitochondrial contamination
# We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a 
# set of features. We use the set of all genes starting with MT- as a set of mitochondrial genes)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC metrics are automatically calculated during CreateSeuratObject. Let's see:
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# --- Visualize QC metrics ----
QC_VlnPlot <- VlnPlot(pbmc, 
                      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                      ncol = 3, cols = c("#e18aff"))
pdf("QC_violinplots.pdf", width = 11, height = 6)
print(QC_VlnPlot)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

FSplot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = c("#7cffe3"))
FSplot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cols = c("#7cb0ff"))
FSplots <- FSplot1 + FSplot2
pdf("FeatureScatterplots.pdf", width = 11, height = 6)
print(FSplots)
dev.off()

# --- Filtering ----
#
# - We filter cells that have unique feature counts over 2,500 or less than 200
# - We filter cells that have >5% mitochondrial counts

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# ---- Normalization ----
#
# After removing unwanted cells from the dataset, the next step is to normalize the data. 
# The default method is "LogNormalize".

pbmc <- NormalizeData(pbmc)

# ---- Feature Selection ----
#
# We next calculate a subset of features that exhibit high cell-to-cell
# variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). 

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10

# Plot variable features with labels
VFplot1 <- VariableFeaturePlot(pbmc)
VFplot2 <- LabelPoints(plot = VFplot1, points = top10, repel = TRUE)
#VFplot2
pdf("VariableFeaturesPlot.pdf", width = 11, height = 6)
print(VFplot2)
dev.off()

# ---- Data Scaling ----
# A linear transformation is applied, which is a standard pre-processing step prior to dimensional reduction.
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#
#####--- LINEAR DIMENSION REDUCTION ----#####
#

# Next we perform PCA on the scaled data. By default, only the previously 
# determined variable features are used as input, but can be defined using features argument.

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# To overcome the extensive technical noise in any single feature for scRNA-seq data, 
# Seurat clusters cells based on their PCA scores, with each PC essentially representing a 'metafeature' 
# that combines information across a correlated feature set. 
# The top principal components therefore represent a robust compression of the dataset.
# An alternative heuristic method generates an ‘Elbow plot’: 
# a ranking of principle components based on the percentage of variance

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
elbw <- ElbowPlot(pbmc)
pdf("ElbowPlot.pdf", width = 11, height = 6)
print(elbw)
dev.off()

#
#####--- CLUSTERING ----#####
#
# Seurat v3 applies a graph-based clustering approach. These methods embed cells in a graph structure 
# for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature 
# expression patterns, and then attempt to partition this graph into highly interconnected 'quasi-cliques' or 'communities'.

pbmc <- FindNeighbors(pbmc, dims = 1:8) # Based on the elbow plot I decided to put 8 components
pbmc <- FindClusters(pbmc, resolution = 0.5)

head(Idents(pbmc), 5)

# ---- Non-linear dimensional reduction (UMAP/tSNE) ----
pbmc <- RunUMAP(pbmc, dims = 1:8)
umaplot <- DimPlot(pbmc, reduction = "umap")
pdf("UMAP.pdf", width = 7, height = 6)
print(umaplot)
dev.off()
# save the object at this point so that it can easily be loaded back in without having to rerun
saveRDS(pbmc, file = "/pbmc_tutorial.rds")

### ---- Finding diferentially expressed features (cluster biomarkers) ----

# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = "avg_logFC")
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
# VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
# Other visualization tools:
RidgePlot(pbmc, features = c("MS4A1", "CD79A"))
DotPlot(pbmc, features = c("MS4A1", "CD79A"))

ftplt <- FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))
pdf("FeaturePlot.pdf", width = 11, height = 6)
print(ftplt)
dev.off()

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
heatmap <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
pdf("HeatmapMarkers.pdf", width = 8, height = 4)
print(heatmap)
dev.off()

# ---- Assigning cell type identity to clusters ----
# Fortunately in the case of this dataset, we can use canonical markers 
# to easily match the unbiased clustering to known cell types:
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
labelclust <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
pdf("UMAP_label.pdf", width = 7, height = 6)
print(labelclust)
dev.off()

# Save the final object:
saveRDS(pbmc, file = "pbmc3k_final.rds")

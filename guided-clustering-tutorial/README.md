# Seurat Workflows - Guided Clustering Tutorial

This repository collects **figures generated while following Seurat’s [Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) workflow** on the PBMC3k dataset. Rather than a step-by-step tutorial, the goal is to **build familiarity with Seurat’s standard clustering pipeline for scRNA-seq** and to document the key QC and clustering outputs with short explanations.

## Repository Structure

- `README.md` - this document  
- `seurat_tutorial.R` - code for the tutorial (scripts / notebooks)
- `filtered_gene_bc_matrices` - input dataset of Peripheral Blood Mononuclear Cells (PBMC) - 2,700 single cells that were sequenced on the Illumina NextSeq 500. (available online [here](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz))
- `.pdfs` - generated figures and outputs

---

## What This Workflow Covers

- **QC**: `nFeature_RNA`, `nCount_RNA`, `% mitochondrial`
- **Filtering**
- **Normalization** (e.g., `NormalizeData`)
- **Variable feature selection** (`FindVariableFeatures`)
- **Scaling** (`ScaleData`)
- **PCA** (`RunPCA`) + diagnostics (ElbowPlot / PC heatmaps)
- **Neighbors + clustering** (`FindNeighbors`, `FindClusters`)
- **UMAP / t-SNE** (`RunUMAP`, optional `RunTSNE`)
- **Markers** (`FindMarkers`) 
- **Common plots**: violin plots, feature plots, dimensional reductions, dot plots

---

## 1) Quality control overview - violin plots

**File:** [QC_violinplots.pdf](QC_violinplots.pdf)

This panel summarizes basic QC metrics commonly inspected before filtering:
- **nFeature_RNA**: number of detected genes per cell
- **nCount_RNA**: total UMI counts per cell
- **percent.mt**: fraction of counts mapping to mitochondrial genes (often used as a proxy for stressed/low-quality cells)

These distributions help decide reasonable thresholds for filtering low-quality cells (e.g., extremely low features/counts or high mitochondrial percentage).

---

## 2) QC relationships - feature scatter plots

**File:** [FeatureScatterplots.pdf](FeatureScatterplots.pdf)

Scatterplots are useful to inspect relationships between QC metrics, typically:
- **nCount_RNA vs nFeature_RNA** (cells with more UMIs usually have more genes detected)
- **nCount_RNA vs percent.mt** (high mitochondrial fraction at a given depth can indicate poor-quality cells)

This helps identify outliers and supports choosing filtering cutoffs.

---

## 3) Highly variable features

**File:** [VariableFeaturesPlot.pdf](VariableFeaturesPlot.pdf)

Seurat selects **highly variable genes** (HVGs) because they contain the strongest biological signal for downstream dimensionality reduction and clustering. This plot shows feature variability across the dataset; the most variable genes are used for PCA and neighbor graph construction.

---

## 4) Choosing the number of PCs - elbow plot

**File:** [ElbowPlot.pdf](ElbowPlot.pdf)

PCA reduces dimensionality by capturing major axes of variation.
The elbow plot shows how much variance each principal component explains; the “elbow” region is often used to pick a reasonable number of PCs to retain for:
- `FindNeighbors()`
- `FindClusters()`
- `RunUMAP()` / `RunTSNE()`

---

## 5) Cluster visualization - UMAP

**File:** [UMAP.pdf](UMAP.pdf)

UMAP is a non-linear dimensionality reduction method that projects cells into 2D while preserving local neighborhood structure. Here, cells are shown in UMAP space and colored by **Seurat clusters** (computed from the neighbor graph). This is typically the “main map” used to interpret cluster structure and separation.

---

## 6) Cluster visualization - labeled UMAP

**File:** [UMAP_label.pdf](UMAP_label.pdf)

Same UMAP embedding as above, but with **cluster labels** added on top. It is useful for presentations and for quickly referencing clusters during marker exploration and annotation.

---

## 7) Canonical marker expression - feature plot

**File:** [FeaturePlot.pdf](FeaturePlot.pdf)

Feature plots overlay gene expression on the UMAP embedding. This is often used to:
- check **canonical markers** for known PBMC populations (e.g., T cells, B cells, NK cells, monocytes)
- support **cluster annotation** by comparing expression patterns across clusters

---

## 8) Marker-based cluster characterization - heatmap

**File:** [HeatmapMarkers.pdf](HeatmapMarkers.pdf)

After identifying cluster markers (e.g., with `FindAllMarkers()`), a heatmap summarizes how **top marker genes** vary across clusters.
This helps:
- confirm clusters are transcriptionally distinct
- compare clusters by their characteristic gene programs
- guide biological interpretation and labeling

---

## Notes

- Exact cluster shapes/labels can vary depending on Seurat version and random seeds used in UMAP.
- These figures are meant as a learning artifact to understand what each stage of the workflow produces and how to interpret it.

---

## Reference

- Seurat Guided Clustering Tutorial (PBMC3k):  
  [https://satijalab.org/seurat/articles/pbmc3k_tutorial.html](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)





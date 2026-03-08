# Seurat Workflows — Guided Clustering Tutorial (Draft)

This repository contains a **reproducible Seurat workflow** based on the official **Guided Clustering** tutorial available [here](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).
The goal of this workflow is to get familiar with Seurat’s standard clustering procedure for scRNA-seq data including **QC**, **normalization**, **feature selection**, **dimensionality reduction**, **clustering**, **marker discovery**, and **visualization**.

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

## Before Getting Started

Make sure you have **R** installed and a working R environment.  
Recommended:

- R (>= 4.1)
- Seurat (v4 or v5)
- ggplot2, patchwork (often pulled in by Seurat)
  
---





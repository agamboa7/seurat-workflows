# Seurat Workflows — Guided Clustering Tutorial (Draft)

This repository contains a **reproducible Seurat workflow** based on the official **Guided Clustering** tutorial available [here](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).
The goal of this workflow is to get familiar with Seurat’s standard clustering procedure for scRNA-seq data including **QC**, **normalization**, **feature selection**, **dimensionality reduction**, **clustering**, **marker discovery**, and **visualization**.

## Repository Structure (Suggested)

Depending on what you uploaded, you can adapt this section to match your repo layout:

- `README.md` — this document  
- `seurat_tutorial.R` — code for the tutorial (scripts / notebooks)
- `filtered_gene_bc_matrices` — input data (available online)
- `.pdfs` — generated figures and outputs

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
- **Markers** (`FindAllMarkers`) + basic interpretation
- **Common plots**: violin plots, feature plots, dimensional reductions, dot plots

---

## Before Getting Started

Make sure you have **R** installed and a working R environment.  
Recommended:

- R (>= 4.1)
- Seurat (v4 or v5)
- ggplot2, patchwork (often pulled in by Seurat)
- Optional: `cowplot`, `dplyr`

---

## Option A: Quick Install (install.packages)

In an R session:


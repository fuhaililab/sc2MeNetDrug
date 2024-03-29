---
layout: default
title: Clustering
permalink: /Clustering/
nav_order: 6
---

# Clustering

## Introduction

This section becomes accessible only after uploading data in the "Upstream Analysis Data" and completing the dimension reduction. Clustering analysis enables the grouping of scRNA-seq data into distinct clusters, potentially representing different cell types.

<p align="center"><img src="../pic/clusteringPanel.png" alt="Clustering panel" style="zoom:80%;" /></p>

1. First, you need to select two parameters. The first is "the number of PCs", which determines the number of PCs used to identify neighboring cells. The second is "resolution", which dictates the granularity of clustering. A lower value typically yields fewer clusters, while a higher value produces more clusters.

2. Next, click the blue "Run clustering" button to initiate the analysis. The results will be displayed in the right panel.

  <p align="center"><img src="../pic/clusteringResult.png" alt="Clustering result" style="zoom:67%;" /></p>


## Data

After clustering, you will obtain one `.RData` files in your working directory:

* `clusters.RData`:  An data.frame variable named `clusters` saves the clustering result for each cell in the dataset.


## Advanced Hyper-parameter Tuning
All main functions used in clustering module can be located in `R/clustering.R`. Users can adjust all hyper-parameters used in clustering in this file.

For KNN algorithm, change the parameters by changing the following code in the file:

```R
seurat_data <- FindNeighbors(seurat_data,
                               reduction="pca",
                               dims = 1:nPC,
                               verbose = FALSE)
```

Please see more information and changeable parameters about this function in [document](https://satijalab.org/seurat/reference/findneighbors).

For Louvain algorithm, change the parameters by changing the following code in the file:

```R
seurat_data <- FindClusters(seurat_data,
                              resolution = resolution,
                              verbose = FALSE)
```

Please see more information and changeable parameters about this function in [document](https://satijalab.org/seurat/reference/findclusters).

**Importance:** After modifying the file, please make sure to restart the application to let modified parameters to be effective.

## Methodology

The clustering analysis in sc2MeNetDrug is conducted using the FindNeighbors and FindClusters functions from the Seurat package. The FindNeighbors function applies KNN on the selected PCs to gather neighboring information for each cell. Subsequently, the FindClusters function leverages this information to execute the Louvain algorithm, generating the clustering results.
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

## Methodology

The clustering analysis in sc2MeNetDrug is conducted using the FindNeighbors and FindClusters functions from the Seurat package. The FindNeighbors function applies KNN on the selected PCs to gather neighboring information for each cell. Subsequently, the FindClusters function leverages this information to execute the Louvain algorithm, generating the clustering results.
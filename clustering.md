---
layout: default
title: Clustering
permalink: /Clustering/
nav_order: 7
---

# Clustering

## Introduction

This part is only avaliable after you upload data in the "Upstream analysis Data" part and finish the dimension reduction part. Clustering analysis allows you to cluster scRNA-seq data into several clusters, which may indicate different cell types. There are two steps in clustering: main clustering and sub-clustering.

First, you can choose two parameters, The first is "eps", which is the size of the epsilon neighborhood. The second is "minPts", which is the number of minimum points in the eps region. Then, after clicking the blue button "Pre-cluster", the algorithm will order the data according to the parameters you selected. 

<p align="center"><img src="../pic/preCluster.png" alt="preCluster" style="zoom:80%;" /></p>

Next, select the parameter "threshold" to cut and obtain the main clustering result. Generally, a smaller threshold will result in more clusters and more outlier points. After main clustering, you can see the results in a plot on the right.

<p align="center"><img src="../pic/mainCluster.png" alt="mainCluster" style="zoom:67%;" /></p>

After main clustering, you can use GMM to further cluster data. Click the blue button "Sub-clustering" to obtain GMM clustering results.

<p align="center"><img src="../pic/subCluster.png" alt="subCluster" style="zoom:80%;" /></p>

## Data

After clustering, you will obtain three `.RData` files in your working directory:

* `preClusterResult.RData`:  An 'optics' object in variable `dfDensity` saves the ordering information. For more information, you can see [dbscan in R](https://cran.r-project.org/web/packages/dbscan/dbscan.pdf).

* `mainClusterResult.RData`: An 'optics' object saves the main clustering results in a variable named `dfClust`, you can obtain the clustering results by:

  ```R
  dfClust$cluster
  ```

* `subClusteringResult.RData`: A data frame that saves sub-clustering results in a variable named as `sub_df_cluster`.  The data frames have three columns. The first two columns are 2-dimensional cell data after dimension reduction, and the third column is cluster information. 

## Methodology

Clustering analysis in sc2MeNetDrug involves the two steps described below.

### Main clustering

sc2MeNetDrug uses the OPTICS<sup>1</sup> algorithm in main clustering. This is done by the`optics` and `extractDBSCAN` functions in the R package `dbscan`. 

### Sub-clustering

sc2MeNetDrug uses the Gaussian Mixture Model (GMM) in sub-clustering. First, the application will discover big clusters in main clustering and run the GMM algorithm in this cluster. This is done by the `Mclust` function in the`mclust` R package with parameters set as `model="EII"`. The number of clusters in the GMM is set based on the size of main cluster.





## References

1. Ankerst, M., Breunig, M., Kriegel, H.-P. & Sander, J. *OPTICS: Ordering Points to Identify the Clustering Structure*. *Sigmod Record* vol. 28 (1999).
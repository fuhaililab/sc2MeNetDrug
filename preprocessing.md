---
layout: default
title: Preprocessing
permalink: /preprocessing/
nav_order: 5
---

# Preprocessing

## Introduction

Once the required data has been successfully uploaded, the preprocessing analysis can be performed. Preprocessing allows you to do quailty control, normalization, and imputation for the read count data set. Imputation is an optional step in preprocessing as it takes much longer than other steps, but it is recommended because it gives more reasonable results in following analyses. Note that **if you have a group or design file, you should upload it before performing preprocessing.** 

<p align="center"><img src="../pic/preprocessing.png" alt="preprocessing" style="zoom:50%;" /></p>

If you want to do imputation, check the imputation checkbox. Then, click the blue button "Preprocessing" to start preprocessing. Once preprocessing finishes, the results are displayed at the top of the application:

<p align="center"><img src="../pic/preprocessingFinish.png" alt="preprocessingFinish" style="zoom:50%;" /></p>

## Data

After preprocessing, you will see at least two `.RData` files in your working directory. If you upload data in the "Upstream Analysis data" part, you will have:

* `rna_df.RData`: Saves the read count data in a Seurat object after preprocessing. The name of variable is `seurat_data`.  Data in matrix format from the Seurat object can be obtained through:

  ```R
  #If you don't do imputation in preprocessing
  seurat_data[["RNA"]]@data
  #If you do imputation in preprocessing
  seurat_data[["ALRA"]]@data
  #if you are not sure which one, you can just call
  seurat_data[[seurat_data@active.assay]]@data
  ```

  For more information, you can see [Seurat Website](https://satijalab.org/seurat/).

* `rnaPreprocessingResult.RData`: This file saves two variables named `keep_cell_index` and `keep_gene_index`, which are vectors that save the cell and gene index in the original dataset that are kept after preprocessing. 

* `rnaGroupInformation.RData`: If you upload group or design information, this file will save this information after preprocessing in a variable named `group_list`.

If you upload data in the "Downstream Analysis Data" part, you will have:

* `network_df.RData` Saves the read count data in a Seurat object after preprocessing. The name of this variable is `seurat_data`. 
* `networkPreprocessingResult.RData`: This file saves two variables named  `keep_cell_index` and `keep_gene_index`, which are vectors that save the cell and gene index in the original dataset that are kept after preprocessing. 
* `networkGroupInformation.RData`: If you upload group or design information, this file will save this information after preprocessing with the variable named `group_list`.

## Methodology

The preprocessing for scRNA-seq read count data in sc2MeNetDrug consists of three steps: quality control, normalization and imputation.

### Quality Control

Quality control is done in several steps. First, cells with total counts less than the threshold will be removed. The threshold value is computed by 0.012\\(\times\\)​(number of genes). Then, cells with expressed genes less than the threshold (at least 1 read) will be removed. The threshold is computed by 0.012​\\(\times\\)(number of genes). Next, cells with abnormally high ratios of counts mapping to 34 mitochondrial genes (relative to the total number of genes) will be removed. To be specific, we have soft and hard thresholds to discover abnormal cells. The soft threshold to compute the total mitochondrial expression ratio is calculated by:
\\[\frac{\text{Total counts in mitochondrially encoded genes}}{\text{Total counts in all genes}}\\]

Then, a K-means clustering algorithm is applied to cluster cells based on this ratio with the number of clusters set to k=2. If one cluster has a number of cells larger than 5 times the number of cells in another cluster, and this cluster has a lower mean mitochondrial expression ratio, the cells in this cluster will be kept and the cells in the other will be removed. Otherwise, we will apply the hard threshold. First, the 98% quantile of mitochondrial expression ratio of the whole data set is computed, and if the ratio is larger than 0.09, the threshold is set as 0.09, otherwise the threshold is set as this ratio. Finally, all the cells that have a mitochondrial expression ratio larger than the threshold are removed. The fourth step of quality control is to remove all mitochondrial encoded genes.

### Normalization

To normalize scRNA-seq read count data, the read count value for gene $X$ in one cell sample must be scaled using the following expression:

\\[\text{scaled expression for gene X}=\frac{\text{Read count for gene X} }{\text{Total count of cell sample}}\times10000\\]

Then, the data is transformed into log space using the natural logarithm. This is done by the `NormalizeData` function in the`Seurat` package<sup>1</sup>.

### Imputation

Imputation is done by the `runALRA` function in the`Seurat` package with default parameters. This method<sup>2</sup> computes the k-rank approximation to A_norm and adjusts it according to the error distribution learned from the negative values.



## References

1. Stuart, T. *et al.* Comprehensive Integration of Single-Cell Data. *Cell* (2019) doi:10.1016/j.cell.2019.05.031.

2. Linderman, G. C., Zhao, J. & Kluger, Y. Zero-preserving imputation of scRNA-seq data using low-rank approximation. *bioRxiv* 397588 (2018) doi:10.1101/397588.




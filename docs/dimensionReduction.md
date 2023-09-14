---
layout: default
title: Dimension Reduction
permalink: /dimensionReduction/
nav_order: 5
---

# Dimension Reduction

## Introduction

This section becomes accessible only after you've uploaded data under the "Upstream Analysis Data" and completed preprocessing. Dimension reduction enables the conversion of high-dimensional scRNA-seq data into a lower-dimensional space, facilitating visualization and subsequent analysis.

<p align="center">
  <img src="../pic/dimensionReductionPanel.png" alt="Dimension reduction panel" style="zoom:50%;" />
  </p>

1. First, choose the number of PCs in the panel to determine how many principal components will be utilized by the UMAP dimension reduction algorithm. (For more details, refer to the Methodology section).

2. Next, click "Run PCA + UMAP" to start dimension reduction analysis.

3. After dimension reduction, we will see visualization results in the right panel:

    <p align="center"><img src="../pic/dimensionReductionResult.png" alt="Dimension reduction result" style="zoom:60%;" /></p>

4. You can adjust the number of PC to obtain the best results. 
## Data

After dimension reduction, you will receive two `.RData` files in your working directory:

* `encoder_result.RData`: Saves the data after the auto encoder in a variable named `encoder_result`. The data is a data frame with 64 columns representing 64 dimensions, and number of rows is equal to the number of samples after preprocessing.
* `tsne_result.RData`: Saves the data after the T-SNE algorithm in a variable named  `tsne_result`.  The data is a data frame with 2 columns representing 2 dimensions, and the number of rows is equal to the number of samples after preprocessing.


## Methodology

Dimension reduction analysis in sc2MeNetDrug involves several steps:

1. First, select the top 3000 variable genes. To identify these genes, local polynomial regression fits the relationship between log variance and log mean. Subsequently, gene expression values are standardized using the observed mean and the expected variance (determined by the fitted line). The variance of gene expression is calculated on the standardized values after clipping. This procedure is automatically executed by the `SCTransform` function in the Seurat package.

2. Next, Principal Components Analysis (PCA) is applied to these 3000 variable genes. These genes are then projected into 50 dimensions in order served as 50 different principal components (PCs). This is implemented using `RunPCA` function in Seurat.

3. Finally, the UMAP<sup>1</sup> method will be used on the first \\(x\\) PCs and further project data into 2 dimensions, where \\(x\\) is the number of PCs selected by the user. This is implemented using `RunUMAP` function in Seurat.

## Advanced Hyper-parameter Tuning
All main functions used in dimension reduction module can be located in `R/dimensionReduction.R`. Users can adjust all hyper-parameters used in dimension reduction in this file.
For PCA, change the parameters by changing the following code in the file:

```R
seurat_data <- RunPCA(seurat_data, verbose = FALSE)
```

Please see more information and changeable parameters about this function in [document](https://satijalab.org/seurat/reference/runpca).

For UMAP algorithm, change the parameters by changing the following code in the file:

```R
seurat_data <- RunUMAP(seurat_data, reduction="pca", dims = 1:nPC, verbose = FALSE)
```
Please see more information and changeable parameters about this function in [document](https://satijalab.org/seurat/reference/runumap).

**Importance:** After modifying the file, please make sure to restart the application to let modified parameters to be effective.

## References

1. McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018





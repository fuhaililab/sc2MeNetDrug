---
layout: default
title: Dimension Reduction
permalink: /dimensionReduction/
nav_order: 6
---

# Dimension Reduction

## Introduction

This part will only be avaliable after you upload data in the "Upstream analysis Data" part and finish preprocessing. Dimension reduction allows you to reduce high dimensional scRNA-seq data to low dimensional space for visualization and further analysis.

<p align="center">
  <img src="../pic/ae.png" alt="ae" style="zoom:50%;" />
  </p>

First, click the blue button "Auto Encoder" to train the auto encoder.Then, select parameters in the T-SNE panel and perform the T-SNE computation. Notice that you can also skip the auto encoder and directly perform T-SNE. If you want to T-SNE directly, you can uncheck the "Use auto encoder result" check box.

"Maximum iteration time" controls the maximum iteration time of the algorithm. "Perplexity number" controls the number of close neighbors each point has in algorithm assumption. 

<p align="center"><img src="../pic/tsne.png" alt="tsne" style="zoom:50%;" /></p>

You can adjust these two parameters to obtain the best results. After dimension reduction, we will see visualization results:

<p align="center"><img src="../pic/dimensionReduction.png" alt="DimensionReduction" style="zoom:60%;" /></p>

## Data

After dimension reduction, you will receive two `.RData` files in your working directory:

* `encoder_result.RData`: Saves the data after the auto encoder in a variable named `encoder_result`. The data is a data frame with 64 columns representing 64 dimensions, and number of rows is equal to the number of samples after preprocessing.
* `tsne_result.RData`: Saves the data after the T-SNE algorithm in a variable named  `tsne_result`.  The data is a data frame with 2 columns representing 2 dimensions, and the number of rows is equal to the number of samples after preprocessing.

## Methodology

Dimension reduction analysis in sc2MeNetDrug involves several steps:

First, select first 2048 variable genes. To select variable genes, use local polynomial regression to fit the relationship of log(variance) and log(mean). Then, standardize the gene expression values using the observed mean and expected variance (given by the fitted line). Gene expression variance is then calculated on the standardized values after clipping. This is done by `FindVariableFeatures` with `selection.method` set as `"vst"` in the`Seurat` package<sup>[1]</sup>.

Next, use the auto encoder to reduce the dimensions from 2048 to 64. First, use min-max normalization based on the genes to normalize the data in these 2048 genes. Then, train the data in the auto encoder. The structure of the auto encoder is described as follows. In the encoder part, we have four dense layers with output dimensions 1024, 512, 128 and 64. After each dense layer, a batch normalization layer is added to speed up convergence. After the second and third dense layer, we add a dropout layer with drop out percentages set as 0.2 and 0.3. In the decoder part, we have four dense layers with output dimensions 128, 512, 1024 and 2048. After each dense layer, a batch normalization layer is added. The activation function of each layer is ReLU. In the training part, we set the loss function as MSE, the optimizer as “Adam”, the epoch to be 15 and batch size to be 128. The auto encoder is done using the`Keras` package for R.

Finally, T-SNE is used to further reduce the data to two dimensions. If the user checks the "Use auto encoder result" checkbox, the application will use the output of the encoder as the input of the T-SNE algorithm, otherwise, the application will directly use the normalized data with the first 2048 variable genes. T-SNE is done using the`Multicore T-sne `R package.



## References

1. Stuart, T. *et al.* Comprehensive Integration of Single-Cell Data. *Cell* (2019) doi:10.1016/j.cell.2019.05.031.





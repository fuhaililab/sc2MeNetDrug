---
layout: default
title: Gene Feature Exploration
permalink: /geneFeatureExploration/
nav_order: 7
---

# Gene Feature Exploration

## Introduction

This section becomes accessible only after uploading data under "Upstream Analysis Data" and finishing the clustering analysis. It enables users to examine the expression patterns of all genes in the dataset, facilitating the identification of potential marker genes for each cluster.

<p align="center"><img src="../pic/geneFeaturePanel.png" alt="Gene feature exploration panel" style="zoom:80%;" /></p>

1. In the left panel, you can choose the gene you want to explore by select it in the panel or by direct type its name to find it. You can select multiple genes at the same time.

2. Next, click the blue "Generate expression plot" button to produce the visualizations. Two distinct plots will be displayed. The first is a violin plot illustrating the expression distribution within each cluster. The second is a density map showcasing the expression distribution across the entire dataset.

   <p align="center"><img src="../pic/geneFeatureResult.png" alt="Gene feature exploration result" style="zoom:80%;" /></p>
   
3. Once the clustering analysis is done, the sc2MeNetDrug will automatically select the top 8 variable genes in the dataset and plot it.

## Data

The result data of two plots will be saved in the working directory:

* `gene_expression_vln_plot.RData`:  An plot variable named `gene_expression_vln_plot` saves the data of gene expression distribution violin plot.
* `gene_expression_scatter_plot.RData`:  An plot variable named `gene_expression_scatter_plot` saves the data of gene expression density plot.

You can see the these plot directly in the Rstudio by entering:
   ```R
   gene_expression_vln_plot
   ```
   ```R
   gene_expression_scatter_plot
   ```

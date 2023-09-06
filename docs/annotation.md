---
layout: default
title: Cell Annotation
permalink: /annotation/
nav_order: 9
---

# Cell Annotation

## Introduction

This section becomes accessible only after uploading data under the "Upstream Analysis Data" and completing the clustering analysis. Cell annotation enables the identification of cell types for each cell sample in the dataset, using clustering outcomes and expression data.

<p align="center"><img src="../pic/cellAnnotationPanel.png" alt="Cell annotation panel" style="zoom:40%;" /></p>

1. First, in the cell annotation panel above, select all the cell types to be included in the cell annotation analysis. Note that the list encompasses all the different cell types present in the marker gene table. If you've modified the table, this list will update automatically.

2. We also offer preset cell type combinations for Alzheimer's disease and Pancreatic cancer. You can apply these directly by clicking the blue buttons labeled "Alzheimer's disease" or "Pancreatic cancer".

3. After you select the cell type list, you can click the blue button "Cell Annotation" in the cell annotation section. 

4. After computation, you can see the results of cell annotation in :

    <p align="center"><img src="../pic/annotationResult.png" alt="Annotation result" style="zoom:80%;" /></p>

    We also provide cell annotation result for each cluster:

    <p align="center"><img src="../pic/annotationResult2.png" alt="Annotation result2" style="zoom:80%;" /></p>

5. Finally, we also provide a manual label correction panel in the bottom left side if you want to modify the result computed by the sc2MeNetDrug. You can modify the result by directly click the table cell and type the new cell type for each cluster.
       
    <p align="center"><img src="../pic/cellAnnotationCorrection.png" alt="Cell annotation correction panel" style="zoom:80%;" /></p>

## Data

After cell annotation, you will get three `.RData` files in your working directory:

* `annotation_es.RData` : This file stores the enrichment score results for each cluster within each cell marker gene set in a data frame variable named `annotation_es`. In this frame, each column represents a cell, while each row corresponds to a cluster.
* `annotation_result.RData`: The classification results for each cluster are stored in this file in a data frame variable named `annotation_result`. The first column represents the cluster, the second column indicates the corresponding cell type, and the third column combines the cluster and cell type for visualization purposes.
* `cell_annotation.RData`: This file saves the cell type for each cell in a data frame variable named as `cell_annotation`. 


## Methodology

Cell annotation in sc2MeNetDrug is computed using Gene Set Enrichment Analysis (GSEA)<sup>1</sup>. After the user selects candidate cell types in the "Biomarker Gene" section and starts computation, the application will compute the log fold change for cluster \(N\) using the following formula:

\\[\text{log fold change for cluster N}=\text{mean expression for cluster N}-{\text{mean expression for other cells}}\\]

Next, we rank the genes based on log fold changes and calculate the enrichment scores of marker gene sets for each prospective cell type. Ultimately, the cell type with the highest enrichment score is chosen as the designated type for that cluster. If none of the cell types yield a positive enrichment score, the cluster is labeled as unknown.


## References

1. Subramanian, A. *et al.* Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *Proc. Natl. Acad. Sci.* **102**, 15545 LP – 15550 (2005).
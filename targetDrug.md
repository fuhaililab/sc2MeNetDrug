---
layout: default
title: Drug Discovering and Clustering Based on Target and Chemical Structure
nav_order: 6
parent: Gene Expression, Communication and Drug Discovering

---

# Drug Discovering and Clustering Based on Target and Chemical Structure

## Introduction

You can find this analysis in the "Communication and Drug" section. Drug discovering and clustering based on targets and chemical structure analysis allows you to find potential drugs for diseases based on targets and chemical structure of drugs.

This part will be computed along with the cell-cell communication portion.More information can be found here: [cell-cell communication analysis](/cell-cellCommunication.md). You don't need to do any additional operation for this analysis, just set all the parameters for the cell-cell communication and Drug Discovering and Clustering Based on Signaling Signatures, then click blue button  "Generate Communication Network and Drug" to start computation.

<p align="center"><img src="pic/downstreamNetwork.png" alt="downstreamNetwork" style="zoom:40%;" /></p>

After computation, the drug discovering and clustering results will be shown after drug clustering based on signaling signatures. Drug discovering is based on the result of cell-cell communication. for each communication network, there will have three result. the first is also a drug network plot. the genes and edges in network is equal to cell-cell communication network. But here the size and transparency of genes indicate wheter specific nodes have drug target on them. a larger gene size and less transparency indicate this gene have drug target on it. You can view the drug list for specific gene by click the node in plot. 

<p align="center"><img src="pic/targetDrugNetwork.png" alt="targetDrugNetwork" style="zoom:67%;" /></p>

Below the network plot is drug summary table. The table contains all the drugs that have target gene on the network above. Table have five columns. The first column is drug name, the second column is target gene set for the drugs. the third is drug bank id for drugs. the fourth is description of drug. The last column is SMILES string for the drug.

<p align="center"><img src="pic/targetDrugTable.png" alt="targetDrugTable" style="zoom:50%;" /></p>

Below drug summary table is drug clustering result.  All the nodes that link by edges are from same cluster. The text in the plot shows the name of the drug.

<p align="center"><img src="pic/targetDrugClustering.png" alt="targetDrugClustering" style="zoom:50%;" /></p>

## Data

All data for drug discovering and clustering based on targets and chemical structure will be saved in "cellCommunication" . For specific cell-cell combination, data is saved in the "CellType1-CellType2" directory inside "cellCommunication". Inside the "CellType1-CellType2" file, you will see:

* "CellType1_CellType2":This directory saves data for drug discovering data from cell type 1 to cell type 2 communication network, where "CellType1" and  "CellType2" are cell types selected by the user. Inside the directory,  you can see:
  * `targetDrug.RData`: Saves all the drug discovering and clustering based on target and chemical structure analysis result in list variable `targetDrug1`. the variable have four member, the first is drug summary table for all the drugs discoverd in cell type1 to cell type2 communication network. The second is drug network data saved in `.json` format. The third is APcluster object saved APclustering result for drugs. For more information about APcluster object, you can see [APCluster](https://cran.r-project.org/web/packages/apcluster/vignettes/apcluster.pdf). The final is drug clustering result saved in `.json` format.
* "CellType2_CellType1":This directory saves data for drug discovering data from cell type 2 to cell type 1 communication network, where "CellType1" and  "CellType2" are cell types selected by the user. Inside the directory,  you can see:
  * `targetDrug.RData`: Saves all the drug discovering and clustering based on target and chemical structure analysis result in list variable `targetDrug2`. the variable have four member, the first is drug summary table for all the drugs discoverd in cell type2 to cell type1 communication network. The second is drug network data saved in `.json` format. The third is APcluster object saved APclustering result for drugs. For more information about APcluster object, you can see [APCluster](https://cran.r-project.org/web/packages/apcluster/vignettes/apcluster.pdf).The final is drug clustering result saved in `.json` format.

## Methodology

### Drug Discovering

 We collect totally 6650 drugs from drug bank database and target genes for all the drugs. After we compute and generate cell-cell communication network, we look through each gene in network and search drug that target on this gene. Finally, we will have drug list for each gene in the cell-cell communication network. Gether all the drugs and we can get drug summary table.

### Drug Clustering

Drug clustering for target drugs is based on the chemical structre of drugs<sup>[1]</sup>. first, we calculate fingerprint of drug by extend method. First, construct molecule object using SMILES information of drugs using `parse.smiles` function in R package `rcdk`. Next, get fingerprint information of drug using `get.fingerprint` function in R package `fingerprint`.  Based on fingerprint of drugs, we can compute similarity between drugs by Tanimoto index. The formulation of Tanimoto index is follow:

\\[S_{A,B}=\frac{c}{a+b-c}\\]

Where \\(S_{A,B}\\) is the similarity between drug \\(A​\\) and drug \\(B​,\\) \\(a​\\) is number of bits in drug \\(A\\)​ and \\(b\\) is number of bits in drug \\(B\\)​, \\(c\\) is number of bits in both two drugs. This is done by `fp.sim.matrix` function in R package `fingerprint` and set parameter `method="tanimoto"`.

Finally, we the similarity score between each two drug, we can construct similarity matrix and use AP clustering to clustering drugs. AP clustering is done using the R package `apcluster` .

## Reference

1. Voicu, A., Duteanu, N., Voicu, M., Vlad, D. & Dumitrascu, V. The rcdk and cluster R packages applied to drug candidate selection. *J. Cheminform.* **12**, 3 (2020).


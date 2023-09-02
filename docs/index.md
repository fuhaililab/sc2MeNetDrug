---
layout: page
title: Introduction
permalink: /
nav_order: 1
---

# sc2MeNetDrug

Single-cell RNA sequencing (scRNA-seq) data analysis has rapidly evolved in recent years, offering various methods and theories to gain insights into the complex biological processes in both animals and humans. However, a comprehensive analysis of scRNA-seq data typically involves multiple steps, and the results from one step can influence the outcomes of the next, making the analysis particularly challenging. Furthermore, most tools and packages designed for scRNA-seq data analysis only cover specific parts of the analysis process, necessitating researchers to download and become familiar with numerous packages.

We introduce sc2MeNetDrug, a tool designed to facilitate efficient, reliable, and convenient analysis of scRNA-seq data. In this tool, we encompass **data quality control, imputation, normalization, data dimension reduction, cell population clustering, cell type annotation, upstream network analysis, cell-cell communication network analysis, drug discovery, Gene Ontology analysis**, and many other useful methods for analyzing scRNA-seq data. We provide a user-friendly interface, detailed instructions, a powerful visualization tool, and streamlined data-saving and retrieval methods. Additionally, researchers can perform data analysis without the need for programming.

<span class="fs-8">
[Download](./downloadRequest.md){: .btn .btn-green .ml-auto }
</span>

## Overall Workflow

<p align="center"><img src="pic/overview.png" alt="overview" style="zoom:67%;" /></p>

The picture above shows the overall workflow of sc2MeNetDrug. All the analyses done in sc2MeNetDrug can be divided into two parts: upstream analysis and downstream analysis. In upstream analysis, the user can upload read count scRNA-seq data and the application will go through preprocessing, dimension reduction, clustering, and obtain cell annotations for each cell sample in the data set. All downstream analyses can be performed like upstream network analyses, GO analysis, cell-cell communication network analysis, and drug discovery. Meanwhile, if you already have cell annotation results, you can upload it with your read count data to perform downstream analysis directly.

sc2MeNetDrug has seven sections:

* **Upload Data**: In this section, data can be uploaded for upstream and downstream analysis and preprocessing.
* **Dimension Reduction**: Dimension reduction analysis is performed in this section.
* **Clustering**: Main clustering and sub-clustering are performed in this section.
* **Biomarker Gene**: The biomarker gene database can be viewed and edited in this section. 
* **Cell Annotation**: Cell annotation analysis is performed in this section.
* **Gene Expression**: Upstream network analysis, EMT-PRO analysis, and GO analysis can be performed in this section. In addition, you can set up the data set and ligand-receptor database at the top of this section.
* **Communication and Drug**: Cell-cell communication network analysis and drug discovering analysis are performed in this section.




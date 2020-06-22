---
layout: default
title: Biomarker Gene
permalink: /biomarkerGene/
nav_order: 8
---

# Biomarker Gene

## Introduction

In order to do cell annotation, you need marker genes for each cell type. Here we collect biomarker genes for more than 50 different cell types.

<p align="center"><img src="../pic/markerGeneSea.png" alt="markerGeneSea" style="zoom:80%;" /></p>

You can see the detailed gene list for cell type by clicking the bubble in the plot, or, you can find them in table:

<p align="center"><img src="../pic/markerGeneTable.png" alt="markerGeneTable" style="zoom:50%;" /></p>

Each column represents cell type and each row represents genes. If the value is 1, this gene is a marker of this cell type. Meanwhile, if you have your own marker gene, you can click the blue button "Add new gene" to add a new marker gene. Once you click the button, we can find the new row in the last page with the gene symbol name as "newGene", you can then edit the name and change the value from 0 to 1 in the corresponding cell type in the table directly. If you want to add a new cell type, you need to first enter the name of cell type in the text input frame and then click the button "Add new cell". Then you can scroll to the right of table and find your new cell type column. Notice that currently we don't support reloading edited marker gene information, which means if you want use your own data, you need add it every time you open the application. You can also delete or retrive original marker gene tables by clicking the "Delete selected gene" or the "Original gene table" button, respectively.


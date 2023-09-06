---
layout: page
title: Installation
permalink: /installation/
nav_order: 2
---

# Installation

Here we provide open source code of sc2MeNetDrug in [GitHub](https://github.com/fuhaililab/sc2MeNetDrug). 

## Requirements
1. Typically, sc2MeNetDrug may require a substantial amount of memory when processing large scRNA-seq datasets. Given the typical size of scRNA-seq data, we recommend using a computing resource with **at least 16GB of RAM** (32GB or more is highly recommended).
2. R environment with `R>= 4.2.1`.
3. RStudio or other applications that can deploy R Shiny apps.
4. All denpendences R packages:

   ```
   shinydashboard
   shiny
   shinyjs
   shinyBS
   plotly
   r2d3
   DT
   stringr
   shinyFiles
   fs
   data.table
   org.Hs.eg.db
   GO.db
   GOSemSim
   jsonlite
   dplyr
   Seurat
   SeuratWrappers
   igraph
   rcdk
   fingerprint
   glmGamPoi
   ```
   
5. If you are in a Mac environment, you may need environment that can compile c++ packages. You can install that by entering the following script in the terminal:
   
   `xcode-select --install`

## Instruction for the installation
1. In your terminal, `cd` to your desired directory and clone the source code by:
   `git clone https://github.com/fuhaililab/sc2MeNetDrug.git`

2. Open the source code in RStudio and run the `dependencies.R` script to install all dependent packages used in sc2MeNetDrug (ignore this step if you already have all packages installed).

3. Go to `R` directory and open `ui.R` or `server.R` in Rstudio, you will see **Run App** in the top right of the interface, click it to open the sc2MeNetDrug.



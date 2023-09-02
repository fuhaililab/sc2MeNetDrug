# sc2MeNetDrug
A user-friendly computational tool for analyzing scRNA-seq data to discover inter-cell signaling communication and drugs. In this tool, we include data quality control, imputation, normalization, data dimension reduction, cell population clustering, cell type annotation, upstream network analysis, cell-cell communication network analysis, drug discovery, Gene Ontology analysis, and many other useful methods for analyzing scRNA-seq data. We provide a user-friendly interface, instructions, a powerful visualization tool, and convenient data-saving and retrieving methods. In addition, researchers can analyze data without any programming. For more information and documentation, please visit [our website](https://fuhaililab.github.io/sc2MeNetDrug/).

### Overview
![overview](docs/pic/overview.png)

### Requirement
```
R >= 4.2.1
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

### Download
First, make sure you have the required R environment (>=4.2.1) and Rstudio IDE. Next, in the local terminal, `cd` to the desired directory and clone the sc2MeNetDrug repository by:
```
git clone https://github.com/fuhaililab/sc2MeNetDrug.git
```
Then, open the repository in the Rstudio and run the `dependencies.R` to install all dependent packages used in sc2MeNetDrug. Finally, go to R directory, open `ui.R` and click `run app` in the top right part of Rstudio to open the sc2MeNetDrug.

### Usage
The tool is built up on R shiny and powered by Seurat for easy use and better visualization. We provide detailed instruction for using sc2MeNetDrug in [documents](https://fuhaililab.github.io/sc2MeNetDrug/).


### Citation
If you find the tool helpful, considering cite our paper:
```
@article {Feng2021sc2menetdrug,
	author = {Jiarui Feng and S. Peter Goedegebuure and Amanda Zeng and Ye Bi and Ting Wang and Philip Payne and Li Ding and David DeNardo and William Hawkins and Ryan C. Fields and Fuhai Li},
	title = {SC2MeNetDrug: A computational tool to uncover inter-cell signaling targets and identify relevant drugs based on single cell RNA-seq data},
	elocation-id = {2021.11.15.468755},
	year = {2021},
	doi = {10.1101/2021.11.15.468755},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2021/11/19/2021.11.15.468755},
	eprint = {https://www.biorxiv.org/content/early/2021/11/19/2021.11.15.468755.full.pdf},
	journal = {bioRxiv}
}

```



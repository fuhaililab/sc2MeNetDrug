install.packages(c("shinydashboard", "shiny", "shinyjs", "shinyBS", "plotly", "r2d3", "DT",
                 "stringr", "shinyFiles", "fs", "data.table", "jsonlite", "dplyr", "Seurat",
                 "igraph", "fingerprint", "rcdk", "remotes", "R.utils"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("org.Hs.eg.db", "GO.db", "GOSemSim", "glmGamPoi"))
remotes::install_github('satijalab/seurat-wrappers')
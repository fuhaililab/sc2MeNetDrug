#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The main UI file of the app.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


library(shinydashboard)
library(shiny)
library(shinyjs)
library(shinyBS)
library(plotly)
library(r2d3)
library(DT)
library(stringr)
library(shinyFiles)
library(fs)
library(data.table)
library(org.Hs.eg.db)
library(GO.db)
library(GOSemSim)
library(dbscan)
library(jsonlite)
library(dplyr)
library(mclust)
library(Seurat)
library(SeuratWrappers)
library(igraph)
library(apcluster)
library(rcdk)
library(fingerprint)
library(glmGamPoi)
dashboardPage(
  title = "sc2MeNetDrug",
  skin = "green",
  
  #HEAD--------------------------------------
  
  dashboardHeader(title = "SC2MeNetDrug", titleWidth = 240),
  
  dashboardSidebar(
    sidebarMenu(
      #menuItem("About",tabName = "about",icon=icon("fas fa-bell")),
      menuItem(
        "Upload Data",
        tabName = "uploadData",
        icon = icon("fas fa-folder-open"),
        selected = TRUE
      ),
      menuItem(
        "Dimension Reduction",
        tabName = "dimensionReduction",
        icon = icon("fas fa-filter")
      ),
      menuItem(
        "Clustering",
        tabName = "cluster",
        icon = icon("fas fa-boxes")
      ),
      menuItem(
        "Gene Feature Exploration",
        tabName = "genePlot",
        icon = icon("fas fa-chart-bar")
      ),
      menuItem(
        "Biomarker Gene",
        tabName = "markerGene",
        icon = icon("fas fa-table")
      ),
      menuItem(
        "Cell Annotation",
        tabName = "classification",
        icon = icon("fas fa-layer-group")
      ),
      menuItem(
        "Gene Expression",
        tabName = "geneExpression",
        icon = icon("far fa-chart-bar")
      ),
      menuItem(
        "Communication and Drug",
        tabName = "network",
        icon = icon("fas fa-project-diagram")
      ),
      fluidRow(
        tipify(
          bsButton(
            inputId = "continueWork",
            label = "Continue Your Work",
            style = "primary",
            size = "default",
            icon = icon("fas fa-box")
          ),
          "You can only reload your last work session in working directory. Make sure you do not move any files generated in the directory, otherwise it may cause some unexpected errors!"
        ),
        id = "tabs"
      )
    ),
    width = 240
  ),
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet",
                type = "text/css",
                href = "page.css"),
      tags$style(
        "font-size:18px",
        ".shiny-output-error { visibility: hidden; }",
        ".shiny-output-error:before { visibility: hidden; }"
      )
      
    ),
    useShinyjs(),
    bsAlert("errorAlert"),
    tags$script(src = "saveSvg.js"),
    tabItems(
      tabItem(
        tabName = "uploadData",
        fluidRow(
          box(
            title = "Working Directory",
            status = "info",
            width = 12,
            solidHeader = TRUE,
            tags$h5(
              "You can set up or change your working directory here.
                Note that resetting the working directory will reload the entire application."
            ),
            tags$div(id = "currentDirBar",
                     verbatimTextOutput("currentDir")),
            fluidRow(
              shinyDirButton(
                "workingDir",
                label = "Select Working Directory",
                title = "Choose your working directory:",
                style =
                  "width:180px"
              ),
              tipify(
                bsButton(
                  inputId = "workingDirSet",
                  "Set Working Directory",
                  size = "default",
                  style = "info"
                ),
                "After selecting the directory, click this button to reload the application."
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "Uploading Upstream Analysis Data",
            status = "primary",
            width = 6,
            solidHeader = TRUE,
            fluidRow(
              id = "uploadRnaData",
              fileInput(
                "rnaDataFile",
                label = "Choose your read count data file",
                multiple = FALSE,
                width = '450px',
                placeholder = "Please select valid file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".rds",
                  "RDS"
                )
              ),
              fileInput(
                "rnaGroupFile",
                label = "Choose your group/design file",
                multiple = FALSE,
                width = '450px',
                placeholder = "Please select valid file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".rds",
                  "RDS"
                )
              )
              
            ),
            column(
              id = "rnaOption",
              width = 12,
              # Input: Checkbox if file has header ----
              checkboxInput("rnaHeader", "Header", TRUE),
              
              # Input: Select separator ----
              radioButtons(
                "rnaSep",
                "Separator",
                choices = c(
                  Comma = ",",
                  Semicolon = ";",
                  Tab = "\t"
                ),
                selected = ","
              ),
              
              # Input: Select quotes ----
              radioButtons(
                "rnaQuote",
                "Quote",
                choices = c(
                  None = "",
                  "Double Quote" = '"',
                  "Single Quote" = "'"
                ),
                selected = '"'
              )
            ),
            column(
              class = "preprocessing",
              width = 12,
              fluidRow(
                tags$h3("Preprocessing"),
                popify(
                  checkboxInput("rnaImputation", "Imputation", F, width = "100px"),
                  "Please be aware that the imputation would take a long time depending on your data size.",
                  placement = "bottom",
                  options = list(container = "body")
                ),
                popify(
                  checkboxInput("rnaMouseToHuman", "Mouse gene convert to human gene", F, width =
                                  "300px"),
                  "If the dataset is collected using mouse cells, please convert the mouse gene symbols to human gene symbols.",
                  placement = "bottom",
                  options = list(container = "body")
                )
              ),
              bsButton(
                inputId = "rnaPreprocessing",
                label = "Do Preprocessing",
                size = "default",
                style = "primary"
              )
            )
          ),
          box(
            title = "Instructions for Uploading Upstream Analysis Data",
            status = "primary",
            width = 6,
            collapsible = TRUE,
            solidHeader = TRUE,
            tags$h4(
              id = "rnaInstructionTitle",
              "To ensure proper functionality of the application,
              please carefully read the upload instructions:"
            ),
            tags$ol(
              class = "content",
              tags$li(
                "For read count data, you can upload a .RDS file contains dgCMatrix or Seurat object or .csv file."
              ),
              tags$li(
                "The read count data should be a sparse matrix or data frame where each row represents
                                          a gene and each column represents a cell sample."
              ),
              tags$li(
                'The row names of the read count data should consist of unique gene symbols.
                If your .csv data frame includes column names, please ensure to check the "Header" checkbox.'
              ),
              tags$li(
                "The group or design file is optional, but we recommend that you upload it to obtain
                meaningful analysis results. The data you upload should be a .csv file or .RDS file
                containing a data frame with only one column. The length of the column should be
                equal to the number of cells in the read count data. Please refrain from including
                row names in your uploaded data."
              )
            ),
            bsButton(
              inputId = "rnaDataExpend",
              label = "",
              icon = icon("search-plus", class = "opt"),
              style = "warning",
              size = "extra-small"
            )
          )
        ),
        fluidRow(
          box(
            title = "Uploading Downstream Analysis Data ",
            status = "warning",
            width = 6,
            solidHeader = TRUE,
            fluidRow(
              id = "uploadNetworkData",
              fileInput(
                "networkDataFile",
                label = "Choose your read count data file",
                multiple = FALSE,
                width = '450px',
                placeholder = "Please select valid file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".rds",
                  "RDS"
                )
              ),
              fileInput(
                "networkGroupFile",
                label = "Choose your cell annotation and group/design file",
                multiple = FALSE,
                width = '450px',
                placeholder = "Please select valid file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".rds",
                  ".RDS"
                )
              )
            ),
            column(
              id = "networkOption",
              width = 12,
              # Input: Checkbox if file has header ----
              checkboxInput("networkHeader", "Header", TRUE),
              
              # Input: Select separator ----
              radioButtons(
                "networkSep",
                "Separator",
                choices = c(
                  Comma = ",",
                  Semicolon = ";",
                  Tab = "\t"
                ),
                selected = ","
              ),
              
              # Input: Select quotes ----
              radioButtons(
                "networkQuote",
                "Quote",
                choices = c(
                  None = "",
                  "Double Quote" = '"',
                  "Single Quote" = "'"
                ),
                selected = '"'
              )
            ),
            
            column(
              class = "preprocessing",
              width = 12,
              fluidRow(
                tags$h3("Preprocessing"),
                popify(
                  checkboxInput("networkImputation", "Imputation", F, width = "100px"),
                  "Please be aware that the imputation would take a long time
                  depending on your data size.",
                  placement = "bottom",
                  options = list(container = "body")
                ),
                popify(
                  checkboxInput(
                    "networkMouseToHuman",
                    "Mouse gene convert to human gene",
                    F,
                    width = "300px"
                  ),
                  "If the dataset is collected using mouse cells,
                  please convert the mouse gene symbols to human gene symbols.",
                  placement = "bottom",
                  options =
                    list(container = "body")
                )
              ),
              bsButton(
                inputId = "networkPreprocessing",
                label = "Do Preprocessing",
                size = "default",
                style = "primary"
              )
            )
          ),
          box(
            title = "Instruction for Uploading Downstream Analysis Data",
            status = "warning",
            width = 6,
            collapsible = TRUE,
            solidHeader = TRUE,
            tags$h4(
              id = "networkInstructionTitle",
              "To ensure proper functionality of the application,
              please carefully read the upload instructions:"
            ),
            tags$ol(
              class = "content",
              tags$li(
                "For read count data, you can upload a .RDS file contains dgCMatrix or Seurat object or .csv file."
              ),
              tags$li(
                "The read count data should be a sparse matrix or data frame where each row represents
                                          a gene and each column represents a cell sample."
              ),
              tags$li(
                'The row names of the read count data should consist of unique gene symbols.
                If your .csv data frame includes column names, please ensure to check the "Header" checkbox.'
              ),
              tags$li(
                "The group or design file is required in the downstream analysis part.
                The data you upload should be a .csv file or .RDS file containing a
                data frame with one or two columns. The first column indicates the cell type for each cell.
                The second column is optional and can be used to specify a group or design for each cell.
                The length of the rows should be equal to the number of cells in the read count data.
                Please avoid including row names in your data."
              )
            ),
            bsButton(
              inputId = "geneExpressionExpend",
              label = "",
              icon = icon("search-plus", class = "opt"),
              style = "warning",
              size = "extra-small"
            )
          )
        ),
        fluidRow(
          box(
            title = "Processing Drug Data",
            status = "success",
            width = 6,
            solidHeader = T,
            tags$h5(
              "Before conducting cell-cell communication and drug discovery analysis,
              you need to follow the instructions to process the drug data file."
            ),
            fluidRow(
              shinyDirButton(
                "drugFileDir",
                label = "Select Drug File Directory",
                title = "Select Drug File Directory",
                style =
                  "width:180px"
              ),
              tipify(
                bsButton(
                  inputId = "drugFileProcessing",
                  "Process Drug Data",
                  size = "default",
                  style = "info"
                ),
                "Click the button to process the drug file. If you have already
                processed it before, there is no need to redo it after restarting
                the application."
              )
            )
          ),
          box(
            title = "Instruction for processing Drug File",
            status = "success",
            width = 6,
            colllapsible = T,
            solidHeader = T,
            tags$h4(
              id = "drugInstructionTitle",
              "To perform cell-cell communication and drug discovery analysis,
              you must download the drug data and process it according to the provided
              instructions."
            ),
            tags$ol(
              class = "content",
              tags$li(
                "Download file from website:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742. You can
                                    find more information about data in https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit."
              ),
              tags$li(
                "The data you need to download includes the following files:
                GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz,
                GSE92742_Broad_LINCS_gene_info.txt.gz,
                GSE92742_Broad_LINCS_sig_info.txt.gz,
                and GSE92742_Broad_LINCS_pert_info.txt.gz."
              ),
              tags$li(
                "After downloading the data, unzip the files and place all four data files in the same directory."
              ),
              tags$li(
                "Click on 'Select Drug File Directory' and choose the directory
                where you have placed all the data files. Afterward, click 'Process Drug Data'
                to initiate the processing."
              ),
              tags$li(
                "After the processing is complete, you should now be able to perform
                cell-cell communication and drug discovery analysis.
                If the 'Communication and Drug' section remains locked, try clicking
                'Load Data' in the 'Gene Expression' section to refresh the page."
              ),
            ),
            bsButton(
              inputId = "drugFileExpend",
              label = "",
              icon = icon("search-plus", class = "opt"),
              style = "warning",
              size = "extra-small"
            )
          )
        )
      ),
      
      
      tabItem(
        tabName = "dimensionReduction",
        fluidRow(
          id = "drNotification",
          box(
            title = "Notification",
            width = 12,
            status = "danger",
            solidHeader = TRUE,
            tags$h1(
              class = "notification",
              "Please upload the read count data in the upstream analysis section
              and perform preprocessing first."
            )
          )
        ),
        fluidRow(
          id = "dimensionReductionPart1",
          box(
            title = "Instruction for Dimension Reduction",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            tags$p(
              "Here, you can conduct dimension reduction analysis.
              The data will initially be projected into 50 dimensions using PCA analysis.
              Then, UMAP analysis will be applied to further project the data into 2 dimensions.
              You have the option to select the top k dimensions computed by PCA for the UMAP analysis."
            )
          )
        ),
        fluidRow(
          id = "dimensionReductionPart2",
          box(
            title = "PCA + UMAP",
            width = 4,
            status = "warning",
            solidHeader = TRUE,
            sliderInput(
              inputId = "nPC",
              label = "Select tge number of PCs:",
              min = 10,
              max = 50,
              value = 30,
              step = 1,
              width = '200px'
            ),
            bsButton("drButton", "Run PCA + UMAP", size = "default", style =
                       "primary")
          ),
          box(
            title = "Results Visualization",
            width = 8,
            status = "warning",
            solidHeader = TRUE,
            plotlyOutput("visualization", width = "auto", height =
                           "auto"),
            downloadButton(
              outputId = "drDownload",
              label = "Download Result Data",
              icon = icon("download")
            )
          ),
        )
      ),
      
      tabItem(
        tabName = "cluster",
        fluidRow(
          id = "clusterNotification",
          box(
            title = "Notification",
            width = 12,
            status = "danger",
            solidHeader = TRUE,
            tags$h1(class = "notification", "Please perform the dimension reduction analysis first.")
          )
        ),
        fluidRow(
          id = "clusterPart1",
          box(
            title = "Instruction for Clustering Analysis",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            tags$p(
              "Here, you can conduct the clustering analysis.
              The data will be clustered using the Louvain algorithm.
              You can choose to use the top k dimensions computed by
              PCA and select the clustering resolution to be used in the algorithm."
            )
          )
        ),
        fluidRow(
          id = "clusterPart2",
          box(
            title = "Clustering",
            status = "success",
            width = 4,
            solidHeader = TRUE,
            sliderInput(
              inputId = "neighbor_nPC",
              label = "Select the number of PCs:",
              min = 10,
              max = 50,
              value = 30,
              step = 1,
              width = '200px'
            ),
            sliderInput(
              inputId = "cluster_resolution",
              label = "Select the resolution:",
              min = 0.1,
              max = 2.0,
              value = 0.8,
              step = 0.1,
              width = '200px'
            ),
            
            bsButton(
              inputId = "ClusterButton",
              label = "Run clustering",
              size = "default",
              style = "primary"
            )
          ),
          box(
            title = "Results Visualization",
            status = "success",
            width = 8,
            solidHeader = TRUE,
            plotlyOutput("ClusterPlot", width = "auto", height =
                           "auto") ,
            downloadButton(
              outputId = "ClusterDownload",
              label = "Download Result Data",
              icon = icon("download")
            )
          )
        ),
        
      ),
      
      tabItem(
        tabName = "genePlot",
        fluidRow(
          id = "genePlotNotification",
          box(
            title = "Notification",
            width = 12,
            status = "danger",
            solidHeader = TRUE,
            tags$h1(class = "notification", "Please perform the clustering analysis first.")
          )
        ),
        fluidRow(
          id = "geneFeaturePlotPart1",
          box(
            title = "Instruction for Gene Feature Exploration",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            tags$p(
              "Here, you can explore the expression patterns of all genes across different
            clusters to identify the potential marker genes for each cluster."
            )
          )
        ),
        fluidRow(
          id = "geneFeaturePlotPart2",
          box(
            title = "Gene Expression Exploration",
            status = "success",
            width = 4,
            solidHeader = TRUE,
            uiOutput("genePlotSelection"),
            bsButton(
              "genePlotButton",
              "Generate expression plot",
              size = "default",
              style = "primary"
            )
          ),
          column(
            width = 8,
            
            box(
              title = "Gene Expression Violin Plot",
              status = "success",
              width = NULL,
              solidHeader = TRUE,
              plotOutput("geneExpressionViolinPlot")
            ),
            box(
              title = "Gene Expression Distribution Plot",
              status = "success",
              width = NULL,
              solidHeader = TRUE,
              plotOutput("geneExpressionScatterPlot")
            )
          )
        )
      ),
      
      tabItem(
        tabName = "markerGene",
        fluidRow(
          box(
            title = "Instruction for Biomarker Gene",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            tags$p(
              "In this section, you can explore the biomarker gene database used for cell annotation.
              You can also add or delete marker genes by editing the marker gene table. "
            )
          )
        ),
        fluidRow(
          box(
            title = "Biomarker Gene Table",
            width = 6,
            status = "success",
            solidHeader = TRUE,
            DTOutput("geneTable")
          ),
          
          box(
            title = "Biomarker Gene Database",
            width = 6,
            status = "warning",
            solidHeader = TRUE,
            d3Output("bubblePlot", width = "auto", height = "600px")
          ),
        ),
        
        fluidRow(
          box(
            title = "Updating Biomarker Gene Table",
            width = 6,
            status = "success",
            solidHeader = TRUE,
            bsButton(
              "deleteRowButton",
              "Delete selected gene",
              size = "default",
              style = "primary"
            ),
            popify(
              bsButton(
                "originalButton",
                "Original marker gene table",
                size = "default",
                style = "primary"
              ),
              "Replace current marker gene table with original marker gene table",
              placement = "bottom",
              options = list(container = "body")
            ),
            
            popify(
              bsButton(
                "saveMarkerButton",
                "Save current marker gene table",
                size = "default",
                style = "primary"
              ),
              "If you modify original marker gene table, you can save it to current working directory. Next time when
                you working in this directory, the sc2MeNetDrug will automatically use saved marker gene table.",
              placement = "bottom",
              options = list(container = "body")
            ),
          ),
          box(
            title = "Updating Biomarker Gene Table",
            width = 6,
            status = "success",
            solidHeader = TRUE,
            fluidRow(tags$div(
              id = "newCellBar",
              textInput(
                inputId = "newCellType",
                label = "Name of new cell type",
                width = "200px"
              ),
              uiOutput("markerGeneTableSelection")
            ), ),
            bsButton(
              "addRowButton",
              "Add new gene",
              size = "default",
              style = "primary"
            ),
          )
        )
      ),
      tabItem(
        tabName = "classification",
        fluidRow(
          id = "classificationNotification",
          box(
            id = "classificationNotification",
            title = "Notification",
            width = 12,
            status = "danger",
            solidHeader = TRUE,
            tags$h1(class = "notification", "Please do clustering analysis first.")
          )
        ),
        
        fluidRow(
          id = "classificationPart1",
          box(
            title = "Instruction for Cell Annotation",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            tags$p(
              "In this section, you can perform cell annotation analysis.
              Before conducting the annotation, you need to select all possible cell
              types to be used for cell annotation. If the cell type you want to include
              is not in our database, you can manually add it to the Biomarker Gene section.
              Start cell annotation by clicking 'Cell Annotation' buttom. The cell annotation
              analysis will initially assign each cluster a cell type based
              on the Gene Set Enrichment Analysis result using the biomarker gene database.
              Users can then manually modify the results based on their domain knowledge."
            ),
          )
        ),
        fluidRow(
          id = "classificationPart2",
          column(width = 4,
                 fluidRow(
                   box(
                     title = "Cell Type Selection",
                     width = 12,
                     status = "success",
                     solidHeader = TRUE,
                     tags$p("Select cell type you want to include in cell annotation."),
                     uiOutput("cellSelection"),
                     fluidRow(
                       bsButton(
                         inputId = "adButton",
                         label = "Alzheimer’s disease",
                         size = "default",
                         style = "primary"
                       ),
                       bsButton(
                         inputId = "pdacButton",
                         label = "Pancreatic Cancer",
                         size = "default",
                         style = "primary"
                       ),
                       bsButton(
                         inputId = "doClassifyButton",
                         label = "Cell Annotation",
                         size = "default",
                         style = "primary"
                       )
                     )
                   )
                 ),
                 fluidRow(
                   box(
                     title = "Mannual label correction",
                     width = 12,
                     status = "success",
                     solidHeader = TRUE,
                     tags$p("You can manually modify/assign cell annotation result here."),
                     DTOutput("annotationTable")
                   ),
                   
                 )),
          column(
            width = 8,
            box(
              width = 12,
              title = "Cell Annotation Results",
              status = "warning",
              solidHeader = TRUE,
              plotlyOutput("classifyPlot", width = "auto", height =
                             "auto")
            ),
            box(
              width = 12,
              title = "Annotation-Cluster Results",
              status = "warning",
              solidHeader = TRUE,
              plotlyOutput("classifyClusterPlot", width = "auto", height =
                             "auto"),
              downloadButton(
                outputId = "classificationDownload",
                label = "Download Result Data",
                icon = icon("download")
              )
            )
            
          )
          
        ),
      ),
      tabItem(
        tabName = "geneExpression",
        fluidRow(
          id = "geneExpressionNotification",
          box(
            title = "Notification",
            width = 12,
            status = "danger",
            solidHeader = TRUE,
            tags$h1(
              class = "notification",
              "Please upload the read count data in the downstream analysis section
              or perform the cell annotation analysis first."
            )
          )
        ),
        fluidRow(
          id = "geneExpressionPart1",
          box(
            title = "Instruction for Downstream Analysis",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            tags$p(
              "In this section, you can conduct upstream network analysis,
              EMT-PRO analysis, and Gene Ontology (GO) analysis."
            ),
            radioButtons(
              "dataUseSelect",
              "Which dataset should be used in all downstream analyses?",
              choices = c(
                "Result data from the upstream analysis" = 1,
                "Uploaded data from the downstream analysis." =
                  0
              ),
              selected = 1
            ),
            bsButton(
              "loadExpressionDataButton",
              "Load Data",
              size = "default",
              style = "primary"
            ),
            
            
            checkboxGroupInput(
              inputId = "ligRecData",
              label = "Choose Ligands-receptors intercations database
              used in all downstream analyses.",
              choices = c("DLRP", "nicheNet", "BaderLab"),
              choiceValues = c("DLRP", "nicheNet", "BaderLab"),
              width = "300px",
              selected = c("DLRP", "nicheNet", "BaderLab")
            )
            
          )
        ),
        fluidRow(
          id = "geneExpressionPart2",
          box(
            title = "Cell Distribution",
            width = 12,
            status = "warning",
            solidHeader = TRUE,
            plotlyOutput("cellDistribution", width = "auto", height =
                           "auto")
          )
        ),
        fluidRow(
          id = "geneExpressionPart3",
          box(
            title = "Upstream Network Analysis",
            width = 3,
            status = "success",
            solidHeader = TRUE,
            sliderInput(
              inputId = "upRegFc",
              label = "Please choose the Log-fold-change threshold:",
              min = 0.0,
              max = 1,
              value = 0.2,
              step = 0.01,
              width = '200px'
            ),
            sliderInput(
              inputId = "upRegPValue",
              label = "Please choose the P-value threshold:",
              min = 0.001,
              max = 1,
              value = 0.05,
              step = 0.001,
              width = '200px'
            ),
            tags$h4(id = "upRegTestTitle", "Test group"),
            uiOutput("upRegTestGroup"),
            tags$h4(id = "upRegNormalTitle", "Control group"),
            uiOutput("upRegNormalGroup"),
            popify(
              tags$div(
                checkboxInput("upRegUseOldResult", "Use saved calculation", T, width = "200px"),
                width = "200px"
              ),
              "If you have previously computed an analysis for this
              group combination and want to adjust the threshold,
              you can utilize the results of the differential genes test
              from the old analysis to save a significant portion of the
              calculations. However, avoid checking it if you have
              changed the group combination.",
              placement = "right",
              options = list(container = "body")
            ),
            bsButton(
              "upRegButton",
              "Upstream Network Analysis",
              size = "default",
              style = "primary"
            )
          ),
          column(
            width = 9,
            box(
              title = "Up-regulated Ligands",
              width = 12,
              status = "success",
              solidHeader = TRUE,
              collapsible = TRUE,
              d3Output("upRegLigandPlot", width = "auto", height =
                         "500px")
            ),
            box(
              title = "Up-regulated Receptors",
              width = 12,
              status = "success",
              solidHeader = TRUE,
              collapsible = TRUE,
              d3Output("upRegReceptorPlot", width =
                         "auto", height = "500px")
            )
          )
        ),
        fluidRow(
          id = "geneExpressionPart4",
          box(
            title = "Up-regulated Ligands to Expressed Receptors Network",
            width = 12,
            status = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            tags$div(class = "networkSupportBar",
                     d3Output(
                       "up_expPlot", width = "auto", height = "1200px"
                     )),
            bsButton(
              inputId = "upStreamPlotDownload1",
              label = "Download Network Chart",
              icon = icon("download"),
              size = "default",
              style = "info"
            )
          ),
          box(
            title = "Expressed Ligands to Up-regulated Receptors Network",
            width = 12,
            status = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            tags$div(class = "networkSupportBar",
                     d3Output(
                       "exp_upPlot", width = "auto", height = "1200px"
                     )),
            bsButton(
              inputId = "upStreamPlotDownload2",
              label = "Download Network Chart",
              icon = icon("download"),
              size = "default",
              style = "info"
            )
          ),
          box(
            title = "Up-regulated/Expressed Ligands to Up-regulated/Expressed Receptors Network",
            width = 12,
            status = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            tags$div(
              class = "networkSupportBar",
              d3Output("combinePlot", width = "auto", height =
                         "1200px")
            ),
            bsButton(
              inputId = "upStreamPlotDownload3",
              label = "Download Network Chart",
              icon = icon("download"),
              size = "default",
              style = "info"
            )
          ),
          box(
            title = "Up-regulated Ligands to Up-regulated Receptors Network",
            width = 12,
            status = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            tags$div(class = "networkSupportBar",
                     d3Output(
                       "up_upPlot", width = "auto", height = "1200px"
                     )),
            bsButton(
              inputId = "upStreamPlotDownload4",
              label = "Download Network Chart",
              icon = icon("download"),
              size = "default",
              style = "info"
            )
          )
        ),
        fluidRow(
          id = "geneExpressionPart5",
          box(
            title = "EMT-PRO Analytics",
            width = 3,
            status = "info",
            solidHeader = TRUE,
            tags$h4("Select Cell Type"),
            uiOutput("EPcellTypeSelect"),
            tags$h4(id = "EPGroupTitle", "Select group"),
            uiOutput("EPGroup"),
            bsButton(
              "epButton",
              "Compute EMT-PRO Score",
              size = "default",
              style = "primary"
            )
          ),
          box(
            title = "EMT-PRO Analysis Result",
            width = 9,
            status = "info",
            solidHeader = TRUE,
            plotlyOutput("emt_pro", height = "auto", width =
                           "auto")
          )
        ),
        fluidRow(
          id = "geneExpressionPart6",
          box(
            title = "GO Analytics",
            width = 3,
            status = "danger",
            solidHeader = T,
            tags$h4("Select Cell Type:"),
            uiOutput("GOCellType"),
            tags$h4(id = "GOTestTitle", "Test group"),
            uiOutput("GOTestGroup"),
            tags$h4(id = "GONormalTitle", "Control group"),
            uiOutput("GONormalGroup"),
            sliderInput(
              inputId = "GOUpFc",
              label = "Please choose the up-regulated log-fold-change threshold:",
              min = 0.0,
              max = 1,
              value = 0.2,
              step = 0.01,
              width = '200px'
            ),
            sliderInput(
              inputId = "GODnFc",
              label = "Please choose the down-regulated log-fold-change threshold:",
              min = -1,
              max = 0.0,
              value = -0.2,
              step = 0.01,
              width = '200px'
            ),
            sliderInput(
              inputId = "GOPValue",
              label = "Please choose the P-value threshold:",
              min = 0.001,
              max = 1,
              value = 0.05,
              step = 0.001,
              width = '200px'
            ),
            bsButton("goButton", "GO analytics", size = "default", style =
                       "primary")
          ),
          box(
            title = "GO Result",
            width = 9,
            status = "danger",
            solidHeader = T,
            selectInput(
              "GOTableSelection",
              label = h4("Please choose a table to display:"),
              choices = c(
                "Up-regulated GO Result" = "Up-regulated GO Result",
                "Down-regulated GO Result" =
                  "Down-regulated GO Result"
              ),
              selected = "Up-regulated GO Result"
            ),
            uiOutput("GOTableDisplay")
          )
        ),
        fluidRow(
          id = "geneExpressionPart7",
          box(
            title = "Select GO",
            width = 3,
            status = "danger",
            solidHeader = T,
            uiOutput("GONetworkSelect")
          ),
          box(
            title = "GO Analysis Result",
            width = 9,
            status = "danger",
            solidHeader = T,
            d3Output("GOPlot", width = "auto", height = "500px")
          )
        )
      ),
      tabItem(
        tabName = "network",
        fluidRow(
          id = "communicationNotification",
          box(
            title = "Notification",
            width = 12,
            status = "danger",
            solidHeader = TRUE,
            tags$h1(
              class = "notification",
              "To proceed with this section, please upload the
              read count data in the downstream analysis section
              or perform cell annotation analysis first.
              Additionally, ensure that you process the drug file
              in the upload data section."
            )
          )
        ),
        fluidRow(
          id = "communicationPart1",
          box(
            title = "Instruction for Downstream analysis",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            tags$p(
              "In this section, you can conduct inter-cell signaling
              communication and drug discovery analysis."
            )
          )
        ),
        
        fluidRow(
          id = "communicationPart2",
          box(
            title = "Communication and Drug Discovering",
            width = 4,
            status = "warning",
            solidHeader = TRUE,
            tags$h4("Cell type 1"),
            uiOutput("cellTypeSelect1"),
            tags$h4("Cell type 2"),
            uiOutput("cellTypeSelect2"),
            tags$h4(id = "networkTestTitle", "Test group"),
            uiOutput("controlGroup"),
            tags$h4(id = "networkNormalTitle", "Normal group"),
            uiOutput("normalGroup"),
            radioButtons(
              "networkResource",
              "Signaling resources",
              choices = c(KEGG = "KEGG",
                          STRING = "STRING"),
              selected = "KEGG"
            ),
            sliderInput(
              inputId = "networkFc",
              label = "Please choose the log-fold-change threshold:",
              min = 0.0,
              max = 1,
              value = 0.2,
              step = 0.01,
              width = '200px'
            ),
            sliderInput(
              inputId = "networkPValue",
              label = "Please choose the P-value threshold:",
              min = 0.001,
              max = 1,
              value = 0.05,
              step = 0.001,
              width = '200px'
            ),
            sliderInput(
              inputId = "drugNumber",
              label = "Please choose the number of top drugs:",
              min = 1,
              max = 200,
              value = 50,
              step = 1,
              width = '200px'
            ),
            checkboxInput("networkUseOldResult", "Use the saved results", T, width =
                            "200px"),
            checkboxInput(
              "networkPadjust",
              "Use the adjusted P-value instead of original P-value.",
              T,
              width =
                "200px"
            ),
            bsTooltip(
              "networkUseOldResult",
              "If you have previously conducted an analysis for this
              cell type combination and wish to adjust the threshold,
              you can employ the results of the differential genes test
              from the previous analysis to save a substantial amount of
              computation. However, avoid checking it if you changed the cell type combination.",
              placement = "right",
              options = list(container = "body")
            ),
            checkboxInput(
              "useFDAOnly",
              "Only use drugs in the DrugBank database for drug discovering.",
              T,
              width = "200px"
            ),
            bsButton(
              "communicationButton",
              "Generate Communication Network and Drug",
              size = "default",
              style = "primary"
            )
          ),
          box(
            title = "Gene Expression Information",
            width = 8,
            status = "warning",
            solidHeader = TRUE,
            selectInput(
              "networkTableSelection",
              label = h4("Please choose a table to display:"),
              choices = c(
                "Ligands of cell type 1 " = "ligType1",
                "Receptors of cell type 1" = "recType1",
                "Ligands of cell type 2" = "ligType2",
                "Receptors of cell type 2" = "recType2"
              ),
              selected = "Ligands of cell type 1"
            ),
            uiOutput("networkTableDisplay")
          )
        ),
        fluidRow(
          id = "communicationPart13",
          box(
            title = "Signaling Network Result Display Selection",
            width = 12,
            status = "danger",
            solidHeader = TRUE,
            collapsible = TRUE,
            selectInput(
              "networkDisplaySelect",
              label = h4("Please choose a network to display:"),
              choices = c(
                "Activated Signaling Pathways of Cell Type 1" = "activated1",
                "Activated Signaling Pathways of Cell Type 2" =
                  "activated2",
                "Inter-Cell Signaling Communication from Cell Type 1 to Cell Type 2" =
                  "downstream1",
                "Inter-Cell Signaling Communication  from Cell Type 2 to Cell Type 1" =
                  "downstream2"
              ),
              selected = "activated1"
            )
          )
        ),
        
        fluidRow(
          id = "communicationPart3",
          box(
            title = "Activated Signaliing Pathways",
            width = 12,
            status = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            uiOutput("activatedNetworkTitle1"),
            tags$div(
              class = "networkSupportBar",
              d3Output("activatedNetworkPlot1", width =
                         "auto", height = "1200px")
              
            ),
            bsButton(
              inputId = "activatedNetworkPlotDownload1",
              label = "Download the Network Chart",
              icon = icon("download"),
              size = "default",
              style = "info"
            )
          )
        ),
        fluidRow(
          id = "communicationPart4",
          box(
            title = "Activated Signaliing Pathways",
            width = 12,
            status = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            uiOutput("activatedNetworkTitle2"),
            tags$div(
              class = "networkSupportBar",
              d3Output("activatedNetworkPlot2", width =
                         "auto", height = "1200px")
            ),
            bsButton(
              inputId = "activatedNetworkPlotDownload2",
              label = "Download the Network Chart",
              icon = icon("download"),
              size = "default",
              style = "info"
            )
          )
        ),
        fluidRow(
          id = "communicationPart5",
          box(
            title = "Inter-Cell Signaling Communication",
            width = 12,
            status = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            uiOutput("networkTitle1"),
            tags$div(
              class = "networkSupportBar",
              d3Output("networkPlot1", width = "auto", height =
                         "1200px")
            ),
            bsButton(
              inputId = "networkPlotDownload1",
              label = "Download the Network Chart",
              icon = icon("download"),
              size = "default",
              style = "info"
            )
          )
        ),
        fluidRow(
          id = "communicationPart6",
          box(
            title = "Inter-Cell Signaling Communication",
            width = 12,
            status = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            uiOutput("networkTitle2"),
            tags$div(
              class = "networkSupportBar",
              d3Output("networkPlot2", width = "auto", height =
                         "1200px")
            ),
            bsButton(
              inputId = "networkPlotDownload2",
              label = "Download the Network Chart",
              icon = icon("download"),
              size = "default",
              style = "info"
            )
            
          )
        ),
        fluidRow(
          id = "communicationPart14",
          box(
            title = "Signaling Drug Result Display Selection",
            width = 12,
            status = "danger",
            solidHeader = TRUE,
            collapsible = TRUE,
            selectInput(
              "targetDrugDisplaySelect",
              label = h4("Select Specific network to display:"),
              choices = c(
                "Drug discovery results for inter-cell signaling communication from cell type 1 to cell type 2" =
                  "downstream1",
                "Drug discovery result for inter-cell signaling communication from cell type 2 to cell Type 1" =
                  "downstream2"
              ),
              selected = "downstream1"
            )
          )
        ),
        
        fluidRow(
          id = "communicationPart7",
          box(
            title = "Drug Discovering Based on Signaling Signatures",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            uiOutput("drugTableTitle1"),
            uiOutput("drugTableDisplay1")
          )
        ),
        fluidRow(
          id = "communicationPart8",
          box(
            title = "Drug Clustering Result for Signaling Signature Drug Discovery",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            d3Output("drugClusterPlot1", width = "auto", height =
                       "800px")
          )
        ),
        fluidRow(
          id = "communicationPart9",
          box(
            title = "Drug Discovery Based on Signaling Signatures",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            uiOutput("drugTableTitle2"),
            uiOutput("drugTableDisplay2")
          )
        ),
        fluidRow(
          id = "communicationPart10",
          box(
            title = "Drug Clustering Result for Signaling Signature Drug Discovery",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            d3Output("drugClusterPlot2", width = "auto", height =
                       "800px")
          )
        ),
        fluidRow(
          id = "communicationPart11",
          box(
            title = "Drug Discovery Based on Targeted Genes",
            width = 12,
            status = "success",
            solidHeader = T,
            uiOutput("targetDrugNetworkTitle1"),
            tags$div(
              class = "networkSupportBar",
              d3Output("drugMappingNetwork1", width = "auto", height =
                         "1200px")
            )
          ),
          box(
            title = "Drug Discovery Based on Targeted Genes",
            width = 12,
            status = "success",
            solidHeader = T,
            uiOutput("targetDrugTableTitle1"),
            uiOutput("drugMappingTable1Display")
          ),
          box(
            title = "Drug clustering Result for Drug Discovery Based on Targeted Genes",
            width = 12,
            status = "success",
            solidHeader = T,
            d3Output(
              "targetDrugClusterPlot1",
              width = "auto",
              height = "1200px"
            )
            
          )
        ),
        fluidRow(
          id = "communicationPart12",
          box(
            title = "Drug Discovering Based on Targeted Genes",
            width = 12,
            status = "success",
            solidHeader = T,
            uiOutput("targetDrugNetworkTitle2"),
            tags$div(
              class = "networkSupportBar",
              d3Output("drugMappingNetwork2", width = "auto", height =
                         "1200px")
            )
          ),
          box(
            title = "Drug Discovering Based on Targeted Genes",
            width = 12,
            status = "success",
            solidHeader = T,
            uiOutput("targetDrugTableTitle2"),
            uiOutput("drugMappingTable2Display")
          ),
          box(
            title = "Drug clustering Result for Drug Discovery Based on Targeted Genes",
            width = 12,
            status = "success",
            solidHeader = T,
            d3Output(
              "targetDrugClusterPlot2",
              width = "auto",
              height = "1200px"
            )
          )
        )
      )
      
    )
  )
)

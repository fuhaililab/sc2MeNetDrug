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
library(keras)
library(mclust)
library(Rtsne)
library(Seurat)
library(igraph)
library(apcluster)
library(rcdk)
library(fingerprint)
dashboardPage(
  title="sc2MeNetDrug",
  skin="green",
  
  #HEAD--------------------------------------
  
  dashboardHeader(title="SC2MeNetDrug",titleWidth=240),
  
  dashboardSidebar(
    sidebarMenu(
      #menuItem("About",tabName = "about",icon=icon("fas fa-bell")),
      menuItem("Upload Data",tabName = "uploadData",icon=icon("fas fa-folder-open"),selected=TRUE),
      menuItem("Dimension Reduction",tabName = "dimensionReduction",icon=icon("fas fa-filter")),    
      menuItem("Clustering",tabName = "cluster",icon=icon("fas fa-boxes")), 
      menuItem("Biomarker Gene",tabName = "markerGene",icon=icon("fas fa-table")), 
      menuItem("Cell Annotation",tabName = "classification",icon=icon("far fa-layer-group")), 
      menuItem("Gene Expression",tabName="geneExpression",icon=icon("far fa-chart-bar")),
      menuItem("Communication and Drug",tabName = "network",icon=icon("fas fa-project-diagram")), 
      fluidRow(
        tipify(bsButton(inputId="continueWork",label="Continue Your Work",style="primary",size="default",icon=icon("fas fa-box")),
               "You can only reload your last work session in working directory. Make sure you do not move any files generated in the directory, otherwise it may cause some unexpected errors!"),
        id="tabs"
      )
      ),
    width=240
  ),
  dashboardBody(
    tags$head(
      tags$link(
        rel = "stylesheet", 
        type = "text/css",
        href = "page.css"),
      tags$style("font-size:18px",         
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }")
      
    ),
    useShinyjs(),
    bsAlert("errorAlert"),
    tags$script(src = "saveSvg.js"),
    tabItems(
      # tabItem(tabName = "about",
      #         tags$div(id="aboutContainer",
      #                  tags$h2(id="about_title","SC2MeNetDrug"),
      #                  tags$p(id="about_introduction","SC2MeNetDrug is a tool that help you to analysis scRNA seq data in 
      #                         an efficient, reliable, and convenient way.
      #                         In this tool, we include data quality control, imputation, normalization, 
      #                         data dimension reduction, cell population clustering, cell type annotation,
      #                         upstream network analysis,cell-cell communication network analysis, drug discovering,
      #                         Gene Ontology analysis and many other useful methods to analyze scSeq-RNA data.")
      #                  ),
      #         tags$div(id="bottonContainer",
      # 
      #                  tags$div(id="contactContainer",
      #                           tags$p(id="contact","Contact"),
      #                           tags$p(id="email","E-mail:fuhai.li@wustl.edu")
      #                           ),
      #                  tags$div(id="linkContainer",
      #                           tags$p(id="link","Links"),
      #                           tags$p(id="article","Article Link:"),
      #                           tags$p(id="github","Github Link:")
      #                           )
      #                  )
      #         ),
      tabItem(tabName = "uploadData",
              fluidRow(
                box(title="Working Directory",status="info",width=12,solidHeader = TRUE,
                    tags$h5("You can set up your working directory here. 
                            Please notice that resetting the working directory will reload the whole application."),
                      tags$div(id="currentDirBar",
                               verbatimTextOutput("currentDir")  
                               ),
                    fluidRow(
                      shinyDirButton("workingDir",label="Select Working Directory",title="Choose your Working Directory",style=
                                       "width:180px"),
                      tipify(bsButton(inputId="workingDirSet","Set Working Directory",size="default",style="info"),
                             "After you select directory, click button to reset app.")
                    )
                    )
              ),
              fluidRow(
                      box(title="Upstream Analysis Data Uploading",status="primary",width=6,solidHeader = TRUE,
                          fluidRow(id="uploadRnaData",
                                   fileInput("rnaDataFile",label="Choose your read count data file",multiple=FALSE,
                                             width='450px',
                                             placeholder = "Please select valid file",
                                             accept = c("text/csv",
                                                        "text/comma-separated-values,text/plain",".csv",".rds","RDS")),
                                   fileInput("rnaGroupFile",label="Choose your group or design file",multiple=FALSE,
                                             width='450px',
                                             placeholder = "Please select valid file",
                                             accept = c("text/csv",
                                                        "text/comma-separated-values,text/plain",".csv",".rds","RDS"))
                                   
                          ),
                            column(id="rnaOption",width=12,
                                   # Input: Checkbox if file has header ----
                                   checkboxInput("rnaHeader", "Header", TRUE),
                                   
                                   # Input: Select separator ----
                                   radioButtons("rnaSep", "Separator",
                                                choices = c(Comma = ",",
                                                            Semicolon = ";",
                                                            Tab = "\t"),
                                                selected = ","),
                                   
                                   # Input: Select quotes ----
                                   radioButtons("rnaQuote", "Quote",
                                                choices = c(None = "",
                                                            "Double Quote" = '"',
                                                            "Single Quote" = "'"),
                                                selected = '"')
                            ),
                          column(class="preprocessing",width=12,
                            fluidRow(
                              tags$h3("Preprocessing"),
                              popify(checkboxInput("rnaImputation", "Imputation", T,width="100px"),                           
                                     "Normally, imputation would take several hours depend on you data
                                        size. Please be careful!",placement ="bottom",
                                     options=list(container="body")
                              ),
                              popify(checkboxInput("rnaMouseToHuman", "Mouse gene convert to human gene", F,width="300px"),                           
                                     "If the dataset is collected using mouse cell, please convert mouse gene to human gene.",placement ="bottom",
                                     options=list(container="body")
                              )
                            ),
                            bsButton(inputId = "rnaPreprocessing",label="Do Preprocessing",size="default",style="primary")
                          )
                          ),
                      box(title="Upstream Analysis Data Uploading Instruction",status="primary",width=6,collapsible = TRUE,solidHeader = TRUE,
                          tags$h4(id="rnaInstructionTitle","In order for the application to run properly,
                                  please read the upload instructions carefully:"),
                          tags$ol(class="content",
                                  tags$li("You can upload a csv or RDS file."),
                                  tags$li("Read count data should be a data frame with every row representing
                                          a gene and every column representing a cell sample."),
                                  tags$li("The row name of the read count data should be a unique gene symbol. 
                                          If your data has a column name, please check the “Header” checkbox."
                                          ),
                                  tags$li("Group or design file is optional, but we recommend you upload it in 
                                          order to receive reasonable analysis results. Data you upload should 
                                          be a data frame with only one column. The length of the row should be 
                                          equal to the number of cells in the read count data. Please do not include row name in your uploaded data.")
                                  ),
                          bsButton(inputId = "rnaDataExpend",label="",icon=icon("search-plus", class = "opt"),style="warning",size="extra-small")
                          )
                      ),
              fluidRow(
                box(title="Downstream Analysis Data Upload",status = "warning",width = 6,solidHeader = TRUE,
                    fluidRow(id="uploadNetworkData",
                             fileInput("networkDataFile",label="Choose your read count data file",multiple=FALSE,
                                       width='450px',
                                       placeholder = "Please select valid file",
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",".csv",".rds","RDS")),
                             fileInput("networkGroupFile",label="Choose your group or design file",multiple=FALSE,
                                       width='450px',
                                       placeholder = "Please select valid file",
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",".csv",".rds",".RDS"))
                    ),
                    column(id="networkOption",width=12,
                           # Input: Checkbox if file has header ----
                           checkboxInput("networkHeader", "Header", TRUE),
                           
                           # Input: Select separator ----
                           radioButtons("networkSep", "Separator",
                                        choices = c(Comma = ",",
                                                    Semicolon = ";",
                                                    Tab = "\t"),
                                        selected = ","),
                           
                           # Input: Select quotes ----
                           radioButtons("networkQuote", "Quote",
                                        choices = c(None = "",
                                                    "Double Quote" = '"',
                                                    "Single Quote" = "'"),
                                        selected = '"')
                    ),
                    
                    column(class="preprocessing",width=12,
                           fluidRow(
                             tags$h3("Preprocessing"),
                             popify(checkboxInput("networkImputation", "Imputation", T,width="100px"),                           
                                    "Normally, imputation would take several hours depend on you data
                                        size. Please be careful!",placement ="bottom",
                                    options=list(container="body")
                             ),                            
                             popify(checkboxInput("networkMouseToHuman", "Mouse gene convert to human gene", F,width="300px"),                           
                                                                   "If the dataset is collected using mouse cell, please convert mouse gene to human gene.",placement ="bottom",
                                                                   options=list(container="body")
                             )
                             ),
                             bsButton(inputId = "networkPreprocessing",label="Do Preprocessing",size="default",style="primary")
                           )
                    ),
                box(title="Downstream Analysis Data Upload Instruction",status = "warning",width=6,collapsible = TRUE,solidHeader = TRUE,
                    tags$h4(id="networkInstructionTitle","In order for the application to run properly,
                                  please read the upload instructions carefully:"),
                    tags$ol(class="content",
                            tags$li("You can upload a csv or RDS file."),
                            tags$li("Read count data should be a data frame with every
                                    row representing a gene and every column representing a cell sample."),
                            tags$li("The row name of the read count data should be a unique gene symbol. 
                                    If your data has a column name, please check the “Header” checkbox."
                            ),
                            tags$li("Group or design file in downstream analysis part is required. the data frame should have one or two columns,
                                  with first column indicate the cell type for each cell. Second column is optional, which is used to specify group or design for 
                                  each cell.The length of the row should be equal to the number of cells in the read count data. Please do not include rowname in your data.
                                    ")
                    ),
                    bsButton(inputId = "geneExpressionExpend",label="",icon=icon("search-plus", class = "opt"),style="warning",size="extra-small")
                )
              ),
              fluidRow(
                box(title="Drug File Processing",status="success",width=6,collapsible = T,solidHeader = T,
                    tags$h5("Before doing cell-cell communication and drug discovering analysis, you need to follow the 
                            instruction and process drug data file."),
                    fluidRow(
                      shinyDirButton("drugFileDir",label="Select Drug File Directory",title="Select Drug File Directory",style=
                                       "width:180px"),
                      tipify(bsButton(inputId="drugFileProcessing","Process Drug Data",size="default",style="info"),
                             "Click button to process drug file.I you have processed it before, no need to redo it after you 
                             restart application.")
                    )
                    ),
                box(title="Drug File Processing Instruction",status="success",width=6,colllapsible=T,solidHeader = T,
                    tags$h4(id="drugInstructionTitle","In order to do cell-cell communication and drug discovering 
                            analysis, you need to download drug data and process it follow the instruction"),
                    tags$ol(class="content",
                            tags$li("Download file from website:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742. You can 
                                    find more information about data in https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit."),
                            tags$li("The data you should download are GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz,
                                     GSE92742_Broad_LINCS_gene_info.txt.gz, GSE92742_Broad_LINCS_sig_info.txt.gz and GSE92742_Broad_LINCS_pert_info.txt.gz. "
                            ),
                            tags$li("After you download data, upzip file and put all four data file in same directory."),
                            tags$li("Click 'Select Drug File Directory' and select directory where you
                             put all the data in, then click 'process drug data' to start processing."),
                            tags$li("After processing, you should able to do cell-cell communication and drug discovering
                                     analysis now. If section 'Communication and Drug' don't unlock, try to click 'load data' in
                                     'Gene Expression' section to refresh."),
                            tags$li("If you use Windows or Mac version, the processed file will be automatically saved so you don't
                                     need to process it twice. However, if you use docker version, you need to commit the change after you
                                     process the file in order to save it. For more information, please visit our website:
                                    https://fuhaililab.github.io/sc2MeNetDrug/data/")
                    ),
                    )
              )
          ),
      
      tabItem(tabName = "dimensionReduction",
              fluidRow(id="drNotification",
                       box(title="Notification",width=12,status="danger",solidHeader = TRUE,
                           tags$h1(class="notification","Please upload read count data in upstream analysis part and do preprocessing first.")
                       )
              ),
              fluidRow(id="dimensionReductionPart1",
                box(title="Dimension Reduction Instruction",width=12,status = "primary",solidHeader = TRUE,
                    tags$p("Here you can do the dimension reduction analysis.There are two steps in dimension reduction part. First,
                          the 2048 most variable genes are found and the auto encoder is used to reduce the dimensions from 2048 to 64. Next, the T-SNE
                            is used to further reduce dimensions from 64 to 2. You can also skip the first step and use T-SNE to directly reduce 
                           the dimensions from 2048 to 64.")
                    )
                ),
              fluidRow(id="dimensionReductionPart2",
                  column(  width=4,
                    box(title="Auto Encoder",width=12,status = "warning",solidHeader = TRUE,
                        bsButton("aeButton","Auto Encoder",size="default",style="primary"),
                        downloadButton(outputId="adDownload",label="Download Result Data",icon=icon("download"))
                        ),
                    box(title="T-SNE",width=12,status ="warning",solidHeader = TRUE,
                     sliderInput(inputId = "maxIterationSlider",label="Maximum iteration time:",min=10,max=1000,value=500,step=10,width='300px'),
                     sliderInput(inputId = "perplexitySlider",label="Perplexity number:",min=1,max=500,value=30,step=1,width='300px'),
                     checkboxInput("useAutoEncoder", "Use auto encoder result", TRUE,width="auto"),
                     bsButton("drButton","T-SNE",size="default",style="primary")
                       )
                    ),
                       box(title="Result Visulization",width=8,status = "warning",solidHeader = TRUE,
                           plotlyOutput("visualization",width="auto",height="auto"),
                           downloadButton(outputId="drDownload",label="Download Result Data",icon=icon("download"))
                       )
              )
          ),
      
      tabItem(tabName = "cluster",
              fluidRow(id="clusterNotification",
                box(title="Notification",width=12,status="danger",solidHeader = TRUE,
                    tags$h1(class="notification","Please do dimension reduction analysis first.")
                )),
              fluidRow(id="clusterPart1",
                       box(title="Clustering Instruction",width=12,status = "primary",solidHeader = TRUE,
                           tags$p("Here you can do the clustering analysis. There are two steps in clustering. First, the OPTICS
                                   algorithm is used to cluster the data. Next, for each cluster, the GMM is used to further cluster data into 
                                  sub-cluster.")
                       )
              ),
              fluidRow(id="clusterPart2",
                       box(title="Pre-clustering",status ="warning",width=4,solidHeader = TRUE,
                           sliderInput(inputId = "eps",label="Select eps:",min=0.0001,max=4,value=1,step=0.05,width='200px'),
                           sliderInput(inputId = "minPts",label="Select minPts:",min=0.0001,max=80,value=5,step=1,width='200px'),
                           bsButton(inputId="preClusterButton",label="Pre-clustering",size="default",style="primary")
                           ),
                       box(title="Reachability Chart",status="warning",width=8,solidHeader = TRUE,
                           plotlyOutput("preClusterPlot",width="auto",height="auto")
                           )
                       ),
              fluidRow(id="clusterPart3",
                       box(title="Main Clustering",status="success",width=4,solidHeader = TRUE,
                           sliderInput(inputId = "eps_cl",label="Select threshold:",min=0.0001,max=2,value=1,step=0.02,width='200px'),
                           bsButton(inputId="mainClusterButton",label="Main Clustering",size="default",style="primary")
                           ),
                       box(title="Main Clustering Result",status = "success",width=8,solidHeader = TRUE,
                           plotlyOutput("mainClusterPlot",width="auto",height="auto") ,
                           downloadButton(outputId="mianClusterDownload",label="Download Result Data",icon=icon("download"))
                           )
                       ),
              fluidRow(id="clusterPart4",
                       box(title="Sub-clustering",status="info",width=4,solidHeader = TRUE,
                           tags$p("Based on the main clustering result, use GMM to further clustering data to sub-clusters."),
                           bsButton(inputId="subClusterButton",label="Sub-clustering",size="default",style="primary")
                       ),
                       box(title="Sub-clustering Result",status = "info",width=8,solidHeader = TRUE,
                           plotlyOutput("subClusterPlot",width="auto",height="auto") ,
                           downloadButton(outputId="subClusterDownload",label="Download Result Data",icon=icon("download"))
                       )
              )
              ),
      tabItem(tabName = "markerGene",
              fluidRow(
                box(title="Biomarker Gene Instruction",width=12,status = "primary",solidHeader = TRUE,
                    tags$p("This section display the marker genes database used for cell annotation. You can explore or edit it in this section.")
                )
                ),
              fluidRow(
                box(title="Biomarker Gene Database",width=12,status = "warning",solidHeader = TRUE,
                    d3Output("bubblePlot",width="auto",height="600px")
                )
              ),
              fluidRow(
                box(title="Biomarker Gene Table",width=8,status = "success",solidHeader = TRUE,
                    DTOutput("geneTable"),
                    fluidRow(
                      bsButton("addRowButton","Add new gene",size="default",style="primary"),
                      bsButton("deleteRowButton","Delete selected gene",size="default",style="primary"),
                      popify(bsButton("originalButton","Original marker gene table",size="default",style="primary"),
                             "Replace current marker gene table with original marker gene table",placement ="bottom",
                             options=list(container="body")),
                      popify(bsButton("saveMarkerButton","Save current marker gene table",size="default",style="primary"),
                            "If you modify original marker gene table, you can save it to current working directory. Next time when
                             you working in this directory, the sc2MeNetDrug will automatically use saved marker gene table.",placement ="bottom",
                            options=list(container="body"))
                             
                    ),
                    fluidRow(
                      tags$div(id="newCellBar",
                               textInput(inputId="newCellType",label="Name of new cell type",width="200px")
                      ),
                      tags$div(id="addColButtonBar",
                               bsButton("addColButton","Add new cell",size="default",style="primary")
                      )
                    )
                    ),
                box(title="Cell Type Selection",width=4,status = "success",solidHeader = TRUE,
                    tags$p("Select cell type you want to include in cell annotation."),
                    uiOutput("cellSelection"),
                    fluidRow(
                      bsButton(inputId="adButton",label="Alzheimer’s disease",size="default",style="primary"),
                      bsButton(inputId="pdacButton",label="Pancreatic Cancer",size="default",style="primary")
                    )
                    
                    )
              )
          ),
      tabItem(tabName = "classification",
              fluidRow(id="classificationNotification",
                box(id="classificationNotification",title="Notification",width=12,status="danger",solidHeader = TRUE,
                    tags$h1(class="notification","Please do clustering analysis first.")
                )
              ),
              fluidRow(id="classificationPart1",
                box(title="Cell Annotation Instruction",width=12,status = "primary",solidHeader = TRUE,
                    tags$p("Here you can do cell annotation analysis. Notice that you need to first select potential cell type
                           in dataset in biomarker gene section first. Otherwise the application will use all cell types in database for computation, which 
                           may not produce accurate results."),
                    checkboxInput("useSubCluster", "Use sub-clustering result", TRUE,width="auto"),
                    bsButton("doClassifyButton","Cell Annotation",size="default",style="primary"))
              ),
              fluidRow(id="classificationPart2",
                box(title="Cell Annotation Results",width=12,status = "warning",solidHeader = TRUE,
                    plotlyOutput("classifyPlot",width="auto",height="auto")
                    )
              ),
              fluidRow(id="classificationPart3",
                       box(title="Annotation-Cluster Results",width=12,status = "success",solidHeader = TRUE,
                           plotlyOutput("classifyClusterPlot",width="auto",height="auto"),
                           downloadButton(outputId="classificationDownload",label="Download Result Data",icon=icon("download"))
                       )
                    )
            ),
      tabItem(tabName = "geneExpression",
              fluidRow(id="geneExpressionNotification",
                       box(title="Notification",width=12,status="danger",solidHeader = TRUE,
                           tags$h1(class="notification","Please upload the read count data in the downstream analysis part,
                                   or perform the cell annotation analysis first.")
                       )
              ),
              fluidRow(id="geneExpressionPart1",
                       box(title="Downstream Analysis Instruction",width=12,status = "primary",solidHeader = TRUE,
                           tags$p("Here you can do upstreaem network analysis, EMT-PRO analysis and GO analysis."),
                           radioButtons("dataUseSelect", "Which one dataset should use in downstream analysis?",
                                        choices = c("Upstream analysis result data" = 1,
                                                    "Downstream analysis uploaded data."=0),
                                        selected = 1),
                           bsButton("loadExpressionDataButton","Load Data",size="default",style="primary"),
                           
                           
                           checkboxGroupInput(inputId="ligRecData",label="Choose Ligands-receptors intercations database.",
                                              choices = c("DLRP","nicheNet","BaderLab"),choiceValues=c("DLRP","nicheNet","BaderLab"),width="300px",
                                              selected=c("DLRP","nicheNet","BaderLab"))
                           
                       )
              ),
              fluidRow(id="geneExpressionPart2",
                       box(title="Cell Distribution",width=12,status="warning",solidHeader = TRUE,
                           plotlyOutput("cellDistribution",width="auto",height="auto")
                       )
              ),
              fluidRow(id="geneExpressionPart3",
                       box(title="Upstream Network Analysis",width=3,status = "success",solidHeader = TRUE,
                           sliderInput(inputId = "upRegFc",label="Select Log Fold Change Threshold:",min=0.01,max=2,value=0.2,step=0.01,width='200px'),
                           sliderInput(inputId = "upRegPValue",label="Select Pvalue Threshold:",min=0.001,max=1,value=0.05,step=0.001,width='200px'),
                           tags$h4(id="upRegTestTitle","Test group"),
                           uiOutput("upRegTestGroup"),
                           tags$h4(id="upRegNormalTitle","Control group"),
                           uiOutput("upRegNormalGroup"),
                           popify(tags$div(checkboxInput("upRegUseOldResult", "Use saved calculation", T,width="200px"),width="200px"),                           
                                  "If you have computed analysis for this group combination before and want to adjust threshold, you can use old different genes test result to save majority of calculation.
                                  But avoid use it if you change group combination.",placement ="right",
                                  options=list(container="body")
                           ),
                           bsButton("upRegButton","Upstream Network Analysis",size="default",style="primary")
                       ),
                       column(width=9,
                              box(title="Up-regulated Ligands",width=12,status = "success",solidHeader = TRUE,collapsible = TRUE,
                                  d3Output("upRegLigandPlot",width="auto",height="500px")
                              ),
                              box(title="Up-regulated Receptors",width=12,status = "success",solidHeader = TRUE,collapsible = TRUE,
                                  d3Output("upRegReceptorPlot",width="auto",height="500px")
                              )
                       )
              ),
              fluidRow(id="geneExpressionPart4",
                       box(title="Up-regulated Ligands to Expressed Receptors Network",width=12,status = "success",
                           solidHeader = TRUE,collapsible = TRUE,
                           tags$div(class="networkSupportBar",
                             d3Output("up_expPlot",width="auto",height="1200px")
                           ),
                           bsButton(inputId = "upStreamPlotDownload1",label="Download Network Chart",icon=icon("download"),size="default",style="info")
                       ),
                       box(title="Expressed Ligands to Up-regulated Receptors Network",width=12,status = "success",
                           solidHeader = TRUE,collapsible = TRUE,
                           tags$div(class="networkSupportBar",
                                    d3Output("exp_upPlot",width="auto",height="1200px")
                                    ),
                           bsButton(inputId = "upStreamPlotDownload2",label="Download Network Chart",icon=icon("download"),size="default",style="info")
                       ),
                       box(title="Up-regulated/Expressed Ligands to Up-regulated/Expressed Receptors Network",width=12,status = "success",
                           solidHeader = TRUE,collapsible = TRUE,
                           tags$div(class="networkSupportBar",
                             d3Output("combinePlot",width="auto",height="1200px")
                           ),
                           bsButton(inputId = "upStreamPlotDownload3",label="Download Network Chart",icon=icon("download"),size="default",style="info")
                       ),
                       box(title="Up-regulated Ligands to Up-regulated Receptors Network",width=12,status = "success",
                           solidHeader = TRUE,collapsible = TRUE,
                           tags$div(class="networkSupportBar",
                             d3Output("up_upPlot",width="auto",height="1200px")
                           ),
                           bsButton(inputId = "upStreamPlotDownload4",label="Download Network Chart",icon=icon("download"),size="default",style="info")
                       )
                       ),
              fluidRow(id="geneExpressionPart5",
                       box(title="EMT-PRO Analytics",width=3,status = "info",solidHeader = TRUE,
                           tags$h4("Select Cell Type"),
                           uiOutput("EPcellTypeSelect"),
                           tags$h4(id="EPGroupTitle","Select group"),
                           uiOutput("EPGroup"),
                           bsButton("epButton","Compute EMT-PRO Score",size="default",style="primary")
                           ),
                       box(title="EMT-PRO Analysis Result",width=9,status = "info",solidHeader = TRUE,
                           plotlyOutput("emt_pro",height="auto",width="auto")
                           )
                       ),
              fluidRow(id="geneExpressionPart6",
                       box(title="GO Analytics",width=3,status="danger",solidHeader = T,
                           tags$h4("Select Cell Type:"),
                           uiOutput("GOCellType"),
                           tags$h4(id="GOTestTitle","Test group"),
                           uiOutput("GOTestGroup"),
                           tags$h4(id="GONormalTitle","Control group"),
                           uiOutput("GONormalGroup"),
                           sliderInput(inputId = "GOUpFc",label="Select Up-regulated Log Fold Change Threshold:",min=0.01,max=2,value=1,step=0.01,width='200px'),
                           sliderInput(inputId = "GODnFc",label="Select Down-regulated Log Fold Change Threshold:",min=-2,max=-0.01,value=-0.2,step=0.01,width='200px'),
                           sliderInput(inputId = "GOPValue",label="Select P value Threshold:",min=0.001,max=1,value=0.05,step=0.001,width='200px'),
                           bsButton("goButton","GO analytics",size="default",style="primary")
                           ),
                       box(title="GO Result",width=9,status="danger",solidHeader = T,
                           selectInput("GOTableSelection",label=h4("Select table to display:"),
                                       choices=c("Up-regulated GO Result"="Up-regulated GO Result",
                                                 "Down-regulated GO Result"="Down-regulated GO Result"),
                                       selected="Up-regulated GO Result"),
                           uiOutput("GOTableDisplay")
                           )
              ),
              fluidRow(id="geneExpressionPart7",
                       box(title="Select GO",width=3,status = "danger",solidHeader = T,
                             uiOutput("GONetworkSelect")
                           ),
                       box(title="GO Analysis Result",width=9,status = "danger",solidHeader = T,
                           d3Output("GOPlot",width="auto",height="500px")
                           )
                       )
              ),
      tabItem(tabName = "network",
              fluidRow(id="communicationNotification",
                  box(title="Notification",width=12,status="danger",solidHeader = TRUE,
                      tags$h1(class="notification","Please upload the read count data in downstream analysis part, or perform cell annotation analysis 
                                   first. Meanwhile, please process the drug file in upload data section first.")
                  )
              ),
              fluidRow(id="communicationPart1",
                box(title="Downstream Analysis Instruction",width=12,status = "primary",solidHeader = TRUE,
                    tags$p("Here you can do inter-cell signaling communication and drug discovering analytics.")
                   )
              ),

              fluidRow(id="communicationPart2",
                box(title="Communication and Drug Discovering",width=4,status = "warning",solidHeader = TRUE,
                    tags$h4("Cell type 1"),
                    uiOutput("cellTypeSelect1"),
                    tags$h4("Cell type 2"),
                    uiOutput("cellTypeSelect2"),
                    tags$h4(id="networkTestTitle","Test group"),
                    uiOutput("controlGroup"),
                    tags$h4(id="networkNormalTitle","Normal group"),
                    uiOutput("normalGroup"),
                    radioButtons("networkResource", "Signaling resources",
                    choices = c(KEGG = "KEGG",
                               STRING="STRING"),
                    selected = "KEGG"),
                    sliderInput(inputId = "networkFc",label="Select Log Fold Change Threshold:",min=0.01,max=2,value=0.2,step=0.01,width='200px'),
                    sliderInput(inputId = "networkPValue",label="Select P value Threshold:",min=0.001,max=1,value=0.05,step=0.001,width='200px'),
                    sliderInput(inputId = "drugNumber",label="Select Number of Top Drugs :",min=1,max=200,value=50,step=1,width='200px'),
                    checkboxInput("networkUseOldResult", "Use saved calculation", T,width="200px"),
                    checkboxInput("networkPadjust", "Use adjusted p-value.", T,width="200px"),
                    bsTooltip( "networkUseOldResult", "If you have computed analysis for this group combination before and want to adjust threshold, you can use old different genes test result to save majority of calculation.
                                  But avoid use it if you change cell type combination.",placement ="right",
                           options=list(container="body")),
                    checkboxInput("useFDAOnly","Only use drug bank drugs in drug discovering.",T,width="200px"),
                    bsButton("communicationButton","Generate Communication Network and Drug",size="default",style="primary")
                ),
                box(title="Gene Expression Information",width=8,status = "warning",solidHeader = TRUE,
                    selectInput("networkTableSelection",label=h4("Select table to display:"),
                                choices=c("Ligands of cell type 1 "="ligType1",
                                          "Receptors of cell type 1"="recType1",
                                          "Ligands of cell type 2"="ligType2",
                                          "Receptors of cell type 2"="recType2"),
                                selected="Ligands of cell type 1"),
                    uiOutput("networkTableDisplay")
              )
              ),
            fluidRow(id="communicationPart13",
                     box(title="Signaling Network Result Display Selection",width=12,status = "danger",solidHeader = TRUE,collapsible = TRUE,
                         selectInput("networkDisplaySelect",label=h4("Select Specific network to display:"),
                                     choices=c("Activated Signaling Pathways of Cell Type 1"="activated1",
                                               "Activated Signaling Pathways of Cell Type 2"="activated2",
                                               "Inter-Cell Signaling Communication from Cell Type 1 to Cell Type 2"="downstream1",
                                               "Inter-Cell Signaling Communication  from Cell Type 2 to Cell Type 1"="downstream2"),
                                     selected="activated1")
                     )
                     ),

              fluidRow(id="communicationPart3",
                       box(title="Activated Signaliing Pathways",width=12,status = "success",solidHeader = TRUE,collapsible = TRUE,
                           uiOutput("activatedNetworkTitle1"),
                           tags$div(class="networkSupportBar",
                                    d3Output("activatedNetworkPlot1",width="auto",height="1200px")

                           ),
                           bsButton(inputId = "activatedNetworkPlotDownload1",label="Download Network Chart",icon=icon("download"),size="default",style="info")
                       )),
              fluidRow(id="communicationPart4",
                       box(title="Activated Signaliing Pathways",width=12,status = "success",solidHeader = TRUE,collapsible = TRUE,
                           uiOutput("activatedNetworkTitle2"),
                           tags$div(class="networkSupportBar",
                                    d3Output("activatedNetworkPlot2",width="auto",height="1200px")
                           ),
                           bsButton(inputId = "activatedNetworkPlotDownload2",label="Download Network Chart",icon=icon("download"),size="default",style="info")
                       )),
              fluidRow(id="communicationPart5",
                box(title="Inter-Cell Signaling Communication",width=12,status = "success",solidHeader = TRUE,collapsible = TRUE,
                             uiOutput("networkTitle1"),
                    tags$div(class="networkSupportBar",
                             d3Output("networkPlot1",width="auto",height="1200px")
                             ),
                    bsButton(inputId = "networkPlotDownload1",label="Download Network Chart",icon=icon("download"),size="default",style="info")
              )),
              fluidRow(id="communicationPart6",
                       box(title="Inter-Cell Signaling Communication",width=12,status = "success",solidHeader = TRUE,collapsible = TRUE,
                           uiOutput("networkTitle2"),
                           tags$div(class="networkSupportBar",
                                    d3Output("networkPlot2",width="auto",height="1200px")
                           ),
                           bsButton(inputId = "networkPlotDownload2",label="Download Network Chart",icon=icon("download"),size="default",style="info")
                           
                       )),    
            fluidRow(id="communicationPart14",
                     box(title="Signaling Drug Result Display Selection",width=12,status = "danger",solidHeader = TRUE,collapsible = TRUE,
                         selectInput("targetDrugDisplaySelect",label=h4("Select Specific network to display:"),
                                     choices=c(
                                       "Drug result of inter-cell signaling communication from cell type 1 to cell type 2"="downstream1",
                                       "Drug result of inter-cell signaling communication from cell type 2 to cell Type 1"="downstream2"),
                                     selected="downstream1")
                     )
            ),
            
              fluidRow(id="communicationPart7",
                           box(title="Drug Discovering Based on Signaling Signatures",width=12,status="info",solidHeader = TRUE,
                               uiOutput("drugTableTitle1"),
                               uiOutput("drugTableDisplay1"))
                       ),
              fluidRow(id="communicationPart8",
                       box(title="Drug Clustering Result for Signaling Signaturee Drug Discovering",width=12,status="info",solidHeader = TRUE,
                           d3Output("drugClusterPlot1",width="auto",height="800px")
                       )
              ),
              fluidRow(id="communicationPart9",
                       box(title="Drug Discovering Based on Signaling Signatures",width=12,status="info",solidHeader = TRUE,
                           uiOutput("drugTableTitle2"),
                           uiOutput("drugTableDisplay2"))
                       ),
             fluidRow(id="communicationPart10",
                      box(title="Drug Clustering Result for Signaling Signature Drug Discovering",width=12,status="info",solidHeader = TRUE,
                          d3Output("drugClusterPlot2",width="auto",height="800px")
                      )
              ),
             fluidRow(id="communicationPart11",
                      box(title="Drug Discovering Based on Targets",width=12,status="success",solidHeader = T,
                          uiOutput("targetDrugNetworkTitle1"),
                          tags$div(class="networkSupportBar",
                          d3Output("drugMappingNetwork1",width="auto",height="1200px")
                          )
                          ),
                      box(title="Drug Discovering Based on Targets",width=12,status="success",solidHeader = T,
                          uiOutput("targetDrugTableTitle1"),
                          uiOutput("drugMappingTable1Display")
                          ),                 
                      box(title="Drug clustering Result for Targets Drug Discovering",width=12,status="success",solidHeader = T,
                                d3Output("targetDrugClusterPlot1",width="auto",height="1200px")

                      )
             ),
             fluidRow(id="communicationPart12",
                      box(title="Drug Discovering Based on Targets",width=12,status="success",solidHeader = T,
                          uiOutput("targetDrugNetworkTitle2"),
                          tags$div(class="networkSupportBar",
                          d3Output("drugMappingNetwork2",width="auto",height="1200px")
                          )            
                          ),
                      box(title="Drug Discovering Based on Targets",width=12,status="success",solidHeader = T,
                          uiOutput("targetDrugTableTitle2"),
                          uiOutput("drugMappingTable2Display")
                      ),                 
                      box(title="Drug clustering Result for Targets Drug Discovering",width=12,status="success",solidHeader = T,
                          d3Output("targetDrugClusterPlot2",width="auto",height="1200px"))
             )
      )

    )
  )
)

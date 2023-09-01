#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# All loading modules used in the app. Mainly contain functions for loading
# different UI given different conditions/input.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Initialization of all variable used during the application session.
#'
#' @return
#' @export
#'
#' @examples
initialization_rv <- function() {
  rv <- reactiveValues(
    keggInfo = NULL,
    original_marker_gene = NULL,
    marker_gene = NULL,
    select_marker_gene = NULL,
    display_marker_gene = NULL,
    bubbleData = NULL,
    df_table = NULL,
    
    genePlotSelectionUI = {
      
    },
    markerGeneTableSelectionUI = {
      
    },
    cellSelectionUI = {
      
    },
    adButtonUI = {
      
    },
    DieaseSpecificTitleUI = {
      
    },
    cellTypeSelect1UI = {
      
    },
    cellTypeSelect2UI = {
      
    },
    networkTableDisplayUI = {
      
    },
    drugTable1UI = {
      
    },
    drugTable2UI = {
      
    },
    networkData = 0,
    drugNetworkData = 0,
    drug_table1 = NULL,
    drug_table2 = NULL,
    rna_df = NULL,
    rna_have_group = 0,
    rna_num_sample = 0,
    rna_group_list = NULL,
    rna_gene_list = NULL,
    rna_type_list = NULL,
    network_df = NULL,
    network_have_group = 0,
    network_group_list = NULL,
    network_gene_list = NULL,
    network_type_list = NULL,
    df_dr = NULL,
    encoder_result = NULL,
    df_cluster = NULL,
    sub_df_cluster = NULL,
    gene_expression_vln_plot = NULL,
    annotation_result = NULL,
    df_classify = NULL,
    rna_data_mean = NULL,
    network_data_mean = NULL,
    useRnaData = 1,
    cell_type1 = NULL,
    cell_type2 = NULL,
    downstream_network_nodes1 = NULL,
    downstream_network_edges1 = NULL,
    downstream_network_json1 = NULL,
    networkData1 = 0,
    haveGenerated1 = 0,
    networkTitle1UI = {
      
    },
    downstream_network_nodes2 = NULL,
    downstream_network_edges2 = NULL,
    downstream_network_json2 = NULL,
    networkData2 = 0,
    haveGenerated2 = 0,
    networkTitle2UI = {
      
    },
    cell_type1_ligands = NULL,
    cell_type1_receptors = NULL,
    cell_type2_ligands = NULL,
    cell_type2_receptors = NULL,
    drug_network_data = NULL,
    drug_table1 = NULL,
    drug_table2 = NULL,
    upRegLigand_network = NULL,
    upRegReceptor_network = NULL,
    upRegGenerated = 0,
    controlGroupUI = {
      
    },
    normalGroupUI = {
      
    },
    networkDownloadUI = {
      
    },
    cell_count = NULL,
    distribution_have_group = 0,
    EPcellTypeSelectUI = {
      
    },
    EPGroupUI = {
      
    },
    outputDir = NULL,
    gSymZs = NULL,
    rankMatrix = NULL,
    drug_mapping = NULL,
    drugData1 = 0,
    drugData2 = 0,
    drugTableTitle1UI = {
      
    },
    drugTableTitle2UI = {
      
    },
    drugGenerated = 0,
    drugClusterData1 = 0,
    drugClusterData2 = 0,
    drug_json1 = NULL,
    drug_json2 = NULL,
    networkCellSelectTitleUI = {
      
    },
    GONetwork1 = NULL,
    GOGenerated = 0,
    GOData = 0,
    GOTestGroupUI = {
      
    },
    GONormalGroupUI = {
      
    },
    GOCellTypeUI = {
      
    },
    GO_up_table = NULL,
    GO_dn_table = NULL,
    GO_table = NULL,
    GOTableDisplayUI = {
      
    },
    GONetworkSelectUI = {
      
    },
    GOList = NULL,
    netGO_up = NULL,
    netGO_dn = NULL,
    upRegTestGroupUI = {
      
    },
    upRegNormalGroupUI = {
      
    },
    upStreamGenerated = 0,
    up_exp_network = NULL,
    exp_up_network = NULL,
    combine_network = NULL,
    up_up_network = NULL,
    DLRP = NULL,
    nicheNet = NULL,
    baderLab = NULL,
    ligRecDatabase = NULL,
    drugNetworkGenerated = 0,
    drugNetworkJson1 = NULL,
    drugNetworkJson2 = NULL,
    drug_mapping_table1 = NULL,
    drug_mapping_table2 = NULL,
    drugNetworkJson1 = NULL,
    drugNetworkJson2 = NULL,
    drugBankInteraction = NULL,
    drugBankInformation = NULL,
    targetDrug_json1 = NULL,
    targetDrug_json2 = NULL,
    activatedNetworkTitle1UI = {
      
    },
    activatedNetworkTitle2UI = {
      
    },
    activatedNetworkGenerated1 = 0,
    activatedNetworkGenerated2 = 0,
    targetDrugTableTitle1UI = {
      
    },
    targetDrugTableTitle2UI = {
      
    },
    targetDrugNetworkTitle1UI = {
      
    },
    targetDrugNetworkTitle2UI = {
      
    }
  )
  
}

#' Initialization of the UI interface.
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
initialization_UI <- function(rv) {
  shinyjs:::hide("dimensionReductionPart1")
  shinyjs:::hide("dimensionReductionPart2")
  shinyjs:::show("drNotification")
  shinyjs:::show("clusterNotification")
  shinyjs:::hide("clusterPart1")
  shinyjs:::hide("clusterPart2")
  shinyjs:::show("genePlotNotification")
  shinyjs:::hide("geneFeaturePlotPart1")
  shinyjs:::hide("geneFeaturePlotPart2")
  shinyjs:::hide("clusterPart4")
  shinyjs:::show("classificationNotification")
  shinyjs:::hide("classificationPart1")
  shinyjs:::hide("classificationPart2")
  #shinyjs:::hide("classificationPart3")
  shinyjs::show("geneExpressionNotification")
  shinyjs::hide("geneExpressionPart1")
  shinyjs::hide("geneExpressionPart2")
  shinyjs::hide("geneExpressionPart3")
  shinyjs::hide("geneExpressionPart4")
  shinyjs::hide("geneExpressionPart5")
  shinyjs::hide("geneExpressionPart6")
  shinyjs::hide("geneExpressionPart7")
  hideEPSelection(rv)
  shinyjs:::show("communicationNotification")
  shinyjs:::hide("communicationPart1")
  shinyjs:::hide("communicationPart2")
  shinyjs:::hide("communicationPart3")
  shinyjs:::hide("communicationPart4")
  shinyjs:::hide("communicationPart5")
  shinyjs:::hide("communicationPart6")
  shinyjs:::hide("communicationPart7")
  shinyjs:::hide("communicationPart8")
  shinyjs:::hide("communicationPart9")
  shinyjs:::hide("communicationPart10")
  shinyjs:::hide("communicationPart11")
  shinyjs:::hide("communicationPart12")
  shinyjs:::hide("communicationPart13")
  shinyjs:::hide("communicationPart14")
  
  shinyjs:::show("drugNotification")
  shinyjs:::hide("drugPart1")
  shinyjs:::hide("drugPart2")
  shinyjs:::hide("drugPart3")
  if (!dir.exists("./cache")) {
    dir.create("./cache")
  }
}


#' load rank matrix data for drug discovery. The data will be auto saved for
#' loading after drug file processing.
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadRankMatrix <- function(rv) {
  if (file.exists("./data/rankMatrix92742a.RData") &&
      file.exists("./data/drug_mapping.RData")) {
    load("./data/rankMatrix92742a.RData")
    rv$gSymZs <- gSymZs
    rv$rankMatrix <- rankMatrix
    load("./data/drug_mapping.RData")
    rv$drug_mapping <- drug_mapping
  }
}

#' Load working directory used last time.
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadlastWorkingDir <- function(rv) {
  if (file.exists("./cache/lastWorkingDir.RData")) {
    load("./cache/lastWorkingDir.RData")
    rv$outputDir <- newDir
  } else{
    rv$outputDir <- NULL
  }
}


#' Load Full ligand-receptor database (including Bader Lab, DLRP, NicheNet).
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadLigRecDatabase <- function(rv) {
  progress <- shiny::Progress$new()
  progress$set(message = "Change database", value = 0)
  on.exit(progress$close())
  progress$set(0.5, detail = "Processing")
  load("./data/baderLab_lr_database.RData")
  load("./data/DLRP_lr_database.RData")
  load("./data/nicheNet_lr_database.RData")
  rv$DLRP <- DLRP
  rv$nicheNet <- nicheNet
  rv$baderLab <- baderLab
  ligRecDatabase <- rbind(DLRP, nicheNet, baderLab)
  ligRecDatabase <- unique(ligRecDatabase)
  ligRecDatabase <- as.matrix(ligRecDatabase)
  rv$ligRecDatabase <- ligRecDatabase
  progress$set(1, detail = "Finish")
  rm(DLRP, baderLab, nicheNet, ligRecDatabase)
  gc()
}


#' Load KEGG signaling interaction database.
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadKeggInfo <- function(rv) {
  load("./data/keggInfo.RData")
  rv$keggInfo <- keggInfo
  rm(keggInfo)
  gc()
}


#' Load transcriptional factor-target interaction database.
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadTfTargetInteraction <- function(rv) {
  load('./data/TfTargetInteraction_human.RData')
  rv$TfTargetInteraction <- TfTargetInteraction
  rm(TfTargetInteraction)
  gc()
}


#' Load STRING signaling interaction database.
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadSTRING <- function(rv) {
  load("./data/STRING.RData")
  rv$STRING = STRING
  rm(STRING)
  gc()
}

#' Load marker gene database.
#'
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadMarkerGene <- function(rv, session) {
  if ((!is.null(rv$outputDir)) &&
      file.exists(paste(rv$outputDir, "marker_gene_table.RData", sep = "/"))) {
    load(paste(rv$outputDir, "marker_gene_table.RData", sep = "/"))
  } else{
    load("./data/marker_gene_table.RData")
  }
  
  #generate bubble plot data
  bubbleData <- c()
  unique_cell_type <- unique(marker_gene_table$cell_type)
  for (i in 1:length(unique_cell_type)) {
    cell_type = unique_cell_type[i]
    gene_list <-
      marker_gene_table[marker_gene_table[, 1] == cell_type, 2]
    cell_marker_gene <- paste(gene_list, collapse = ",")
    name_length <- nchar(cell_type) + sample(4:10, 1)
    row_content <- c(cell_type, name_length, cell_marker_gene)
    if (i == 1) {
      bubbleData <- row_content
    }
    else{
      bubbleData <- rbind(bubbleData, row_content)
    }
  }
  
  colnames(bubbleData) <- c("cellType", "length", "markerGene")
  
  rv$bubbleData <- as.data.frame(bubbleData)
  rownames(marker_gene_table) <- c(1:nrow(marker_gene_table))
  rv$original_marker_gene <- marker_gene_table
  rv$select_marker_gene <- marker_gene_table
  rv$display_marker_gene <- marker_gene_table
  
  session$sendCustomMessage("refreshBubble", "refresh")
  rm(
    bubbleData,
    marker_gene_table,
    gene_list,
    cell_marker_gene,
    name_length,
    row_content
  )
  gc()
}



#' Load drug bank database.
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadDrugBankData <- function(rv) {
  load("data/DrugBank_drug-target_new.RData")
  load("data/DrugBank_drugInformation.RData")
  rv$drugBankInteraction <- drugBank_drug_target
  rv$drugBankInformation <- drugBank_drugInformation
  rm(drugBank_drug_target, drugBank_drugInformation)
  gc()
}


#' Load gene set selection input in the Gene feature plot section
#' (need to dynamically load selection list from app variable).
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadGenePlotSelection <- function(rv) {
  rv$genePlotSelectionUI <- {
    selectInput(
      "genePlotSelection",
      "Select gene to explore:",
      rv$rna_gene_list,
      multiple = T,
      selectize = T
    )
  }
  
}

#' Load gene set selection input in the marker gene section
#' (need to dynamically load selection list from app variable).
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadMarkerGeneTableSelection <- function(rv) {
  rv$markerGeneTableSelectionUI <- {
    selectInput(
      "markerGeneTableSelection",
      "Select gene to add as marker gene:",
      rv$rna_gene_list,
      multiple = T,
      selectize = T,
      width = '200px'
    )
  }
  
}

#' Load cell type set selection input in the cell annotation section
#' (need to dynamically load selection list from app variable).
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadCellSelection <- function(rv) {
  cell_list <- unique(rv$display_marker_gene$cell_type)
  rv$cellSelectionUI <- {
    tags$div(id = "supportBar",
             tags$div(
               id = "cellSelectionFlow",
               checkboxGroupInput(
                 inputId = "cellSelection",
                 label = NULL,
                 choices = cell_list,
                 choiceValues = cell_list,
                 width = "100%"
               )
             ))
  }
  rm(cell_list)
  gc()
}


#' Load raw seq-RNA data uploaded by user.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadRnaFile <- function(input, rv, session) {
  progress <- shiny::Progress$new()
  progress$set(message = "Load file", value = 0)
  progress$set(0.1, detail = "read data")
  on.exit(progress$close())
  tryCatch({
    suffix <-
      strsplit(input$rnaDataFile$datapath, split = "\\.")[[1]][2]
    if (tolower(suffix) == "csv") {
      data <- read.csv(
        input$rnaDataFile$datapath,
        header = input$rnaHeader,
        sep = input$rnaSep,
        quote = input$rnaQuote,
        row.names = 1,
        check.names = FALSE,
        as.is = FALSE
      )
      data <- Matrix(data, sparse = TRUE)
    } else if (tolower(suffix) == "rds") {
      data <- readRDS(input$rnaDataFile$datapath)
      if (class(data)[1] == "Seurat") {
        # Due to compatibility issue,
        # we will recreate a new Seurat object using raw count data.
        data <- data[["RNA"]]@counts
      } else if (class(data)[1] != "dgCMatrix") {
        stop(safeError("invalid data"))
      }
    } else{
      stop(safeError("invalid data"))
    }
    rv$rna_num_sample <- ncol(data)
    rv$rna_gene_list <- as.character(rownames(data))
    rv$rna_df <- data
    progress$set(1, detail = "finish")
    rm(data)
    gc()
  },
  error = function(e) {
    #    return a safeError if a parsing error occurs
    if (e$message == "invalid data") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Upload Error",
        content = "Only accept .csv or .rds file of Seurat object or dgCMatrix.",
        style = "danger"
      )
    } else{
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Upload Error",
        content = paste0(
          "The file you uploaded is invalid, Please read upload instruction carefully.",
          e
        ),
        style = "danger"
      )
    }
    
    
  },
  finally = {
    progress$inc(1, detail = "transposition finish")
    on.exit(progress$close())
  })
}



#' Load group/design information corresponding to the uploaded raw seq-RNA data.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadRnaGroupFile <- function(input, rv, session) {
  tryCatch({
    if (is.null(rv$rna_df)) {
      stop(safeError("no RNA data"))
    }
    suffix <-
      strsplit(input$rnaGroupFile$datapath, split = "\\.")[[1]][2]
    if (tolower(suffix) == "csv") {
      df <- read.csv(
        input$rnaGroupFile$datapath,
        header = input$rnaHeader,
        sep = input$rnaSep,
        quote = input$rnaQuote,
        check.names = FALSE,
        as.is = FALSE
      )
    } else if (tolower(suffix) == "rds") {
      df <- readRDS(input$rnaGroupFile$datapath)
    } else{
      stop(safeError("invalid data"))
    }
    group_list <- as.character(df[, 1])
    unique_group <- unique(group_list)
    if (length(unique_group) <= 1 ||
        length(group_list) != rv$rna_num_sample) {
      stop(safeError("incorrect data"))
    }
    rv$rna_group_list <- group_list
    save(group_list,
         file = paste0(rv$outputDir, "/rnaGroupInformation.RData"))
    rv$rna_have_group = 1
    rm(group_list, unique_group, df, suffix)
    gc()
  },
  error = function(e) {
    #    return a safeError if a parsing error occurs
    if (e$message == "no RNA data") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Upload Error",
        content = "Please unpload seq-RNA data first!.",
        style = "danger"
      )
    } else if (e$message == "invalid data") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Upload Error",
        content = "Only accept .csv or .rds file of data.frame.",
        style = "danger"
      )
    } else if (e$message == "incorrect data") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Upload Error",
        content = "the length of group file not match count data or the number of unique group is only 1.",
        style = "danger"
      )
    } else{
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Upload Error",
        content = paste0(
          "The file you uploaded is invalid, Please read upload instruction carefully.",
          e
        ),
        style = "danger"
      )
    }
  })
}


#' Load scRNA-seq data for downstream analysis uploaded by user.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadCommunicationFile <- function(input, rv, session) {
  progress <- shiny::Progress$new()
  progress$set(message = "Data processing", value = 0)
  on.exit(progress$close())
  progress$set(0.1, detail = "read data")
  tryCatch({
    suffix <-
      strsplit(input$networkDataFile$datapath, split = "\\.")[[1]][2]
    if (tolower(suffix) == "csv") {
      network_df <- read.csv(
        input$networkDataFile$datapath,
        header = input$networkHeader,
        sep = input$networkSep,
        quote = input$networkQuote,
        row.names = 1,
        check.names = FALSE,
        as.is = FALSE
      )
      network_df <- Matrix(network_df, sparse = TRUE)
    } else if (tolower(suffix) == "rds") {
      network_df <- readRDS(input$networkDataFile$datapath)
      if (class(network_df)[1] == "Seurat" &
          class(network_df)[2] == "SeuratObject") {
        network_df <- network_df[["RNA"]]@counts
      } else if (class(network_df)[1] != "dgCMatrix") {
        stop(safeError("invalid data"))
      }
    } else{
      stop(safeError("invalid data"))
    }
    rv$network_gene_list <- as.character(rownames(network_df))
    rv$network_num_sample <- ncol(network_df)
    rv$network_df <- network_df
    rm(suffix, network_df)
    progress$inc(1, detail = "Finish")
  },
  error = function(e) {
    #    return a safeError if a parsing error occurs
    if (e$message == "invalid data") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Upload Error",
        content = "Only accept .csv or .rds file.",
        style = "danger"
      )
    } else{
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Upload Error",
        content = "The file you uploaded is invalid, Please read upload instruction carefully.",
        style = "danger"
      )
    }
    
  },
  finally = {
    progress$inc(1, detail = "Finish")
    on.exit(progress$close())
  })
}


#' Load group/design information corresponding to the uploaded seq-RNA data for downstream analysis.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadNetworkGroupFile <- function(input, rv, session) {
  tryCatch({
    if (is.null(rv$network_df)) {
      stop(safeError("no network data"))
    }
    suffix <-
      strsplit(input$networkGroupFile$datapath, split = "\\.")[[1]][2]
    if (tolower(suffix) == "csv") {
      df <- read.csv(
        input$networkGroupFile$datapath,
        header = input$networkHeader,
        sep = input$networkSep,
        quote = input$networkQuote,
        check.names = FALSE,
        as.is = FALSE
      )
    } else if (tolower(suffix) == "rds") {
      df <- readRDS(input$networkGroupFile$datapath)
    } else{
      stop(safeError("invalid data"))
    }
    
    num_row <- nrow(df)
    if (num_row != rv$network_num_sample) {
      stop(safeError("Information not correct."))
    }
    num_col <- ncol(df)
    if (num_col == 1) {
      type_list <- as.character(df[, 1])
      rv$network_type_list <- type_list
      save(type_list,
           file = paste0(rv$outputDir, "/networkTypeInformation.RData"))
      rv$network_have_group = 0
      rm(type_list)
    } else if (num_col == 2) {
      type_list <- as.character(df[, 1])
      rv$network_type_list <- type_list
      save(type_list,
           file = paste0(rv$outputDir, "/networkTypeInformation.RData"))
      group_list <- as.character(df[, 2])
      unique_group <- unique(group_list)
      if (length(unique_group) <= 1 ||
          length(group_list) != rv$network_num_sample) {
        stop(safeError("Group information not correct."))
      }
      rv$network_group_list <- group_list
      save(group_list,
           file = paste0(rv$outputDir, "/networkGroupInformation.RData"))
      rv$network_have_group = 1
      rm(type_list, group_list, unique_group)
    } else{
      stop(safeError("Group information not correct."))
    }
    rm(suffix, df)
    gc()
  },
  error = function(e) {
    #    return a safeError if a parsing error occurs
    if (e == "no network data") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Upload Error",
        content = "Please unpload gene communication data first!.",
        style = "danger"
      )
    } else if (e$message == "invalid data") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Upload Error",
        content = "Only accept .csv or .rds file.",
        style = "danger"
      )
    } else{
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Upload Error",
        content = "The file you uploaded is invalid, Please read upload instruction carefully.",
        style = "danger"
      )
    }
  })
}



#' Load up-regulated Ligands and Receptors group selection
#' (need to dynamically load selection list from app variable).
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadUpReg <- function(rv) {
  if (rv$useRnaData == 0) {
    if (rv$network_have_group == 1) {
      group_list <- unique(rv$network_group_list)
      shinyjs::show("upRegNormalTitle")
      shinyjs::show("upRegTestTitle")
      #test group
      rv$upRegTestGroupUI <- {
        tags$div(class = "GroupSupportBar",
                 tags$div(
                   class = "groupSelectionFlow",
                   checkboxGroupInput(
                     inputId = "upRegTestGroup",
                     label = NULL,
                     choices = group_list,
                     choiceValues = group_list,
                     width = "100%"
                   )
                 ))
      }
      #normal group
      rv$upRegNormalGroupUI <- {
        tags$div(class = "GroupSupportBar",
                 tags$div(
                   class = "groupSelectionFlow",
                   checkboxGroupInput(
                     inputId = "upRegNormalGroup",
                     label = NULL,
                     choices = group_list,
                     choiceValues = group_list,
                     width = "100%"
                   )
                 ))
      }
      rm(group_list)
    } else{
      shinyjs::hide("upRegNormalTitle")
      shinyjs::hide("upRegTestTitle")
      #test group
      rv$upRegTestGroupUI <- {
        
      }
      #normal group
      rv$upRegNormalGroupUI <- {
        
      }
    }
  } else{
    if (rv$rna_have_group == 1) {
      shinyjs::show("upRegNormalTitle")
      shinyjs::show("upRegTestTitle")
      group_list <- unique(rv$rna_group_list)
      #test group
      rv$upRegTestGroupUI <- {
        tags$div(class = "GroupSupportBar",
                 tags$div(
                   class = "groupSelectionFlow",
                   checkboxGroupInput(
                     inputId = "upRegTestGroup",
                     label = NULL,
                     choices = group_list,
                     choiceValues = group_list,
                     width = "100%"
                   )
                 ))
      }
      #normal group
      rv$upRegNormalGroupUI <- {
        tags$div(class = "GroupSupportBar",
                 tags$div(
                   class = "groupSelectionFlow",
                   checkboxGroupInput(
                     inputId = "upRegNormalGroup",
                     label = NULL,
                     choices = group_list,
                     choiceValues = group_list,
                     width = "100%"
                   )
                 ))
      }
      rm(group_list)
    } else{
      shinyjs::hide("upRegNormalTitle")
      shinyjs::hide("upRegTestTitle")
      #test group
      rv$upRegTestGroupUI <- {
        
      }
      #normal group
      rv$upRegNormalGroupUI <- {
        
      }
    }
  }
  gc()
}


#' Load cell-cell communication selection UI from upstream analysis data.
#' (need to dynamically load selection list from app variable).
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadRnaDataSelection <- function(rv) {
  type_list <- unique(rv$rna_type_list)
  #cell type 1
  rv$cellTypeSelect1UI <- {
    tags$div(class = "GroupSupportBar",
             tags$div(
               class = "groupSelectionFlow",
               checkboxGroupInput(
                 inputId = "cellTypeSelect1",
                 label = NULL,
                 choices = type_list,
                 choiceValues = type_list,
                 width = "100%"
               )
             ))
  }
  #cell type 2
  rv$cellTypeSelect2UI <- {
    # tags$h4("Stroma Cell")
    tags$div(class = "GroupSupportBar",
             tags$div(
               class = "groupSelectionFlow",
               checkboxGroupInput(
                 inputId = "cellTypeSelect2",
                 label = NULL,
                 choices = type_list,
                 choiceValues = type_list,
                 width = "100%"
               )
             ))
  }
  if (rv$rna_have_group == 1) {
    shinyjs::show("networkNormalTitle")
    shinyjs::show("networkTestTitle")
    group_list <- unique(rv$rna_group_list)
    #control group
    rv$controlGroupUI <- {
      tags$div(class = "GroupSupportBar",
               tags$div(
                 class = "groupSelectionFlow",
                 checkboxGroupInput(
                   inputId = "controlGroup",
                   label = NULL,
                   choices = group_list,
                   choiceValues = group_list,
                   width = "100%"
                 )
               ))
    }
    #normal group
    rv$normalGroupUI <- {
      tags$div(class = "GroupSupportBar",
               tags$div(
                 class = "groupSelectionFlow",
                 checkboxGroupInput(
                   inputId = "normalGroup",
                   label = NULL,
                   choices = group_list,
                   choiceValues = group_list,
                   width = "100%"
                 )
               ))
    }
  } else{
    shinyjs::hide("networkNormalTitle")
    shinyjs::hide("networkTestTitle")
    #control group
    rv$controlGroupUI <- {
      
    }
    #normal group
    rv$normalGroupUI <- {
      
    }
  }
}

#' Load cell-cell communication selection UI from downstream anaylsis data.
#' (need to dynamically load selection list from app variable).
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadNetworkDataSelection <- function(rv) {
  type_list <- unique(rv$network_type_list)
  #cell type 1
  rv$cellTypeSelect1UI <- {
    tags$div(class = "GroupSupportBar",
             tags$div(
               class = "groupSelectionFlow",
               checkboxGroupInput(
                 inputId = "cellTypeSelect1",
                 label = NULL,
                 choices = type_list,
                 choiceValues = type_list,
                 width = "100%"
               )
             ))
  }
  #cell type 2
  rv$cellTypeSelect2UI <- {
    tags$div(class = "GroupSupportBar",
             tags$div(
               class = "groupSelectionFlow",
               checkboxGroupInput(
                 inputId = "cellTypeSelect2",
                 label = NULL,
                 choices = type_list,
                 choiceValues = type_list,
                 width = "100%"
               )
             ))
  }
  if (rv$network_have_group == 1) {
    group_list <- unique(rv$network_group_list)
    shinyjs::show("networkNormalTitle")
    shinyjs::show("networkTestTitle")
    #control group
    rv$controlGroupUI <- {
      tags$div(class = "GroupSupportBar",
               tags$div(
                 class = "groupSelectionFlow",
                 checkboxGroupInput(
                   inputId = "controlGroup",
                   label = NULL,
                   choices = group_list,
                   choiceValues = group_list,
                   width = "100%"
                 )
               ))
    }
    #normal group
    rv$normalGroupUI <- {
      tags$div(class = "GroupSupportBar",
               tags$div(
                 class = "groupSelectionFlow",
                 checkboxGroupInput(
                   inputId = "normalGroup",
                   label = NULL,
                   choices = group_list,
                   choiceValues = group_list,
                   width = "100%"
                 )
               ))
    }
  } else{
    shinyjs::hide("networkNormalTitle")
    shinyjs::hide("networkTestTitle")
    #control group
    rv$controlGroupUI <- {
      
    }
    #normal group
    rv$normalGroupUI <- {
      
    }
  }
}


#' Update the ligand-receptor database given the selection result from the user.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
ligRecDatabaseSelect <- function(input, rv) {
  ligRecDatabase <- c()
  if (length(input$ligRecData) == 0 ||
      length(input$ligRecData) == 3) {
    database <- rbind(rv$DLRP, rv$nicheNet, rv$baderLab)
    ligRecDatabase <- database
  } else{
    for (i in 1:length(input$ligRecData)) {
      database <- input$ligRecData[i]
      if (database == "DLRP") {
        if (i == 1) {
          ligRecDatabase <- rv$DLRP
        } else{
          ligRecDatabase <- rbind(ligRecDatabase, rv$DLRP)
        }
      } else if (database == "nicheNet") {
        if (i == 1) {
          ligRecDatabase <- rv$nicheNet
        } else{
          ligRecDatabase <- rbind(ligRecDatabase, rv$nicheNet)
        }
      } else{
        if (i == 1) {
          ligRecDatabase <- rv$baderLab
        } else{
          ligRecDatabase <- rbind(ligRecDatabase, rv$baderLab)
        }
      }
    }
  }
  ligRecDatabase <- unique(ligRecDatabase)
  ligRecDatabase <- as.matrix(ligRecDatabase)
  rv$ligRecDatabase <- ligRecDatabase
  rm(ligRecDatabase, database)
  gc()
}

#' Load gene ontology selection UI.
#' (need to dynamically load selection list from app variable).
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadGO <- function(rv) {
  if (rv$useRnaData == 0) {
    type_list <- unique(rv$network_type_list)
    #cell type
    rv$GOCellTypeUI <- {
      tags$div(class = "GroupSupportBar",
               tags$div(
                 class = "groupSelectionFlow",
                 checkboxGroupInput(
                   inputId = "GOCellType",
                   label = NULL,
                   choices = type_list,
                   choiceValues = type_list,
                   width = "100%"
                 )
               ))
    }
    if (rv$network_have_group == 1) {
      shinyjs::show("GONormalTitle")
      shinyjs::show("GOTestTitle")
      group_list <- unique(rv$network_group_list)
      #control group
      rv$GOTestGroupUI <- {
        tags$div(class = "GroupSupportBar",
                 tags$div(
                   class = "groupSelectionFlow",
                   checkboxGroupInput(
                     inputId = "GOTestGroup",
                     label = NULL,
                     choices = group_list,
                     choiceValues = group_list,
                     width = "100%"
                   )
                 ))
      }
      #normal group
      rv$GONormalGroupUI <- {
        tags$div(class = "GroupSupportBar",
                 tags$div(
                   class = "groupSelectionFlow",
                   checkboxGroupInput(
                     inputId = "GONormalGroup",
                     label = NULL,
                     choices = group_list,
                     choiceValues = group_list,
                     width = "100%"
                   )
                 ))
      }
    } else{
      shinyjs::hide("GONormalTitle")
      shinyjs::hide("GOTestTitle")
      #control group
      rv$GOTestGroupUI <- {
        
      }
      #normal group
      rv$GONormalGroupUI <- {
        
      }
    }
  } else{
    type_list <- unique(rv$rna_type_list)
    #cell type
    rv$GOCellTypeUI <- {
      tags$div(class = "GroupSupportBar",
               tags$div(
                 class = "groupSelectionFlow",
                 checkboxGroupInput(
                   inputId = "GOCellType",
                   label = NULL,
                   choices = type_list,
                   choiceValues = type_list,
                   width = "100%"
                 )
               ))
    }
    if (rv$rna_have_group == 1) {
      shinyjs::show("GONormalTitle")
      shinyjs::show("GOTestTitle")
      group_list <- unique(rv$rna_group_list)
      #control group
      rv$GOTestGroupUI <- {
        tags$div(class = "GroupSupportBar",
                 tags$div(
                   class = "groupSelectionFlow",
                   checkboxGroupInput(
                     inputId = "GOTestGroup",
                     label = NULL,
                     choices = group_list,
                     choiceValues = group_list,
                     width = "100%"
                   )
                 ))
      }
      #normal group
      rv$GONormalGroupUI <- {
        tags$div(class = "GroupSupportBar",
                 tags$div(
                   class = "groupSelectionFlow",
                   checkboxGroupInput(
                     inputId = "GONormalGroup",
                     label = NULL,
                     choices = group_list,
                     choiceValues = group_list,
                     width = "100%"
                   )
                 ))
      }
    } else{
      shinyjs::hide("GONormalTitle")
      shinyjs::hide("GOTestTitle")
      #control group
      rv$GOTestGroupUI <- {
        
      }
      #normal group
      rv$GONormalGroupUI <- {
        
      }
    }
  }
}

#' Load gene ontology network selection UI.
#' (need to dynamically load selection list from app variable).
#'
#' @param rv Variables list in the Shiny app.
#' @param GO_list All GO terms computed from the GO analysis.
#'
#' @return
#' @export
#'
#' @examples
loadGONetworkSelect <- function(rv, GO_list) {
  rv$GONetworkSelectUI <- {
    selectInput(
      "GOList",
      label = h4("Select GO"),
      choices = GO_list,
      selected = GO_list[1]
    )
  }
}

#' Load dimension reduction section.
#'
#' @return
#' @export
#'
#' @examples
loadDimensionReduction <- function() {
  shinyjs:::show("dimensionReductionPart1")
  shinyjs:::show("dimensionReductionPart2")
  shinyjs:::hide("drNotification")
}


#' Load clustering section.
#'
#' @return
#' @export
#'
#' @examples
loadCluster <- function() {
  shinyjs:::hide("clusterNotification")
  shinyjs:::show("clusterPart1")
  shinyjs:::show("clusterPart2")
}

#' Load gene feature exploration section.
#'
#' @return
#' @export
#'
#' @examples
loadGeneExploration <- function() {
  shinyjs::hide("genePlotNotification")
  shinyjs::show("geneFeaturePlotPart1")
  shinyjs::show("geneFeaturePlotPart2")
}


#' Load cell annotation section.
#'
#' @return
#' @export
#'
#' @examples
loadClassification <- function() {
  shinyjs:::hide("classificationNotification")
  shinyjs:::show("classificationPart1")
  shinyjs:::show("classificationPart2")
  #shinyjs:::show("classificationPart3")
}


#' Load gene expression part 1.
#'
#' @return
#' @export
#'
#' @examples
loadGeneExpression1 <- function() {
  shinyjs::hide("geneExpressionNotification")
  shinyjs::show("geneExpressionPart1")
}

#' Load gene expression part 2.
#'
#' @return
#' @export
#'
#' @examples
loadGeneExpression2 <- function() {
  shinyjs::show("geneExpressionPart2")
  shinyjs::show("geneExpressionPart3")
  shinyjs::show("geneExpressionPart4")
  shinyjs::show("geneExpressionPart5")
  shinyjs::show("geneExpressionPart6")
  shinyjs::show("geneExpressionPart7")
}

#' Hide EMT-PRO UI
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
hideEPSelection <- function(rv) {
  rv$EPcellTypeSelectUI <- {
    
  }
  rv$EPGroupUI <- {
    
  }
  shinyjs::hide("epButton")
}

#' Load EMT-PRO UI
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadEPSelection <- function(rv) {
  if (rv$useRnaData == 1) {
    type_list <- unique(rv$rna_type_list)
    rv$EPcellTypeSelectUI <- {
      tags$div(class = "GroupSupportBar",
               tags$div(
                 class = "groupSelectionFlow",
                 checkboxGroupInput(
                   inputId = "epType",
                   label = NULL,
                   choices = type_list,
                   choiceValues = type_list,
                   width = "100%"
                 )
               ))
    }
    
    if (rv$rna_have_group == 1) {
      shinyjs::show("EPGroupTitle")
      group_list <- unique(rv$rna_group_list)
      rv$EPGroupUI <- {
        tags$div(class = "GroupSupportBar",
                 tags$div(
                   class = "groupSelectionFlow",
                   checkboxGroupInput(
                     "epGroup",
                     label = NULL,
                     choices = group_list,
                     choiceValues = group_list,
                     width = "100%"
                   )
                 ))
      }
    } else{
      shinyjs::hide("EPGroupTitle")
    }
  } else{
    type_list <- unique(rv$network_type_list)
    rv$EPcellTypeSelectUI <- {
      tags$div(class = "GroupSupportBar",
               tags$div(
                 class = "groupSelectionFlow",
                 checkboxGroupInput(
                   inputId = "epType",
                   label = NULL,
                   choices = type_list,
                   choiceValues = type_list,
                   width = "100%"
                 )
               ))
    }
    if (rv$network_have_group == 1) {
      shinyjs::show("EPGroupTitle")
      group_list <- unique(rv$network_group_list)
      rv$EPGroupUI <- {
        tags$div(class = "GroupSupportBar",
                 tags$div(
                   class = "groupSelectionFlow",
                   checkboxGroupInput(
                     "epGroup",
                     label = NULL,
                     choices = group_list,
                     choiceValues = group_list,
                     width = "100%"
                   )
                 ))
      }
    } else{
      shinyjs::hide("EPGroupTitle")
    }
  }
  shinyjs::show("epButton")
}

#' Load inter-cell communication section.
#'
#' @return
#' @export
#'
#' @examples
loadCommunication <- function() {
  shinyjs:::hide("communicationNotification")
  shinyjs:::show("communicationPart1")
  shinyjs:::show("communicationPart2")
  shinyjs:::show("communicationPart3")
  shinyjs:::hide("communicationPart4")
  shinyjs:::hide("communicationPart5")
  shinyjs:::hide("communicationPart6")
  shinyjs:::show("communicationPart7")
  shinyjs:::show("communicationPart8")
  shinyjs:::hide("communicationPart9")
  shinyjs:::hide("communicationPart10")
  shinyjs:::show("communicationPart11")
  shinyjs:::hide("communicationPart12")
  shinyjs:::show("communicationPart13")
  shinyjs:::show("communicationPart14")
}

#' Load inter-cell signaling communication network 1.
#'
#' @param rv Variables list in the Shiny app.
#' @param cell_type1 The source cell type.
#' @param cell_type2 The target cell type.
#'
#' @return
#' @export
#'
#' @examples
loadNetwork1 <- function(rv, cell_type1, cell_type2) {
  rv$networkTitle1UI = {
    tags$h3(paste(
      "Inter-Cell Signaling Communication from",
      paste(cell_type1, collapse = "+"),
      "to",
      paste(cell_type2, collapse = "+"),
      sep = " "
    ))
  }
}

#' Load inter-cell signaling communication network 1.
#'
#' @param rv Variables list in the Shiny app.
#' @param cell_type1 The target cell type.
#' @param cell_type2 The source cell type.
#'
#' @return
#' @export
#'
#' @examples
loadNetwork2 <- function(rv, cell_type1, cell_type2) {
  rv$networkTitle2UI = {
    tags$h3(paste(
      "Inter-Cell Signaling Communication from",
      paste(cell_type2, collapse = "+"),
      "to",
      paste(cell_type1, collapse = "+"),
      sep = " "
    ))
  }
}


#' Load activated signaling pathway network 1.
#'
#' @param rv Variables list in the Shiny app.
#' @param cell_type1 The selected cell type.
#'
#' @return
#' @export
#'
#' @examples
loadActivatedNetwork1 <- function(rv, cell_type1) {
  rv$activatedNetworkTitle1UI = {
    tags$h3(paste0(
      "Activated Signaling Pathways for ",
      paste(cell_type1, collapse = "+")
    ))
  }
}

#' Load activated signaling pathway network 2.
#'
#' @param rv Variables list in the Shiny app.
#' @param cell_type2 The selected cell type.
#'
#' @return
#' @export
#'
#' @examples
loadActivatedNetwork2 <- function(rv, cell_type2) {
  rv$activatedNetworkTitle2UI = {
    tags$h3(paste0(
      "Activated Signaling Pathways for ",
      paste(cell_type2, collapse = "+")
    ))
  }
}

#' Load drug discovery result network 1.
#'
#' @param rv Variables list in the Shiny app.
#' @param cell_type1 The source cell type.
#' @param cell_type1 The target cell type.
#'
#' @return
#' @export
#'
#' @examples
loadTargetDrugNetwork1 <- function(rv, cell_type1, cell_type2) {
  rv$targetDrugNetworkTitle1UI = {
    tags$h3(
      paste(
        "Targets Drug Discovery Result for Downstream Network from",
        paste(rv$cell_type1, collapse = "+"),
        "to"
        ,
        paste(rv$cell_type2, collapse = "+"),
        sep = " "
      )
    )
  }
}

#' Load drug discovery result network 2.
#'
#' @param rv Variables list in the Shiny app.
#' @param cell_type1 The target cell type.
#' @param cell_type1 The source cell type.
#'
#' @return
#' @export
#'
#' @examples
loadTargetDrugNetwork2 <- function(rv, cell_type1, cell_type2) {
  rv$targetDrugNetworkTitle2UI = {
    tags$h3(
      paste(
        "Targets Drug Discovery Result for Downstream Network from",
        paste(rv$cell_type2, collapse = "+"),
        "to"
        ,
        paste(rv$cell_type1, collapse = "+"),
        sep = " "
      )
    )
  }
}



# Load raw scRNA-seq data uploading instruction expend
#'
#' @return
#' @export
#'
#' @examples
rnaDataExpendLoading <- function() {
  showModal(
    modalDialog(
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
          "The group or design file is optional, but we recommend that you upload it to obtain
                meaningful analysis results. The data you upload should be a .csv file or .RDS file
                containing a data frame with only one column. The length of the column should be
                equal to the number of cells in the read count data. Please refrain from including
                row names in your uploaded data."
        ),
      ),
      title = "Instructions for Uploading Upstream Analysis Data",
      easyClose = TRUE
    )
  )
}

#' Load downstream analysis data uploading instruction expend
#'
#' @return
#' @export
#'
#' @examples
geneExpressionExpendLoading <- function() {
  showModal(
    modalDialog(
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
          "The group or design file is required in the downstream analysis part.
                The data you upload should be a .csv file or .RDS file containing a
                data frame with one or two columns. The first column indicates the cell type for each cell.
                The second column is optional and can be used to specify a group or design for each cell.
                The length of the rows should be equal to the number of cells in the read count data.
                Please avoid including row names in your data."
        )
      ),
      title = "Instruction for Uploading Downstream Analysis Data",
      easyClose = TRUE
    )
  )
}



#' Load drug file uploading insturction expand.
#'
#' @return
#' @export
#'
#' @examples
drugFileExpendLoading <- function() {
  showModal(
    modalDialog(
      tags$h4(
        id = "drugProcessingInstruction",
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
      title = "Instruction for processing Drug File",
      easyClose = TRUE
    )
  )
}
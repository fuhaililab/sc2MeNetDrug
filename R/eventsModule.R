#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# All event modules used in the app. Mainly contain functions for processing 
# events conducted by the user.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#' Application event function for setting working directory.
#' The function will first change and save the new working directory selected
#' by user and then reload all data saved in the new directory to the application.
#'
#' @param defaultDir The root directory.
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
workingDirSetEvent <- function(defaultDir, input, rv, session) {
  tryCatch({
    if (is.integer(input$workingDir)) {
      stop(safeError("Directory does not exist."))
    }
    
    newDir <- parseDirPath(defaultDir, input$workingDir)
    initialization_rv()
    rv$outputDir <- newDir
    initialization_UI(rv)
    loadMarkerGene(rv, session)
    save(newDir, file = "./cache/lastWorkingDir.RData")
  },
  error = function(e) {
    if (e$message == "Directory does not exist.") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Load Error",
        content = "Please select a directory first.",
        style = "warning"
      )
    } else{
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Load Error",
        content = paste0(
          "Something wrong when setting working directory:",
          e$message
        ),
        style = "warning"
      )
    }
    
  })
  
}

#' Application event function for preprocessing scRNA-seq data.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#' @param use_rna_data If true, the function will preprocess the data uploaded from the upstream analysis.
#' @param imp If true, do the imputation step additionally.
#'
#' @return
#' @export
#'
#' @examples
preprocessingEvent <-
  function(input, rv, session, use_rna_data, imp) {
    progress <- shiny::Progress$new()
    progress$set(message = "Preprocessing:", value = 0)
    on.exit(progress$close())
    tryCatch({
      if (use_rna_data) {
        data <- rv$rna_df
        group_list <- rv$rna_group_list
        gene_list <- rv$rna_gene_list
        type_list <- NULL
      }
      else{
        if (is.null(rv$network_type_list)) {
          stop(safeError("no cell type data"))
        }
        data <- rv$network_df
        type_list <- rv$network_type_list
        group_list <- rv$network_group_list
        gene_list <- rv$network_gene_list
      }
      
      raw_num_cell <- ncol(data)
      raw_num_gene <- nrow(data)
      gene_num_after_convert <- raw_num_gene
      
      
      if (input$rnaMouseToHuman == T) {
        progress$set(0, detail = "Convert mouse gene symbol to human gene symbol")
        data <- mouse_to_human_convert(data, gene_list)
        gene_list <- rownames(data)
        gene_num_after_convert <- length(gene_list)
      }
      
      # convert data to Seurat data object
      if (!is.null(type_list)) {
        additional_type_data = data.frame(cell_type = type_list)
        rownames(additional_type_data) <- colnames(data)
        data <-
          CreateSeuratObject(counts = data, meta.data = additional_type_data)
        
      } else{
        data <- CreateSeuratObject(counts = data)
      }
      
      if (!is.null(group_list)) {
        data <- SetIdent(data, value = group_list)
      }
      
      #quality control
      progress$set(0.2, detail = "Remove low quality samples and genes")
      data <- quality_control(data)
      
      
      #normalization
      progress$set(0.4, detail = "Normalization")
      data <- sctransform_seurat(data)
      
      
      progress$set(0.5, detail = "Imputation")
      if (imp) {
        data <- imputation(data)
      }
      
      progress$set(0.95, detail = "Finish")
      final_num_cell <- dim(data)[2]
      if (use_rna_data) {
        save(data, file = paste0(rv$outputDir, "/rna_df.RData"))
        
        if (!is.null(group_list)) {
          group_list <- as.character(Idents(data))
          save(group_list,
               file = paste0(rv$outputDir, "/rnaGroupInformation.RData"))
        }
        
        rv$rna_df <- data
        rv$rna_gene_list <- rownames(data)
        rv$rna_group_list <- group_list
        
        drop_num_cell <- raw_num_cell - final_num_cell
        rm(data, group_list, gene_list)
        gc()
        #load dimension reduction part
        loadDimensionReduction()
      }
      else{
        
        type_list <- data@meta.data$cell_type
        
        if (!is.null(group_list)) {
          group_list <- as.character(Idents(data))
          save(group_list,
               file = paste0(rv$outputDir, "/networkGroupInformation.RData"))
          group_type_list <- paste(group_list, type_list, sep = "_")
          data <-
            SetIdent(data, value = group_type_list)
        } else{
          data <- SetIdent(data, value = type_list)
        }
        
        save(data, file = paste0(rv$outputDir, "/network_df.RData"))
        rv$network_df <- data
        

        
        drop_num_cell <- raw_num_cell - final_num_cell
        
        rm(data, group_list, type_list)
        gc()
        #load ExpressionPart
        loadGeneExpression1()
      }
      
      progress$set(1, detail = "Finish")
      if (input$rnaMouseToHuman == F) {
        createAlert(
          session = session,
          anchorId = "errorAlert",
          title = "Preprocessing Completed",
          content = paste0(
            "The data have a total of ",
            raw_num_cell,
            " cells and ",
            raw_num_gene,
            " genes. After preprocessing, ",
            drop_num_cell,
            " cells are removed."
          ),
          style = "success"
        )
      } else{
        createAlert(
          session = session,
          anchorId = "errorAlert",
          title = "Preprocessing finish",
          content = paste0(
            "The data have a total of ",
            raw_num_cell,
            " cells and ",
            raw_num_gene,
            " genes. After converting mouse to human gene,
                  The data have a total of ",
            gene_num_after_convert,
            " gene. After preprocessing, ",
            drop_num_cell,
            " cells are removed."
          ),
          style = "success"
        )
      }
      
    },
    error = function(e) {
      if (e$message == "no cell type data") {
        createAlert(
          session = session,
          anchorId = "errorAlert",
          title = "Computation Error",
          content = "Please upload cell type data first.",
          style = "danger"
        )
      } else{
        createAlert(
          session = session,
          anchorId = "errorAlert",
          title = "Computation Error",
          content = paste(
            "There appears to be an issue when conducting the preprocessing analysis.",
            e,
            sep =
              ":"
          ),
          style = "danger"
        )
      }
      
    })
    
    
    
  }




#' Application event function for dimension reduction analysis.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
drEvent <- function(input, rv, session) {
  tryCatch({
    progress <- shiny::Progress$new()
    progress$set(message = "Dimension reduction:", value = 0)
    on.exit(progress$close())
    
    data <- rv$rna_df
    progress$set(0.2, detail = "Run PCA")
    #PCA dimension reduction
    data <- PCA(data)
    pca_result = data.frame(data[["pca"]]@cell.embeddings)
    progress$set(0.5, detail = "Run UMAP")
    data <- UMAP(data, input$nPC)
    umap_result <- data.frame(data[["umap"]]@cell.embeddings)
    colnames(umap_result) <- c("UMAP_1", "UMAP_2")
    rv$df_dr <- umap_result
    rv$rna_df <- data
    #save results
    progress$set(0.7, detail = "Save results")
    save(pca_result, file = paste0(rv$outputDir, "/pca_result.RData"))
    save(umap_result,
         file = paste0(rv$outputDir, "/umap_result.RData"))
    save(data, file = paste0(rv$outputDir, "/rna_df.RData"))
    rm(umap_result, pca_result)
    gc()
    #load clustering part
    loadCluster()
    progress$set(1, detail = "Finish")
  },
  error = function(e) {
    createAlert(
      session = session,
      anchorId = "errorAlert",
      title = "Computation Error",
      content = paste(
        "There appears to be an issue when conducting the dimension reduction analysis.",
        e,
        sep = ":"
      ),
      style = "danger"
    )
  })
}


#' Application event function for clustering analysis.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
ClusterEvent <- function(input, rv, session) {
  tryCatch({
    progress <- shiny::Progress$new()
    progress$set(message = "Clustering:", value = 0)
    on.exit(progress$close())
    
    data <- rv$rna_df
    progress$set(0.2, detail = "Find neighbors.")
    data <- find_neighbors(data, input$neighbor_nPC)
    progress$set(0.4, detail = "Do clustering.")
    data <- clustering(data, input$cluster_resolution)
    clusters <- data@meta.data["seurat_clusters"]
    progress$set(0.7, detail = "Save results.")
    save(clusters, file = paste0(rv$outputDir, "/clusters.RData"))
    rv$df_cluster <- cbind(rv$df_dr, clusters)
    colnames(rv$df_cluster)[3] <- "cluster"
    save(data, file = paste0(rv$outputDir, "/rna_df.RData"))
    rv$rna_df <- data
    selected_gene_list <-
      data[[data@active.assay]]@var.features[c(1:8)]
    plot_results <- geneExpressionPlot(data, selected_gene_list)
    gene_expression_vln_plot <- plot_results[[1]]
    gene_expression_scatter_plot <- plot_results[[2]]
    
    rv$gene_expression_vln_plot <- gene_expression_vln_plot
    rv$gene_expression_scatter_plot <- gene_expression_scatter_plot
    save(
      gene_expression_vln_plot,
      file = paste0(rv$outputDir, "/gene_expression_vln_plot.RData")
    )
    save(
      gene_expression_scatter_plot,
      file = paste0(rv$outputDir, "/gene_expression_scatter_plot.RData")
    )
    
    rm(clusters, data, plot_results)
    gc()
    
    loadClassification()
    progress$set(1, detail = "Finish")
  },
  error = function(e) {
    createAlert(
      session = session,
      anchorId = "errorAlert",
      title = "Computation Error",
      content = paste(
        "There appears to be an issue when conducting the clustering analysis.",
        e,
        sep = ":"
      ),
      style = "danger"
    )
  })
}


#' Plot gene expression and distribution.
#'
#' @param seurat_data Seurat object.
#' @param selected_gene_list Selected gene list for the plot.
#'
#' @return
#' @export
#'
#' @examples
geneExpressionPlot <- function(seurat_data, selected_gene_list) {
  gene_expression_vln_plot <-
    VlnPlot(
      seurat_data,
      features = selected_gene_list,
      pt.size = 0.2,
      ncol = 4,
      group.by = "seurat_clusters"
    )
  
  gene_expression_scatter_plot <-
    FeaturePlot(seurat_data,
                features = selected_gene_list,
                pt.size = 0.2,
                ncol = 4)
  list(gene_expression_vln_plot, gene_expression_scatter_plot)
}


#' plot expression distribution violin plot and scatter plot
#' for user selected gene set.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
geneExpressionPlotEvent <- function(input, rv) {
  selected_gene_list <- input$genePlotSelection
  plot_results <- geneExpressionPlot(rv$rna_df, selected_gene_list)
  gene_expression_vln_plot <- plot_results[[1]]
  gene_expression_scatter_plot <- plot_results[[2]]
  
  rv$gene_expression_vln_plot <- gene_expression_vln_plot
  rv$gene_expression_scatter_plot <- gene_expression_scatter_plot
  save(
    gene_expression_vln_plot,
    file = paste0(rv$outputDir, "/gene_expression_vln_plot.RData")
  )
  save(
    gene_expression_scatter_plot,
    file = paste0(rv$outputDir, "/gene_expression_scatter_plot.RData")
  )
  rm(plot_results)
  gc()
}


#' Application event function for cell annotation analysis.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
classificationEvent <- function(input, rv, session) {
  tryCatch({
    clProgress <- shiny::Progress$new()
    clProgress$set(message = "Cell annotation", value = 0)
    on.exit(clProgress$close())
    df <- rv$rna_df[[rv$rna_df@active.assay]]@data
    annotation_result <-
      gseaClassification(
        df,
        rv$df_cluster$cluster,
        rv$rna_gene_list,
        rv$select_marker_gene,
        clProgress,
        rv$outputDir
      )
    clProgress$inc(0.95, detail = "Intgreate result")
    # prevent auto sort in merge
    df_cluster <- rv$df_cluster
    df_cluster$id <- c(1:nrow(df_cluster))
    #join classification result to every sample
    cell_annotation <-
      merge(df_cluster,
            annotation_result,
            by = "cluster",
            all.x = TRUE)
    cell_annotation <-
      cell_annotation[order(cell_annotation$id), ]
    cell_annotation <- cell_annotation[, -4]
    
    rownames(annotation_result) <- annotation_result$cluster
    annotation_result <- annotation_result[, c(1, 2)]
    save(annotation_result,
         file = paste0(rv$outputDir, "/annotation_result.RData"))
    save(cell_annotation,
         file = paste0(rv$outputDir, "/cell_annotation.RData"))
    
    rv$annotation_result <- annotation_result
    rv$df_classify <- cell_annotation
    rv$rna_type_list <- cell_annotation$type
    rm(df, df_cluster, annotation_result, cell_annotation)
    gc()
    loadGeneExpression1()
    clProgress$inc(1, detail = "Finish")
  },
  error = function(e) {
    createAlert(
      session = session,
      anchorId = "errorAlert",
      title = "Computation Error",
      content = paste(
        "There appears to be an issue when conducting the cell annotation analysis.",
        e,
        sep = ":"
      ),
      style = "danger"
    )
  })
}


#' Application event function for loading data in the downstream analysis part.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
geneExpressionDataLoadingEvent <- function(input, rv, session) {
  #load corresponding data user choose,
  tryCatch({
    loadingProgress <- shiny::Progress$new()
    loadingProgress$set(message = "Loading data", value = 0)
    on.exit(loadingProgress$close())
    
    if (input$dataUseSelect == 1) {
      if (is.null(rv$df_classify) || is.null(rv$rna_df)) {
        stop(safeError("noRnaData"))
      }
      rv$useRnaData = 1
      loadRnaDataSelection(rv)
      loadingProgress$inc(0.2, detail = "Set and save identity information for Seurat object.")
      setRnaDataIdents(rv)
      loadingProgress$inc(0.8, detail = "Compute cell distribution.")
      
      if (rv$rna_have_group == 0) {
        cell_count <- cellDistribution(rv$rna_type_list, NULL)
        rv$cell_count <- cell_count
        rv$distribution_have_group = 0
        
      } else{
        cell_count <-
          cellDistribution(rv$rna_type_list, rv$rna_group_list)
        rv$cell_count <- cell_count
        rv$distribution_have_group = 1
      }
    } else{
      if (is.null(rv$network_df)) {
        stop(safeError("noNetworkData"))
      }
      rv$useRnaData = 0
      loadNetworkDataSelection(rv)
      loadingProgress$inc(0.8, detail = "Compute cell distribution.")
      if (rv$network_have_group == 0) {
        cell_count <- cellDistribution(rv$network_type_list, NULL)
        rv$cell_count <- cell_count
        rv$distribution_have_group = 0
        
        
      } else{
        cell_count <-
          cellDistribution(rv$network_type_list, rv$network_group_list)
        rv$cell_count <- cell_count
        rv$distribution_have_group = 1
      }
    }
    distribution_have_group <- rv$distribution_have_group
    save(
      distribution_have_group,
      cell_count,
      file = paste0(rv$outputDir, "/cellDistributionPlot.RData")
    )
    rm(cell_count, distribution_have_group)
    gc()
    loadGeneExpression2()
    loadEPSelection(rv)
    loadGO(rv)
    loadUpReg(rv)
    if (file.exists("./data/rankMatrix92742a.RData") &&
        file.exists("./data/drug_mapping.RData")) {
      loadCommunication()
    }
    useRnaData = rv$useRnaData
    save(useRnaData,
         file = paste0(rv$outputDir, "/lastAnalysisDataIndex.RData"))
    loadingProgress$inc(1, detail = "Finish")
    
    
  },
  error = function(e) {
    if (e$message == "noRnaData") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Data Error",
        content = "There are no results available from raw data analysis.",
        style = "danger"
      )
    } else if (e$message == "noNetworkData") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Data Error",
        content = "There are no data available from downstream data uploading.",
        style = "danger"
      )
    } else{
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Loading Error",
        content = paste("Something wrong", e, sep = ":"),
        style = "danger"
      )
    }
  })
}



#' Application event function for cell-cell communication analysis.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
communicationEvent <- function(input, rv, session) {
  coProgress <- shiny::Progress$new()
  coProgress$set(message = "Discover Communication", value = 0)
  on.exit(coProgress$close())
  tryCatch({
    cell_type1 <- input$cellTypeSelect1
    cell_type2 <- input$cellTypeSelect2
    if (is.null(cell_type1) || is.null(cell_type2)) {
      stop(safeError("empty input"))
    }
    rv$cell_type1 <- cell_type1
    rv$cell_type2 <- cell_type2
    fc_thres <- input$networkFc
    pv_thres <- input$networkPValue
    fc_downThres <- input$networkDFc
    
    fcDir <- paste0(rv$outputDir, "/cellCommunication")
    if (!dir.exists(fcDir)) {
      dir.create(fcDir)
    }
    
    if (rv$useRnaData == 1) {
      data = rv$rna_df
      gene_symbol = rv$rna_gene_list
      type_list = rv$rna_type_list
      if (rv$rna_have_group == 0) {
        testGroup <- NULL
        normalGroup <- NULL
        group_list <- NULL
        
        networkDir <-
          paste0("/", paste(
            paste(rv$cell_type1, collapse = "+"),
            paste(rv$cell_type2, collapse = "+"),
            sep = "-"
          ))
      }
      else {
        testGroup <- input$controlGroup
        normalGroup <- input$normalGroup
        if (is.null(testGroup) || is.null(normalGroup)) {
          stop(safeError("empty input"))
        }
        networkDir <-
          paste0("/", paste(
            paste(rv$cell_type1, collapse = "+"),
            paste(rv$cell_type2, collapse = "+"),
            paste(normalGroup, collapse = "+"),
            paste(testGroup, collapse = "+"),
            sep = "-"
          ))
        
        group_list <- rv$rna_group_list
      }
    }
    else{
      data = rv$network_df
      gene_symbol = rv$network_gene_list
      type_list = rv$network_type_list
      if (rv$network_have_group == 0) {
        testGroup <- NULL
        normalGroup <- NULL
        group_list <- NULL
        networkDir <-
          paste0("/", paste(
            paste(rv$cell_type1, collapse = "+"),
            paste(rv$cell_type2, collapse = "+"),
            sep = "-"
          ))
        
      }
      else {
        testGroup <- input$controlGroup
        normalGroup <- input$normalGroup
        if (is.null(testGroup) || is.null(normalGroup)) {
          stop(safeError("empty input"))
        }
        group_list <- rv$network_group_list
        
        networkDir <-
          paste0("/", paste(
            paste(rv$cell_type1, collapse = "+"),
            paste(rv$cell_type2, collapse = "+"),
            paste(normalGroup, collapse = "+"),
            paste(testGroup, collapse = "+"),
            sep = "-"
          ))
        
      }
    }
    
    netDir <- paste0(fcDir, networkDir)
    if (!dir.exists(netDir)) {
      dir.create(netDir)
    }
    
    
    result <-
      communication(
        data = data,
        gene_symbol = gene_symbol,
        group_list = group_list,
        type_list = type_list,
        normalGroup = normalGroup,
        testGroup = testGroup,
        fc_thres = fc_thres,
        pv_thres = pv_thres,
        cell_type1 = rv$cell_type1,
        cell_type2 = rv$cell_type2,
        netDir = netDir,
        ligRecDatabase = rv$ligRecDatabase,
        TfTargetInteraction = rv$TfTargetInteraction,
        keggInfo = rv$keggInfo,
        useOld = input$networkUseOldResult,
        STRING = rv$STRING,
        resource = input$networkResource,
        padjust = input$networkPadjust,
        progress = coProgress
      )
    
    
    #load gene expression information table
    ligRecInformation <- result[[1]]
    rv$cell_type1_ligands <- ligRecInformation[[1]]
    rv$cell_type1_receptors <- ligRecInformation[[2]]
    rv$cell_type2_ligands <- ligRecInformation[[3]]
    rv$cell_type2_receptors <- ligRecInformation[[4]]
    rm(ligRecInformation)
    
    
    #set directory for file saving
    netDir <- paste0(fcDir, networkDir)
    network1Dir <-
      paste0(netDir, "/", paste(
        paste(rv$cell_type1, collapse = "+"),
        paste(rv$cell_type2, collapse = "+"),
        sep = "_"
      ))
    
    if (!dir.exists(network1Dir)) {
      dir.create(network1Dir)
    }
    
    network2Dir <-
      paste0(netDir, "/", paste(
        paste(rv$cell_type2, collapse = "+"),
        paste(rv$cell_type1, collapse = "+"),
        sep = "_"
      ))
    if (!dir.exists(network2Dir)) {
      dir.create(network2Dir)
    }
    
    
    
    type1_to_type2_result <- result[[2]]
    network_infor1 <- type1_to_type2_result[[1]]
    activated_network_nodes2 <- type1_to_type2_result[[2]]
    activated_network_edges2 <- type1_to_type2_result[[3]]
    downstream_network_nodes1 <- type1_to_type2_result[[4]]
    downstream_network_edges1 <- type1_to_type2_result[[5]]
    
    
    type2_to_type1_result <- result[[3]]
    network_infor2 <- type2_to_type1_result[[1]]
    activated_network_nodes1 <- type2_to_type1_result[[2]]
    activated_network_edges1 <- type2_to_type1_result[[3]]
    downstream_network_nodes2 <- type2_to_type1_result[[4]]
    downstream_network_edges2 <- type2_to_type1_result[[5]]
    
    
    rm(
      type1_to_type2_result,
      type2_to_type1_result,
      result,
      data,
      gene_symbol,
      group_list,
      type_list,
      testGroup,
      normalGroup
    )
    gc()
    
    # drug discovery part
    coProgress$set(0.8, detail = "Discover drug")
    if (is.null(rv$rankMatrix)) {
      loadRankMatrix(rv)
    }
    
    
    if (is.null(downstream_network_nodes1)) {
      targetDrug1 <- NULL
      signalingDrug1 <- NULL
    } else{
      #targets drug discoveriing
      targetDrug1 <-
        findDrug(
          downstream_network_nodes1,
          downstream_network_edges1,
          rv$drugBankInteraction,
          rv$drugBankInformation
        )
      #signaling signature drug discovering
      signalingDrug1 <-
        findDrug2(
          downstream_network_nodes1,
          rv$gSymZs,
          rv$rankMatrix,
          rv$drug_mapping,
          input$drugNumber,
          input$useFDAOnly,
          rv$drugBankInformation
        )
      
    }
    if (is.null(downstream_network_nodes2)) {
      targetDrug2 <- NULL
      signalingDrug2 <- NULL
    } else{
      #targets drug discoveriing
      targetDrug2 <-
        findDrug(
          downstream_network_nodes2,
          downstream_network_edges2,
          rv$drugBankInteraction,
          rv$drugBankInformation
        )
      #signaling signature drug discovering
      signalingDrug2 <-
        findDrug2(
          downstream_network_nodes2,
          rv$gSymZs,
          rv$rankMatrix,
          rv$drug_mapping,
          input$drugNumber,
          input$useFDAOnly,
          rv$drugBankInformation
        )
      
    }
    
    #load all the results to front and save corresponding files
    #only save file if network exist
    #if at least one of four network don't exist, display notification
    network_exist_flag <- 1
    
    if (!is.null(activated_network_nodes1)) {
      save(
        network_infor1,
        activated_network_nodes1,
        activated_network_edges1,
        cell_type1,
        file = paste0(netDir, "/activatedNetworkData1.RData")
      )
      display_data <-
        processActivateNetworkDisplay(activated_network_nodes1, activated_network_edges1)
      display_json <-
        networkDataToJson(display_data[[1]], display_data[[2]])
      rv$activated_network_json1 <- display_json
      rv$activated_network_nodes1 <- activated_network_nodes1
      rv$activated_network_edges1 <- activated_network_edges1
      shinyjs:::show("activatedNetworkPlotDownload1")
      rm(display_data, display_json)
    } else{
      if (input$networkResource == "KEGG") {
        network_exist_flag <- 0
      }
      rv$activated_network_json1 <- NULL
    }
    
    if (!is.null(downstream_network_nodes1)) {
      save(
        downstream_network_nodes1,
        downstream_network_edges1,
        cell_type1,
        cell_type2,
        file = paste0(network1Dir, "/downstreamNetworkData.RData")
      )
      display_json <-
        networkDataToJson(downstream_network_nodes1, downstream_network_edges1)
      rv$downstream_network_json1 <- display_json
      rv$downstream_network_nodes1 <- downstream_network_nodes1
      rv$downstream_network_edges1 <- downstream_network_edges1
      shinyjs:::show("networkPlotDownload1")
      rm(display_json)
      
    } else{
      network_exist_flag <- 0
      rv$downstream_network_json1 <- NULL
    }
    
    if (!is.null(activated_network_nodes2)) {
      save(
        network_infor2,
        activated_network_nodes2,
        activated_network_edges2,
        cell_type2,
        file = paste0(netDir, "/activatedNetworkData2.RData")
      )
      
      display_data <-
        processActivateNetworkDisplay(activated_network_nodes2, activated_network_edges2)
      display_json <-
        networkDataToJson(display_data[[1]], display_data[[2]])
      rv$activated_network_json2 <- display_json
      rv$activated_network_nodes2 <- activated_network_nodes2
      rv$activated_network_edges2 <- activated_network_edges2
      shinyjs:::show("activatedNetworkPlotDownload2")
      rm(display_data, display_json)
    } else{
      if (input$networkResource == "KEGG") {
        network_exist_flag <- 0
      }
      rv$activated_network_json2 <- NULL
    }
    
    if (!is.null(downstream_network_nodes2)) {
      save(
        downstream_network_nodes2,
        downstream_network_edges2,
        cell_type1,
        cell_type2,
        file = paste0(network2Dir, "/downstreamNetworkData.RData")
      )
      display_json <-
        networkDataToJson(downstream_network_nodes2, downstream_network_edges2)
      rv$downstream_network_json2 <- display_json
      rv$downstream_network_nodes2 <- downstream_network_nodes2
      rv$downstream_network_edges2 <- downstream_network_edges2
      shinyjs:::show("networkPlotDownload2")
      rm(display_json)
      
    } else{
      network_exist_flag <- 0
      rv$downstream_network_json2 <- NULL
    }
    
    if (network_exist_flag == 0) {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = "There is at least one network that doesn't exist.
        You might consider adjusting the thresholds.",
        style = "warning"
      )
    }
    #load title for each network
    loadNetwork1(rv, rv$cell_type1, rv$cell_type2)
    loadNetwork2(rv, rv$cell_type1, rv$cell_type2)
    loadActivatedNetwork1(rv, rv$cell_type1)
    loadActivatedNetwork2(rv, rv$cell_type2)
    
    #load drug discovering result
    
    if (is.null(signalingDrug1)) {
      rv$drug_table1 <- NULL
      rv$drug_json1 <- NULL
    } else{
      rv$drug_table1 <- signalingDrug1[[1]]
      display_data <-
        processDrugClusteringDisplay(signalingDrug1[[3]], signalingDrug1[[4]])
      display_json <-
        networkDataToJson(display_data[[1]], display_data[[2]])
      rv$drug_json1 <- display_json
      rm(display_data, display_json)
      save(signalingDrug1,
           file = paste0(network1Dir, "/signalingDrug.RData"))
    }
    
    if (is.null(signalingDrug2)) {
      rv$drug_table2 <- NULL
      rv$drug_json2 <- NULL
    } else{
      rv$drug_table2 <- signalingDrug2[[1]]
      display_data <-
        processDrugClusteringDisplay(signalingDrug2[[3]], signalingDrug2[[4]])
      display_json <-
        networkDataToJson(display_data[[1]], display_data[[2]])
      rv$drug_json2 <- display_json
      rm(display_data, display_json)
      save(signalingDrug2,
           file = paste0(network2Dir, "/signalingDrug.RData"))
    }
    
    
    
    if (is.null(targetDrug1)) {
      rv$drugNetworkJson1 <- NULL
      rv$drug_mapping_table1 <- NULL
      rv$targetDrug_json1 <- NULL
    } else{
      rv$drugNetworkJson1 <- targetDrug1[[2]]
      rv$drug_mapping_table1 <- targetDrug1[[1]]
      display_data <-
        processDrugClusteringDisplay(targetDrug1[[4]], targetDrug1[[5]])
      display_json <-
        networkDataToJson(display_data[[1]], display_data[[2]])
      rv$targetDrug_json1 <- display_json
      rm(display_data, display_json)
      save(targetDrug1, file = paste0(network1Dir, "/targetDrug.RData"))
    }
    
    if (is.null(targetDrug2)) {
      rv$drugNetworkJson2 <- NULL
      rv$drug_mapping_table2 <- NULL
      rv$targetDrug_json2 <- NULL
    } else{
      rv$drugNetworkJson2 <- targetDrug2[[2]]
      rv$drug_mapping_table2 <- targetDrug2[[1]]
      display_data <-
        processDrugClusteringDisplay(targetDrug2[[4]], targetDrug2[[5]])
      display_json <-
        networkDataToJson(display_data[[1]], display_data[[2]])
      rv$targetDrug_json2 <- display_json
      rm(display_data, display_json)
      save(targetDrug2, file = paste0(network2Dir, "/targetDrug.RData"))
    }
    #load target drug network title
    loadTargetDrugNetwork1(rv, rv$cell_type1, rv$cell_type2)
    loadTargetDrugNetwork2(rv, rv$cell_type1, rv$cell_type2)
    
    #refresh all the D3 plot
    #remove previous nodes when user load network plot more than once
    if (rv$haveGenerated1 == 1) {
      session$sendCustomMessage("refreshNetwork1", "refresh")
    }
    rv$haveGenerated1 = 1
    rm(downstream_network_edges1, downstream_network_nodes1)
    
    #remove previous nodes when user load network plot more than once
    if (rv$haveGenerated2 == 1) {
      session$sendCustomMessage("refreshNetwork2", "refresh")
    }
    rv$haveGenerated2 = 1
    rm(downstream_network_edges2, downstream_network_nodes2)
    
    #remove previous nodes when user load network plot more than once
    if (rv$activatedNetworkGenerated1 == 1) {
      session$sendCustomMessage("refreshActivatedNetwork1", "refresh")
    }
    rv$activatedNetworkGenerated1 = 1
    rm(activated_network_edges1, activated_network_nodes1)
    
    #remove previous nodes when user load network plot more than once
    if (rv$activatedNetworkGenerated2 == 1) {
      session$sendCustomMessage("refreshActivatedNetwork2", "refresh")
    }
    rv$activatedNetworkGenerated2 = 1
    rm(activated_network_edges2, activated_network_nodes2)
    
    
    
    #remove previous nodes when user load network plot more than once
    if (rv$drugGenerated == 1) {
      session$sendCustomMessage("refreshCluster1", "refresh")
      session$sendCustomMessage("refreshCluster2", "refresh")
    }
    rv$drugGenerated = 1
    #remove previous nodes when user load network plot more than once
    if (rv$drugNetworkGenerated == 1) {
      session$sendCustomMessage("refreshDrugNetwork1", "refresh")
      session$sendCustomMessage("refreshDrugNetwork2", "refresh")
      session$sendCustomMessage("refreshTargetCluster1", "refresh")
      session$sendCustomMessage("refreshTargetCluster2", "refresh")
    }
    rv$drugNetworkGenerated = 1
    
    save(cell_type1,
         cell_type2,
         netDir,
         file = paste0(rv$outputDir, "/lastNetworkIndexing.RData"))
    rm(
      targetDrug1,
      targetDrug2,
      signalingDrug1,
      signalingDrug2,
      cell_type1,
      cell_type2
    )
    gc()
    coProgress$set(1, detail = "Finish")
    
  },
  error = function(e) {
    if (e$message == "empty input") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = "Empty input.",
        style = "danger"
      )
    } else if (e$message == "no old result") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = paste("You don't have any previously saved results."),
        style = "danger"
      )
    }
    else{
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = paste(
          "There seems to be an issue when performing the downstream network analysis.
                        You might want to consider adjusting the thresholds.",
          e,
          sep = ":"
        ),
        style = "danger"
      )
    }
  })
}

#' Application event function for finding up-regulated ligands
#' and receptors using Seurat DEG analysis.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
upRegulateEvent <- function(input, rv, session) {
  tryCatch({
    fc_thres <- input$upRegFc
    pv_thres <- input$upRegPValue
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculate fold change and P-value", value = 0)
    upRegDir <- paste0(rv$outputDir, "/upRegulatedLigandsReceptors")
    if (!dir.exists(upRegDir)) {
      dir.create(upRegDir)
    }
    if (rv$useRnaData == 1) {
      data <- rv$rna_df
      cell_type_list <- rv$rna_type_list
      gene_list <- rv$rna_gene_list
      if (rv$rna_have_group == 0) {
        testGroup <- NULL
        normalGroup <- NULL
        groupDir <- paste0(upRegDir, "/noGroup")
        if (!dir.exists(groupDir)) {
          dir.create(groupDir)
        }
      }
      else {
        if (is.null(input$upRegTestGroup) ||
            is.null(input$upRegNormalGroup)) {
          stop(safeError("empty input"))
        }
        
        testGroup <- input$upRegTestGroup
        normalGroup <- input$upRegNormalGroup
        cell_group_list <- rv$rna_group_list
        
        groupDir <-
          paste0(upRegDir, paste0("/", paste(
            paste(normalGroup, collapse = "+"),
            paste(testGroup, collapse = "+"),
            sep = "-"
          )))
        if (!dir.exists(groupDir)) {
          dir.create(groupDir)
        }
      }
    }
    else {
      data <- rv$network_df
      cell_type_list <- rv$network_type_list
      gene_list <- rv$network_gene_list
      if (rv$network_have_group == 0) {
        testGroup <- NULL
        normalGroup <- NULL
        groupDir <- paste0(upRegDir, "/noGroup")
        if (!dir.exists(groupDir)) {
          dir.create(groupDir)
        }
      }
      else {
        if (is.null(input$upRegTestGroup) ||
            is.null(input$upRegNormalGroup)) {
          stop(safeError("empty input"))
        }
        testGroup <- input$upRegTestGroup
        normalGroup <- input$upRegNormalGroup
        cell_group_list <- rv$network_group_list
        groupDir <-
          paste0(upRegDir, paste0("/", paste(
            paste(normalGroup, collapse = "+"),
            paste(testGroup, collapse = "+"),
            sep = "-"
          )))
        if (!dir.exists(groupDir)) {
          dir.create(groupDir)
        }
      }
    }
    
    if (input$upRegUseOldResult == 0) {
      genesInformation <-
        upRegTest(
          data = data,
          cell_type_list = cell_type_list,
          normalGroup = normalGroup,
          testGroup = testGroup,
          cell_group_list = cell_group_list,
          progress = progress
        )
    }
    else{
      if (file.exists(paste0(groupDir, "/genesInformation.RData"))) {
        load(paste0(groupDir, "/genesInformation.RData"))
      } else{
        stop(safeError("no old result"))
      }
    }
    
    save(genesInformation,
         file = paste0(groupDir, "/genesInformation.RData"))
    
    
    upReg_network <-
      up_regulated_network(
        genesInformation,
        gene_list,
        fc_thres,
        fc_thres,
        pv_thres,
        ligRecDatabase = rv$ligRecDatabase
      )
    upStreamData <-
      upStreamNetwork(
        genesInformation,
        gene_list,
        fc_thres,
        pv_thres,
        ligRecDatabase = rv$ligRecDatabase
      )
    
    progress$set(0.9, detail = "Integration")
    rv$upRegLigand_network <- upReg_network[[1]][[1]]
    rv$upRegReceptor_network <- upReg_network[[2]][[1]]
    if (rv$upRegGenerated == 1) {
      session$sendCustomMessage("refreshUpRegReceptor", "refresh")
      session$sendCustomMessage("refreshUpRegLigand", "refresh")
    }
    rv$upRegGenerated = 1
    if (!is.null(upStreamData[[1]])) {
      rv$up_exp_network <- upStreamData[[1]][[1]]
    } else{
      rv$up_exp_network <- NULL
    }
    if (!is.null(upStreamData[[2]])) {
      rv$exp_up_network <- upStreamData[[2]][[1]]
      
    } else{
      rv$exp_up_network <- NULL
    }
    if (!is.null(upStreamData[[3]])) {
      rv$combine_network <- upStreamData[[3]][[1]]
    } else{
      rv$combine_network <- NULL
    }
    if (!is.null(upStreamData[[4]])) {
      rv$up_up_network <- upStreamData[[4]][[1]]
      
    } else{
      rv$up_up_network <- NULL
    }
    if (rv$upStreamGenerated == 1)
    {
      session$sendCustomMessage("refreshUpExp", "refresh")
      session$sendCustomMessage("refreshExpUp", "refresh")
      session$sendCustomMessage("refreshCombine", "refresh")
      session$sendCustomMessage("refreshUpUp", "refresh")
    }
    rv$upStreamGenerated = 1
    save(upReg_network, file = paste0(groupDir, "/upRegNetwork.RData"))
    save(upStreamData, file = paste0(groupDir, "/upStreamNetwork.RData"))
    save(groupDir,
         file = paste0(rv$outputDir, "/lastUpRegIndexing.RData"))
    if (is.null(rv$up_up_network) ||
        is.null(rv$combine_network) ||
        is.null(rv$exp_up_network) || is.null(rv$up_exp_network)) {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = "At least one of the upstream networks doesn't exist.
                  You might want to consider adjusting the thresholds.",
        style = "warning"
      )
    }
    rm(upReg_network, upStreamData, genesInformation)
    gc()
    progress$set(1, detail = "Finish")
    
  },
  error = function(e) {
    if (e$message == "no old result") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = paste("You do not have any previously saved results."),
        style = "danger"
      )
    } else if (e$message == "empty input") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = "Empty input.",
        style = "danger"
      )
    }
    else{
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = paste(
          "There appears to be an issue when conducting the up-regulated gene
                        analysis. You might want to consider adjusting the thresholds.",
          e
        ),
        style = "danger"
      )
    }
  })
  
}

#' Application event function for computing EMT-PRO score.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
emt_proEvent <- function(input, rv, session) {
  tryCatch({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Compute EMT-PRO score", value = 0)
    if (rv$useRnaData == 1) {
      data <- rv$rna_df
      gene_list <- rv$rna_gene_list
      type_list <- rv$rna_type_list
      group_list <- rv$rna_group_list
      if (rv$rna_have_group == 1) {
        if (is.null(input$epGroup)) {
          stop(safeError("empty input"))
        }
        epGroup <- input$epGroup
        
      }
      else {
        epGroup <- NULL
      }
    }
    else{
      data <- rv$network_df
      gene_list <- rv$network_gene_list
      type_list <- rv$network_type_list
      group_list <- rv$network_group_list
      if (rv$network_have_group == 1) {
        if (is.null(input$epGroup)) {
          stop(safeError("empty input"))
        }
        epGroup <- input$epGroup
        
      }
      else {
        epGroup <- NULL
      }
    }
    
    if (is.null(input$epType)) {
      stop(safeError("empty input"))
    }
    else{
      epType <- input$epType
    }
    
    df <- as.matrix(data[[data@active.assay]]@data)
    progress$set(message = "Compute EMT-PRO score", value = 0.4)
    result <-
      emt_pro(df,
              gene_list,
              group_list,
              type_list,
              epGroup,
              epType)
    emt_score <- result[[1]]
    pro_score <- result[[2]]
    
    
    progress$set(message = "Compute EMT-PRO score", value = 0.9)
    emt_pro_score <-
      data.frame(
        EMT = emt_score,
        PRO = pro_score,
        type = paste(epType, collapse = "-")
      )
    
    rv$emt_pro_score <- emt_pro_score
    save(emt_pro_score,
         file = paste0(rv$outputDir, "/emt_pro_score.RData"))
    progress$set(message = "Finish", value = 1)
    rm(emt_pro_score, result, df, emt_score, pro_score)
    gc()
  },
  error = function(e) {
    if (e$message == "no enough gene") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = "There are not enough marker genes available for computation.",
        style = "danger"
      )
    } else if (e$message == "empty input") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = "Empty input",
        style = "danger"
      )
    }
    else{
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = paste(
          "There seems to be an issue when conducting EMT-PRO analysis.",
          e,
          sep =
            ":"
        ),
        style = "danger"
      )
    }
  })
  
}


#' Add a new marker gene to the marker gene table.
#' If the input is not specified by user, name the cell as "new cell type"
#' and gene as "new gene".
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
addRowButtonEvent <- function(input, rv) {
  cell_type <- input$newCellType
  if (cell_type == "" | is.null(cell_type)) {
    cell_type = "new cell type"
  }
  gene_list <- input$markerGeneTableSelection
  if (is.null(gene_list)) {
    gene_list <- c("new gene")
  }
  
  for (gene in gene_list) {
    new_row <- c(cell_type, gene)
    rv$display_marker_gene <- rbind(rv$display_marker_gene, new_row)
  }
  
}




#' Delete the selected marker gene from the marker gene table.
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
deleteRowButtonEvent <- function(input, rv) {
  if (!is.null(input$geneTable_rows_selected)) {
    rv$display_marker_gene <-
      rv$display_marker_gene[-as.numeric(input$geneTable_rows_selected), ]
  }
}

#' Select all cell type related to the Alzheimer's disease for cell annotation.
#'
#'
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
adButtonEvent <- function(rv, session) {
  cell_list <- unique(rv$display_marker_gene$cell_type)
  updateCheckboxGroupInput(
    session = session,
    inputId = "cellSelection",
    label = NULL,
    choices = cell_list,
    choiceValues = cell_list,
    selected = c("Endothelial", "Per", "In", "Ex", "Ast", "OPC", "Oli", "Mic")
  )
}

#' Select all cell type related to the pancreatic cancer for cell annotation.
#'
#'
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
pdacButtonEvent <- function(rv, session) {
  cell_list <- unique(rv$display_marker_gene$cell_type)
  updateCheckboxGroupInput(
    session = session,
    inputId = "cellSelection",
    label = NULL,
    choices = cell_list,
    choiceValues = cell_list,
    selected = c(
      "Ductal1",
      "Ductal2",
      "Acinar",
      "Endocrine",
      "Endothelial",
      "Fibroblast",
      "Stellate",
      "Macrophage",
      "Tcell",
      "Bcell"
    )
  )
}


#' Only display the marker gene table that the cell types are selected by the user.
#'
#' @param rv Variables list in the Shiny app.
#' @param input Input variables of the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
cellSelectionEvent <- function(rv, input) {
  #if user do not choose any type, display original table
  if (is.null(input$cellSelection)) {
    rv$select_marker_gene <- rv$display_marker_gene
  }
  else{
    rv$select_marker_gene <-
      rv$display_marker_gene[rv$display_marker_gene$cell_type
                             %in% input$cellSelection,]
  }
}


#' Display the original marker gene table.
#'
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
originalButtonEvent <- function(rv) {
  rv$display_marker_gene <- rv$original_marker_gene
}

#download cell-cell communication plot
downloadNetworkEvent <- function(session) {
  session$sendCustomMessage("downloadNetwork", "dir")
}


#' Application event function for Gene Ontology analysis given the selecetd
#' cell type and fold-change threshold.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
GOEvent <- function(input, rv, session) {
  tryCatch({
    progress <- shiny::Progress$new()
    progress$set(message = "GO analytics", value = 0)
    on.exit(progress$close())
    
    cell_type <- input$GOCellType
    upFc_thres <- input$GOUpFc
    dnFc_thres <- input$GODnFc
    pv_thres <- input$GOPValue
    
    fcDir <- paste0(rv$outputDir, "/GeneOntology")
    if (!dir.exists(fcDir)) {
      dir.create(fcDir)
    }
    
    if (is.null(cell_type)) {
      stop(safeError("empty input"))
    }
    if (rv$useRnaData == 1) {
      gene_list <- rv$rna_gene_list
      df = rv$rna_df[[rv$rna_df@active.assay]]@data
      
      if (rv$rna_have_group == 0) {
        t_keep_index <- rv$rna_type_list %in% cell_type
        n_keep_index <- !rv$rna_type_list %in% cell_type
        GoDir <- paste0(fcDir, paste0("/", paste(cell_type, collapse = "+")))

      }
      else{
        testGroup <- input$GOTestGroup
        normalGroup <- input$GONormalGroup
        
        
        t_keep_index <- rv$rna_type_list == cell_type &
          rv$rna_group_list %in% testGroup
        n_keep_index <- rv$rna_type_list == cell_type &
          rv$rna_group_list %in% normalGroup
        
        GoDir <-
          paste0(fcDir, paste0("/", paste(
            paste(cell_type, collapse = "+"),
            paste(normalGroup, collapse = "+"),
            paste(testGroup, collapse = "+"),
            sep = "-"
          )))
      }
    }
    else {
      gene_list <- rv$network_gene_list
      df = rv$network_df[[rv$network_df@active.assay]]@data
      
      if (rv$network_have_group == 0) {
        t_keep_index <- rv$network_type_list %in% cell_type
        n_keep_index <- !rv$network_type_list %in% cell_type
        GoDir <- paste0(fcDir, paste0("/", paste(cell_type, collapse = "+")))
        
      }
      else{
        testGroup <- input$GOTestGroup
        normalGroup <- input$GONormalGroup
        t_keep_index <- rv$network_type_list == cell_type &
          rv$network_group_list %in% testGroup
        n_keep_index <- rv$network_type_list == cell_type &
          rv$network_group_list %in% normalGroup

        
        GoDir <-
          paste0(fcDir, paste0("/", paste(
            paste(cell_type, collapse = "+"),
            paste(normalGroup, collapse = "+"),
            paste(testGroup, collapse = "+"),
            sep = "-"
          )))
      }
    }
    
    if (!dir.exists(GoDir)) {
      dir.create(GoDir)
    }
    df_n <- df[, n_keep_index]
    df_t <- df[, t_keep_index]
    
    
    GO_result <-
      GOAnalytics(df_t,
                  df_n,
                  gene_list,
                  upFc_thres,
                  dnFc_thres,
                  pv_thres,
                  rv$outputDir,
                  progress)
    

    save(GO_result, file = paste0(GoDir, "/GO_result.RData"))
    rv$GO_up_table <- GO_result[[1]]
    rv$GO_dn_table <- GO_result[[2]]
    rv$netGO_up <- GO_result[[3]]
    rv$netGO_dn <- GO_result[[4]]
    GO_list <-
      unique(c(
        as.character(rv$GO_up_table[, 1]),
        as.character(rv$GO_dn_table[, 1])
      ))
    GO_name <-
      unique(c(as.character(
        paste(rv$GO_up_table[, 1], rv$GO_up_table[, 2], sep = "  ")
      ),
      as.character(
        paste(rv$GO_dn_table[, 1], rv$GO_dn_table[, 2], sep = "  ")
      )))
    names(GO_list) <- GO_name
    rv$GO_list <- GO_list
    selectGONetwork(input, rv, session)
    loadGONetworkSelect(rv, rv$GO_list)
    
    
    save(GoDir, file = paste0(rv$outputDir, "/lastGOIndexing.RData"))
    
    
    rm(df, df_n, df_t, GO_result, GO_list, GO_name)
    gc()
    progress$set(1, detail = "Finish")
  },
  error = function(e) {
    if (e$message == "empty input") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = "Empty input",
        style = "danger"
      )
    }
    createAlert(
      session = session,
      anchorId = "errorAlert",
      title = "Computation Error",
      content = paste("Some thing wrong when carry out GO analytics", e, sep =
                        ":"),
      style = "danger"
    )
  })
  
  
}


#' Display the Gene Ontology network selected by the user.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
selectGONetwork <- function(input, rv, session) {
  GO <- input$GOList
  GO <- rv$GO_list[rv$GO_list == GO]
  GO_name <- names(GO)
  rv$GONetwork <-
    GONetworkGenerating(rv$netGO_up, rv$netGO_dn, GO, GO_name)
  if (rv$GOGenerated == 1) {
    session$sendCustomMessage("refreshGONetwork1", "refresh")
  }
  rv$GOGenerated = 1
  rv$GOData = 1
  rm(GO, GO_name)
  gc()
}


#' Application event function to process drug file for drug-related analysis.
#'
#' @param defaultDir Directory for loading drug files.
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
processDrugFileEvent <- function(defaultDir, input, rv, session) {
  progress <- shiny::Progress$new()
  progress$set(message = "Processing Drug File", value = 0)
  on.exit(progress$close())
  tryCatch({
    if (is.integer(input$drugFileDir)) {
      stop(safeError("no directory"))
    }
    dataDir <- parseDirPath(defaultDir, input$drugFileDir)
    if (file.exists(
      paste(
        dataDir,
        'GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx',
        sep = "/"
      )
    ) &&
    file.exists(paste(dataDir, 'GSE92742_Broad_LINCS_sig_info.txt', sep =
                      '/')) &&
    file.exists(paste(dataDir, 'GSE92742_Broad_LINCS_gene_info.txt', sep =
                      '/')) &&
    file.exists(paste(dataDir, "GSE92742_Broad_LINCS_pert_info.txt", sep =
                      "/"))) {
      dataList <-
        processDrugFile(dataDir, rv$drugBankInformation, progress)
      rv$rankMatrix <- dataList[[1]]
      rv$drug_mapping <- dataList[[2]]
      rv$gSymZs <- dataList[[3]]
      loadRankMatrix(rv)
      rm(dataList)
      gc()
      progress$set(1, detail = "Finish")
    } else{
      stop(safeError("no data"))
    }
  },
  error = function(e) {
    if (e$message == "no directory") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = "You select wrong directory",
        style = "danger"
      )
    } else if (e$message == "no data") {
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = "Some data don't exist in directory, check instruction carefully.",
        style = "danger"
      )
    } else{
      createAlert(
        session = session,
        anchorId = "errorAlert",
        title = "Computation Error",
        content = paste("Some thing wrong:", e, sep = ":"),
        style = "danger"
      )
    }
  })
}

#' Display user selected cell-cell communication networks.
#'
#' @param input Input variables of the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
networkDisplayEvent <- function(input) {
  if (input$networkDisplaySelect == "activated1") {
    shinyjs:::show("communicationPart3")
    shinyjs:::hide("communicationPart4")
    shinyjs:::hide("communicationPart5")
    shinyjs:::hide("communicationPart6")
  } else if (input$networkDisplaySelect == "activated2") {
    shinyjs:::show("communicationPart4")
    shinyjs:::hide("communicationPart3")
    shinyjs:::hide("communicationPart5")
    shinyjs:::hide("communicationPart6")
  } else if (input$networkDisplaySelect == "downstream1") {
    shinyjs:::show("communicationPart5")
    shinyjs:::hide("communicationPart3")
    shinyjs:::hide("communicationPart4")
    shinyjs:::hide("communicationPart6")
  } else{
    shinyjs:::show("communicationPart6")
    shinyjs:::hide("communicationPart3")
    shinyjs:::hide("communicationPart4")
    shinyjs:::hide("communicationPart5")
  }
}


#' Display user-selected drug networks.
#'
#' @param input Input variables of the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
drugDisplayEvent <- function(input) {
  if (input$targetDrugDisplaySelect == "downstream1") {
    shinyjs:::show("communicationPart7")
    shinyjs:::show("communicationPart8")
    shinyjs:::show("communicationPart11")
    shinyjs:::hide("communicationPart9")
    shinyjs:::hide("communicationPart10")
    shinyjs:::hide("communicationPart12")
  } else{
    shinyjs:::hide("communicationPart7")
    shinyjs:::hide("communicationPart8")
    shinyjs:::hide("communicationPart11")
    shinyjs:::show("communicationPart9")
    shinyjs:::show("communicationPart10")
    shinyjs:::show("communicationPart12")
  }
}


#' Save the modified marker gene table to the working directory.
#' Next time, load the saved marker gene table instead of the original one as long as the working
#' directory not changed.
#'
#' @param rv Variables list in the Shiny app.
#' @param session Running session variables in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
saveMarkerGeneEvent <- function(rv, session) {
  progress <- shiny::Progress$new()
  progress$set(message = "Save new marker gene tablee", value = 0)
  on.exit(progress$close())
  tryCatch({
    if (!is.null(rv$display_marker_gene)) {
      if (!is.null(rv$outputDir)) {
        marker_gene_table <- rv$display_marker_gene
        save(
          marker_gene_table,
          file = paste(rv$outputDir, "marker_gene_table.RData", sep = "/")
        )
      } else{
        stop(safeError("Please set working directory first."))
      }
      
    } else{
      stop(safeError("Marker gene tale do not exist."))
    }
  },
  error = function(e) {
    createAlert(
      session = session,
      anchorId = "errorAlert",
      title = "Save Error",
      content = e$message,
      style = "warning"
    )
  })
  
  
  progress$set(1, detail = "Finish.")
}

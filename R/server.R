# Set The maximum upload size to 8GB, which should be suitable in the most of the cases.  
options(shiny.maxRequestSize=8*1024^3)

# Load all modules.
source("loadingModule.R")
source("preprocessing.R")
source("utils.R")
source("reactiveModule.R")
source("eventsModule.R")
source("autoEncoder.R")
source("dimensionReduction.R")
source("continueWorkModule.R")
source("clustering.R")
source("subClustering.R")
source("gseaClassification.R")
source("iCSC.R")
source("findDrug.R")
source("upRegulatedGenes.R")
source("GOAnalytics.R")
source("EMTPRO.R")
source("processDrugFile.R")

# Main APP function.
function(input, output,session) {
  # Initialize all global variables.
  rv<-reactiveValues(
    keggInfo=NULL,
    original_marker_gene=NULL,
    marker_gene=NULL,
    select_marker_gene=NULL,
    display_marker_gene=NULL,
    bubbleData=NULL,
    df_table=NULL,
    cellSelectionUI={},
    adButtonUI={},
    DieaseSpecificTitleUI={},
    cellTypeSelect1UI={},
    cellTypeSelect1UI={},
    networkTableDisplayUI={},
    drugTable1UI={},
    drugTable2UI={},
    networkData=0,
    drugNetworkData=0,
    drug_table1=NULL,
    drug_table2=NULL,
    rna_df=NULL,
    rna_have_group=0,
    rna_num_sample=0,
    rna_group_list=NULL,
    rna_gene_list=NULL,
    rna_type_list=NULL,
    network_df=NULL,
    network_have_group=0,
    network_group_list=NULL,
    network_gene_list=NULL,
    network_type_list=NULL,
    df_dr=NULL,
    encoder_result=NULL,
    df_cluster=NULL,
    sub_df_cluster=NULL,
    df_classify=NULL,
    rna_data_mean=NULL,
    network_data_mean=NULL,
    useRnaData=1,
    cell_type1=NULL,
    cell_type2=NULL,
    downstream_network_nodes1=NULL,
    downstream_network_edges1=NULL,
    downstream_network_json1=NULL,
    networkData1=0,
    haveGenerated1=0,
    networkTitle1UI={},
    downstream_network_nodes2=NULL,
    downstream_network_edges2=NULL,
    downstream_network_json2=NULL,
    networkData2=0,
    haveGenerated2=0,
    networkTitle2UI={},
    cell_type1_ligands=NULL,
    cell_type1_receptors=NULL,
    cell_type2_ligands=NULL,
    cell_type2_receptors=NULL,
    drug_network_data=NULL,
    drug_table1=NULL,
    drug_table2=NULL, 
    upRegLigand_network=NULL,
    upRegReceptor_network=NULL,
    upRegGenerated=0,
    controlGroupUI={},
    normalGroupUI={},
    networkDownloadUI={},
    cell_count=NULL,
    distribution_have_group=0,
    EPcellTypeSelectUI={},
    EPGroupUI={},
    outputDir =NULL,
    gSymZs=NULL,
    rankMatrix=NULL,
    drug_mapping=NULL,
    drugData1=0,
    drugData2=0,
    drugTableTitle1UI={},
    drugTableTitle2UI={},
    drugGenerated=0,
    drugClusterData1=0,
    drugClusterData2=0,
    drug_json1=NULL,
    drug_json2=NULL,
    networkCellSelectTitleUI={},
    GONetwork1=NULL,
    GOGenerated=0,
    GOData=0,
    GOTestGroupUI={},
    GONormalGroupUI={},
    GOCellTypeUI={},
    GO_up_table=NULL,
    GO_dn_table=NULL,
    GO_table=NULL,
    GOTableDisplayUI={},
    GONetworkSelectUI={},
    GOList=NULL,
    netGO_up=NULL,
    netGO_dn=NULL,
    upRegTestGroupUI={},
    upRegNormalGroupUI={},
    upStreamGenerated=0,
    up_exp_network=NULL,
    exp_up_network=NULL,
    combine_network=NULL,
    up_up_network=NULL,
    DLRP=NULL,
    nicheNet=NULL,
    baderLab=NULL,
    ligRecDatabase=NULL,
    drugNetworkGenerated=0,
    drugNetworkJson1=NULL,
    drugNetworkJson2=NULL,
    drug_mapping_table1=NULL,
    drug_mapping_table2=NULL,
    drugNetworkJson1=NULL,
    drugNetworkJson2=NULL,
    drugBankInteraction=NULL,
    drugBankInformation=NULL,
    targetDrug_json1=NULL,
    targetDrug_json2=NULL,
    activatedNetworkTitle1UI={},
    activatedNetworkTitle2UI={},
    activatedNetworkGenerated1=0,
    activatedNetworkGenerated2=0,    
    targetDrugTableTitle1UI={},
    targetDrugTableTitle2UI={},
    targetDrugNetworkTitle1UI={},
    targetDrugNetworkTitle2UI={}
  )
  

  # Dynamic UI render functions.
  output$cellSelection<-renderUI(rv$cellSelectionUI)
  output$DieaseSpecificTitle<-renderUI(rv$DieaseSpecificTitleUI)
  output$cellTypeSelect1<-renderUI({rv$cellTypeSelect1UI})
  output$cellTypeSelect2<-renderUI(rv$cellTypeSelect2UI)
  output$networkTableDisplay<-renderUI(load_table())
  output$drugTableDisplay1<-renderUI(load_drug_table1())
  output$drugTableDisplay2<-renderUI(load_drug_table2())
  output$normalGroup<-renderUI(rv$normalGroupUI)
  output$controlGroup<-renderUI(rv$controlGroupUI)
  output$activatedNetworkTitle1<-renderUI(rv$activatedNetworkTitle1UI)
  output$activatedNetworkTitle2<-renderUI(rv$activatedNetworkTitle2UI)
  output$networkTitle1<-renderUI(rv$networkTitle1UI)
  output$networkTitle2<-renderUI(rv$networkTitle2UI)
  output$EPcellTypeSelect<-renderUI(rv$EPcellTypeSelectUI)
  output$EPGroup<-renderUI(rv$EPGroupUI)
  output$drugTableTitle1<-renderUI(rv$drugTableTitle1UI)
  output$drugTableTitle2<-renderUI(rv$drugTableTitle2UI)
  output$GOCellType<-renderUI(rv$GOCellTypeUI)
  output$GOTestGroup<-renderUI(rv$GOTestGroupUI)
  output$GONormalGroup<-renderUI(rv$GONormalGroupUI)
  output$GOTableDisplay<-renderUI(load_GO_table())
  output$drugMappingTable1Display<-renderUI(load_drug_mapping_table1())
  output$drugMappingTable2Display<-renderUI(load_drug_mapping_table2())
  output$GONetworkSelect<-renderUI(rv$GONetworkSelectUI)
  output$upRegTestGroup<-renderUI(rv$upRegTestGroupUI)
  output$upRegNormalGroup<-renderUI(rv$upRegNormalGroupUI)
  output$targetDrugTableTitle1<-renderUI(rv$targetDrugTableTitle1UI)
  output$targetDrugTableTitle2<-renderUI(rv$targetDrugTableTitle2UI)
  output$targetDrugNetworkTitle1<-renderUI(rv$targetDrugNetworkTitle1UI)
  output$targetDrugNetworkTitle2<-renderUI(rv$targetDrugNetworkTitle2UI)  
  
  # setting and display current work directory
  volumes <- c(getVolumes()())
  defaultDir<-getVolumes()()
  shinyDirChoose(input,"workingDir",roots=volumes, filetypes=c('', 'txt'))
  output$currentDir<-renderPrint({
    if (is.null(rv$outputDir)) {
      cat("No directory selected")
    } else {
      cat("Current working Directory: ",rv$outputDir)
    }
  })
  
  shinyDirChoose(input,"drugFileDir",roots=volumes, filetypes=c('', 'txt'))

  
  # Download data button
  output$adDownload <- downloadHandler(
    filename = function() {
      "autoencoderResult.csv"
    },
    content = function(file) {
      tryCatch(
        {
          if(is.null(rv$encoder_result)){
            stop(safeError("No data exist."))
          }
          write.csv(rv$encoder_result, file)
        },
        error=function(e){
          createAlert(session=session,anchorId = "errorAlert",title="Download Error",
                      content = "No data exist.",
                      style="warning")
        }
      )
    }
  )
  
  
  # Download data button event handle.
  output$drDownload <- downloadHandler(
    filename = function() {
      "DimensionReductionResult.csv"
    },
    content = function(file) {
      tryCatch(
        {
          if(is.null(rv$df_dr)){
            stop(safeError("No data exist."))
          }
          write.csv(rv$df_dr, file)
        },
        error=function(e){
          createAlert(session=session,anchorId = "errorAlert",title="Download Error",
                      content = "No data exist.",
                      style="warning")
        }
      )
    }
  )
  output$mianClusterDownload<-downloadHandler(
    filename = function() {
      "mainClusteringResult.csv"
    },
    content = function(file) {
      tryCatch(
        {
          if(is.null(rv$df_cluster)){
            stop(safeError("No data exist."))
          }
          write.csv(rv$df_cluster, file)
        },
        error=function(e){
          createAlert(session=session,anchorId = "errorAlert",title="Download Error",
                      content = "No data exist.",
                      style="warning")
        }
      )
    }
  )
  output$subClusterDownload<-downloadHandler(
    filename = function() {
      "subClusteringResult.csv"
    },
    content = function(file) {
      tryCatch(
        {
          if(is.null(rv$sub_df_cluster)){
            stop(safeError("No data exist."))
          }
          write.csv(rv$sub_df_cluster, file)
        },
        error=function(e){
          createAlert(session=session,anchorId = "errorAlert",title="Download Error",
                      content = "No data exist.",
                      style="warning")
        }
      )
    }
  )
  
  output$classificationDownload<-downloadHandler(
    filename = function() {
      "cellAnnotationResult.csv"
    },
    content = function(file) {
      tryCatch(
        {
          if(is.null(rv$df_classify)){
            stop(safeError("No data exist."))
          }
          write.csv(rv$df_classify, file)
        },
        error=function(e){
          createAlert(session=session,anchorId = "errorAlert",title="Download Error",
                      content = "No data exist.",
                      style="warning")
        }
      )
    }
  )
  
  output$networkDataDownload<-downloadHandler(
    filename = function() {
      "communicationResult.zip"
    },
    content = function(file) {
      tryCatch(
        {
          if(is.null(rv$networkData)){
            stop(safeError("No data exist."))
          }
          zip(file,"./outputData/t2sRes")
        },
        error=function(e){
          createAlert(session=session,anchorId = "errorAlert",title="Download Error",
                      content = "No data exist.",
                      style="warning")
        }
      )
    },
    contentType = "application/zip"
  )
  
  
  # Button event handle 
  observeEvent("",loadlastWorkingDir(rv))
  observeEvent("",  loadMarkerGene(rv,session))
  observeEvent("",  loadCellSelection(rv))
  observeEvent("",{initialization(rv)})
  observeEvent("",{loadLigRecDatabase(rv)})
  observeEvent("",loadDrugBankData(rv))
  observeEvent("",loadKeggInfo(rv))
  observeEvent("",loadTfTargetInteraction(rv))
  observeEvent("",loadSTRING(rv))
  
  observeEvent(input$networkDisplaySelect,{networkDisplayEvent(input)})
  observeEvent(input$targetDrugDisplaySelect,drugDisplayEvent(input))
  observeEvent(input$workingDirSet,workingDirSetEvent(defaultDir,input,rv,session))
  observeEvent(input$rnaDataFile,loadRnaFile(input,rv,session))
  observeEvent(input$rnaGroupFile,loadRnaGroupFile(input,rv,session))
  observeEvent(input$networkDataFile,loadCommunicationFile(input,rv,session))
  observeEvent(input$networkGroupFile,loadNetworkGroupFile(input,rv,session))
  observeEvent(input$rnaPreprocessing,preprocessingEvent(input,rv,session,use_rna_data = T,
                                                         imp=input$rnaImputation))
  observeEvent(input$networkPreprocessing,preprocessingEvent(input,rv,session,use_rna_data = F,
                                                             imp=input$networkImputation))
  observeEvent(input$aeButton,aeEvent(input,rv,session))
  observeEvent(input$drButton,tsneEvent(input,rv,session))
  observeEvent(input$continueWork,continueWork(input,rv,session))
  observeEvent(input$demoData,demoData(input,rv,session))
  observeEvent(input$preClusterButton,preClusterEvent(input,rv,session))
  observeEvent(input$mainClusterButton,mainClusterEvent(input,rv,session))
  observeEvent(input$subClusterButton,subClusterEvent(input,rv,session))
  observeEvent(input$doClassifyButton,classificationEvent(input,rv,session))
  observeEvent(input$loadExpressionDataButton,geneExpressionDataLoadingEvent(input,rv,session))
  observeEvent(input$addRowButton,addRowButtonEvent(rv))
  observeEvent(input$deleteRowButton,deleteRowButtonEvent(input,rv))
  observeEvent(input$geneTable_cell_edit, markerGeneEditEvent(input,rv))
  observeEvent(input$adButton,adButtonEvent(rv,session))
  observeEvent(input$cellSelection,cellSelectionEvent(rv,input))
  observeEvent(input$originalButton,originalButtonEvent(rv))
  observeEvent(input$geneExpressionExpend,geneExpressionExpendLoading())
  observeEvent(input$rnaDataExpend,rnaDataExpendLoading())
  observeEvent(input$communicationButton,communicationEvent(input,rv,session))
  observeEvent(input$networkPlotDownload1,{session$sendCustomMessage("downloadNetwork1","dir")})
  observeEvent(input$networkPlotDownload2,{session$sendCustomMessage("downloadNetwork2","dir")})
  observeEvent(input$activatedNetworkPlotDownload1,{session$sendCustomMessage("downloadActivatedNetwork1","dir")})
  observeEvent(input$activatedNetworkPlotDownload2,{session$sendCustomMessage("downloadActivatedNetwork2","dir")})
  observeEvent(input$upStreamPlotDownload1,{session$sendCustomMessage("downloadUpStreamNetwork1","dir")})
  observeEvent(input$upStreamPlotDownload2,{session$sendCustomMessage("downloadUpStreamNetwork2","dir")})
  observeEvent(input$upStreamPlotDownload3,{session$sendCustomMessage("downloadUpStreamNetwork3","dir")})
  observeEvent(input$upStreamPlotDownload4,{session$sendCustomMessage("downloadUpStreamNetwork4","dir")})
  observeEvent(input$upRegButton,upRegulateEvent(input,rv,session))
  observeEvent(input$pdacButton,pdacButtonEvent(rv,session))
  observeEvent(input$epButton,emt_proEvent(input,rv,session))
  observeEvent(input$goButton,GOEvent(input,rv,session))
  observeEvent(input$GOList,selectGONetwork(input,rv,session))
  observeEvent(input$addColButton,addColButtonEvent(input,rv))
  observeEvent(input$ligRecData,ligRecDatabaseSelect(input,rv))
  observeEvent(input$drugFileProcessing,processDrugFileEvent(defaultDir,input,rv,session))
  observeEvent(rv$display_marker_gene,loadCellSelection(rv))
  observeEvent(input$saveMarkerButton,saveMarkerGeneEvent(rv,session))
  
  
  # reactive function to show table dynamicly
  load_table<-reactive({loadTableRactive(input,rv)})
  load_drug_table1<-reactive({loadDrugTable1Ractive(rv)})
  load_drug_table2<-reactive({loadDrugTable2Ractive(rv)})
  load_GO_table<-reactive({loadGOTableRactive(input,rv)})
  load_drug_mapping_table1<-reactive({loadDrugMappingTable1(input,rv)})
  load_drug_mapping_table2<-reactive({loadDrugMappingTable2(input,rv)})
  
  #edit marker gene button event(would cause bug if put in another script )
  markerGeneEditEvent<-function(input,rv){
    info = input$geneTable_cell_edit
    str(info)
    i = info$row
    j = info$col
    v = info$value
    rv$display_marker_gene[i, j] <<- DT::coerceValue(v, rv$display_marker_gene[i, j])
  }
  
  
  #render Data table
  output$networkTable<-renderDT(
        rv$df_table,server=FALSE,width = 600, height = 700,
        selection = "single",
        extensions = c('FixedColumns',"Buttons","FixedHeader"),
        options = list(
          pageLength=10,
          dom = 'Bfrtip',
          scrollX = TRUE,
          fixedColumns = FALSE,
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
        )
  )
  output$geneTable<-renderDT(
        rv$display_marker_gene,server=FALSE,width = 600, height = 700,
        selection = "single",
        editable=TRUE,
        extensions = c('FixedColumns',"Buttons","FixedHeader"),
        options = list(
          pageLength=13,
          dom = 'Bfrtip',
          scrollX = TRUE,
          fixedColumns = FALSE,
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
        )
  )
  
  output$drugTable1<-renderDT(
    rv$drug_table1,server=FALSE,width = 1000, height = 800,
    selection = "single",
    editable=FALSE,
    extensions = c('FixedColumns',"Buttons","FixedHeader"),
    options = list(
      pageLength=15,
      dom = 'Bfrtip',
      scrollX = TRUE,
      fixedColumns = FALSE,
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
    )
  )
  
  output$drugTable2<-renderDT(
    rv$drug_table2,server=FALSE,width = 1000, height = 800,
    selection = "single",
    editable=FALSE,
    extensions = c('FixedColumns',"Buttons","FixedHeader"),
    options = list(
      pageLength=15,
      dom = 'Bfrtip',
      scrollX = TRUE,
      fixedColumns = FALSE,
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
    )
  )
  output$GOTable<-renderDT(
    rv$GO_table,server=FALSE,width = 1000, height = 600,
    selection = "single",
    editable=FALSE,
    extensions = c('FixedColumns',"Buttons","FixedHeader"),
    options = list(
      pageLength=10,
      dom = 'Bfrtip',
      scrollX = TRUE,
      fixedColumns = FALSE,
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
    )
  )
  
  output$drug_mapping_table1<-renderDT(
    rv$drug_mapping_table1,server=FALSE,width = 1000, height = 800,
    selection = "single",
    editable=FALSE,
    extensions = c('FixedColumns',"Buttons","FixedHeader"),
    options = list(
      pageLength=10,
      dom = 'Bfrtip',
      scrollX = TRUE,
      fixedColumns = FALSE,
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
    )
  )
  output$drug_mapping_table2<-renderDT(
    rv$drug_mapping_table2,server=FALSE,width = 1000, height = 800,
    selection = "single",
    editable=FALSE,
    extensions = c('FixedColumns',"Buttons","FixedHeader"),
    options = list(
      pageLength=10,
      dom = 'Bfrtip',
      scrollX = TRUE,
      fixedColumns = FALSE,
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
    )
  )
  #render plot
  output$bubblePlot <- renderD3({
    if(is.null(rv$bubbleData)){
      return(NULL)
    }
    bubble_json=bubbleDataToJson(rv$bubbleData)
    rv$drugNetworkPlot<-r2d3(
      data=bubble_json,d3_version = 4,
      script = "./www/bubble.js",
      width=900,
      height=798,viewer = "internal"
    )
    
  })
  
  
  output$visualization<-renderPlotly({
    if(is.null(rv$df_dr))
      return(NULL)
    #plot the data after dimension reduction(the first plot)
    plot_ly(rv$df_dr,x=~V1,y=~V2,type="scatter",mode="markers",opacity=0.8,width=600,height=600)%>%toWebGL()
  })
  
  
  #render pre cluster plot
  output$preClusterPlot<-renderPlotly({
    if(is.null(rv$dfDensity)){
      return(NULL)
    }
    max_index<-length(rv$dfDensity$reachdist)
    plot_ly(x=c(1:max_index),y=rv$dfDensity$reachdist,type="bar",width=700,height=450)%>%toWebGL()
  })
  
  #render main cluster plot
  output$mainClusterPlot<-renderPlotly({
    if(is.null(rv$df_cluster)){
      return (NULL)
    }
    plot_ly(rv$df_cluster,x=~V1,y=~V2,type="scatter",mode="markers"
            ,color=~factor(cluster),colors="Accent",width=700,height=600,
            text=~factor(cluster),hovertemplate = paste('<b>Cluster</b>: %{text}'))%>%
      layout(legend=c(itemdoubleclick="toggleothers"))%>%toWebGL()
  }
  )
  #render sub cluster plot
  output$subClusterPlot<-renderPlotly({
    if(is.null(rv$sub_df_cluster)){
      return (NULL)
    }
    plot_ly(rv$sub_df_cluster,x=~V1,y=~V2,type="scatter",mode="markers"
            ,color=~factor(cluster),colors="Accent",width=700,height=600,
            text=~factor(cluster),hovertemplate = paste('<b>Cluster</b>: %{text}'))%>%
      layout(legend=c(itemdoubleclick="toggleothers"))%>%toWebGL()
  }
  )
  
  #render classification plot1
  output$classifyPlot<-renderPlotly({
    if(is.null(rv$df_classify)){
      return (NULL)
    }
    plot_ly(rv$df_classify,x=~V1,y=~V2,type="scatter",mode="markers"
            ,color=~factor(type),colors="Accent",width=700,height=600,
            text=~factor(type),hovertemplate = paste('<b>Cluster</b>: %{text}'))%>%
      layout(legend=c(itemdoubleclick="toggleothers"))%>%toWebGL()
  }
  )
  #render classification plot2
  output$classifyClusterPlot<-renderPlotly({
    if(is.null(rv$df_classify)){
      return (NULL)
    }
    plot_ly(rv$df_classify,x=~V1,y=~V2,type="scatter",mode="markers"
            ,color=~factor(cluster_type),colors="Accent",width=700,height=600,
            text=~factor(cluster_type),hovertemplate = paste('<b>Cluster</b>: %{text}'))%>%
      layout(legend=c(itemdoubleclick="toggleothers"))%>%toWebGL()
  }
  )
  
  
  output$cellDistribution<-renderPlotly({
    if(is.null(rv$cell_count)){
      return (NULL)
    }
    if(rv$distribution_have_group==1){
      plot_ly(rv$cell_count,x=~group,y=~value,type="bar",color=~variable,alpha=0.9,colors="Dark2",width=1000,height=600)%>%
        layout(title="Cell Distribution for Each Cell Type",yaxis = list(title = 'percentage of Cell'),xaxis=list(title="Type of Cell"), barmode = 'stack')%>%toWebGL()
    }else{
      plot_ly(y=rv$cell_count$count,x=rv$cell_count$type,type="bar",width=1000,height=600)%>%
        layout(title="Cell Distribution for Each Cell Type",yaxis = list(title = 'percentage of Cell'),xaxis=list(title="Type of Cell"), barmode = 'stack')%>%toWebGL()
    }
  })
  
  output$emt_pro<-renderPlotly({
    if(is.null(rv$emt_pro_score)){
      return(NULL)
    }
    p<-plot_ly(rv$emt_pro_score,x=~PRO,y=~EMT,width=800,height=600,showlegend=FALSE)
    subplot(p%>%
              add_markers(color=~type)%>%
              layout(xaxis = list(range = c(0,0.8),title="PRO"),yaxis = list(range = c(0,0.8),title="EMT")),
            p%>%
              add_trace(type='histogram2dcontour',colorscale="Viridis",autocontour=FALSE,contours=list(start=0,end=200,size=10),showscale = T)%>%
              layout(xaxis = list(range = c(0,0.8)),yaxis = list(range = c(0,0.8))),
            titleY = TRUE,
            titleX= TRUE,
            shareY = TRUE)
    
  })
  
  output$networkPlot1 <- renderD3({
    if(!is.null(rv$downstream_network_json1)){

      r2d3(
        data=rv$downstream_network_json1,d3_version = 4,
        script = "./www/networkPlot.js",
        viewer="internal" 
      )
    }
  })
  
  output$networkPlot2 <- renderD3({
    if(!is.null(rv$downstream_network_json2)){
      r2d3(
        data=rv$downstream_network_json2,d3_version = 4,
        script = "./www/networkPlot2.js",
        viewer="internal" 
      )
    }
  })
  
  output$activatedNetworkPlot1 <- renderD3({
    if(!is.null(rv$activated_network_json1)){
      
      r2d3(
        data=rv$activated_network_json1,d3_version = 4,
        script = "./www/activatedNetworkPlot1.js",
        viewer="internal" 
      )
    }
  })
  
  output$activatedNetworkPlot2 <- renderD3({
    if(!is.null(rv$activated_network_json2)){
      r2d3(
        data=rv$activated_network_json2,d3_version = 4,
        script = "./www/activatedNetworkPlot2.js",
        viewer="internal" 
      )
    }
  })
  
  output$upRegLigandPlot <- renderD3({
    if(!is.null(rv$upRegLigand_network)){
      r2d3(
        data=rv$upRegLigand_network,d3_version = 4,
        script = "./www/cell_ligand_network.js",
        width=500,
        height=500,viewer = "internal"
      )
    }
  })
  output$upRegReceptorPlot <- renderD3({
    if(!is.null(rv$upRegReceptor_network)){
      r2d3(
        data=rv$upRegReceptor_network,d3_version = 4,
        script = "./www/cell_receptor_network.js",
        width=500,
        height=500,viewer = "internal"
      )
    }
  })
  
  output$GOPlot<- renderD3({
    if(rv$GOData!=0){
      r2d3(
        data=rv$GONetwork,d3_version = 4,
        script = "./www/GONetwork1.js",
        viewer = "internal"
      )
    }
  })
  
  output$drugClusterPlot1 <- renderD3({
    if(!is.null(rv$drug_json1)){
      r2d3(
        data=rv$drug_json1,d3_version = 4,
        script = "./www/drugClusterPlot1.js",
        viewer = "internal"
      )
    }
  })
  
  output$drugClusterPlot2 <- renderD3({
    if(!is.null(rv$drug_json2)){
      r2d3(
        data=rv$drug_json2,d3_version = 4,
        script = "./www/drugClusterPlot2.js",
        viewer = "internal"
      )
    }
  })
  
  output$up_expPlot <- renderD3({
    if(!is.null(rv$up_exp_network)){
      r2d3(
        data=rv$up_exp_network,d3_version = 4,
        script = "./www/up_expNetwork.js",
        viewer = "internal"
      )
    }
  })
  output$exp_upPlot <- renderD3({
    if(!is.null(rv$exp_up_network)){
      r2d3(
        data=rv$exp_up_network,d3_version = 4,
        script = "./www/exp_upNetwork.js",
        viewer = "internal"
      )
    }
  })
  output$combinePlot <- renderD3({
    if(!is.null(rv$combine_network)){
      r2d3(
        data=rv$combine_network,d3_version = 4,
        script = "./www/combineNetwork.js",
        viewer = "internal"
      )
    }
  })
  output$up_upPlot <- renderD3({
    if(!is.null(rv$up_up_network)){
      r2d3(
        data=rv$up_up_network,d3_version = 4,
        script = "./www/up_upNetwork.js",
        viewer = "internal"
      )
    }
  })
  output$drugMappingNetwork1 <- renderD3({
    if(!is.null(rv$drugNetworkJson1)){
      r2d3(
        data=rv$drugNetworkJson1,d3_version = 4,
        script = "./www/drugNetworkPlot1.js",
        viewer = "internal"
      )
    }
  })
  output$drugMappingNetwork2 <- renderD3({
    if(!is.null(rv$drugNetworkJson2)){
      r2d3(
        data=rv$drugNetworkJson2,d3_version = 4,
        script = "./www/drugNetworkPlot2.js",
        viewer = "internal"
      )
    }
  })
  
  output$targetDrugClusterPlot1 <- renderD3({
    if(!is.null(rv$targetDrug_json1)){
      r2d3(
        data=rv$targetDrug_json1,d3_version = 4,
        script = "./www/targetDrugClusterPlot1.js",
        viewer = "internal"
      )
    }
  })
  
  output$targetDrugClusterPlot2 <- renderD3({
    if(!is.null(rv$targetDrug_json2)){
      r2d3(
        data=rv$targetDrug_json2,d3_version = 4,
        script = "./www/targetDrugClusterPlot2.js",
        viewer = "internal"
      )
    }
  })
}


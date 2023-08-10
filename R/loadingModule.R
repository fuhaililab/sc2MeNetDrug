#initialization, load data and display
#default ui
initialization<-function(rv){
  shinyjs:::hide("dimensionReductionPart1")
  shinyjs:::hide("dimensionReductionPart2")
  shinyjs:::show("drNotification")
  shinyjs:::show("clusterNotification")
  shinyjs:::hide("clusterPart1")
  shinyjs:::hide("clusterPart2")
  shinyjs:::hide("clusterPart3")
  shinyjs:::hide("clusterPart4")
  shinyjs:::show("classificationNotification")
  shinyjs:::hide("classificationPart1")
  shinyjs:::hide("classificationPart2")
  shinyjs:::hide("classificationPart3")
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
  if(!dir.exists("./cache")){
    dir.create("./cache")
  }
}

#load rank matrix data for drug discovering
loadRankMatrix<-function(rv){
  if(file.exists("./data/rankMatrix92742a.RData")&&file.exists("./data/drug_mapping.RData")){
    load("./data/rankMatrix92742a.RData")
    rv$gSymZs<-gSymZs
    rv$rankMatrix<-rankMatrix
    load("./data/drug_mapping.RData")
    rv$drug_mapping<-drug_mapping
  }
}

#load last working directory
loadlastWorkingDir<-function(rv){
  if(file.exists("./cache/lastWorkingDir.RData")){
    load("./cache/lastWorkingDir.RData")
    rv$outputDir<-newDir
  }else{
    rv$outputDir<-NULL
  }
}

loadLigRecDatabase<-function(rv){
  progress <- shiny::Progress$new()
  progress$set(message = "Change database", value = 0)
  on.exit(progress$close())
  progress$set(0.5,detail="Processing")
  load("./data/baderLab_lr_database.RData")
  load("./data/DLRP_lr_database.RData")
  load("./data/nicheNet_lr_database.RData")
  rv$DLRP<-DLRP
  rv$nicheNet<-nicheNet
  rv$baderLab<-baderLab
  ligRecDatabase<-rbind(DLRP,nicheNet,baderLab)
  ligRecDatabase<-unique(ligRecDatabase)
  ligRecDatabase<-as.matrix(ligRecDatabase)
  rv$ligRecDatabase<-ligRecDatabase
  progress$set(1,detail="Finish")
  rm(DLRP,baderLab,nicheNet,ligRecDatabase)
  gc()
}
loadKeggInfo<-function(rv){
  load("./data/keggInfo.RData")
  rv$keggInfo<-keggInfo
  rm(keggInfo)
  gc()
}

loadTfTargetInteraction<-function(rv){
  load('./data/TfTargetInteraction_human.RData')
  rv$TfTargetInteraction<-TfTargetInteraction
  rm(TfTargetInteraction)
  gc()
}
loadSTRING<-function(rv){
  load("./data/STRING.RData")
  rv$STRING=STRING
  rm(STRING)
  gc()
}
#set new working directory
#after seting directory, reload whole application
workingDirSetEvent<-function(defaultDir,input,rv,session){
  tryCatch(
    {
      if(is.integer(input$workingDir)){
        stop(safeError("no directory"))
      }
      
      newDir<-parseDirPath(defaultDir, input$workingDir)
      #reload all application
      rv$df_table=NULL
      rv$adButtonUI={}
      rv$DieaseSpecificTitleUI={}
      rv$cellTypeSelect1UI={}
      rv$cellTypeSelect1UI={}
      rv$networkTableDisplayUI={}
      rv$drugTable1UI={}
      rv$drugTable2UI={}
      rv$networkData=0
      rv$drugNetworkData=0
      rv$drug_table1=NULL
      rv$drug_table2=NULL
      rv$rna_df=NULL
      rv$rna_have_group=0
      rv$rna_num_sample=0
      rv$rna_group_list=NULL
      rv$rna_gene_list=NULL
      rv$rna_type_list=NULL
      rv$network_df=NULL
      rv$network_have_group=0
      rv$network_group_list=NULL
      rv$network_gene_list=NULL
      rv$network_type_list=NULL
      rv$df_dr=NULL
      rv$encoder_result=NULL
      rv$df_cluster=NULL
      rv$sub_df_cluster=NULL
      rv$df_classify=NULL
      rv$rna_data_mean=NULL
      rv$network_data_mean=NULL
      rv$useRnaData=1
      rv$cell_type1=NULL
      rv$cell_type2=NULL
      rv$downstream_network_nodes1=NULL
      rv$downstream_network_edges1=NULL
      rv$downstream_network_json1=NULL
      rv$networkData1=0
      rv$networkTitle1UI={}
      rv$downstream_network_nodes2=NULL
      rv$downstream_network_edges2=NULL
      rv$downstream_network_json2=NULL
      rv$networkData2=0
      rv$networkTitle2UI={}
      rv$cell_type1_ligands=NULL
      rv$cell_type1_receptors=NULL
      rv$cell_type2_ligands=NULL
      rv$cell_type2_receptors=NULL
      rv$drug_network_data=NULL
      rv$drug_table1=NULL
      rv$drug_table2=NULL
      rv$upRegLigand_network=NULL
      rv$upRegReceptor_network=NULL
      rv$controlGroupUI={}
      rv$normalGroupUI={}
      rv$networkDownloadUI={}
      rv$cell_count=NULL
      rv$distribution_have_group=0
      rv$EPcellTypeSelectUI={}
      rv$EPGroupUI={}
      rv$drugData1=0
      rv$drugData2=0
      rv$drugTableTitle1UI={}
      rv$drugTableTitle2UI={}
      rv$drugClusterData1=0
      rv$drugClusterData2=0
      rv$drug_json1=NULL
      rv$drug_json2=NULL
      rv$networkCellSelectTitleUI={}
      rv$GONetwork1=NULL
      rv$GOData=0
      rv$GOTestGroupUI={}
      rv$GONormalGroupUI={}
      rv$GOCellTypeUI={}
      rv$GO_up_table=NULL
      rv$GO_dn_table=NULL
      rv$GO_table=NULL
      rv$GOTableDisplayUI={}
      rv$GONetworkSelectUI={}
      rv$GOList=NULL
      rv$netGO_up=NULL
      rv$netGO_dn=NULL
      rv$upRegTestGroupUI={}
      rv$upRegNormalGroupUI={}
      rv$up_exp_network=NULL
      rv$exp_up_network=NULL
      rv$combine_network=NULL
      rv$up_up_network=NULL
      rv$drugNetworkJson1=NULL
      rv$drugNetworkJson2=NULL
      rv$drug_mapping_table1=NULL
      rv$drug_mapping_table2=NULL
      rv$drugNetworkJson1=NULL
      rv$drugNetworkJson2=NULL
      rv$targetDrug_json1=NULL
      rv$targetDrug_json2=NULL
      rv$activatedNetworkTitle1UI={}
      rv$activatedNetworkTitle2UI={}
      rv$outputDir =newDir
      initialization(rv)
      loadMarkerGene(rv,session)
      save(newDir,file="./cache/lastWorkingDir.RData")
    },
    error=function(e){
      if(e$message=="no directory"){
        createAlert(session=session,anchorId = "errorAlert",title="Load Error",
                    content ="You haven't select directory",
                    style="warning")
      }else{
        createAlert(session=session,anchorId = "errorAlert",title="Load Error",
                    content =paste0("Something wrong when setting working directory:",e$message),
                    style="warning")
      }

    }
  )

}

#load marker gene data
loadMarkerGene<-function(rv,session){
  if((!is.null(rv$outputDir))&&file.exists(paste(rv$outputDir,"markerGeneTable.RData",sep="/"))){
    load(paste(rv$outputDir,"markerGeneTable.RData",sep="/"))
    marker_gene<-markerGeneTable
    rm(markrGeneTable)
  }else{
    marker_gene<-read.csv('./data/marker_gene_result.csv',stringsAsFactors=FALSE,check.names = FALSE)
  }
  num_col<-ncol(marker_gene)
  cell_list<-colnames(marker_gene)
  gene_symbol<-as.character(marker_gene[,1])
  #generate bubble plot data
  bubbleData<-c()
  for (i in 2:num_col){
    gene_list<-as.numeric(marker_gene[,i])
    cell_marker_gene_index<-sapply(gene_list,function(x){if(x==0){F}else{T}})
    cell_marker_gene<-paste(gene_symbol[cell_marker_gene_index],collapse=",")
    name_length<-nchar(cell_list[i])+sample(4:10,1)
    row_content<-c(cell_list[i],name_length,cell_marker_gene)
    if(i==2){
      bubbleData<-row_content
    }
    else{
      bubbleData<-rbind(bubbleData,row_content)
    }
  }
  
  colnames(bubbleData)<-c("cellType","length","markerGene")
  rv$bubbleData<-as.data.frame(bubbleData)
  colnames(marker_gene)[1]<-"Gene Symbol"
  rownames(marker_gene)<-c(1:nrow(marker_gene))
  marker_gene[,1]<-toupper(marker_gene[,1])
  rv$original_marker_gene<-marker_gene
  rv$select_marker_gene<-marker_gene
  rv$display_marker_gene<-marker_gene
  
  session$sendCustomMessage("refreshBubble","refresh")
  rm(bubbleData,marker_gene,gene_list,cell_marker_gene,name_length)
  gc()
}

loadDrugBankData<-function(rv){
  load("data/DrugBank_drug-target_new.RData")
  load("data/DrugBank_drugInformation.RData")
  rv$drugBankInteraction<-drugBank_drug_target
  rv$drugBankInformation<-drugBank_drugInformation
  rm(drugBank_drug_target,drugBank_drugInformation)
  gc()
}

loadCellSelection<-function(rv){
  cellList<-colnames(rv$display_marker_gene)[-1]
  rv$cellSelectionUI<-{
    tags$div(id="supportBar",
             tags$div(id="cellSelectionFlow",
                      checkboxGroupInput(inputId="cellSelection",label=NULL,
                                         choices = cellList,choiceValues=cellList,width="100%"))   
    )
  }
  rm(cellList)
  gc()
}


#load raw seq-RNA data uploaded by user
loadRnaFile<-function(input,rv,session){
      progress <- shiny::Progress$new()
      progress$set(message = "Load file", value = 0)
      progress$set(0.1, detail ="read data")
      on.exit(progress$close())
      tryCatch(
        {
        suffix<-strsplit(input$rnaDataFile$datapath,split="\\.")[[1]][2]
        if(tolower(suffix)=="csv"){
          data <- read.csv(input$rnaDataFile$datapath,
                           header = input$rnaHeader,
                           sep = input$rnaSep,
                           quote = input$rnaQuote,row.names = 1,check.names = FALSE,as.is=FALSE)
        }else if(tolower(suffix)=="rds"){
          data<-readRDS(input$rnaDataFile$datapath)
        }else{
          stop(safeError("invalid data"))
        }
         rv$rna_num_sample<-ncol(data)
         rv$rna_gene_list<-as.character(rownames(data))
         rv$rna_df<-data
         progress$set(1, detail ="finish")
         rm(data)
         gc()
        },
        error = function(e) {
          #    return a safeError if a parsing error occurs
          if(e$message=="invalid data"){
            createAlert(session=session,anchorId = "errorAlert",title="Upload Error",
                        content = "Only accept .csv or .rds file.",
                        style="danger")
          }else{
            createAlert(session=session,anchorId = "errorAlert",title="Upload Error",
                        content = paste0("The file you uploaded is invalid, Please read upload instruction carefully.",e),
                        style="danger")
          }


        },
        finally={
          progress$inc(1, detail ="transposition finish")
          on.exit(progress$close())
        }
      )
}


#load gourp information correspond to raw seq-RNA data
loadRnaGroupFile<-function(input,rv,session){
  tryCatch(
    {
      if(is.null(rv$rna_df)){
        stop(safeError("No upstream read count data uploaded."))
      }
      suffix<-strsplit(input$rnaGroupFile$datapath,split="\\.")[[1]][2]
      if(tolower(suffix)=="csv"){
        df <- read.csv(input$rnaGroupFile$datapath,
                       header = input$rnaHeader,
                       sep = input$rnaSep,
                       quote = input$rnaQuote,check.names = FALSE,as.is=FALSE)
      }else if(tolower(suffix)=="rds"){
        df<-readRDS(input$rnaGroupFile$datapath)
      }else{
        stop(safeError("invalid data"))
      }
        group_list<-as.character(df[,1])
        unique_group<-unique(group_list)
        if(length(unique_group)<=1||length(group_list)!=rv$rna_num_sample){
          stop(safeError("Group information not correct."))
        }
        rv$rna_group_list<-group_list
        save(group_list,file=paste0(rv$outputDir,"/rnaGroupInformation.RData"))
        rv$rna_have_group=1
        rm(group_list,unique_group,df,suffix)
        gc()
    },
    error = function(e) {
      #    return a safeError if a parsing error occurs
      if(e$message=="no RNA data"){
        createAlert(session=session,anchorId = "errorAlert",title="Upload Error",
                    content = "Please unpload seq-RNA data first!.",
                    style="danger")
      }else if(e$message=="invalid data"){
        createAlert(session=session,anchorId = "errorAlert",title="Upload Error",
                    content = "Only accept .csv or .rds file.",
                    style="danger")
      }else{
        createAlert(session=session,anchorId = "errorAlert",title="Upload Error",
                    content ="The file you uploaded is invalid, Please read upload instruction carefully.",
                    style="danger")
      }
    }
  )
}


#load gene expression analytics data uploaded by user
loadCommunicationFile<-function(input,rv,session){
  progress <- shiny::Progress$new()
  progress$set(message = "Data processing", value = 0)
  on.exit(progress$close())
  progress$set(0.1, detail ="read data")
  tryCatch(
    {
      suffix<-strsplit(input$networkDataFile$datapath,split="\\.")[[1]][2]
      if(tolower(suffix)=="csv"){
        network_df <- read.csv(input$networkDataFile$datapath,
                               header = input$networkHeader,
                               sep = input$networkSep,
                               quote = input$networkQuote,row.names = 1,check.names = FALSE,as.is=FALSE)
      }else if(tolower(suffix)=="rds"){
        network_df<-readRDS(input$networkDataFile$datapath)
      }else{
        stop(safeError("invalid data"))
      }
      rv$network_gene_list<-as.character(rownames(network_df))
      rv$network_num_sample<-ncol(network_df)
      rv$network_df<-network_df
      rm(suffix,network_df)
      progress$inc(1, detail ="Finish")
    },
    error = function(e) {
      #    return a safeError if a parsing error occurs
      if(e$message=="invalid data"){
        createAlert(session=session,anchorId = "errorAlert",title="Upload Error",
                    content = "Only accept .csv or .rds file.",
                    style="danger")
      }else{
        createAlert(session=session,anchorId = "errorAlert",title="Upload Error",
                    content = "The file you uploaded is invalid, Please read upload instruction carefully.",
                    style="danger")
      }

    },
    finally={
      progress$inc(1, detail ="Finish")
      on.exit(progress$close())
    }
  )
}

#load group information corresponding to gene expression data
loadNetworkGroupFile<-function(input,rv,session){
  tryCatch(
    {
      if(is.null(rv$network_df)){
        stop(safeError("No downstream read count data uploaded"))
      }
      suffix<-strsplit(input$networkGroupFile$datapath,split="\\.")[[1]][2]
      if(tolower(suffix)=="csv"){
        df <- read.csv(input$networkGroupFile$datapath,
                       header = input$networkHeader,
                       sep = input$networkSep,
                       quote = input$networkQuote,check.names = FALSE,as.is=FALSE)
      }else if(tolower(suffix)=="rds"){
        df<-readRDS(input$networkGroupFile$datapath)
      }else{
        stop(safeError("invalid data"))
      }
      
      num_row<-nrow(df)
      if(num_row!=rv$network_num_sample){
        stop(safeError("Information not correct."))
      }
      num_col<-ncol(df)
      if(num_col==1){
        type_list<-as.character(df[,1])
        rv$network_type_list<-type_list
        save(type_list,file=paste0(rv$outputDir,"/networkTypeInformation.RData"))
        rv$network_have_group=0
        rm(type_list)
      }else if(num_col==2){
        type_list<-as.character(df[,1])
        rv$network_type_list<-type_list
        save(type_list,file=paste0(rv$outputDir,"/networkTypeInformation.RData"))
        group_list<-as.character(df[,2])
        unique_group<-unique(group_list)
        if(length(unique_group)<=1||length(group_list)!=rv$network_num_sample){
          stop(safeError("Group information not correct."))
        }
        rv$network_group_list<-group_list
        save(group_list,file=paste0(rv$outputDir,"/networkGroupInformation.RData"))
        rv$network_have_group=1
        rm(type_list,group_list,unique_group)
      }else{
        stop(safeError("Group information not correct."))
      }
      rm(suffix,df)
      gc()
    },
    error = function(e) {
      #    return a safeError if a parsing error occurs
      if(e=="no network data"){
        createAlert(session=session,anchorId = "errorAlert",title="Upload Error",
                    content = "Please unpload gene communication data first!.",
                    style="danger")
      }else if(e$message=="invalid data"){
        createAlert(session=session,anchorId = "errorAlert",title="Upload Error",
                    content = "Only accept .csv or .rds file.",
                    style="danger")
      }else{
        createAlert(session=session,anchorId = "errorAlert",title="Upload Error",
                    content = "The file you uploaded is invalid, Please read upload instruction carefully.",
                    style="danger")
      }
    }
  )
}

#load up-regulated Ligands and Receptors group selection
loadUpReg<-function(rv){
  
  if(rv$useRnaData==0){
    if(rv$network_have_group==1){
      group_list<-unique(rv$network_group_list)
      shinyjs::show("upRegNormalTitle")
      shinyjs::show("upRegTestTitle")
      #test group
      rv$upRegTestGroupUI<-{
        tags$div(class="GroupSupportBar",
                 tags$div(class="groupSelectionFlow",
                          checkboxGroupInput(inputId="upRegTestGroup",label=NULL,
                                             choices = group_list,choiceValues=group_list,width="100%"))   
        )
      }
      #normal group
      rv$upRegNormalGroupUI<-{
        tags$div(class="GroupSupportBar",
                 tags$div(class="groupSelectionFlow",
                          checkboxGroupInput(inputId="upRegNormalGroup",label=NULL,
                                             choices =group_list,choiceValues=group_list,width="100%"))   
        )
      }
      rm(group_list)
    }else{
      shinyjs::hide("upRegNormalTitle")
      shinyjs::hide("upRegTestTitle")
      #test group
      rv$upRegTestGroupUI<-{
      }
      #normal group
      rv$upRegNormalGroupUI<-{
      }
    }
  }else{
    if(rv$rna_have_group==1){
      shinyjs::show("upRegNormalTitle")
      shinyjs::show("upRegTestTitle")
      group_list<-unique(rv$rna_group_list)
      #test group
      rv$upRegTestGroupUI<-{
        tags$div(class="GroupSupportBar",
                 tags$div(class="groupSelectionFlow",
                          checkboxGroupInput(inputId="upRegTestGroup",label=NULL,
                                             choices = group_list,choiceValues=group_list,width="100%"))   
        )
      }
      #normal group
      rv$upRegNormalGroupUI<-{
        tags$div(class="GroupSupportBar",
                 tags$div(class="groupSelectionFlow",
                          checkboxGroupInput(inputId="upRegNormalGroup",label=NULL,
                                             choices =group_list,choiceValues=group_list,width="100%"))   
        )
      }
      rm(group_list)
    }else{
      shinyjs::hide("upRegNormalTitle")
      shinyjs::hide("upRegTestTitle")
      #test group
      rv$upRegTestGroupUI<-{
        
      }
      #normal group
      rv$upRegNormalGroupUI<-{
      }
    }
  }
  gc()
}


#load cell-cell communication UI corresponding to
#seq-RNA data
loadRnaDataSelection<-function(rv){
  type_list<-unique(rv$rna_type_list)
  #cell type 1
  rv$cellTypeSelect1UI<-{
    tags$div(class="GroupSupportBar",
             tags$div(class="groupSelectionFlow",
                      checkboxGroupInput(inputId="cellTypeSelect1",label=NULL,
                                         choices = type_list,choiceValues=type_list,width="100%"))   
    )
  }
  #cell type 2 
  rv$cellTypeSelect2UI<-{
    tags$h4("Stroma Cell")
    tags$div(class="GroupSupportBar",
             tags$div(class="groupSelectionFlow",
                      checkboxGroupInput(inputId="cellTypeSelect2",label=NULL,
                                         choices = type_list,choiceValues=type_list,width="100%"))   
    )
  }
  if(rv$rna_have_group==1){
    shinyjs::show("networkNormalTitle")
    shinyjs::show("networkTestTitle")
    group_list<-unique(rv$rna_group_list)
    #control group
    rv$controlGroupUI<-{
      tags$div(class="GroupSupportBar",
               tags$div(class="groupSelectionFlow",
                        checkboxGroupInput(inputId="controlGroup",label=NULL,
                                           choices = group_list,choiceValues=group_list,width="100%"))   
      )
    }
    #normal group
    rv$normalGroupUI<-{
      tags$div(class="GroupSupportBar",
               tags$div(class="groupSelectionFlow",
                        checkboxGroupInput(inputId="normalGroup",label= NULL,
                                           choices = group_list,choiceValues=group_list,width="100%"))   
      )
    }
  }else{
    shinyjs::hide("networkNormalTitle")
    shinyjs::hide("networkTestTitle")
    #control group
    rv$controlGroupUI<-{
    }
    #normal group
    rv$normalGroupUI<-{
    }
  }
}

#load cell-cell communication UI corresponding to
#gene expression data
loadNetworkDataSelection<-function(rv){
  type_list<-unique(rv$network_type_list)
  #cell type 1
  rv$cellTypeSelect1UI<-{
    tags$div(class="GroupSupportBar",
             tags$div(class="groupSelectionFlow",
                      checkboxGroupInput(inputId="cellTypeSelect1",label=NULL,
                                         choices = type_list,choiceValues=type_list,width="100%"))   
    )
  }
  #cell type 2
  rv$cellTypeSelect2UI<-{
    tags$div(class="GroupSupportBar",
             tags$div(class="groupSelectionFlow",
                      checkboxGroupInput(inputId="cellTypeSelect2",label=NULL,
                                         choices = type_list,choiceValues=type_list,width="100%"))   
    )
  }
  if(rv$network_have_group==1){
    group_list<-unique(rv$network_group_list)
    shinyjs::show("networkNormalTitle")
    shinyjs::show("networkTestTitle")
    #control group
    rv$controlGroupUI<-{
      tags$div(class="GroupSupportBar",
               tags$div(class="groupSelectionFlow",
                        checkboxGroupInput(inputId="controlGroup",label=NULL,
                                           choices = group_list,choiceValues=group_list,width="100%"))   
      )
    }
    #normal group
    rv$normalGroupUI<-{
      tags$div(class="GroupSupportBar",
               tags$div(class="groupSelectionFlow",
                        checkboxGroupInput(inputId="normalGroup",label=NULL,
                                           choices =group_list,choiceValues=group_list,width="100%"))   
      )
    }
  }else{
    shinyjs::hide("networkNormalTitle")
    shinyjs::hide("networkTestTitle")
    #control group
    rv$controlGroupUI<-{
    }
    #normal group
    rv$normalGroupUI<-{
    }
  }
}

ligRecDatabaseSelect<-function(input,rv){
  ligRecDatabase<-c()
  if(length(input$ligRecData)==0||length(input$ligRecData)==3){
    ligRecDatabase<-rbind(rv$DLRP,rv$nicheNet,rv$baderLab)
  }else{
    for(i in 1:length(input$ligRecData)){
      database<-input$ligRecData[i]
      if(database=="DLRP"){
        if(i==1){
          ligRecDatabase<-rv$DLRP
        }else{
          ligRecDatabase<-rbind(ligRecDatabase,rv$DLRP)
        }
      }else if(database=="nicheNet"){
        if(i==1){
          ligRecDatabase<-rv$nicheNet
        }else{
          ligRecDatabase<-rbind(ligRecDatabase,rv$nicheNet)
        }
      }else{
        if(i==1){
          ligRecDatabase<-rv$baderLab
        }else{
          ligRecDatabase<-rbind(ligRecDatabase,rv$baderLab)
        }
      }
    }
  }
  ligRecDatabase<-unique(ligRecDatabase)
  ligRecDatabase<-as.matrix(ligRecDatabase)
  rv$ligRecDatabase<-ligRecDatabase
  rm(ligRecDatabase,database)
  gc()
}

loadGO<-function(rv){
  
  if(rv$useRnaData==0){
    type_list<-unique(rv$network_type_list)
    #cell type
    rv$GOCellTypeUI<-{
      tags$div(class="GroupSupportBar",
               tags$div(class="groupSelectionFlow",
                        checkboxGroupInput(inputId="GOCellType",label=NULL,
                                           choices = type_list,choiceValues=type_list,width="100%"))   
      )
    }
    if(rv$network_have_group==1){
      shinyjs::show("GONormalTitle")
      shinyjs::show("GOTestTitle")
      group_list<-unique(rv$network_group_list)
      #control group
      rv$GOTestGroupUI<-{
        tags$div(class="GroupSupportBar",
                 tags$div(class="groupSelectionFlow",
                          checkboxGroupInput(inputId="GOTestGroup",label=NULL,
                                             choices = group_list,choiceValues=group_list,width="100%"))   
        )
      }
      #normal group
      rv$GONormalGroupUI<-{
        tags$div(class="GroupSupportBar",
                 tags$div(class="groupSelectionFlow",
                          checkboxGroupInput(inputId="GONormalGroup",label=NULL,
                                             choices =group_list,choiceValues=group_list,width="100%"))   
        )
      }
    }else{
      shinyjs::hide("GONormalTitle")
      shinyjs::hide("GOTestTitle")
      #control group
      rv$GOTestGroupUI<-{
      }
      #normal group
      rv$GONormalGroupUI<-{
      }
    }
  }else{
    type_list<-unique(rv$rna_type_list)
    #cell type
    rv$GOCellTypeUI<-{
      tags$div(class="GroupSupportBar",
               tags$div(class="groupSelectionFlow",
                        checkboxGroupInput(inputId="GOCellType",label=NULL,
                                           choices = type_list,choiceValues=type_list,width="100%"))   
      )
    }
    if(rv$rna_have_group==1){
      shinyjs::show("GONormalTitle")
      shinyjs::show("GOTestTitle")
      group_list<-unique(rv$rna_group_list)
      #control group
      rv$GOTestGroupUI<-{
        tags$div(class="GroupSupportBar",
                 tags$div(class="groupSelectionFlow",
                          checkboxGroupInput(inputId="GOTestGroup",label=NULL,
                                             choices = group_list,choiceValues=group_list,width="100%"))   
        )
      }
      #normal group
      rv$GONormalGroupUI<-{
        tags$div(class="GroupSupportBar",
                 tags$div(class="groupSelectionFlow",
                          checkboxGroupInput(inputId="GONormalGroup",label=NULL,
                                             choices =group_list,choiceValues=group_list,width="100%"))   
        )
      }
    }else{
      shinyjs::hide("GONormalTitle")
      shinyjs::hide("GOTestTitle")
      #control group
      rv$GOTestGroupUI<-{
      }
      #normal group
      rv$GONormalGroupUI<-{
      }
  }
}
}

loadGONetworkSelect<-function(rv,GO_list){
  rv$GONetworkSelectUI<-{
    selectInput("GOList",label=h4("Select GO"),
                choices=GO_list,selected=GO_list[1])
  }
}

loadDimensionReduction<-function(){
  shinyjs:::show("dimensionReductionPart1")
  shinyjs:::show("dimensionReductionPart2")
  shinyjs:::hide("drNotification")
}


loadPreCluster<-function(){
  shinyjs:::hide("clusterNotification")
  shinyjs:::show("clusterPart1")
  shinyjs:::show("clusterPart2")
}

loadMainCluster<-function(){
  shinyjs:::show("clusterPart3")
}

loadSubCluster<-function(){
  shinyjs:::show("clusterPart4")
}

loadClassification<-function(){
  shinyjs:::hide("classificationNotification")
  shinyjs:::show("classificationPart1")
  shinyjs:::show("classificationPart2")
  shinyjs:::show("classificationPart3")
}
loadGeneExpression1<-function(){
  shinyjs::hide("geneExpressionNotification")
  shinyjs::show("geneExpressionPart1")
}
loadGeneExpression2<-function(){
  shinyjs::show("geneExpressionPart2")
  shinyjs::show("geneExpressionPart3")
  shinyjs::show("geneExpressionPart4")
  shinyjs::show("geneExpressionPart5")
  shinyjs::show("geneExpressionPart6")
  shinyjs::show("geneExpressionPart7")
}

#hide EMT-PRO UI
hideEPSelection<-function(rv){
  rv$EPcellTypeSelectUI<-{}
  rv$EPGroupUI<-{}
  shinyjs::hide("epButton")
}

#load EMT-PRO UI
loadEPSelection<-function(rv){
  if(rv$useRnaData==1){
    type_list<-unique(rv$rna_type_list)
    rv$EPcellTypeSelectUI<-{ 
      tags$div(class="GroupSupportBar",
      tags$div(class="groupSelectionFlow",
      checkboxGroupInput(inputId="epType",label=NULL,
                       choices = type_list,choiceValues=type_list,width="100%")))
      }

    if(rv$rna_have_group==1){
      shinyjs::show("EPGroupTitle")
      group_list<-unique(rv$rna_group_list)
      rv$EPGroupUI<-{
        tags$div(class="GroupSupportBar",
        tags$div(class="groupSelectionFlow",
        checkboxGroupInput("epGroup",label=NULL,
                           choices=group_list,choiceValues=group_list,width="100%")))
      }
    }else{
      shinyjs::hide("EPGroupTitle")
    }
  }else{
    type_list<-unique(rv$network_type_list)
    rv$EPcellTypeSelectUI<-{      
      tags$div(class="GroupSupportBar",
      tags$div(class="groupSelectionFlow",
      checkboxGroupInput(inputId="epType",label=NULL,
                         choices = type_list,choiceValues=type_list,width="100%")))
      }
    if(rv$network_have_group==1){
      shinyjs::show("EPGroupTitle")
      group_list<-unique(rv$network_group_list)
      rv$EPGroupUI<-{
        tags$div(class="GroupSupportBar",
        tags$div(class="groupSelectionFlow",
        checkboxGroupInput("epGroup",label=NULL,
                           choices=group_list,choiceValues=group_list,width="100%")))
      }
    }else{
      shinyjs::hide("EPGroupTitle")
    }
  }
  shinyjs::show("epButton")
}

loadCommunication<-function(){
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

# loaddrug<-function(){
#   shinyjs:::hide("drugNotification")
#   shinyjs:::show("drugPart1")
#   shinyjs:::show("drugPart2")
#   shinyjs:::show("drugPart3")
# }
loadNetwork1<-function(rv,cell_type1,cell_type2){
  rv$networkTitle1UI={
    tags$h3(paste("Inter-Cell Signaling Communication from",paste(cell_type1,collapse = "+"),"to",paste(cell_type2,collapse = "+"),sep=" "))
  }
}
loadNetwork2<-function(rv,cell_type1,cell_type2){
  rv$networkTitle2UI={
    tags$h3(paste("Inter-Cell Signaling Communication from",paste(cell_type2,collapse = "+"),"to",paste(cell_type1,collapse = "+"),sep=" "))
  }
}

loadActivatedNetwork1<-function(rv,cell_type1){
  rv$activatedNetworkTitle1UI={
    tags$h3(paste0("Activated Signaling Pathways for ",paste(cell_type1,collapse = "+")))
  }
}

loadActivatedNetwork2<-function(rv,cell_type2){

  rv$activatedNetworkTitle2UI={
    tags$h3(paste0("Activated Signaling Pathways for ",paste(cell_type2,collapse = "+")))
  }
}

loadTargetDrugNetwork1<-function(rv,cell_type1,cell_type2){
  rv$targetDrugNetworkTitle1UI={
    tags$h3(paste("Targets Drug Discovering Result for Downstream Network from",paste(rv$cell_type1,collapse = "+"),"to" 
                  ,paste(rv$cell_type2,collapse = "+"), sep=" "))  }
}

loadTargetDrugNetwork2<-function(rv,cell_type1,cell_type2){
  rv$targetDrugNetworkTitle2UI={
    tags$h3(paste("Targets Drug Discovering Result for Downstream Network from",paste(rv$cell_type2,collapse = "+"),"to" 
                  ,paste(rv$cell_type1,collapse = "+"), sep=" "))  }
}



#load gene expression data uploading instruction expend 
geneExpressionExpendLoading<-function(){
  showModal(modalDialog(
    tags$h4(id="networkInstructionTitle","In order application could run properly,Please read the 
                                 upload instruction carefully:"),
    tags$ol(class="content",
            tags$li("You can upload csv or RDS file."),
            tags$li("Read count data should be a data frame, with every row represent gene and every column represent
                                           cell sample."),
            tags$li("the row name of read count data should be unique gene symbol.If your data have column name, please 
                                          check the 'header' checkbox."
            ),
            tags$li("Group or design file is an optional choice, but we recommend you to upload it in order to get
                                           reasonable analysis result. Data you upload should be data frame with only one column. The length
                                           of row should be equal to the number of cell in read count data. Please don't include rowname in you data."),
    ),
    title="Gene Expression Analytics Data Upload Instruction",
    easyClose = TRUE
  ))
}

#load raw-seq RNA data uploading instruction expend 
rnaDataExpendLoading<-function(){
  showModal(modalDialog(
    tags$h4(id="rnaInstructionTitle","In order application could run properly,Please read the 
                                 upload instruction carefully:"),
    tags$ol(class="content",
            tags$li("You can upload csv or RDS file."),
            tags$li("Read count data should be a data frame, with every row represent gene and every column represent
                                           cell sample."),
            tags$li("the row name of read count data should be unique gene symbol.If your data have column name, please 
                                          check the 'header' checkbox."
            ),
            tags$li("Group or design file in downstream analysis part is required. the data show have one or two columns,
                                  with first column indicate the cell type for each cell. Second column is optional, which is used to specify group or design for 
                                  each cell.The length of row should be equal to the number of cell in read count data. Please don't include rowname in you data.
                                    ")
    ),
    title="RNA Data Upload Instruction",
    easyClose = TRUE
  ))
}
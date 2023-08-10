#continue work button event
#check all the data user have in current working directory
#load corresponding data
#raw seq-RNA data analytics and gene expression part 
#would be loaded independently
continueWork<-function(input,rv,session){
  #raw data reload
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Data reloading", value = 0)
  #gene expression data
  result1<-loadRawAnalyticsData(rv,progress)
  result2<-loadNetworkData(rv,progress)
  if(result1==-1&&result2==-1){
    createAlert(session=session,anchorId = "errorAlert",title="Warning",
                content = "You do not have any previous work.",
                style="warning")
    progress$inc(1, detail ="Finish")
  }else if(result1==1||result2==1){
    if(file.exists(paste0(rv$outputDir,"/lastAnalysisDataIndex.RData"))){
      load(paste0(rv$outputDir,"/lastAnalysisDataIndex.RData"))
      if(useRnaData==1){
        rv$useRnaData=1
        loadRnaDataSelection(rv)
      }else{
        rv$useRnaData=0
        loadNetworkDataSelection(rv)
      }
      load(paste0(rv$outputDir,"/cellDistributionPlot.RData"))
      rv$cell_count<-cell_count
      rv$distribution_have_group<-distribution_have_group
      loadGeneExpression2()
      loadEPSelection(rv)
      rm(cell_count,distribution_have_group)
      gc()
      if(file.exists("./data/rankMatrix92742a.RData")&&file.exists("./data/drug_mapping.RData")){
        loadCommunication()
      }
      loadGO(rv)
      loadUpReg(rv)
    }else{
      progress$inc(1, detail ="Finish")
      return()
    }

    progress$set(0.8,detail="reload up-regulate gene data")
    if(file.exists(paste0(rv$outputDir,"/lastUpRegIndexing.RData"))){
      load(paste0(rv$outputDir,"/lastUpRegIndexing.RData"))
      if(file.exists(paste0(groupDir,"/upRegNetwork.RData"))&&file.exists(paste0(groupDir,"/upStreamNetwork.RData"))){
        load(paste0(groupDir,"/upRegNetwork.RData"))
        load(paste0(groupDir,"/upStreamNetwork.RData"))
        # rv$foldChangeandPValue<-foldChangeandPValue
        rv$upRegLigand_network<-upReg_network[[1]][[1]]
        rv$upRegReceptor_network<-upReg_network[[2]][[1]]
        if(rv$upRegGenerated==1){
          session$sendCustomMessage("refreshUpRegReceptor", "refresh")
          session$sendCustomMessage("refreshUpRegLigand", "refresh")
        }
        rv$upRegGenerated=1
        
        rv$up_exp_network<-upStreamData[[1]][[1]]
        rv$exp_up_network<-upStreamData[[2]][[1]]
        rv$combine_network<-upStreamData[[3]][[1]]
        rv$up_up_network<-upStreamData[[4]][[1]]
        if(rv$upStreamGenerated==1)
        {
          session$sendCustomMessage("refreshUpExp", "refresh")
          session$sendCustomMessage("refreshExpUp", "refresh")
          session$sendCustomMessage("refreshCombine", "refresh")
          session$sendCustomMessage("refreshUpUp", "refresh")
        }      
        rv$upStreamGenerated=1
        rm(upReg_network,upStreamData)
        gc()
      }
    }

    
    progress$set(0.82,detail="reload PRO-EMT data")
    if(file.exists(paste0(rv$outputDir,"/emt_pro_score.RData"))){
      load(paste0(rv$outputDir,"/emt_pro_score.RData"))
      rv$emt_pro_score<-emt_pro_score
      rm(emt_pro_score)
      gc()
    }
    progress$set(0.83,detail="reload GO analytics data")
    if(file.exists(paste0(rv$outputDir,"/GO_result.RData"))){
      load(paste0(rv$outputDir,"/GO_result.RData"))
      rv$GO_up_table<-GO_result[[1]]
      rv$GO_dn_table<-GO_result[[2]]
      rv$netGO_up<-GO_result[[3]]
      rv$netGO_dn<-GO_result[[4]]
      GO_list<-c(as.character(rv$GO_up_table[,1]),as.character(rv$GO_dn_table[,1]))
      GO_name<-c(as.character(paste(rv$GO_up_table[,1],rv$GO_up_table[,2],sep="  ")),as.character(paste(rv$GO_dn_table[,1],rv$GO_dn_table[,2],sep="  ")))
      names(GO_list)<-GO_name
      rv$GO_list<-GO_list
      loadGONetworkSelect(rv,rv$GO_list)
      selectGONetwork(input,rv,session)
      # GO<-rv$GO_list[1]
      # GO_name<-names(rv$GO_list)[1]
      # rv$GONetwork<-GONetworkGenerating(rv$netGO_up,rv$netGO_dn,GO,GO_name)
      # if(rv$GOGenerated==1){
      #   session$sendCustomMessage("refreshGONetwork1","refresh")
      # }
      # rv$GOGenerated=1
      # rv$GOData=1
      rm(GO_result,GO_list,GO_name,GO)
      gc()
    }

    progress$set(0.84,detail="reload Network data")
    
    haveNetwork1=0
    haveNetwork2=0
    haveActivatedNetwork1=0
    haveActivatedNetwork2=0
    networkComputed=0
    if(file.exists(paste0(rv$outputDir,"/lastNetworkIndexing.RData"))){
      load(paste0(rv$outputDir,"/lastNetworkIndexing.RData"))
      networkComputed=1
      rv$cell_type1=cell_type1
      rv$cell_type2=cell_type2
    }
    
    if(networkComputed==1){
      fcDir<-paste0(rv$outputDir,"/cellCommunication")
      networkDir<-paste0("/",paste(paste(rv$cell_type1,collapse = "+"),paste(rv$cell_type2,collapse = "+"),sep="-"))
      netDir<-paste0(fcDir,networkDir)
      if(file.exists(paste0(netDir,"/genesInformation.RData"))){
        load(paste0(netDir,"/ligRecInformation.RData"))
        rv$cell_type1_ligands<-ligRecInformation[[1]]
        rv$cell_type1_receptors<-ligRecInformation[[2]]
        rv$cell_type2_ligands<-ligRecInformation[[3]]
        rv$cell_type2_receptors<-ligRecInformation[[4]]
        rm(ligRecInformation)
      }

      if(file.exists(paste0(netDir,"/activatedNetworkData1.RData"))){
        load(paste0(netDir,"/activatedNetworkData1.RData"))
        if(!is.null(activated_network_nodes1)){
          haveActivatedNetwork1=1
        }
      }

      if(file.exists(paste0(netDir,"/activatedNetworkData2.RData"))){
        load(paste0(netDir,"/activatedNetworkData2.RData"))
        if(!is.null(activated_network_nodes2)){
          haveActivatedNetwork2=1
        }
      }

      
      network1Dir<-paste0(netDir,"/",paste(paste(rv$cell_type1,collapse = "+"),paste(rv$cell_type2,collapse = "+"),sep="_"))
      
      if(file.exists(paste0(network1Dir,"/downstreamNetworkData.RData"))){
        load(paste0(network1Dir,"/downstreamNetworkData.RData"))
        if(!is.null(downstream_network_nodes1)){
          haveNetwork1=1
        }
      }
      
      network2Dir<-paste0(netDir,"/",paste(paste(rv$cell_type2,collapse = "+"),paste(rv$cell_type1,collapse = "+"),sep="_"))
      if(file.exists(paste0(network2Dir,"/downstreamNetworkData.RData"))){
        load(paste0(network2Dir,"/downstreamNetworkData.RData"))
        if(!is.null(downstream_network_nodes2)){
          haveNetwork2=1
        }
      }
      
      if(haveActivatedNetwork1==1){
        rv$activated_network_nodes1<-activated_network_nodes1
        rv$activated_network_edges1<-activated_network_edges1
        display_data<-processActivateNetworkDisplay(activated_network_nodes1,activated_network_edges1)
        display_json<-networkDataToJson(display_data[[1]],display_data[[2]])
        rv$activated_network_json1<-display_json
        if(rv$activatedNetworkGenerated1==1){
          session$sendCustomMessage("refreshActivatedNetwork1", "refresh")
        }
        rv$activatedNetworkGenerated1=1
        loadActivatedNetwork1(rv,rv$cell_type1)
        rm(activated_network_nodes1,activated_network_edges1,display_data,display_json)
      }
      
      if(haveActivatedNetwork2==1){
        rv$activated_network_nodes2<-activated_network_nodes2
        rv$activated_network_edges2<-activated_network_edges2
        display_data<-processActivateNetworkDisplay(activated_network_nodes2,activated_network_edges2)
        display_json<-networkDataToJson(display_data[[1]],display_data[[2]])
        rv$activated_network_json2<-display_json
        if(rv$activatedNetworkGenerated2==1){
          session$sendCustomMessage("refreshActivatedNetwork2", "refresh")
        }
        rv$activatedNetworkGenerated2=1
        loadActivatedNetwork2(rv,rv$cell_type2)
        rm(activated_network_edges2,activated_network_nodes2,display_json,display_data)
      }
      
      if(haveNetwork1==1){
        rv$downstream_network_nodes1<-downstream_network_nodes1
        rv$downstream_network_edges1<-downstream_network_edges1
        display_json<-networkDataToJson(downstream_network_nodes1,downstream_network_edges1)
        rv$downstream_network_json1<-display_json
        if(rv$haveGenerated1==1){
          session$sendCustomMessage("refreshNetwork1", "refresh")
        }
        rv$haveGenerated1=1
        loadNetwork1(rv,rv$cell_type1,rv$cell_type2)
        rm(downstream_network_edges1,downstream_network_nodes1,display_json)
      }
      
      if(haveNetwork2==1){
        rv$downstream_network_nodes2<-downstream_network_nodes2
        rv$downstream_network_edges2<-downstream_network_edges2
        display_json<-networkDataToJson(downstream_network_nodes2,downstream_network_edges2)
        rv$downstream_network_json2<-display_json
        if(rv$haveGenerated2==1){
          session$sendCustomMessage("refreshNetwork2", "refresh")
        }
        rv$haveGenerated2=1
        loadNetwork2(rv,rv$cell_type1,rv$cell_type2)
        rm(downstream_network_edges2,downstream_network_nodes2,display_json)
      }
      
      progress$set(0.9,detail="reload drug data")
        if(file.exists(paste0(network1Dir,"/targetDrug.RData"))){
          load(paste0(network1Dir,"/targetDrug.RData"))
          rv$drugNetworkJson1<-targetDrug1[[2]]
          rv$drug_mapping_table1<-targetDrug1[[1]]
          display_data<-processDrugClusteringDisplay(targetDrug1[[4]],targetDrug1[[5]])
          display_json<-networkDataToJson(display_data[[1]],display_data[[2]])
          rv$targetDrug_json1<-display_json
          loadTargetDrugNetwork1(rv,rv$cell_type1,rv$cell_type2)
          rm(targetDrug1,display_data,display_json)
        }
        if(file.exists(paste0(network2Dir,"/targetDrug.RData"))){
          load(paste0(network2Dir,"/targetDrug.RData"))
          rv$drugNetworkJson2<-targetDrug2[[2]]
          rv$drug_mapping_table2<-targetDrug2[[1]]
          display_data<-processDrugClusteringDisplay(targetDrug2[[4]],targetDrug2[[5]])
          display_json<-networkDataToJson(display_data[[1]],display_data[[2]])
          rv$targetDrug_json2<-display_json
          loadTargetDrugNetwork2(rv,rv$cell_type1,rv$cell_type2)
          rm(targetDrug2,display_data,display_json)
        }
      
      if(file.exists(paste0(network1Dir,"/signalingDrug.RData"))){
        load(paste0(network1Dir,"/signalingDrug.RData"))
        rv$drug_table1<-signalingDrug1[[1]]
        display_data<-processDrugClusteringDisplay(signalingDrug1[[3]],signalingDrug1[[4]])
        display_json<-networkDataToJson(display_data[[1]],display_data[[2]])
        rv$drug_json1<-display_json
        rm(display_data,display_json)
      }
      if(file.exists(paste0(network2Dir,"/signalingDrug.RData"))){
        load(paste0(network2Dir,"/signalingDrug.RData"))
        rv$drug_table2<-signalingDrug2[[1]]
        display_data<-processDrugClusteringDisplay(signalingDrug2[[3]],signalingDrug2[[4]])
        display_json<-networkDataToJson(display_data[[1]],display_data[[2]])
        rv$drug_json2<-display_json
        rm(display_data,display_json)
      }
      
      
 
      
        if(rv$drugNetworkGenerated==1){
          session$sendCustomMessage("refreshDrugNetwork1", "refresh")
          session$sendCustomMessage("refreshDrugNetwork2", "refresh")
          session$sendCustomMessage("refreshTargetCluster1", "refresh")
          session$sendCustomMessage("refreshTargetCluster2", "refresh")
        }
        rv$drugNetworkGenerated=1
        
        if(rv$drugGenerated==1){
          session$sendCustomMessage("refreshCluster1", "refresh")
          session$sendCustomMessage("refreshCluster2", "refresh")
        }
        rv$drugGenerated=1
        gc()
      }
  }
  progress$set(1,detail="Finish")
}



loadNetworkData<-function(rv,progress){
  progress$set(0.6, detail ="Reload communication data")
  if(file.exists(paste0(rv$outputDir,"/network_df.RData"))){
    load(paste0(rv$outputDir,"/network_df.RData"))
    rv$network_df<-seurat_data
    rv$network_gene_list<-rownames(rv$network_df)
    rv$network_num_sample<-ncol(rv$network_df)
    load(paste0(rv$outputDir,"/networkTypeInformation.RData"))
    rv$network_type_list<-type_list
    loadGeneExpression1()
    rm(seurat_data,type_list)
    gc()
  }else{
    return (-1)
  }
  if(file.exists(paste0(rv$outputDir,"/networkGroupInformation.RData"))){
    load(paste0(rv$outputDir,"/networkGroupInformation.RData"))
    rv$network_have_group=1
    rv$network_group_list<-group_list
    rm(group_list)
    gc()
  }else{
    rv$network_have_group=0
  }
  return(1)
}



loadRawAnalyticsData<-function(rv,progress){
  progress$set(0.1, detail ="Reload raw data")
  if(file.exists(paste0(rv$outputDir,"/rna_df.RData"))){
    #   if(input$rowName==0){
    #    rv$raw_df <- read.csv("./outputData/raw_data.csv",
    #                     header = input$header,
    #                      sep = input$sep,
    #                      quote = input$quote,row.names = NULL,check.names = FALSE)
    #  }else{
    #    rv$raw_df <- read.csv("./outputData/raw_data.csv",
    #                      header = input$header,
    #                      sep = input$sep,
    #                      quote = input$quote,row.names = 1,check.names = FALSE)
    # }
    load(paste0(rv$outputDir,"/rna_df.RData"))
    rv$rna_df<-seurat_data
    rv$rna_gene_list<-rownames(seurat_data)
    rv$rna_data_mean<-rowMeans(seurat_data)
    rv$rna_num_sample<-ncol(seurat_data)
    #load dimension reduction part
    loadDimensionReduction()
    rm(seurat_data)
    gc()
  }else{
    return(-1)
  }
  if(file.exists(paste0(rv$outputDir,"/rnaGroupInformation.RData"))){
    load(paste0(rv$outputDir,"/rnaGroupInformation.RData"))
    if(!is.null(group_list)){
      rv$rna_have_group=1
      rv$rna_group_list<-group_list
      rm(group_list)
      gc()
    }else{
      rv$rna_have_group=0
    }
  }else{
    rv$rna_have_group=0
  }
  
  progress$set(0.4, detail ="Reload dimension reduction data")
  if(file.exists(paste0(rv$outputDir,"/encoder_result.RData"))){
    load(paste0(rv$outputDir,"/encoder_result.RData"))
    rv$encoder_result<-encoder_result
    rm(encoder_result)
    gc()
  }else{
    return(0)
  }
  
  #dimension reduction result reload
  if(file.exists(paste0(rv$outputDir,"/tsne_result.RData"))){
    load(paste0(rv$outputDir,"/tsne_result.RData"))
    rv$df_dr<-tsne_result
    #load pre cluster part
    loadPreCluster()
    rm(tsne_result)
    gc()
  }else{
    return(0)
  }
  
  progress$set(0.45, detail ="Reload pre cluster data")
  #pre-cluster result reload
  if(file.exists(paste0(rv$outputDir,"/preClusterResult.RData"))){
    load(paste0(rv$outputDir,"/preClusterResult.RData"))
    rv$dfDensity<-dfDensity
    #load pre cluster part
    loadMainCluster()
    rm(dfDensity)
    gc()
  }else{
    return(0)
  }
  
  progress$set(0.48, detail ="Reload main cluster data")
  #main cluster result reload
  if(file.exists(paste0(rv$outputDir,"/mainClusterResult.RData"))){
    load(paste0(rv$outputDir,"/mainClusterResult.RData"))
    rv$dfClust<-dfClust
    clusters<-dfClust$cluster
    #load pre cluster part
    loadSubCluster()
    rv$df_cluster<-cbind(rv$df_dr,clusters)
    colnames(rv$df_cluster)[3]<-"cluster"
    rm(dfClust)
    gc()
  }else{
    return(0)
  }
  
  progress$set(0.51, detail ="Reload sub cluster data")
  #sub cluster result reload
  if(file.exists(paste0(rv$outputDir,"/subClusterResult.RData"))){
    load(paste0(rv$outputDir,"/subClusterResult.RData"))
    rv$sub_df_cluster<-sub_df_cluster
    loadClassification()
    rm(sub_df_cluster)
    gc()
  }else{
    return(0)
  }
  progress$set(0.55, detail ="Reload Cell annotation data")
  #classification result reload
  if(file.exists(paste0(rv$outputDir,"/cell_annotation.RData"))){
    load(paste0(rv$outputDir,"/cell_annotation.RData"))
    rv$df_classify<-cell_annotation
    rv$rna_type_list<-cell_annotation$type
    loadGeneExpression1()
    rm(cell_annotation)
    gc()
  }else{
    progress$inc(1,detail="Finish")
    return(1)
  }
  return(1)
}

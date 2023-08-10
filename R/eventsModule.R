#preproessing
preprocessingEvent<-function(input,rv,session,use_rna_data,imp){
  progress<-shiny::Progress$new()
  progress$set(message="Preprocessing",value=0)
  on.exit(progress$close())
  tryCatch({
    if(use_rna_data){
      data<-rv$rna_df
      group_list<-rv$rna_group_list
      gene_list<-rv$rna_gene_list
      raw_num_cell<-ncol(data)
      raw_num_gene<-nrow(data)
      gene_num_after_convert<-raw_num_gene
      #quality control
      progress$set(0, detail ="Remove low quality samples and genes")
      if(input$rnaMouseToHuman==T){
        data<-mouse_to_human_convert(data,gene_list)
        gene_list<-rownames(data)
        gene_num_after_convert<-length(gene_list)
      }
      if(is.null(group_list)){
        data_list<-qualityControl(data=data,gene_list=gene_list)
        data<-data_list[[1]]
        gene_list<-data_list[[2]]
        keep_cell_index<-data_list[[3]]
        keep_gene_index<-data_list[[4]]
      }else{
        data_list<-qualityControl(data=data,gene_list=gene_list,group_list=group_list)
        data<-data_list[[1]]
        gene_list<-data_list[[2]]
        group_list<-data_list[[3]]
        keep_cell_index<-data_list[[4]]
        keep_gene_index<-data_list[[5]]
      }
    }else{
      if(is.null(rv$network_type_list)){
        stop(safeError("no cell type data"))
      }
      data<-rv$network_df
      type_list<-rv$network_type_list
      group_list<-rv$network_group_list
      gene_list<-rv$network_gene_list
      raw_num_cell<-ncol(data)
      raw_num_gene<-nrow(data)
      gene_num_after_convert<-raw_num_gene
      #quality control
      progress$set(0, detail ="Remove low quality samples and genes")
      if(input$networkMouseToHuman==T){
        data<-mouse_to_human_convert(data,gene_list)
        gene_list<-rownames(data)
        gene_num_after_convert<-length(gene_list)
      }
      if(is.null(group_list)){
        data_list<-qualityControl(data=data,gene_list=gene_list,type_list=type_list)
        data<-data_list[[1]]
        gene_list<-data_list[[2]]
        type_list<-data_list[[3]]
        keep_cell_index<-data_list[[4]]
        keep_gene_index<-data_list[[5]]
      }else{
        data_list<-qualityControl(data=data,gene_list=gene_list,group_list=group_list,type_list=type_list)
        data<-data_list[[1]]
        gene_list<-data_list[[2]]
        type_list<-data_list[[3]]
        group_list<-data_list[[4]]
        keep_cell_index<-data_list[[5]]
        keep_gene_index<-data_list[[6]]
        
      }
    }

    print(dim(data))
    #normalization
    progress$set(0.4, detail ="Normalization")
    seurat_data<-CreateSeuratObject(counts=data)
    seurat_data<-normalization_v2(seurat_data)
    
    progress$set(0.5, detail ="Imputation")
    if(imp){
      seurat_data<-imputation(seurat_data)
    }
    
    # progress$set(0.5, detail ="Regress out cell cycle effect and rescale")
    # if(input$regress_out){
    #   seurat_data<-rescale(seurat_data,regress=T)
    # }else{
    #   seurat_data<-rescale(seurat_data,regress=F)
    # }
    # data<-seurat_data[["alra"]]@scale.data
    progress$set(0.9, detail ="Fina Variable genes")   
    seurat_data<-FindVariableFeatures(seurat_data,nfeatures = 2048,selection.method = "disp")
    
    progress$set(0.95, detail ="Finish")
    if(use_rna_data){
      rv$rna_df<-seurat_data
      rv$rna_group_list<-group_list
      rv$rna_gene_list<-gene_list
      save(seurat_data,file=paste0(rv$outputDir,"/rna_df.RData"))
      if(!is.null(group_list)){
        save(group_list,file=paste0(rv$outputDir,"/rnaGroupInformation.RData"))
      }
      drop_num_cell<-raw_num_cell-length(keep_cell_index)
      drop_num_gene<-gene_num_after_convert-length(keep_gene_index)
      save(keep_cell_index,keep_gene_index,file=paste0(rv$outputDir,"/rnaPreproceessingResult.RData"))
      
      rm(seurat_data,group_list,gene_list,data_list)
      gc()
      #load dimension reduction part
      loadDimensionReduction()
    }else{
      if(!is.null(group_list)){
        group_type_list<-paste(group_list,type_list,sep="_")
        seurat_data<-SetIdent(seurat_data,value=group_type_list)
      }else{
        seurat_data<-SetIdent(seurat_data,value=type_list)
      }
      rv$network_df<-seurat_data
      rv$network_group_list<-group_list
      rv$network_gene_list<-gene_list
      rv$network_type_list<-type_list
      save(seurat_data,file=paste0(rv$outputDir,"/network_df.RData"))
      save(type_list,file=paste0(rv$outputDir,"/networkTypeInformation.RData"))
      if(!is.null(group_list)){
        save(group_list,file=paste0(rv$outputDir,"/networkGroupInformation.RData"))
      }
      drop_num_cell<-raw_num_cell-length(keep_cell_index)
      drop_num_gene<-gene_num_after_convert-length(keep_gene_index)
      save(keep_cell_index,keep_gene_index,file=paste0(rv$outputDir,"/networkPreproceessingResult.RData"))
      rm(seurat_data,group_list,gene_list,type_list,data_list)
      gc()
      #load ExpressionPart
      loadGeneExpression1()
    }
    rm(keep_cell_index,keep_gene_index)
    gc()
    progress$set(1, detail ="Finish")
    if(input$rnaMouseToHuman==F){
      createAlert(session=session,anchorId = "errorAlert",title="Preprocessing Completed",
                  content = paste0("The data have total of ",raw_num_cell," cells and ",raw_num_gene," genes. After preprocessing, ", drop_num_cell,
                                   " of cells and ",drop_num_gene," genes were removed."),
                  style="success")
    }else{
      createAlert(session=session,anchorId = "errorAlert",title="Preprocessing finish",
                  content = paste0("The data have total of ",raw_num_cell," cells and ",raw_num_gene," genes. After converting mouse to human gene, 
                  The data have total of ",gene_num_after_convert," gene. After preprocessing, ", drop_num_cell," cells and ",
                                   drop_num_gene," genes were removed."),
                  style="success")
    }

  },
  error=function(e){
    if(e$message=="no cell type data"){
      createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                  content ="Please upload cell type data first.",
                  style="danger")
    }else{
      createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                  content = paste("Some thing wrong when carry out preprocessing",e,sep=":"),
                  style="danger")
    }

  }
  )
  
  
  
}

#auto encoder event
aeEvent<-function(input,rv,session){
  aeProgress <- shiny::Progress$new()
  aeProgress$set(message = "Auto Encoder", value = 0)
  on.exit(aeProgress$close())
  tryCatch({
    aeProgress$set(0, detail ="Find high variance gene in data set")
    #find first 2048 largest variance gene in data
    df_sd<-findSD(rv$rna_df)
    df_matrix<-transposition(df_sd)
    aeProgress$set(0.15, detail ="normalization")
    df_normal<-normalization(df_matrix)
    df_normal<-as.matrix(df_normal)
    aeProgress$set(0.37, detail ="training auto encoder")
    encoder_result<-autoEncoder(df_normal,2048,aeProgress)
    rv$encoder_result<-encoder_result
    save(encoder_result,file=paste0(rv$outputDir,"/encoder_result.RData"))
    rm(df_sd,df_matrix,df_normal,encoder_result)
    gc()
    aeProgress$set(1, detail ="Finish")
  },
  error=function(e){
    print(e)
    createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                content = paste("Some thing wrong when carry out dimension reduction analytics",e,sep=":"),
                style="danger")
  }
  
  )
}

#T-sne event
tsneEvent<-function(input,rv,session){
  tryCatch({
    progress<-shiny::Progress$new()
    progress$set(message="T-sne",value=0)
    on.exit(progress$close())
    if(input$useAutoEncoder==T){
      if(!is.null(rv$encoder_result)){
        progress$set(0.2, detail ="run T-sne")
        #tsne dimension reduction
        tsne_result<-runTsne(rv$encoder_result,input$maxIterationSlider,input$perplexitySlider)
        rv$df_dr<-tsne_result
        #save tsne result
        save(tsne_result,file=paste0(rv$outputDir,"/tsne_result.RData"))
        rm(tsne_result)
        gc()
        #load pre clustering part
        loadPreCluster()
      }else{
        stop(safeError("no aeData"))
      }
    }else{
      progress$set(0.05, detail ="Find high variance gene in data set")
      df_sd<-findSD(rv$rna_df)
      df_normal<-normalization(df_sd)
      df_normal<-as.matrix(df_normal)
      progress$set(0.3, detail ="run T-sne")
      tsne_result<-runTsne(df_normal,input$maxIterationSlider,input$perplexitySlider)
      rv$df_dr<-tsne_result
      #save tsne result
      save(tsne_result,file=paste0(rv$outputDir,"/tsne_result.RData"))
      #load pre clustering part
      rm(df_sd,df_normal,tsne_result)
      gc()
      loadPreCluster()
    }
    progress$set(1, detail ="Finish")
  },
  error=function(e){
    print(e)
    if(e$message=="no aeData"){
      createAlert(session=session,anchorId = "errorAlert",title="Data Error",
                  content ="Please carry out auto encoder analytics first!",
                  style="danger")
    }else{
      createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                  content = paste("Some thing wrong when carry out dimension reduction analytics",e,sep=":"),
                  style="danger")
    }
  }
  
  )
}

#pre-clustering event
preClusterEvent<-function(input,rv,session){
  tryCatch({
    dfDensity<-findOrder(rv$df_dr)
    save(dfDensity,file=paste0(rv$outputDir,"/preClusterResult.RData"))
    rv$dfDensity<-dfDensity
    rm(dfDensity)
    gc()
    loadMainCluster()
  },
  error=function(e){
    createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                content = paste("Some thing wrong when carry out pre clustering analytics",e,sep=":"),
                style="danger")
  }
  )
}

#main cluster event
mainClusterEvent<-function(input,rv,session){
  tryCatch({
    dfClust<-extract(rv$dfDensity,input$eps_cl)
    save(dfClust,file=paste0(rv$outputDir,"/mainClusterResult.RData"))
    rv$dfClust<-dfClust
    clusters=rv$dfClust$cluster
    rv$df_cluster<-cbind(rv$df_dr,clusters)
    colnames(rv$df_cluster)[3]<-"cluster"
    rm(clusters,dfClust)
    gc()
    loadSubCluster()
  },
  error=function(e){
    createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                content = paste("Some thing wrong when carry out main clustering analytics",e,sep=":"),
                style="danger")
  }
  )
}

#sub cluster event
subClusterEvent<-function(input,rv,session){
  tryCatch({
    progress <- shiny::Progress$new()
    progress$set(message = "Sub-Clustering", value = 0)
    on.exit(progress$close())
    sub_df_cluster<-sub_cluster(rv$df_cluster,progress)
    save(sub_df_cluster,file=paste0(rv$outputDir,"/subClusterResult.RData"))
    rv$sub_df_cluster<-sub_df_cluster
    rm(sub_df_cluster)
    gc()
    loadClassification()
    progress$set(1,detail="Finish")
  },
  error=function(e){
    createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                content = paste("Some thing wrong when carry out sub clustering analytics",e,sep=":"),
                style="danger")
  }
  )
}
#classification event
classificationEvent<-function(input,rv,session){
  tryCatch({
    clProgress <- shiny::Progress$new()
    clProgress$set(message = "Cell annotation", value = 0)
    on.exit(clProgress$close())
    df<-rv$rna_df[[rv$rna_df@active.assay]]@data
    if(input$useSubCluster==T){
      annotation_result<-gseaClassification(df,rv$sub_df_cluster$cluster,rv$rna_gene_list, rv$select_marker_gene,clProgress,rv$outputDir)
      clProgress$inc(0.95, detail ="Integrate result")
      save(annotation_result,file=paste0(rv$outputDir,"/annotation_result.RData"))
      # prevent auto sort in merge
      sub_df_cluster<-rv$sub_df_cluster
      sub_df_cluster$id<-c(1:nrow(sub_df_cluster))
      #join classification result to every sample
      cell_annotation<-merge(sub_df_cluster,annotation_result,by="cluster",all.x = TRUE)
      cell_annotation<-cell_annotation[order(cell_annotation$id),]
      cell_annotation<-cell_annotation[,-4]
      rm(sub_df_cluster)
    }else{
      annotation_result<-gseaClassification(df,rv$df_cluster$cluster,rv$rna_gene_list,rv$select_marker_gene,clProgress,rv$outputDir)
      clProgress$inc(0.95, detail ="Intgreate result")
      save(annotation_result,file=paste0(rv$outputDir,"/annotation_result.RData"))
      # prevent auto sort in merge
      df_cluster<-rv$df_cluster
      df_cluster$id<-c(1:nrow(df_cluster))
      #join classification result to every sample
      cell_annotation<-merge(df_cluster,annotation_result,by="cluster",all.x = TRUE)
      cell_annotation<-cell_annotation[order(cell_annotation$id),]
      cell_annotation<-cell_annotation[,-4]
      rm(df_cluster)
    }
    save(cell_annotation,file=paste0(rv$outputDir,"/cell_annotation.RData"))
    rv$df_classify<-cell_annotation
    rv$rna_type_list<-cell_annotation$type
    
    if(!is.null(rv$rna_group_list)){
      group_type_list<-paste(rv$rna_group_list,rv$rna_type_list,sep="_")
      seurat_data<-SetIdent(rv$rna_df,value=group_type_list)
      rm(group_type_list)
    }else{
      seurat_data<-SetIdent(rv$rna_df,value=rv$rna_type_list)
    }
    rv$rna_df<-seurat_data
    save(seurat_data,file=paste0(rv$outputDir,"/rna_df.RData"))
    rm(df,annotation_result,cell_annotation,seurat_data)
    gc()
    loadGeneExpression1()
    clProgress$inc(1, detail ="Finish")
  },
  error=function(e){
    createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                content = paste("Some thing wrong when carry out classification analytics",e,sep=":"),
                style="danger")
  }
  )
}

#cell-cell communciation event
communicationEvent<-function(input,rv,session){
  coProgress <- shiny::Progress$new()
  coProgress$set(message = "Discover Communication", value = 0)
  on.exit(coProgress$close())
  tryCatch(
    {
  cell_type1<-input$cellTypeSelect1
  cell_type2<-input$cellTypeSelect2
  if(is.null(cell_type1)||is.null(cell_type2)){
    stop(safeError("empty input"))
  }
  rv$cell_type1<-cell_type1
  rv$cell_type2<-cell_type2
  fc_thres<-input$networkFc
  pv_thres<-input$networkPValue
  fc_downThres<-input$networkDFc
  #check which data user use and whether have group information
  #if no group, use total mean result as normal group


      fcDir<-paste0(rv$outputDir,"/cellCommunication")
      networkDir<-paste0("/",paste(paste(rv$cell_type1,collapse = "+"),paste(rv$cell_type2,collapse = "+"),sep="-"))

      if(rv$useRnaData==1){
        if(rv$rna_have_group==0){
          result<-communication(data=rv$rna_df,gene_symbol=rv$rna_gene_list,type_list=rv$rna_type_list,
                                fc_thres=fc_thres,pv_thres=pv_thres,cell_type1=rv$cell_type1,
                                cell_type2=rv$cell_type2,outputDir=rv$outputDir,ligRecDatabase=rv$ligRecDatabase,
                                TfTargetInteraction=rv$TfTargetInteraction,keggInfo=rv$keggInfo,STRING=rv$STRING,resource=input$networkResource,
                                useOld=input$networkUseOldResult,padjust=input$networkPadjust,
                                progress=coProgress)
        }else{
          testGroup<-input$controlGroup
          normalGroup<-input$normalGroup
          if(is.null(testGroup)||is.null(normalGroup)){
            stop(safeError("empty input"))
          }
          result<-communication(data=rv$rna_df,gene_symbol=rv$rna_gene_list,group_list=rv$rna_group_list,type_list=rv$rna_type_list,
                                normalGroup=normalGroup,testGroup=testGroup,fc_thres=fc_thres,pv_thres=pv_thres,cell_type1=rv$cell_type1,
                                cell_type2=rv$cell_type2,outputDir=rv$outputDir,ligRecDatabase=rv$ligRecDatabase,
                                TfTargetInteraction=rv$TfTargetInteraction,keggInfo=rv$keggInfo,STRING=rv$STRING,resource=input$networkResource,
                                useOld=input$networkUseOldResult,padjust=input$networkPadjust,progress=coProgress)
        }
      }else{
        if(rv$network_have_group==0){
          result<-communication(data=rv$network_df,gene_symbol=rv$network_gene_list,type_list=rv$network_type_list,
                                fc_thres=fc_thres,pv_thres=pv_thres,cell_type1=rv$cell_type1,cell_type2=rv$cell_type2,
                                outputDir=rv$outputDir,ligRecDatabase=rv$ligRecDatabase,
                                TfTargetInteraction=rv$TfTargetInteraction,keggInfo=rv$keggInfo,STRING=rv$STRING,
                                resource=input$networkResource,useOld=input$networkUseOldResult,padjust=input$networkPadjust,
                                progress=coProgress)
        }else{
          testGroup<-input$controlGroup
          normalGroup<-input$normalGroup
          if(is.null(testGroup)||is.null(normalGroup)){
            stop(safeError("empty input"))
          }
          result<-communication(data=rv$network_df,gene_symbol=rv$network_gene_list,group_list=rv$network_group_list,type_list=rv$network_type_list,
                                normalGroup=normalGroup,testGroup=testGroup, fc_thres=fc_thres, pv_thres=pv_thres,cell_type1=rv$cell_type1,
                                cell_type2=rv$cell_type2,outputDir=rv$outputDir,ligRecDatabase=rv$ligRecDatabase,
                                TfTargetInteraction=rv$TfTargetInteraction,keggInfo=rv$keggInfo,useOld=input$networkUseOldResult,
                                STRING=rv$STRING,resource=input$networkResource,padjust=input$networkPadjust,progress=coProgress)
        }
      }
      
      #load gene expression information table
      ligRecInformation<-result[[1]]
      rv$cell_type1_ligands<-ligRecInformation[[1]]
      rv$cell_type1_receptors<-ligRecInformation[[2]]
      rv$cell_type2_ligands<-ligRecInformation[[3]]
      rv$cell_type2_receptors<-ligRecInformation[[4]]
      rm(ligRecInformation)
      
      
      #Drug discovering 
      #set directory for file saving
      netDir<-paste0(fcDir,networkDir)
      network1Dir<-paste0(netDir,"/",paste(paste(rv$cell_type1,collapse = "+"),paste(rv$cell_type2,collapse = "+"),sep="_"))
      network2Dir<-paste0(netDir,"/",paste(paste(rv$cell_type2,collapse = "+"),paste(rv$cell_type1,collapse = "+"),sep="_"))
     
      type1_to_type2_result<-result[[2]]
      activated_network_nodes2<-type1_to_type2_result[[1]]
      activated_network_edges2<-type1_to_type2_result[[2]]
      downstream_network_nodes1<-type1_to_type2_result[[3]]
      downstream_network_edges1<-type1_to_type2_result[[4]]

      
      type2_to_type1_result<-result[[3]]  
      activated_network_nodes1<-type2_to_type1_result[[1]]
      activated_network_edges1<-type2_to_type1_result[[2]]
      downstream_network_nodes2<-type2_to_type1_result[[3]]
      downstream_network_edges2<-type2_to_type1_result[[4]]
      
      
      rm(type1_to_type2_result,type2_to_type1_result,result)
      
      coProgress$set(0.8,detail="Discover drug")
      
      
      if(is.null(rv$rankMatrix)){
        loadRankMatrix(rv)
      }
      


      if(is.null(downstream_network_nodes1)){
        targetDrug1<-NULL
        signalingDrug1<-NULL
      }else{
        #targets drug discoveriing
        targetDrug1<-findDrug(downstream_network_nodes1,downstream_network_edges1,rv$drugBankInteraction,rv$drugBankInformation)
        #signaling signature drug discovering
        signalingDrug1<-findDrug2(downstream_network_nodes1,rv$gSymZs,rv$rankMatrix,rv$drug_mapping,input$drugNumber,input$useFDAOnly,rv$drugBankInformation)
        
      }
      if(is.null(downstream_network_nodes2)){
        targetDrug2<-NULL
        signalingDrug2<-NULL
      }else{
        #targets drug discoveriing
        targetDrug2<-findDrug(downstream_network_nodes2,downstream_network_edges2,rv$drugBankInteraction,rv$drugBankInformation)
        #signaling signature drug discovering
        signalingDrug2<-findDrug2(downstream_network_nodes2,rv$gSymZs,rv$rankMatrix,rv$drug_mapping,input$drugNumber,input$useFDAOnly,rv$drugBankInformation)
        
      }
      
      #load all the results to front and save corresponding files
      #only save file if network exist
      #if at least one of four network don't exist, display notification
      network_exist_flag<-1
      
      if(!is.null(activated_network_nodes1)){
        save(activated_network_nodes1,activated_network_edges1,cell_type1,file=paste0(netDir,"/activatedNetworkData1.RData"))
        display_data<-processActivateNetworkDisplay(activated_network_nodes1,activated_network_edges1)
        display_json<-networkDataToJson(display_data[[1]],display_data[[2]])
        rv$activated_network_json1<-display_json
        rv$activated_network_nodes1<-activated_network_nodes1
        rv$activated_network_edges1<-activated_network_edges1
        shinyjs:::show("activatedNetworkPlotDownload1")
        rm(display_data,display_json)
      }else{
        if (input$networkResource=="KEGG"){
          network_exist_flag<-0
        }          
        rv$activated_network_json1<-NULL
      }
      
      if(!is.null(downstream_network_nodes1)){
        save(downstream_network_nodes1,downstream_network_edges1,cell_type1,cell_type2,file=paste0(network1Dir,"/downstreamNetworkData.RData"))
        display_json<-networkDataToJson(downstream_network_nodes1,downstream_network_edges1)
        rv$downstream_network_json1<-display_json
        rv$downstream_network_nodes1<-downstream_network_nodes1
        rv$downstream_network_edges1<-downstream_network_edges1
        shinyjs:::show("networkPlotDownload1")
        rm(display_json)

      }else{
        network_exist_flag<-0
        rv$downstream_network_json1<-NULL
      }

      if(!is.null(activated_network_nodes2)){
        save(activated_network_nodes2,activated_network_edges2,cell_type2,file=paste0(netDir,"/activatedNetworkData2.RData"))
        
        display_data<-processActivateNetworkDisplay(activated_network_nodes2,activated_network_edges2)
        display_json<-networkDataToJson(display_data[[1]],display_data[[2]])
        rv$activated_network_json2<-display_json
        rv$activated_network_nodes2<-activated_network_nodes2
        rv$activated_network_edges2<-activated_network_edges2
        shinyjs:::show("activatedNetworkPlotDownload2")
        print(rv$activated_network_json2)
        rm(display_data,display_json)
      }else{
        if (input$networkResource=="KEGG"){
        network_exist_flag<-0
        }
        rv$activated_network_json2<-NULL
      }
      
      if(!is.null(downstream_network_nodes2)){
        save(downstream_network_nodes2,downstream_network_edges2,cell_type1,cell_type2,file=paste0(network2Dir,"/downstreamNetworkData.RData"))
        display_json<-networkDataToJson(downstream_network_nodes2,downstream_network_edges2)
        rv$downstream_network_json2<-display_json
        rv$downstream_network_nodes2<-downstream_network_nodes2
        rv$downstream_network_edges2<-downstream_network_edges2
        shinyjs:::show("networkPlotDownload2")
        rm(display_json)
        
      }else{
        network_exist_flag<-0
        rv$downstream_network_json2<-NULL
      }
      
      if(network_exist_flag==0){
        createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                    content ="At least one network don't exist. You may consider to adjust the
                                  thresholds.",
                    style="warning")
      }
      #load title for each network 
      loadNetwork1(rv,rv$cell_type1,rv$cell_type2)
      loadNetwork2(rv,rv$cell_type1,rv$cell_type2)
      loadActivatedNetwork1(rv,rv$cell_type1)
      loadActivatedNetwork2(rv,rv$cell_type2)
      
      #load drug discovering result
      
      if(is.null(signalingDrug1)){
        rv$drug_table1<-NULL
        rv$drug_json1<-NULL
      }else{
        rv$drug_table1<-signalingDrug1[[1]]
        display_data<-processDrugClusteringDisplay(signalingDrug1[[3]],signalingDrug1[[4]])
        display_json<-networkDataToJson(display_data[[1]],display_data[[2]])
        rv$drug_json1<-display_json
        rm(display_data,display_json)
        save(signalingDrug1,file=paste0(network1Dir,"/signalingDrug.RData"))
      }

      if(is.null(signalingDrug2)){
        rv$drug_table2<-NULL
        rv$drug_json2<-NULL
      }else{
        rv$drug_table2<-signalingDrug2[[1]]
        display_data<-processDrugClusteringDisplay(signalingDrug2[[3]],signalingDrug2[[4]])
        display_json<-networkDataToJson(display_data[[1]],display_data[[2]])
        rv$drug_json2<-display_json
        rm(display_data,display_json)
        save(signalingDrug2,file=paste0(network2Dir,"/signalingDrug.RData"))
      }
      
    

      if(is.null(targetDrug1)){
        rv$drugNetworkJson1<-NULL
        rv$drug_mapping_table1<-NULL
        rv$targetDrug_json1<-NULL
      }else{
        rv$drugNetworkJson1<-targetDrug1[[2]]
        rv$drug_mapping_table1<-targetDrug1[[1]]
        display_data<-processDrugClusteringDisplay(targetDrug1[[4]],targetDrug1[[5]])
        display_json<-networkDataToJson(display_data[[1]],display_data[[2]])
        rv$targetDrug_json1<-display_json
        rm(display_data,display_json)
        save(targetDrug1,file=paste0(network1Dir,"/targetDrug.RData"))
      }
      
      if(is.null(targetDrug2)){
        rv$drugNetworkJson2<-NULL
        rv$drug_mapping_table2<-NULL
        rv$targetDrug_json2<-NULL
      }else{
        rv$drugNetworkJson2<-targetDrug2[[2]]
        rv$drug_mapping_table2<-targetDrug2[[1]]
        display_data<-processDrugClusteringDisplay(targetDrug2[[4]],targetDrug2[[5]])
        display_json<-networkDataToJson(display_data[[1]],display_data[[2]])
        rv$targetDrug_json2<-display_json
        rm(display_data,display_json)
        save(targetDrug2,file=paste0(network2Dir,"/targetDrug.RData"))
      }
      #load target drug network title
      loadTargetDrugNetwork1(rv,rv$cell_type1,rv$cell_type2)
      loadTargetDrugNetwork2(rv,rv$cell_type1,rv$cell_type2)
      
      #refresh all the D3 plot
      #remove previous nodes when user load network plot more than once
      if(rv$haveGenerated1==1){
        session$sendCustomMessage("refreshNetwork1", "refresh")
      }
      rv$haveGenerated1=1
      rm(downstream_network_edges1,downstream_network_nodes1)
      
      #remove previous nodes when user load network plot more than once
      if(rv$haveGenerated2==1){
        session$sendCustomMessage("refreshNetwork2", "refresh")
      }
      rv$haveGenerated2=1
      rm(downstream_network_edges2,downstream_network_nodes2)
      
      #remove previous nodes when user load network plot more than once
      if(rv$activatedNetworkGenerated1==1){
        session$sendCustomMessage("refreshActivatedNetwork1", "refresh")
      }
      rv$activatedNetworkGenerated1=1
      rm(activated_network_edges1,activated_network_nodes1)
      
      #remove previous nodes when user load network plot more than once
      if(rv$activatedNetworkGenerated2==1){
        session$sendCustomMessage("refreshActivatedNetwork2", "refresh")
      }
      rv$activatedNetworkGenerated2=1
      rm(activated_network_edges2,activated_network_nodes2)
      

      
      #remove previous nodes when user load network plot more than once
      if(rv$drugGenerated==1){
        session$sendCustomMessage("refreshCluster1", "refresh")
        session$sendCustomMessage("refreshCluster2", "refresh")
      }
      rv$drugGenerated=1
      #remove previous nodes when user load network plot more than once
      if(rv$drugNetworkGenerated==1){
        session$sendCustomMessage("refreshDrugNetwork1", "refresh")
        session$sendCustomMessage("refreshDrugNetwork2", "refresh")
        session$sendCustomMessage("refreshTargetCluster1", "refresh")
        session$sendCustomMessage("refreshTargetCluster2", "refresh")
      }
      rv$drugNetworkGenerated=1
      
      save(cell_type1,cell_type2,file=paste0(rv$outputDir,"/lastNetworkIndexing.RData"))
      rm(targetDrug1,targetDrug2,signalingDrug1,signalingDrug2,cell_type1,cell_type2)
      gc()
      coProgress$set(1,detail="Finish")
      
    },
    error=function(e){
      if(e$message=="empty input"){
        createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                    content = "Empty input",
                    style="danger")
      }else if(e$message=="no old result"){
        createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                    content = paste("You have no old result saved."),
                    style="danger")
      }
      else{
        createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                    content = paste("Error occur",e,sep=":"),
                    style="danger")
      }
    }
  )
}

#load data in the second part analyics
geneExpressionDataLoadingEvent<-function(input,rv,session){
  #load corresponding data user choose, 
  tryCatch(
    {
      if(input$dataUseSelect==1){
        if(is.null(rv$df_classify)||is.null(rv$rna_df)){
          stop(safeError("noRnaData"))
        }
        rv$useRnaData=1
        loadRnaDataSelection(rv)
        if(rv$rna_have_group==0){
          cell_count<-cellDistribution(rv$rna_type_list,NULL)
          rv$cell_count<-cell_count
          rv$distribution_have_group=0
          
        }else{
          cell_count<-cellDistribution(rv$rna_type_list,rv$rna_group_list)
          rv$cell_count<-cell_count
          rv$distribution_have_group=1
        }
      }else{
        if(is.null(rv$network_df)){
          stop(safeError("noNetworkData"))
        }
        rv$useRnaData=0
        loadNetworkDataSelection(rv)
        if(rv$network_have_group==0){
          cell_count<-cellDistribution(rv$network_type_list,NULL)
          rv$cell_count<-cell_count
          rv$distribution_have_group=0
          
          
        }else{
          cell_count<-cellDistribution(rv$network_type_list,rv$network_group_list)
          rv$cell_count<-cell_count
          rv$distribution_have_group=1
        }
      }
      distribution_have_group<-rv$distribution_have_group
      save(distribution_have_group,cell_count,file=paste0(rv$outputDir,"/cellDistributionPlot.RData"))
      rm(cell_count,distribution_have_group  )
      gc()
      loadGeneExpression2()
      loadEPSelection(rv)
      loadGO(rv)
      loadUpReg(rv)
      if(file.exists("./data/rankMatrix92742a.RData")&&file.exists("./data/drug_mapping.RData")){
        loadCommunication()
      }
      useRnaData=rv$useRnaData
      save(useRnaData,file=paste0(rv$outputDir,"/lastAnalysisDataIndex.RData"))
      
    },
    error=function(e){
      if(e$message=="noRnaData"){
        createAlert(session=session,anchorId = "errorAlert",title="Data Error",
                    content = "There do not have raw data analysis result",
                    style="danger")
      }else if(e$message=="noNetworkData"){
        createAlert(session=session,anchorId = "errorAlert",title="Data Error",
                    content = "There do not have gene expression data",
                    style="danger")
      }else{
        createAlert(session=session,anchorId = "errorAlert",title="Loading Error",
                    content = paste("Something wrong",e,sep=":"),
                    style="danger")
      }
    }
  )
}

#find up-regulated ligands and receptors 
#use total mean as normal group
upRegulateEvent<-function(input,rv,session){
  tryCatch(
    {
      fc_thres<-input$upRegFc
      pv_thres<-input$upRegPValue
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Calculate fold change and P-value", value = 0)
      upRegDir<-paste0(rv$outputDir,"/upRegulatedLigandsReceptors")
      if(!dir.exists(upRegDir)){
        dir.create(upRegDir)
      }
      if(rv$useRnaData==1){
        if(rv$rna_have_group==0){
          groupDir<-paste0(upRegDir,"/noGroup")
          if(!dir.exists(groupDir)){
            dir.create(groupDir)
          }
          if(input$upRegUseOldResult==0){
            genesInformation<-upRegTest(data=rv$rna_df,cell_type_list=rv$rna_type_list,
                                        progress=progress)
          }else{
            if(file.exists(paste0(groupDir,"/genesInformation.RData"))){
              load(paste0(groupDir,"/genesInformation.RData"))
            }else{
              stop(safeError("no old result"))
            }
          }

        }else{
          if(is.null(input$upRegTestGroup)||is.null(input$upRegNormalGroup)){
            stop(safeError("empty input"))
          }
          
          testGroup<-input$upRegTestGroup
          normalGroup<-input$upRegNormalGroup
          groupDir<-paste0(upRegDir,paste0("/",paste(paste(normalGroup,collapse = "+"),paste(testGroup,collapse = "+"),sep="-")))
          if(!dir.exists(groupDir)){
            dir.create(groupDir)
          }
          if(input$upRegUseOldResult==0){
            genesInformation<-upRegTest(data=rv$rna_df,normalGroup=normalGroup,testGroup=testGroup,cell_type_list = rv$rna_type_list,
                                        cell_group_list=rv$rna_group_list,progress = progress)
          }else{
            if(file.exists(paste0(groupDir,"/genesInformation.RData"))){
              load(paste0(groupDir,"/genesInformation.RData"))
            }else{
              stop(safeError("no old result"))
            }
          }

        }
      }else{
        if(rv$network_have_group==0){
          groupDir<-paste0(upRegDir,"/noGroup")
          if(!dir.exists(groupDir)){
            dir.create(groupDir)
          }
          if(input$upRegUseOldResult==0){
          genesInformation<-upRegTest(data=rv$network_df,cell_type_list=rv$network_type_list,
                                progress=progress)
          }else{
            if(file.exists(paste0(groupDir,"/genesInformation.RData"))){
              load(paste0(groupDir,"/genesInformation.RData"))
            }else{
              stop(safeError("no old result"))
            }
          }
        }else{
          if(is.null(input$upRegTestGroup)||is.null(input$upRegNormalGroup)){
            stop(safeError("empty input"))
          }
          testGroup<-input$upRegTestGroup
          normalGroup<-input$upRegNormalGroup
          groupDir<-paste0(upRegDir,paste0("/",paste(paste(normalGroup,collapse = "+"),paste(testGroup,collapse = "+"),sep="-")))
          if(!dir.exists(groupDir)){
            dir.create(groupDir)
          }
          if(input$upRegUseOldResult==0){
          genesInformation<-upRegTest(data=rv$network_df,normalGroup=normalGroup,testGroup=testGroup,cell_type_list = rv$network_type_list,
                                cell_group_list=rv$network_group_list,progress = progress)    
          }else{
            if(file.exists(paste0(groupDir,"/genesInformation.RData"))){
              load(paste0(groupDir,"/genesInformation.RData"))
            }else{
              stop(safeError("no old result"))
            }
          }
        }
      }
      save(genesInformation,file=paste0(groupDir,"/genesInformation.RData"))
      
      
      
      if(rv$useRnaData==1){
        upReg_network<-up_regulated_network(genesInformation,rv$rna_gene_list,fc_thres,fc_thres,pv_thres,ligRecDatabase=rv$ligRecDatabase)
        upStreamData<-upStreamNetwork(genesInformation,rv$rna_gene_list,fc_thres,pv_thres,ligRecDatabase=rv$ligRecDatabase)
      }else{

        upReg_network<-up_regulated_network(genesInformation,rv$network_gene_list,fc_thres,fc_thres,pv_thres,ligRecDatabase=rv$ligRecDatabase)
        upStreamData<-upStreamNetwork(genesInformation,rv$network_gene_list,fc_thres,pv_thres,ligRecDatabase=rv$ligRecDatabase)
      }
      progress$set(0.9,detail="Integration")
      rv$upRegLigand_network<-upReg_network[[1]][[1]]
      rv$upRegReceptor_network<-upReg_network[[2]][[1]]
      if(rv$upRegGenerated==1){
        session$sendCustomMessage("refreshUpRegReceptor", "refresh")
        session$sendCustomMessage("refreshUpRegLigand", "refresh")
      }
      rv$upRegGenerated=1
      if(!is.null(upStreamData[[1]])){
        rv$up_exp_network<-upStreamData[[1]][[1]]
      }else{
        rv$up_exp_network<-NULL
      }
      if(!is.null(upStreamData[[2]])){
        rv$exp_up_network<-upStreamData[[2]][[1]]
        
      }else{
        rv$exp_up_network<-NULL
        print(1)
      }
      if(!is.null(upStreamData[[3]])){
        rv$combine_network<-upStreamData[[3]][[1]]
      }else{
        rv$combine_network<-NULL
      }
      if(!is.null(upStreamData[[4]])){
        rv$up_up_network<-upStreamData[[4]][[1]]
        
      }else{
        rv$up_up_network<-NULL
      }
      if(rv$upStreamGenerated==1)
      {
        session$sendCustomMessage("refreshUpExp", "refresh")
        session$sendCustomMessage("refreshExpUp", "refresh")
        session$sendCustomMessage("refreshCombine", "refresh")
        session$sendCustomMessage("refreshUpUp", "refresh")
      }      
      rv$upStreamGenerated=1
      save(upReg_network,file=paste0(groupDir,"/upRegNetwork.RData"))
      save(upStreamData,file=paste0(groupDir,"/upStreamNetwork.RData"))
      save(groupDir,file=paste0(rv$outputDir,"/lastUpRegIndexing.RData"))
      if(is.null(rv$up_up_network)||is.null(rv$combine_network)||is.null(rv$exp_up_network)||is.null(rv$up_exp_network)){
        createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                    content = "At least one of upstream network don't exist. You may consider to adjust the
                                  thresholds.",
                    style="warning")
      }
      rm(upReg_network,upStreamData,genesInformation)
      gc()
      progress$set(1,detail="Finish")
      
    },
    error=function(e){
      if(e$message=="no old result"){
        createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                    content = paste("You have no old result saved."),
                    style="danger")
      }else if(e$message=="empty input"){
        createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                    content = "Empty input",
                    style="danger")
      }
      else{
        createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                    content = paste("Some thing wrong when carry out up regulate gene analytics.",e),
                    style="danger")
      }
    }
  )
  
}

#calculate EMT-PRO score event
emt_proEvent<-function(input,rv,session){
  tryCatch({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Compute EMT-PRO score", value = 0)
    if(rv$useRnaData==1){
      if(rv$rna_have_group==1){
        if(is.null(input$epGroup)||is.null(input$epType)){
          stop(safeError("empty input"))
        }
        df<-as.matrix(rv$rna_df[[rv$rna_df@active.assay]]@data)
        progress$set(message = "Compute EMT-PRO score", value = 0.4)
        result<-emt_pro(df,rv$rna_gene_list,rv$rna_group_list,rv$rna_type_list,input$epGroup,input$epType)
        emt_score<-result[[1]]
        pro_score<-result[[2]]
      }else{
        if(is.null(input$epType)){
          stop(safeError("empty input"))
        }
        df<-as.matrix(rv$rna_df[[rv$rna_df@active.assay]]@data)
        print(dim(df))
        progress$set(message = "Compute EMT-PRO score", value = 0.4)
        result<-emt_pro(df,rv$rna_gene_list,rv$rna_group_list,rv$rna_type_list,input$epGroup,input$epType)
        emt_score<-result[[1]]
        pro_score<-result[[2]]
      }
    }else{
      if(rv$network_have_group==1){
        if(is.null(input$epGroup)||is.null(input$epType)){
          stop(safeError("empty input"))
        }
        df<-as.matrix(rv$network_df[[rv$network_df@active.assay]]@data)
        progress$set(message = "Compute EMT-PRO score", value = 0.4)
        result<-emt_pro(df,rv$network_gene_list,rv$network_group_list,rv$network_type_list,input$epGroup,input$epType)
        emt_score<-result[[1]]
        pro_score<-result[[2]]
      }else{
        if(is.null(input$epType)){
          stop(safeError("empty input"))
        }
        df<-as.matrix(rv$network_df[[rv$network_df@active.assay]]@data)
        progress$set(message = "Compute EMT-PRO score", value = 0.4)
        result<-emt_pro(df,rv$network_gene_list,rv$network_group_list,rv$network_type_list,input$epGroup,input$epType)
        emt_score<-result[[1]]
        pro_score<-result[[2]]
      }
    }
    progress$set(message = "Compute EMT-PRO score", value = 0.9)
    emt_pro_score<-data.frame(EMT=emt_score,PRO=pro_score,type=paste(input$epType,collapse = "-"))
    rv$emt_pro_score<-emt_pro_score
    save(emt_pro_score,file=paste0(rv$outputDir,"/emt_pro_score.RData"))
    progress$set(message = "Finish", value = 1)
    rm(emt_pro_score,result,df,emt_score,pro_score)
    gc()
  },
  error=function(e){
    if(e$message=="no enough gene"){
      createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                  content ="There do not have enough marker genes use for computation ",
                  style="danger")
    }else if(e$message=="empty input"){
      createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                  content = "Empty input",
                  style="danger")
    }
    else{
      createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                  content = paste("Some thing wrong when carry out EMT-PRO analytics",e,sep=":"),
                  style="danger")
    }
  }
  
  )
  
}


#add new gene for marker gene table
addRowButtonEvent<-function(rv){
  colNum<-ncol(rv$display_marker_gene)
  newRow<-rep(0,times=colNum)
  newRow<-c("newGene",newRow)
  rv$display_marker_gene<-rbind(rv$display_marker_gene,newRow)
}

addColButtonEvent<-function(input,rv){
  rowNum<-nrow(rv$display_marker_gene)
  newCol<-rep(0,times=rowNum)
  rv$display_marker_gene<-cbind(rv$display_marker_gene,newCol)
  print(input$addColButton)
  colNum<-ncol(rv$display_marker_gene)
  colnames(rv$display_marker_gene)[colNum]<-input$newCellType 
}



#delete gene for marker gene table
deleteRowButtonEvent<-function(input,rv){
  if (!is.null(input$geneTable_rows_selected)) {
    rv$display_marker_gene <- rv$display_marker_gene[-as.numeric(input$geneTable_rows_selected),]
  }
}

#select cell type related to ad
adButtonEvent<-function(rv,session){
  cellList<-colnames(rv$display_marker_gene)[-1]
  updateCheckboxGroupInput(session=session, inputId="cellSelection", label = NULL,
                           choices = cellList,choiceValues=cellList,
                           selected=c("Endothelial","Per","In","Ex","Ast","OPC","Oli","Mic"))
}
#select cell type related to pdac
pdacButtonEvent<-function(rv,session){
  cellList<-colnames(rv$display_marker_gene)[-1]
  updateCheckboxGroupInput(session=session, inputId="cellSelection", label = NULL,
                           choices = cellList,choiceValues=cellList,
                           selected=c("Ductal1","Ductal2","Acinar","Endocrine","Endothelial","Fibroblast","Stellate",
                                      "Macrophage","Tcell","Bcell"))
}

#select corresponding gene cell type and refresh marker gene table
cellSelectionEvent<-function(rv,input){
  #if user do not choose any type, display original table
  if(is.null(input$cellSelection)){
    rv$select_marker_gene<-rv$display_marker_gene
  }
  else{
    list<-c("Gene Symbol",input$cellSelection)
    rv$select_marker_gene<-rv$display_marker_gene[,list]
  }
}

#display original table
originalButtonEvent<-function(rv){
  rv$display_marker_gene<-rv$original_marker_gene
}

#download cell-cell communication plot
downloadNetworkEvent<-function(session){
  session$sendCustomMessage("downloadNetwork","dir")
}

GOEvent<-function(input,rv,session){
  tryCatch(
    {
      progress <- shiny::Progress$new()
      progress$set(message = "GO analytics", value = 0)
      on.exit(progress$close())
      
      cell_type<-input$GOCellType
      upFc_thres<-input$GOUpFc
      dnFc_thres<-input$GODnFc
      pv_thres<-input$GOPValue
      if(is.null(cell_type)){
        stop(safeError("empty input"))
      }
      if(rv$useRnaData==1){
        if(rv$rna_have_group==0){
          df<-rv$rna_df[[rv$rna_df@active.assay]]@data
          df_t<-df[,rv$rna_type_list%in%cell_type]
          df_n<-df[,!rv$rna_type_list%in%cell_type]
          GO_result<-GOAnalytics(df_t,df_n,rv$rna_gene_list,upFc_thres,dnFc_thres,pv_thres,0,rv$outputDir,progress)
        }else{
          testGroup<-input$GOTestGroup
          normalGroup<-input$GONormalGroup
          df<-rv$rna_df[[rv$rna_df@active.assay]]@data
          df_n<-df[,rv$rna_type_list==cell_type&rv$rna_group_list%in%normalGroup]
          df_t<-df[,rv$rna_type_list==cell_type&rv$rna_group_list%in%testGroup]
          GO_result<-GOAnalytics(df_t,df_n,rv$rna_gene_list,upFc_thres,dnFc_thres,pv_thres,1,rv$outputDir,progress)
        }
      }else{
        if(rv$network_have_group==0){
          df<-rv$network_df[[rv$network_df@active.assay]]@data
          df_t<-df[,rv$network_type_list%in%cell_type]
          df_n<-df[,!rv$network_type_list%in%cell_type]
          GO_result<-GOAnalytics(df_t,df_n,rv$network_gene_list,upFc_thres,dnFc_thres,pv_thres,0,rv$outputDir,progress)
        }else{
          testGroup<-input$GOTestGroup
          normalGroup<-input$GONormalGroup
          df<-rv$network_df[[rv$network_df@active.assay]]@data
          df_n<-df[rv$network_type_list==cell_type&rv$network_group_list%in%normalGroup,]
          df_t<-df[rv$network_type_list==cell_type&rv$network_group_list%in%testGroup,]
          GO_result<-GOAnalytics(df_t,df_n,rv$network_gene_list,upFc_thres,dnFc_thres,pv_thres,1,rv$outputDir,progress)
        }
      }
      save(GO_result,file=paste0(rv$outputDir,"/GO_result.RData"))
      rv$GO_up_table<-GO_result[[1]]
      rv$GO_dn_table<-GO_result[[2]]
      rv$netGO_up<-GO_result[[3]]
      rv$netGO_dn<-GO_result[[4]]
      GO_list<-unique(c(as.character(rv$GO_up_table[,1]),as.character(rv$GO_dn_table[,1])))
      GO_name<-unique(c(as.character(paste(rv$GO_up_table[,1],rv$GO_up_table[,2],sep="  ")),
                 as.character(paste(rv$GO_dn_table[,1],rv$GO_dn_table[,2],sep="  "))))
      names(GO_list)<-GO_name
      rv$GO_list<-GO_list
      selectGONetwork(input,rv,session)
      loadGONetworkSelect(rv,rv$GO_list)
      # GO<-rv$GO_list[1]
      # GO_name<-names(rv$GO_list[1])
      # rv$GONetwork<-GONetworkGenerating(rv$netGO_up,rv$netGO_dn,GO,GO_name)
      # if(rv$GOGenerated==1){
      #   session$sendCustomMessage("refreshGONetwork1","refresh")
      # }
      # rv$GOGenerated=1
      # rv$GOData=1
      rm(df,df_n,df_t,GO_result,GO,GO_list,GO_name)
      gc()
      progress$set(1,detail="Finish")
    },
    error=function(e){
      if(e$message=="empty input"){
          createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                      content = "Empty input",
                      style="danger")
      }
      createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                  content = paste("Some thing wrong when carry out GO analytics",e,sep=":"),
                  style="danger")
    }
  )
  
  
}
selectGONetwork<-function(input,rv,session){
  GO<-input$GOList
  GO<-rv$GO_list[rv$GO_list==GO]
  GO_name<-names(GO)
  rv$GONetwork<-GONetworkGenerating(rv$netGO_up,rv$netGO_dn,GO,GO_name)
  if(rv$GOGenerated==1){
    session$sendCustomMessage("refreshGONetwork1","refresh")
  }
  rv$GOGenerated=1
  rv$GOData=1
  rm(GO,GO_name)
  gc()
}


processDrugFileEvent<-function(defaultDir,input,rv,session){
  progress <- shiny::Progress$new()
  progress$set(message = "Processing Drug File", value = 0)
  on.exit(progress$close())
  tryCatch({
    if(is.integer(input$drugFileDir)){
      stop(safeError("no directory"))
    }
    dataDir<-parseDirPath(defaultDir, input$drugFileDir)
    if(file.exists(paste(dataDir,'GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx',sep="/"))&&
       file.exists(paste(dataDir, 'GSE92742_Broad_LINCS_sig_info.txt', sep='/'))&&
       file.exists(paste(dataDir, 'GSE92742_Broad_LINCS_gene_info.txt', sep='/'))&&
       file.exists(paste(dataDir,"GSE92742_Broad_LINCS_pert_info.txt",sep="/"))
    ){
      dataList<-processDrugFile(dataDir,rv$drugBankInformation,progress)
      rv$rankMatrix<-dataList[[1]]
      rv$drug_mapping<-dataList[[2]]
      rv$gSymZs<-dataList[[3]]
      loadRankMatrix(rv)
      rm(dataList)
      gc()
      progress$set(1,detail="Finish")
    }else{
      stop(safeError("no data"))
    }
  },
  error=function(e){
    if(e$message=="no directory"){
      createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                  content ="You select wrong directory",
                  style="danger")
    }else if(e$message=="no data"){
      createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                  content ="Some data don't exist in directory, check instruction carefully.",
                  style="danger")
    }else{
      createAlert(session=session,anchorId = "errorAlert",title="Computation Error",
                  content = paste("Some thing wrong:",e,sep=":"),
                  style="danger")
    }
  }
  )
}

networkDisplayEvent<-function(input){
  if(input$networkDisplaySelect=="activated1"){
    shinyjs:::show("communicationPart3")
    shinyjs:::hide("communicationPart4")
    shinyjs:::hide("communicationPart5")
    shinyjs:::hide("communicationPart6")
  }else if(input$networkDisplaySelect=="activated2"){
    shinyjs:::show("communicationPart4")
    shinyjs:::hide("communicationPart3")
    shinyjs:::hide("communicationPart5")
    shinyjs:::hide("communicationPart6")
  }else if(input$networkDisplaySelect=="downstream1"){
    shinyjs:::show("communicationPart5")
    shinyjs:::hide("communicationPart3")
    shinyjs:::hide("communicationPart4")
    shinyjs:::hide("communicationPart6")
  }else{
    shinyjs:::show("communicationPart6")
    shinyjs:::hide("communicationPart3")
    shinyjs:::hide("communicationPart4")
    shinyjs:::hide("communicationPart5")
  }
}

drugDisplayEvent<-function(input){
if(input$targetDrugDisplaySelect=="downstream1"){
    shinyjs:::show("communicationPart7")
    shinyjs:::show("communicationPart8")
    shinyjs:::show("communicationPart11")
    shinyjs:::hide("communicationPart9")
    shinyjs:::hide("communicationPart10")
    shinyjs:::hide("communicationPart12")
  }else{
    shinyjs:::hide("communicationPart7")
    shinyjs:::hide("communicationPart8")
    shinyjs:::hide("communicationPart11")
    shinyjs:::show("communicationPart9")
    shinyjs:::show("communicationPart10")
    shinyjs:::show("communicationPart12")
  }
}


saveMarkerGeneEvent<-function(rv,session){
  progress <- shiny::Progress$new()
  progress$set(message = "Save new marker gene tablee", value = 0)
  on.exit(progress$close())
  tryCatch({
    if(!is.null(rv$display_marker_gene)){
      if(!is.null(rv$outputDir)){
        markerGeneTable<-rv$display_marker_gene
        save(markerGeneTable,file=paste(rv$outputDir,"markerGeneTable.RData",sep="/"))
      }else{
        stop(safeError("Please set working directory first."))
      }

    }else{
      stop(safeError("Marker gene tale do not exist."))
    }
  },
  error=function(e){
    createAlert(session=session,anchorId = "errorAlert",title="Save Error",
                content =e$message,
                style="warning")
  }
  )
  

  progress$set(1,detail="Finish.")
}









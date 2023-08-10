#cell-cell communication function using KEGG
communication<-function(data,gene_symbol,group_list=NULL,type_list,normalGroup=NULL,testGroup=NULL,
                        fc_thres,pv_thres,cell_type1,cell_type2,outputDir,ligRecDatabase,
                        TfTargetInteraction,keggInfo,STRING,resource,useOld,padjust,progress){
  
  fcDir<-paste0(outputDir,"/cellCommunication")
  if(!file.exists(fcDir)){
    dir.create(fcDir)
  }
  
  networkDir<-paste0("/",paste(paste(cell_type1,collapse = "+"),paste(cell_type2,collapse = "+"),sep="-"))
  netDir<-paste0(fcDir,networkDir)
  if(!file.exists(netDir)){
    dir.create(netDir)
  }
  
  
  # calculate mean 
  progress$set(0.1,detail="Find Up-regulated Genes")
  if(!useOld){
    if(!is.null(group_list)){
      
      num_normal_group<-length(normalGroup)
      num_test_group<-length(testGroup)
      
      normal_group_type1<-paste(normalGroup,rep(cell_type1,num_normal_group),sep ="_")
      normal_group_type2<-paste(normalGroup,rep(cell_type2,num_normal_group),sep="_")
      test_group_type1<-paste(testGroup,rep(cell_type1,num_test_group),sep="_")
      test_group_type2<-paste(testGroup,rep(cell_type2,num_test_group),sep="_")
      each_condition<-unique(Idents(data))
      normal_group_type1<-intersect(normal_group_type1,each_condition)
      normal_group_type2<-intersect(normal_group_type2,each_condition)
      test_group_type1<-intersect(test_group_type1,each_condition)
      test_group_type2<-intersect(test_group_type2,each_condition)
      if(length(normal_group_type1)==0||length(normal_group_type2)==0||length(test_group_type1)==0||length(test_group_type2)==0){
        stop(safeError("no enough cells for at least one group."))
      }
      
      cell_type1_result1<-FindMarkers(data,ident.1=test_group_type1,ident.2=normal_group_type1,logfc.threshold =0,only.pos = T )
      cell_type1_result2<-FindMarkers(data,ident.1=test_group_type1,ident.2=normal_group_type1,logfc.threshold =0,only.pos = T,test.use = "bimod" )
      intersect_gene<-intersect(rownames(cell_type1_result1),rownames(cell_type1_result2))
      cell_type1_result<-data.frame(logFC=cell_type1_result1[intersect_gene,2],pct1=cell_type1_result1[intersect_gene,3],
                                    pct2=cell_type1_result1[intersect_gene,4],wlicox_p_value=cell_type1_result1[intersect_gene,1],wlicox_p_adjust=cell_type1_result1[intersect_gene,5],
                                    bimod_p_value=cell_type1_result2[intersect_gene,1],bimod_p_adjust=cell_type1_result2[intersect_gene,5])
      rownames(cell_type1_result)<-intersect_gene
      progress$set(0.3,detail="Find Up-regulated Genes")
      cell_type2_result1<-FindMarkers(data,ident.1=test_group_type2,ident.2=normal_group_type2,logfc.threshold =0,only.pos = T )
      cell_type2_result2<-FindMarkers(data,ident.1=test_group_type2,ident.2=normal_group_type2,logfc.threshold =0,only.pos = T ,test.use="bimod")
      intersect_gene<-intersect(rownames(cell_type2_result1),rownames(cell_type2_result2))
      cell_type2_result<-data.frame(logFC=cell_type2_result1[intersect_gene,2],pct1=cell_type2_result1[intersect_gene,3],
                                    pct2=cell_type2_result1[intersect_gene,4],wlicox_p_value=cell_type2_result1[intersect_gene,1],
                                    wlicox_p_adjust=cell_type2_result1[intersect_gene,5],bimod_p_value=cell_type2_result2[intersect_gene,1],
                                    bimod_p_adjust=cell_type2_result2[intersect_gene,5])
      rownames(cell_type2_result)<-intersect_gene
      
      
    }else{
      test_group_type1<-cell_type1
      test_group_type2<-cell_type2
      
      cell_type1_result1<-FindMarkers(data,ident.1=test_group_type1,ident.2=NULL,logfc.threshold =0,only.pos = T )
      cell_type1_result2<-FindMarkers(data,ident.1=test_group_type1,ident.2=NULL,logfc.threshold =0,only.pos = T,test.use = "bimod" )
      intersect_gene<-intersect(rownames(cell_type1_result1),rownames(cell_type1_result2))
      cell_type1_result<-data.frame(logFC=cell_type1_result1[intersect_gene,2],pct1=cell_type1_result1[intersect_gene,3],
                                    pct2=cell_type1_result1[intersect_gene,4],wlicox_p_value=cell_type1_result1[intersect_gene,1],wlicox_p_adjust=cell_type1_result1[intersect_gene,5],
                                    bimod_p_value=cell_type1_result2[intersect_gene,1],bimod_p_adjust=cell_type1_result2[intersect_gene,5])
      rownames(cell_type1_result)<-intersect_gene
      progress$set(0.3,detail="Find Up-regulated Genes")
      cell_type2_result1<-FindMarkers(data,ident.1=test_group_type2,ident.2=NULL,logfc.threshold =0,only.pos = T )
      cell_type2_result2<-FindMarkers(data,ident.1=test_group_type2,ident.2=NULL,logfc.threshold =0,only.pos = T ,test.use="bimod")
      intersect_gene<-intersect(rownames(cell_type2_result1),rownames(cell_type2_result2))
      cell_type2_result<-data.frame(logFC=cell_type2_result1[intersect_gene,2],pct1=cell_type2_result1[intersect_gene,3],
                                    pct2=cell_type2_result1[intersect_gene,4],wlicox_p_value=cell_type2_result1[intersect_gene,1],
                                    wlicox_p_adjust=cell_type2_result1[intersect_gene,5],bimod_p_value=cell_type2_result2[intersect_gene,1],
                                    bimod_p_adjust=cell_type2_result2[intersect_gene,5])
      rownames(cell_type2_result)<-intersect_gene
      
    }
    
    cell_type1_result<-cell_type1_result[order(cell_type1_result$logFC,decreasing = T),]
    cell_type2_result<-cell_type2_result[order(cell_type2_result$logFC,decreasing = T),]
    
    
    
    
    genesInformation<-list(cell_type1_result,cell_type2_result,cell_type1,cell_type2)
    save(genesInformation,file=paste0(netDir,"/genesInformation.RData"))
  }else{
    if(file.exists(paste0(netDir,"/genesInformation.RData"))){
      progress$set(0.5,detail="Find Up-regulated Genes")
      load(paste0(netDir,"/genesInformation.RData"))
      cell_type1_result<-genesInformation[[1]]
      cell_type2_result<-genesInformation[[2]]
    }else{
      stop(safeError("no old result"))
    }
    
    
  }
  
  progress$set(0.5,detail="Find Up-regulated Genes")
  fcCellType1<-cell_type1_result[,1]
  fcCellType2<-cell_type2_result[,1]
  if(padjust){
    pvCellType1<-cell_type1_result[,c(5,7)]
    pvCellType2<-cell_type2_result[,c(5,7)]
  }else{
    pvCellType1<-cell_type1_result[,c(4,6)]
    pvCellType2<-cell_type2_result[,c(4,6)]
  }

  
  
  type1_genes<-rownames(cell_type1_result)
  type2_genes<-rownames(cell_type2_result)
  
  
  #  gene symbol
  progress$set(0.65,detail="Integration")
  
  
  # generate variables
  vLigand <- ligRecDatabase[,1]
  vReceptor <- ligRecDatabase[,2]
  Ligand <- unique(vLigand)
  Receptor <- unique(vReceptor)
  
  cell_type1_ligands<-cell_type1_result[rownames(cell_type1_result)%in%Ligand,]
  cell_type1_receptors<-cell_type1_result[rownames(cell_type1_result)%in%Receptor,]
  
  cell_type2_ligands<-cell_type2_result[rownames(cell_type2_result)%in%Ligand,]
  cell_type2_receptors<-cell_type2_result[rownames(cell_type2_result)%in%Receptor,]
  
  
  ligRecInformation<-list(cell_type1_ligands,cell_type1_receptors,cell_type2_ligands,cell_type2_receptors)
  save(ligRecInformation,file=paste0(netDir,"/ligRecInformation.RData"))
  progress$set(0.7,detail="Generate network")
  
  #discover network 
  #cell type2 to cell type1
  networkDir<-paste0(netDir,"/",paste(paste(cell_type1,collapse = "+"),paste(cell_type2,collapse = "+"),sep="_"))
  if(!file.exists(networkDir)){
    dir.create(networkDir)
  }
  if (resource=="KEGG"){
    type1_to_type2_result<-network_kegg(gene_symbol,type1_genes,type2_genes,fcCellType1,pvCellType1,fcCellType2,pvCellType2,
                                        fc_thres,pv_thres,ligRecDatabase,TfTargetInteraction,keggInfo)
  }else{
    type1_to_type2_result<-network_string(gene_symbol,type1_genes,type2_genes,fcCellType1,pvCellType1,fcCellType2,pvCellType2,
                                        fc_thres,pv_thres,ligRecDatabase,TfTargetInteraction,STRING)
  }

  
  #cell type1 to cell type2
  networkDir<-paste0(netDir,"/",paste(paste(cell_type2,collapse = "+"),paste(cell_type1,collapse = "+"),sep="_"))
  if(!file.exists(networkDir)){
    dir.create(networkDir)
  }
  if (resource=="KEGG"){
    type2_to_type1_result<-network_kegg(gene_symbol,type2_genes,type1_genes,fcCellType2,pvCellType2,fcCellType1,pvCellType1,
                                        fc_thres,pv_thres,ligRecDatabase,TfTargetInteraction,keggInfo)
  }else{
    type2_to_type1_result<-network_string(gene_symbol,type2_genes,type1_genes,fcCellType2,pvCellType2,fcCellType1,pvCellType1,
                                        fc_thres,pv_thres,ligRecDatabase,TfTargetInteraction,STRING)
  }

  return(list(ligRecInformation,type1_to_type2_result,type2_to_type1_result))
}



network_kegg<-function(gSyms,type1_genes,type2_genes,fcCellType1,pvCellType1,fcCellType2,pvCellType2,
                  fc_thres,pv_thres,ligRecDatabase,TfTargetInteraction,keggInfo){
  # generate variables
  vLigand <- ligRecDatabase[,1]
  vReceptor <- ligRecDatabase[,2]
  Ligand <- as.character(unique(vLigand))
  Receptor <- as.character(unique(vReceptor))
  vTF <- TfTargetInteraction[,1]
  vTarget <- TfTargetInteraction[,2]
  TFs <- unique(vTF)
  Targets <- unique(vTarget)
  
  names(fcCellType1)<-type1_genes
  names(fcCellType2)<-type2_genes
  
  
  
  # Set Parameters
  thresStromaLigUpFc <- fc_thres  # Threshold to identify activated Lig from Stroma Cells (fold change and p-value)
  thresStromaLigUpPv <- pv_thres
  
  thresStromaLigExpressedFc <- 0  # Threshold to identify expressed Lig from Stroma Cells (fold change and p-value)
  thresStromaLigExpressedPv <- 1.0
  
  thresEpiRecUpFc <- fc_thres  # Threshold to identify activated Receptor from Epi Cells (fold change and p-value)
  thresEpiRecUpPv <- pv_thres  # for identify the up-regulated (actived) receptors
  
  thresEpiGeneUpFc <- fc_thres# Threshold to identify Up-regulated genes from Epi Cells (fold change and p-value)
  thresEpiGeneUpPv <- pv_thres	 # for the fisher's test p-value of pathy encrichment
  
  thresEpiRecExpressedFc <- 0  # Threshold to identify expressed Receptor from Epi Cells (fold change and p-value)
  thresEpiRecExpressedPv <- 1.0  # for identify the expressed receptors
  
  thresEpiGeneExpressedFc <- 0  # Threshold to identify expressed gene from Epi Cells (fold change and p-value)
  thresEpiGeneExpressedPv <- 1.0  # for identify the expressed genes
  

  
  
  
  
  
  #Expressed Ligands from cell type1
  idx1 <- which(type1_genes %in% Ligand)
  gLig1 <- type1_genes[idx1]
  fcLig1 <- fcCellType1[idx1]
  pvLig1 <- pvCellType1[idx1,]
  ligExp1<-gLig1
  
  # up-regulated ligands from cell type1
  idx1a <- which(fcLig1 >= thresStromaLigUpFc & pvLig1[,1] <= thresStromaLigUpPv & pvLig1[,2] <= thresStromaLigUpPv)
  ligUp1 <- gLig1[idx1a]  
  
  # Expressed Receptor from cell type2
  idx1 <- which(type2_genes %in% Receptor)
  gRec1 <- type2_genes[idx1]
  fcRec1 <- fcCellType2[idx1]
  pvRec1 <- pvCellType2[idx1,]
  recExp1 <- gRec1
  
  # Up-regulated Receptors from cell type1
  idx1a <- which(fcRec1 >= thresEpiRecUpFc & pvRec1[,1] <= thresEpiRecUpPv& pvRec1[,2] <= thresEpiRecUpPv)
  recUp1 <- gRec1[idx1a]  # Up-Expressed Receptor in Epi

  if(length(unique(ligUp1))==0&&length(unique(recUp1))==0){
    return(list(NULL,NULL,NULL,NULL))
  }

  
  # unique genes (considering the multiple probes)
  ligUp1u <- unique(ligUp1)
  ligExp1u <- unique(ligExp1)
  recExp1u <- unique(recExp1)
  recUp1u <- unique(recUp1)
  
  
  #generate ligand-receptor network
  #up-regulated ligand->expressed receptors and expressed receptors->up-regulated receptors
  
  upLig_expRec <- ligRecDatabase[(vLigand %in%ligUp1u) & (vReceptor %in% recExp1u), ]
  expLig_upRec<-ligRecDatabase[(vLigand %in%ligExp1u) & (vReceptor %in% recUp1u),]
  lig_rec_network<-unique(rbind(upLig_expRec,expLig_upRec))
  #downstram network
  #cell type 2 network information
  netInfo<-generateKeggNetInfoV1(type2_genes,fcCellType2,keggInfo)
  #activate signaling network among cell type2
  activated_network<-generateKeggActivationNetworkV1(netInfo,fc_thres)
  #downstream network for cell type2
  if(dim(lig_rec_network)[1]!=0){
    downstream_network<-generateKeggReceptorNetworkV1(netInfo,unique(as.character(lig_rec_network[,2])),fc_thres)
  }else{
    downstream_network<-matrix(ncol=3)
    colnames(downstream_network)<-c('source', 'target', 'PathwayName')
  }

  #check if the network exist
  if(dim(activated_network)[1]==1){
    activated_network<-NULL
  }else{
    activated_network<-activated_network[-1,]
    activated_network<-unique(activated_network)
  }
  
  if(dim(downstream_network)[1]==1){
    downstream_network<-NULL
  }else{
    downstream_network<-downstream_network[-1,c(1,2)]
    downstream_network<-unique(downstream_network)
  }
  
  #Integrate activated network
  if(!is.null(activated_network)){
    activated_network[,1]<-paste(activated_network[,1],activated_network[,3],sep="_")
    activated_network[,2]<-paste(activated_network[,2],activated_network[,3],sep="_")
    activated_edges<-data.frame(activated_network)
    nodes<-unique(c(as.character(activated_edges[,1]),as.character(activated_edges[,2])))
    nodeSymbol<-as.character(sapply(nodes,function(x){strsplit(x,split = "_")[[1]][1]}))
    nodeAtt<-as.character(sapply(nodes,function(x){strsplit(x,split = "_")[[1]][2]}))
    nodeSize<-sapply(nodeSymbol,
                       function(x){
                         if(x%in%names(fcCellType2))
                         {return(fcCellType2[x]*1.5+1)
                         }else{
                           return (1)
                         }
                       }
    )
  activated_nodes<-data.frame(id=nodes,name=nodeSymbol,parent=nodeAtt,size=nodeSize)
  activated_nodes<-activated_nodes[order(as.character(activated_nodes$parent)),]
  rownames(activated_edges)<-NULL
  rownames(activated_nodes)<-NULL
  }else{
    activated_edges<-NULL
    activated_nodes<-NULL
  }
  
  #integrate downstream network
  if(!is.null(downstream_network)){
    colnames(lig_rec_network)<-c("source","target")
    downstream_gene<-unique(c(downstream_network[,1],downstream_network[,2]))
    lig_rec_network<-lig_rec_network[lig_rec_network[,2]%in%downstream_gene,]
    lig_rec_network<-rbind(c("source","target"),lig_rec_network)
    
    ligList<-unique(as.character(lig_rec_network[,1]))
    recList<-unique(as.character(lig_rec_network[,2]))
    lig_rec_network[,1]<-paste(lig_rec_network[,1],"Ligand",sep ="_")
    lig_rec_network[,2]<-paste(lig_rec_network[,2],"Receptor",sep ="_")
    
    recIdx1<-which(downstream_network[,1]%in%recList)
    downstream_network[recIdx1,1]<-paste(downstream_network[recIdx1,1],"Receptor",sep ="_")
    recIdx2<-which(downstream_network[,2]%in%recList)
    downstream_network[recIdx2,2]<-paste(downstream_network[recIdx2,2],"Receptor",sep ="_")
    
    tfIdx1<-which(downstream_network[,1]%in%TFs)
    downstream_network[tfIdx1,1]<-paste(downstream_network[tfIdx1,1],"TranscriptionFactor",sep ="_")
    tfIdx2<-which(downstream_network[,2]%in%TFs)
    downstream_network[tfIdx2,2]<-paste(downstream_network[tfIdx2,2],"TranscriptionFactor",sep ="_")
    
    targetIdx1<-which(downstream_network[,1]%in%Targets)
    downstream_network[targetIdx1,1]<-paste(downstream_network[targetIdx1,1],"Target",sep ="_")
    targetIdx2<-which(downstream_network[,2]%in%Targets)
    downstream_network[targetIdx2,2]<-paste(downstream_network[targetIdx2,2],"Target",sep ="_")
    
    num_edge<-dim(downstream_network)[1]
    linkNodeIdx1<-which(!c(1:num_edge)%in%c(recIdx1,tfIdx1,targetIdx1))
    downstream_network[linkNodeIdx1,1]<-paste(downstream_network[linkNodeIdx1,1],"LinkNode",sep ="_")
    linkNodeIdx2<-which(!c(1:num_edge)%in%c(recIdx2,tfIdx2,targetIdx2))
    downstream_network[linkNodeIdx2,2]<-paste(downstream_network[linkNodeIdx2,2],"LinkNode",sep ="_")
    
    
    communication_edges<-data.frame(rbind(lig_rec_network,downstream_network))
    communication_edges<-communication_edges[-1,]
    
    nodes<-unique(c(as.character(communication_edges[,1]),as.character(communication_edges[,2])))
    nodeSymbol<-as.character(sapply(nodes,function(x){strsplit(x,split = "_")[[1]][1]}))
    nodeAtt<-as.character(sapply(nodes,function(x){strsplit(x,split = "_")[[1]][2]}))
    nodeSize<-rep(1,length(nodes))
    ligandSize<-sapply(nodeSymbol[nodeAtt=="Ligand"],
    function(x){
      if(x%in%names(fcCellType1))
      {return(fcCellType1[x]*1.5+1)
      }else{
        return (1)
      }
    }
  )
    nodeSize[nodeAtt=="Ligand"]=ligandSize
    
    otherSize<-sapply(nodeSymbol[nodeAtt!="Ligand"],
        function(x){
          if(x%in%names(fcCellType2))
          {return(fcCellType2[x]*1.5+1)
          }else{
            return (1)
          }
        }
    )
    nodeSize[nodeAtt!="Ligand"]=otherSize
    communication_nodes<-data.frame(id=nodes,name=nodeSymbol,parent=nodeAtt,size=nodeSize)
    communication_nodes<-communication_nodes[order(as.character(communication_nodes$parent)),]
    rownames(communication_edges)<-NULL
    rownames(communication_nodes)<-NULL
  }else{
      communication_nodes<-NULL
      communication_edges<-NULL
    }
  
  return(list(activated_nodes,activated_edges,communication_nodes,communication_edges))
}

network_string<-function(gene_symbol,type1_genes,type2_genes,fcCellType1,pvCellType1,fcCellType2,pvCellType2,
                         fc_thres,pv_thres,ligRecDatabase,TfTargetInteraction,STRING){
  
  # generate variables
  vLigand <- ligRecDatabase[,1]
  vReceptor <- ligRecDatabase[,2]
  Ligand <- as.character(unique(vLigand))
  Receptor <- as.character(unique(vReceptor))
  vTF <- TfTargetInteraction[,1]
  vTarget <- TfTargetInteraction[,2]
  TFs <- unique(vTF)
  Targets <- unique(vTarget)
  
  names(fcCellType1)<-type1_genes
  names(fcCellType2)<-type2_genes
  
  
  
  # Set Parameters
  thresStromaLigUpFc <- fc_thres  # Threshold to identify activated Lig from Stroma Cells (fold change and p-value)
  thresStromaLigUpPv <- pv_thres
  
  thresStromaLigExpressedFc <- 0  # Threshold to identify expressed Lig from Stroma Cells (fold change and p-value)
  thresStromaLigExpressedPv <- 1.0
  
  thresEpiRecUpFc <- fc_thres  # Threshold to identify activated Receptor from Epi Cells (fold change and p-value)
  thresEpiRecUpPv <- pv_thres  # for identify the up-regulated (actived) receptors
  
  thresEpiGeneUpFc <- fc_thres# Threshold to identify Up-regulated genes from Epi Cells (fold change and p-value)
  thresEpiGeneUpPv <- pv_thres	 # for the fisher's test p-value of pathy encrichment
  
  thresEpiRecExpressedFc <- 0  # Threshold to identify expressed Receptor from Epi Cells (fold change and p-value)
  thresEpiRecExpressedPv <- 1.0  # for identify the expressed receptors
  
  thresEpiGeneExpressedFc <- 0  # Threshold to identify expressed gene from Epi Cells (fold change and p-value)
  thresEpiGeneExpressedPv <- 1.0  # for identify the expressed genes
  
  
  
  
  
  
  
  #Expressed Ligands from cell type1
  idx1 <- which(type1_genes %in% Ligand)
  gLig1 <- type1_genes[idx1]
  fcLig1 <- fcCellType1[idx1]
  pvLig1 <- pvCellType1[idx1,]
  ligExp1<-gLig1
  
  # up-regulated ligands from cell type1
  idx1a <- which(fcLig1 >= thresStromaLigUpFc & pvLig1[,1] <= thresStromaLigUpPv & pvLig1[,2] <= thresStromaLigUpPv)
  ligUp1 <- gLig1[idx1a]  
  
  # Expressed Receptor from cell type2
  idx1 <- which(type2_genes %in% Receptor)
  gRec1 <- type2_genes[idx1]
  fcRec1 <- fcCellType2[idx1]
  pvRec1 <- pvCellType2[idx1,]
  recExp1 <- gRec1
  
  # Up-regulated Receptors from cell type1
  idx1a <- which(fcRec1 >= thresEpiRecUpFc & pvRec1[,1] <= thresEpiRecUpPv& pvRec1[,2] <= thresEpiRecUpPv)
  recUp1 <- gRec1[idx1a]  # Up-Expressed Receptor in Epi
  
  if(length(unique(ligUp1))==0&&length(unique(recUp1))==0){
    return(list(NULL,NULL,NULL,NULL))
  }
  
  
  # unique genes (considering the multiple probes)
  ligUp1u <- unique(ligUp1)
  ligExp1u <- unique(ligExp1)
  recExp1u <- unique(recExp1)
  recUp1u <- unique(recUp1)
  
  
  #generate ligand-receptor network
  #up-regulated ligand->expressed receptors and expressed receptors->up-regulated receptors
  
  upLig_expRec <- ligRecDatabase[(vLigand %in%ligUp1u) & (vReceptor %in% recExp1u), ]
  expLig_upRec<-ligRecDatabase[(vLigand %in%ligExp1u) & (vReceptor %in% recUp1u),]
  lig_rec_network<-unique(rbind(upLig_expRec,expLig_upRec))
  receptor_list=as.character(unique(lig_rec_network[,2]))
  #downstram network
  node_String <- union(STRING[,1],STRING[,2])
  G_String <- generateSigNetwork_setNetWeight(STRING, type2_genes,fcCellType2 )
  gs_1a <- intersect(node_String, receptor_list)
  gs_select1<-type2_genes[fcCellType2>fc_thres&pvCellType2[,1]<pv_thres&pvCellType2[,2]<pv_thres]
  gs_2a <- setdiff(as.character(intersect(gs_select1, node_String)), gs_1a) # gs_select1 is the selected genes (for example, top 100 up-regulated genes)
  
  # generate the network
  nt <- length(gs_1a)
  downstream_network <- c('source', 'target')
  for (i in 1:nt){
    print(i)
    net_tmp <- generateSigNetwork_v3a(G_String, gs_1a[i], gs_2a)
    if(length(net_tmp)<3){
      next
    }
    downstream_network <- rbind(downstream_network, net_tmp)
  }
  if (length(downstream_network)<3){
    downstream_network<-NULL
  }else{
    downstream_network <- unique(downstream_network)
    idx_t <- which(downstream_network[,1] %in% 'source')
    downstream_network <- downstream_network[-idx_t,] # remove the 'Source'/'Target' row
    colnames(downstream_network) <- c('source', 'target')
  }

  
  if(!is.null(downstream_network)){
    colnames(lig_rec_network)<-c("source","target")
    downstream_gene<-unique(c(downstream_network[,1],downstream_network[,2]))
    lig_rec_network<-lig_rec_network[lig_rec_network[,2]%in%downstream_gene,]
    lig_rec_network<-rbind(c("source","target"),lig_rec_network)
    
    ligList<-unique(as.character(lig_rec_network[,1]))
    recList<-unique(as.character(lig_rec_network[,2]))
    lig_rec_network[,1]<-paste(lig_rec_network[,1],"Ligand",sep ="_")
    lig_rec_network[,2]<-paste(lig_rec_network[,2],"Receptor",sep ="_")
    
    recIdx1<-which(downstream_network[,1]%in%recList)
    downstream_network[recIdx1,1]<-paste(downstream_network[recIdx1,1],"Receptor",sep ="_")
    recIdx2<-which(downstream_network[,2]%in%recList)
    downstream_network[recIdx2,2]<-paste(downstream_network[recIdx2,2],"Receptor",sep ="_")
    
    tfIdx1<-which(downstream_network[,1]%in%TFs)
    downstream_network[tfIdx1,1]<-paste(downstream_network[tfIdx1,1],"TranscriptionFactor",sep ="_")
    tfIdx2<-which(downstream_network[,2]%in%TFs)
    downstream_network[tfIdx2,2]<-paste(downstream_network[tfIdx2,2],"TranscriptionFactor",sep ="_")
    
    targetIdx1<-which(downstream_network[,1]%in%Targets)
    downstream_network[targetIdx1,1]<-paste(downstream_network[targetIdx1,1],"Target",sep ="_")
    targetIdx2<-which(downstream_network[,2]%in%Targets)
    downstream_network[targetIdx2,2]<-paste(downstream_network[targetIdx2,2],"Target",sep ="_")
    
    num_edge<-dim(downstream_network)[1]
    linkNodeIdx1<-which(!c(1:num_edge)%in%c(recIdx1,tfIdx1,targetIdx1))
    downstream_network[linkNodeIdx1,1]<-paste(downstream_network[linkNodeIdx1,1],"LinkNode",sep ="_")
    linkNodeIdx2<-which(!c(1:num_edge)%in%c(recIdx2,tfIdx2,targetIdx2))
    downstream_network[linkNodeIdx2,2]<-paste(downstream_network[linkNodeIdx2,2],"LinkNode",sep ="_")
    
    
    communication_edges<-data.frame(rbind(lig_rec_network,downstream_network))
    communication_edges<-communication_edges[-1,]
    
    nodes<-unique(c(as.character(communication_edges[,1]),as.character(communication_edges[,2])))
    nodeSymbol<-as.character(sapply(nodes,function(x){strsplit(x,split = "_")[[1]][1]}))
    nodeAtt<-as.character(sapply(nodes,function(x){strsplit(x,split = "_")[[1]][2]}))
    nodeSize<-rep(1,length(nodes))
    ligandSize<-sapply(nodeSymbol[nodeAtt=="Ligand"],
                       function(x){
                         if(x%in%names(fcCellType1))
                         {return(fcCellType1[x]*1.5+1)
                         }else{
                           return (1)
                         }
                       }
    )
    nodeSize[nodeAtt=="Ligand"]=ligandSize
    
    otherSize<-sapply(nodeSymbol[nodeAtt!="Ligand"],
                      function(x){
                        if(x%in%names(fcCellType2))
                        {return(fcCellType2[x]*1.5+1)
                        }else{
                          return (1)
                        }
                      }
    )
    nodeSize[nodeAtt!="Ligand"]=otherSize
    communication_nodes<-data.frame(id=nodes,name=nodeSymbol,parent=nodeAtt,size=nodeSize)
    communication_nodes<-communication_nodes[order(as.character(communication_nodes$parent)),]
    rownames(communication_edges)<-NULL
    rownames(communication_nodes)<-NULL
  }else{
    communication_nodes<-NULL
    communication_edges<-NULL
  }
  
  return(list(NULL,NULL,communication_nodes,communication_edges))
}


generateKeggReceptorNetworkV1 <- function(netInfo, Rs, T0){
  # library(igraph)
  # Rs: receptors, like: Rs <- c("EGFR")
  # T0: threshold,e.g., log2(1.25)
  # T0 <- log2(1.25) # threshold with mean-node score >= T0
  
  if (1 == 1){ # network reconstruction 
    
    xt <- netInfo[[1]]
    v_t_a <- netInfo[[2]]
    v_t_ai <- netInfo[[3]]
    
    # T0 <- log2(1.25) # threshold with mean-node score >= T0
    T1 <- log2(1.25) # threshold with mean-node score <= -T1
    path_t <-matrix(ncol=3)
    colnames(path_t)<-c('source', 'target', 'PathwayName')
    path_ti <- matrix(ncol=3)
    colnames(path_ti)<-c('source', 'target', 'PathwayName')
    # ---
    nt <- length(xt)
    vt <- rep(0.0, nt)
    for (j in 1:nt){
      yt <- xt[[j]]
      zt1a <- yt[[1]][1]
      zt1b <- yt[[1]][2]
      zt2 <- yt[[2]]
      p_n <- zt2[1,3]  # pathway name
      #
      lt <- dim(zt2)[1]
      if (lt <=2){
        next
      }
      #
      if (zt1a >= T0 & zt1b > v_t_a[lt]){  # activation
        if (length(which(zt2[,1] %in% Rs | zt2[,2] %in% Rs))>0){ # include the receptors
            path_t <- rbind(path_t, zt2)
        }
      }
      if (zt1a >= T1 & zt1b < v_t_ai[lt]){  # inhibition
          path_ti <- rbind(path_ti, zt2)
      }
      
    } # end of network construction
  }
  #path_t <- path_t[,c(1:2)]
  return(path_t)
  
}
# ------------------
# generateKeggNetworkV1
# ------------------
generateKeggActivationNetworkV1 <- function(netInfo, T0){
  # library(igraph)
  # Rs: receptors, like: Rs <- c("EGFR")
  # T0: threshold,e.g., log2(1.25)
  # T0 <- log2(1.25) # threshold with mean-node score >= T0
  
  if (1 == 1){ # network reconstruction 
    
    xt <- netInfo[[1]]
    v_t_a <- netInfo[[2]]
    v_t_ai <- netInfo[[3]]
    
    # T0 <- log2(1.25) # threshold with mean-node score >= T0
    T1 <- log2(1.25) # threshold with mean-node score <= -T1
    path_t <-matrix(ncol=3)
    colnames(path_t)<-c('source', 'target', 'PathwayName')
    path_ti <- matrix(ncol=3)
    colnames(path_ti)<-c('source', 'target', 'PathwayName')
    # ---
    nt <- length(xt)
    vt <- rep(0.0, nt)
    for (j in 1:nt){
      yt <- xt[[j]]
      zt1a <- yt[[1]][1]
      zt1b <- yt[[1]][2]
      zt2 <- yt[[2]]
      p_n <- zt2[1,3]  # pathway name
      #
      lt <- dim(zt2)[1]
      if (lt <=2){
        next
      }
      #
      if (zt1a >= T0 & zt1b > v_t_a[lt]){  # activation
        path_t <- rbind(path_t, zt2)
      }
      if (zt1a >= T1 & zt1b < v_t_ai[lt]){  # inhibition
        path_ti <- rbind(path_ti, zt2)
      }
      
    }
    #
    # print('# of Edge in path_t(activation):'); print(dim(unique(path_t)))
    # print('# of Edge in path_ti(inhibition):'); print(dim(unique(path_ti)))
  } # end of network construction
  #path_t <- unique(path_t[,c(1:2)])
  return(path_t)
  
}

# ------------------
# generateKeggReceptorNetworkV1
# ------------------
generateKeggNetInfoV1 <- function(gSyms, fc_tmp, keggInfo){
  # library(igraph)
  # gSyms: gene symbol
  # fc_tmp: fold change (log2-scale) / or significance score (with mean = 0)
  #
  KeggGenes <- keggInfo[[1]]
  KeggPathways <- keggInfo[[2]]
  KeggGraphs <- keggInfo[[3]]
  
  Ta <- 3.0
  fc_tmp[fc_tmp >Ta] = Ta  #
  fc_tmp[fc_tmp < -1*Ta] = -1*Ta  #
  
  # get the genes in KEGG only
  idx_gSym <- which(gSyms %in% KeggGenes)
  gSym1 <- gSyms[idx_gSym]
  fc_tmp1 <- fc_tmp[idx_gSym]
  gSyms <- gSym1
  # take the maximum or minimum value of multiple genes
  st1 <- unique(gSyms)
  nt <- length(st1)
  fc_tmp <- rep(0.0, nt)
  for (i in 1:nt){
    idxt <- which(gSyms %in% st1[i])
    vt <- fc_tmp1[idxt]
    if (length(vt)>1){
      idxt1 <- which(abs(vt) == max(abs(vt)))
      fc_tmp[i] <- vt[idxt1[1]]
    } else {
      fc_tmp[i] <- vt
    }
  }
  gSyms <- st1
  #
  length(gSyms)
  length(fc_tmp)
  max(fc_tmp)
  min(fc_tmp)
  #
  
  if (1 == 1){ # network reconstruction 
    lPath <- c()
    sigPath <- c()
    sigPath <- getSubPathwayKeggV2a1(fc_tmp, gSyms, KeggPathways, KeggGraphs, lPath)
    pt1 <- sigPath$pathTmp
    lPath <- sigPath$lPath
    
    # extract the paths
    path_score_1 <- c()
    path_score_2 <- c()
    path_score_v <- list()
    path_score_a <- list()
    for (i in 1:50){
      path_score_v[[i]] <- c(0)
      path_score_a[[i]] <- c(0)
    }
    #
    path_t <- c('Source', 'Target', 'PathwayName')
    path_ti <- c('Source', 'Target', 'PathwayName')
    xt <- pt1
    #
    nt <- length(xt)
    vt <- rep(0.0, nt)
    for (j in 1:nt){
      yt <- xt[[j]]
      zt1a <- yt[[1]][1]
      zt1b <- yt[[1]][2]
      zt2 <- yt[[2]]
      p_n <- zt2[1,3]  # pathway name
      #
      lt <- dim(zt2)[1]
      path_score_v[[lt]] <- c(path_score_v[[lt]], zt1a)
      path_score_a[[lt]] <- c(path_score_a[[lt]], zt1b)
      #
      path_score_1 <- c(path_score_1, zt1a)
      path_score_2 <- c(path_score_2, zt1b)
    }
    #
    v_t_a <- rep(0.25, 50)  # activation threshold for the signaling paths with different lengths
    v_t_ai <- rep(-0.25, 50) # inhibition score threshold 
    #
    s_p_a <- matrix(0, 50,3) # to calculate the Real # of signaling pathways, mean-activation, mean-inhibition
    s_p_v <- matrix(0, 50,2) # vertex scores-mean 
    n_t <- matrix(0, 50, 2)
    # for (i in 1:50){
    C1 <- 1.25
    for (i in 1:50){
      # activation score
      vt1 <- path_score_a[[i]]
      idxt <- which(vt1 >0)
      if (length(idxt) < 2){next}
      vt1a <- vt1[vt1>0]
      # v_t_a[i] <- max(0.1, mean(vt1a)+1.0*sd(vt1a));
      v_t_a[i] <- mean(vt1a) + C1*sd(vt1a)
      n_t[i,1] <- (sum(vt1a >= v_t_a[i]))
      idxt <- which(vt1 <0)
      if (length(idxt) < 2){next}
      vt1 <- path_score_a[[i]]
      vt1a <- vt1[vt1<0]  # inhibition
      # v_t_ai[i] <- min(-0.1, mean(vt1a)-1.0*sd(vt1a)); 
      v_t_ai[i] <- mean(vt1a) - C1*sd(vt1a);
      n_t[i,2] <- (sum(vt1a <= v_t_ai[i]))
      # print(sum(vt1 < -v_t_ai[i]))
      s_p_a[i,] <- c(length(vt1), mean(vt1[vt1 > 0]), mean(vt1[vt1 < 0]))
      # vertex score
      vt1 <- path_score_v[[i]]
      s_p_v[i,] <- c(length(vt1), mean(vt1[vt1 > 0]))
    }
    v_t_a[c(11:50)] <- 5.0  # will not consider the paths with more than 11 nodes
    v_t_ai[c(11:50)] <- -5.0
    #
    netInfo <- list()
    netInfo[[1]] <- pt1;
    netInfo[[2]] <- v_t_a;
    netInfo[[3]] <- v_t_ai;
    
  } # end of network construction
  return(netInfo)
  
}
# ------------------
# getSubPathwayKeggV2a1
# ------------------
getSubPathwayKeggV2a1 <- function(fc1, gSyms, KeggPathways, KeggGraphs, lPath){
  #
  # print('test')
  pathway_names <- names(KeggPathways)
  n_pathway <- length(KeggPathways)
  # print('npathway')
  # print(n_pathway)
  pathTmp <- list() # save all the paths individually
  n_pt <- 0  # number of paths
  
  n_pathway <- length(KeggGraphs)
  for (j1 in 1:n_pathway){  # for all signaling pathways
    p_name <- pathway_names[j1]
    gTmp <- KeggGraphs[[j1]]
    
    if (is.null(gTmp)){print('no pathway'); next;} # no network/pathway 
    # check if there are TFs and Receptors in this signaling pathway
    nodes_tmp <- V(gTmp)$name
    # z-scores of nodes_tmp
    idxt <- which(gSyms %in% nodes_tmp)
    zScore <- fc1[idxt]
    names(zScore) <- gSyms[idxt]  # for each sample
    # get the source/sink nodes
    # xt1 <- KeggPathways[[j1]]
    # n_From <- xt1[,2]
    # n_To <- xt1[,4]
    
    xt1 <- get.edgelist(gTmp)
    n_From <- xt1[,1]
    n_To <- xt1[,2]
    n_Source <- setdiff(n_From, n_To)
    n_Sink <- setdiff(n_To, n_From)
    #
    rec_set1 <- n_Source
    tf_set1 <- n_Sink
    n_rec1 <- length(rec_set1)
    n_tf1 <- length(tf_set1)
    #
    if (min(n_rec1, n_tf1) < 1){  # no possible Receptor - TF path
      next
    }
    # print(j1)
    for (j2 in 1:n_rec1){  # each source point
      # print("j2"); print(j2)
      pf1 <- rec_set1[j2]  # from point
      paths <- get.shortest.paths(gTmp, pf1, tf_set1, mode='out')
      # paths <- all_simple_paths(gTmp, pf1, tf_set1, mode='out')
      paths <- paths$vpath
      nPath <- length(paths)
      #
      if (nPath > 0){  # connected
        for (k in 1:nPath){  # all paths
          pt <- paths[[k]]
          pt <- pt$name
          nPt <- length(pt)
          if (nPt < 2){  # no path (dis == Inf)
            next
          }
          # for each signaling path
          pathTmp1 <- c('Source', 'Target', 'pathwayName')
          for (l in 1:(nPt-1)){
            vt <- c(pt[l], pt[l+1], p_name)
            pathTmp1 <- rbind(pathTmp1, vt)
          }
          # activation evaluation
          if (nrow(pathTmp1) == 2){
            eTmp1 <- pathTmp1[-1,]
            dim(eTmp1) <- c(1,3)
          } else {
            eTmp1 <- pathTmp1[-1,]
          }
          #
          # st <- getScorePathwayKeggV1a(eTmp1, KeggPathways[[j1]], zScore, 0.0)
          st <- getScorePathwayKeggV1b(eTmp1, KeggPathways[[j1]], zScore, 0.0)
          xt <- list()
          xt[[1]] <- st
          xt[[2]]	<- eTmp1
          n_pt <- (n_pt + 1)
          pathTmp[[n_pt]] <- xt
        }
      }
    }
  }	# end of each signaling pathways?
  # pathways[[j]] <- pathTmp  # pathways of the j-th patient
  res1 <- list(pathTmp=pathTmp, lPath=lPath)
  return(res1)
}

# --------------------
# getScorePathwayKeggV1b
# --------------------
getScorePathwayKeggV1b <- function(edgePath, eType, zScore, alpha1){
  # edgePath: nx2 dimension; each row is an edge
  # eType: edge type, e.g., activation, inhibition,
  # zScore: gene expression of genes, e.g., fold change, z-score, protein expression
  # alpha1: weight of the contribution of the parent node
  vScore <- 0.0  # pathway score initialization
  aScore <- 0.0  # activation score
  # pre-processing
  gSym <- names(zScore)
  nEdge <- nrow(edgePath)
  eType1 <- rep(1.0, nEdge)  # 1.0: activation; -1.0: inhibition
  eScore <- matrix(0.0, nEdge, 2)  # scores of nodes on each edge
  #
  eS <- as.character(eType[,2])  # source nodes
  eT <- as.character(eType[,4])  # target nodes
  for (i in 1:nEdge){
    et <- edgePath[i,]
    # construct the eType1
    idxt <- which((eS %in% et[1]) & (eT %in% et[2]))
    if (length(idxt)>0){
      et1 <- as.character(eType[idxt[1], 6])
      # print(et1)
      if (et1 == 'Process(inhibition)'){
        eType1[i] <- -1.0
      }
    }
    # construct the eScore
    idxt <- which(gSym %in% et[1])
    if (length(idxt)>0){
      eScore[i,1] <- zScore[idxt]
    }
    #
    idxt <- which(gSym %in% et[2])
    if (length(idxt)>0){
      eScore[i,2] <- zScore[idxt]
    }
  }
  # the first node has no parent node
  eType1 <- c(1, eType1)
  # calculate the score of the signaling pathway
  for (i in 1:nEdge){
    vt <- 0.5*abs((eScore[i,1]+eType1[i+1]*eScore[i,2]))  # abs
    vScore <- vScore + vt
    #
    vt <- 0.5*(eScore[i,1]+eScore[i,2])  #
    aScore <- aScore + vt
  }
  vScore <- vScore + 0.5*eType1[i+1]*eScore[i,2] + 0.5*eScore[1,1]  # the last node and first node score
  vScore <- vScore/(nEdge+1)
  #
  aScore <- aScore + 0.5*eScore[i,2] + 0.5*eScore[1,1]  # the last node and first node score
  aScore <- aScore/(nEdge+1)
  # print('testing')
  # if (aScore != vScore) {print(c(aScore, vScore))}
  return(c(vScore,aScore))
}




# sub-functions:
generateSigNetwork_v3a <- function(G_t, gs1a, gs2a){ # 
  #
  library(igraph)
  
  ig_t0 <- G_t[[2]]
  ig_t1 <- G_t[[1]]
  gene_g_t <- V(ig_t0)$name
  
  gs1b <- intersect(gs1a, gene_g_t) # root genes; remove genes are not in the STRING databse
  gs2b <- intersect(gs2a, gene_g_t) # remove genes are not in the
  
  g_net2A <- c('Source', 'Target') # initilization
  V2 <- gs1b # initilization nodes included 
  n_gs1 <- 2
  
  gs2c <- setdiff(gs2b, gs1b)
  V0 <- gs2c
  
  if (n_gs1 > 0){
    for (i_gs in 1:1){ # 
      #
      g_net <- c('Source', 'Target')
      
      flag1 <- 1
      V1 <- V0
      while (flag1 > 0){
        gs_1 <- V1  
        paths_1v <- shortest.paths(ig_t1, V2, gs_1)
        v_t1 <- paths_1v
        
        if (min(v_t1) == Inf){
          print('no reachable targets')
          break; # no reachable targets
        }
        
        idx_t1 <- which(v_t1 == min(v_t1))
        
        if (length(V2) == 1){			
          paths_0 <- get.shortest.paths(ig_t0, V2, gs_1[idx_t1[1]])
        } else {
          idx_t1a <- ceiling(idx_t1[1]/length(V2))
          idx_t1b <- idx_t1[1] %% length(V2)
          if (idx_t1b == 0) {idx_t1b = length(V2)}
          paths_0 <- get.shortest.paths(ig_t0, V2[idx_t1b], gs_1[idx_t1a])
          
        }
        paths_0 <- paths_0$vpath
        # choose the best path and add it into the network 
        pt <- names(paths_0[[1]])
        nPt <- length(pt)
        # convert to 2-column table/network 
        g_net_tmp <- c("tmp", "tmp")
        for (l in 1:(nPt-1)){
          vt <- c(pt[l], pt[l+1])
          g_net_tmp <- rbind(g_net_tmp, vt)
        }
        V_tmp <- union(g_net_tmp[,1], g_net_tmp[,2])
        V_tmp <- V_tmp[-which(V_tmp == "tmp")]
        
        g_net <- rbind(g_net, g_net_tmp[-1,])
        V2 <- union(V2, V_tmp)
        idx_t <- which(V0 %in% V2)
        
        V1 <- V0[-idx_t]
        g_net <- unique(g_net)
        
        flag1 <- length(V1)
      }
      #
      if (length(g_net) > 3){
        g_net1 <- g_net[-1,]
      }else{
        return(g_net2A)
      }
      g_net1 <- unique(g_net1)
      # remove (a,b)-(b,a) edge pairs 
      n_e <- dim(g_net1)[1]
      g_net2 <- g_net1[1,]
      for (i in 2:n_e){
        et <- g_net1[i,]
        net_ta <- rbind(g_net2, et)
        net_tb <- rbind(g_net2, c(et[2], et[1]))
        net_ta <- unique(net_ta)
        net_tb <- unique(net_tb)
        x_ta <- dim(net_ta)[1]
        x_tb <- dim(net_tb)[1]
        if (x_ta <= x_tb){
          g_net2 <- net_ta
        } else {
          g_net2 <- net_tb
        }
      }
      # merge network starting from different nodes
      g_net2A <- rbind(g_net2A, g_net2)
      g_net2A <- unique(g_net2A)
      
    } # end of loop gs1 nodes
  }
  
  g_net2A <- g_net2A[-1,]
  return(g_net2A)
  
}
#
###
#
#
generateSigNetwork_setNetWeight <- function(g_t, gene_symbol_0, fc_t0){ # 
  #
  # gene_symbol_0 <- symbol_2
  # fc_t0 <- 2.0^fc_2
  # gs1a <- gs_t1c
  # gs2a <- gs_t1d
  fc_t0<-exp(fc_t0)
  library(igraph)
  #
  g_t <- cbind(as.character(g_t[,1]), as.character(g_t[,2]))
  g_t <- g_t[-which(g_t[,1] == g_t[,2]),]  # remove the self-connection
  # gene_g_t <- union(g_t[,1], g_t[,2])  # 
  
  ig_t0 <- graph.edgelist(g_t, directed=F)
  ig_t1 <- graph.edgelist(g_t, directed=F)
  #
  gene_g_t <- V(ig_t0)$name
  
  # set the edge weights 
  # case - 1: prefer/biase to up-regulated genes 
  idx_t <- which(gene_symbol_0 %in% gene_g_t)
  gene_symbol_2a <- gene_symbol_0[idx_t]
  fc_t2a <- fc_t0[idx_t]
  
  V0 <- gene_symbol_2a  # common genes 
  n_V0 <- length(V0)
  
  fc_t0a <- fc_t2a  # fold change values
  idx_t <- order(fc_t0a, decreasing=T)
  fc_t0a <- fc_t0a[idx_t]
  V0 <- V0[idx_t]
  
  fc_t0a[fc_t0a <= 0.5] <- 0.5  # ignore the down-regulated genes 
  fc_t0a[fc_t0a >=5.0] <- 5.0  # trucate the extrem values  
  
  e_all <- E(ig_t1)
  e_head <- names(head_of(ig_t1, e_all))
  e_tail <- names(tail_of(ig_t1, e_all))
  #
  e_w1 <- rep(1.0, length(e_all))  # initialization 
  vht <- rep(1.0, length(e_all))
  vtt <- rep(1.0, length(e_all))
  #
  
  for (i in 1:n_V0){  # update the edge; assign the fold change to the nodes
    vt <- fc_t0a[i]
    idx_t <- which(e_head %in% V0[i])
    vht[idx_t] <- vt
    idx_t <- which(e_tail %in% V0[i])
    vtt[idx_t] <- vt
  }
  #
  e_w1 <- 2.0/(vht + vtt)
  # set_edge_attr("weight", value = 1:10)
  ig_t1 <- set.edge.attribute(ig_t1, "weight", value=e_w1)
  # gs1 <- gs1[which(gs1 %in% gene_g_t)]; #print(length(gs1))
  
  net <- list()
  net[[1]] <- ig_t1
  net[[2]] <- ig_t0
  
  return(net)
}
#
###
#
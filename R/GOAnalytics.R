GOAnalytics<-function(testData,normalData,gene_symbol,upFc_thres,dnFc_thres,pv_thres,group=1,outputDir,progress){
  meanTest <-rowMeans(testData)
  meanNormal<-rowMeans(normalData)

  progress$set(0.2,detail="GO analytics")
  fcTN <- meanTest-meanNormal
  gs_up <- gene_symbol[which(fcTN >= upFc_thres)]
  gs_dn <- gene_symbol[which(fcTN <= dnFc_thres&fcTN)]
  # get the enriched GOs
  GO_up <- goEnrichAnalysisV1a(gs_up, gene_symbol)
  progress$set(0.45,detail="GO analytics")
  GO_dn <- goEnrichAnalysisV1a(gs_dn, gene_symbol)
  progress$set(0.7,detail="GO analytics")

  
  
  #
  idx_up <- which(GO_up[,3] == "BP" & as.numeric(GO_up[,4]) <= 1000 & as.numeric(GO_up[,4]) >= 10 & GO_up[,5] <= pv_thres)
  idx_dn<- which(GO_dn[,3] == "BP" & as.numeric(GO_dn[,4]) <= 1000 & as.numeric(GO_dn[,4]) >= 10 & GO_dn[,5] <= pv_thres)
  
  GO_up <- GO_up[idx_up,]
  GO_dn <- GO_dn[idx_dn,]
  
  netGeneGO_up <- getGeneGoNetworkV1a(gs_up, GO_up[,1])
  netGeneGO_dn <- getGeneGoNetworkV1a(gs_dn, GO_dn[,1])
  
  list(GO_up,GO_dn,netGeneGO_up,netGeneGO_dn)
}


# sub-functions ---
goEnrichAnalysisV1a <- function(gs1, all_gene){
  # convert gene symbol to entrezID
  all_gene <- mapIds(org.Hs.eg.db, all_gene, 'ENTREZID', 'SYMBOL')
  all_gene <- all_gene[-which(is.na(all_gene))]
  
  gs1 <- mapIds(org.Hs.eg.db, gs1, 'ENTREZID', 'SYMBOL')
  gs1 <- gs1[-which(is.na(gs1))]
  
  gs2 <- setdiff(all_gene, gs1)
  # check all the GO terms related to these genes
  go_t1 <- mget(gs1, org.Hs.egGO)
  ids_go_1<-unique(as.vector(unlist(lapply(go_t1,names))))
  #nt <- length(go_t1)
  #  ids_go_1 <- c()
  # for (i in 1:nt){
  #    xt1 <- go_t1[[i]]
  #   ids_go_1 <- union(ids_go_1, names(xt1))
  #  }
  
  #
  N_gs1 <- length(gs1)
  N_gs2 <- length(gs2)
  xt <- matrix(1.0, 2, 2)
  nt <- length(ids_go_1)
  pv_go_1 <- matrix(1.0, nt, 4)
  for (i in 1:nt){
    go_idt <- ids_go_1[i]
    gs_t1 <- get(go_idt, org.Hs.egGO2ALLEGS)
    # gs_t2 = unlist(mget(gs_t1,org.Hs.egSYMBOL))
    # fisher test to identify the enriched GOs
    xt[1,1] <- length(intersect(gs1, gs_t1))
    xt[1,2] <- length(intersect(gs2, gs_t1))
    xt[2,1] <- N_gs1 - xt[1,1]
    xt[2,2] <- N_gs2 - xt[1,2]
    pv_go_1[i,3] <- length(gs_t1)
    pv_go_1[i,1] <- Term(GOTERM[[go_idt]])
    pv_go_1[i,2] <- Ontology(GOTERM[[go_idt]])
    pv_go_1[i,4] <- fisher.test(xt)$p.value
  }
  
  pv_go_1 <- cbind(ids_go_1, pv_go_1)
  #
  colnames(pv_go_1) <- c("go_id", "go_term_name", 'go_ontology', '#ofgenes_go', 'p-value')
  return(pv_go_1)
}
#
###
#
getGeneGoNetworkV1a <- function(gs0, GO1){
  #
  library(org.Hs.eg.db)
  library(GO.db)
  
  gs1 <- mapIds(org.Hs.eg.db, gs0, 'ENTREZID', 'SYMBOL')
  idxt <- which(is.na(gs1))
  gs1 <- gs1[-idxt]
  gs0 <-gs0[-idxt]
  
  nt <- length(GO1)
  geneGoNet <- c("GO", "Gene")
  for (i in 1:nt){
    go_idt <- GO1[i]
    gs_t1 <- get(go_idt, org.Hs.egGO2ALLEGS)
    #
    gs_t2 <- intersect(gs1, gs_t1)
    idxt1 <- which(gs1 %in% gs_t2)
    gs_t3 <- gs0[idxt1]
    
    nt1 <- length(gs_t2)
    if (nt1 > 0){
      xt1 <- rep(go_idt, nt1)
      xt2 <- gs_t3
      geneGoNet <- rbind(geneGoNet, cbind(xt1,xt2))
    }
    
  }
  return(geneGoNet[-1,])
  #
}



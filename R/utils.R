
# 
# rescale<-function(seurat_data,regress=T){
#   if(regress){
#     seurat_data<-CellCycleScoring(seurat_data,s.features=cc.genes$s.genes,g2m.features =cc.genes$g2m.genes, set.ident = T)
#     seurat_data<-ScaleData(seurat_data,vars.to.regress = c("S.Score","G2M.Score"),features = rownames(seurat_data))
#   }else{
#     seurat_data<-ScaleData(seurat_data,features=rownames(seurat_data))
#   }
#   seurat_data
# }

#max-min normalization
normalization<-function(data,base=2){
  gene_max<-apply(data,base,max)
  center <- sweep(data, base, apply(data,base,min)) 
  R <- gene_max - apply(data,base,min)
  #prevent Nan
  R[R==0]=0.00001
  df_marker_normal<- sweep(center, base, R, "/") 
  df_marker_normal
}

findSD<-function(seurat_data){
  df<-as.matrix(seurat_data[[seurat_data@active.assay]]@data[VariableFeatures(seurat_data),])
  df
  
}

#data table transposition(row to column)
transposition<-function(df){
  df_trans<-as.data.frame(t(as.matrix(df)))
}


# -Please download (add the function here) and cite the code from: https://rdrr.io/bioc/EnrichmentBrowser/src/R/GSEA.1.0.R 
# GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = 0, correl.vector = NULL) {
#
GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = 0, correl.vector = NULL){  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. It is the same 
  # calculation as in GSEA.EnrichmentScore but faster (x8) without producing the 
  # RES, arg.RES and tag.indicator outputs.
  # This call is intended to be used to asses the enrichment of random 
  # permutations rather than the observed one.
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). 
  # When the score type is 1 or 2 it is necessary to input the correlation vector
  # with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list 
  #       (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set 
  #       (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: 
  #       weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #   correl.vector: A vector with the coorelations (e.g. signal to noise scores)
  #       corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  #
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  peak.res.vector <- valley.res.vector <- 
    tag.diff.vector <- vector(length=Nh, mode="numeric")
  
  loc.vector <- vector(length=N, mode="numeric")
  loc.vector[gene.list] <- seq_len(N)
  tag.loc.vector <- loc.vector[gene.set]
  tag.loc.vector <- sort(tag.loc.vector, decreasing =FALSE)
  
  if (weighted.score.type == 0) tag.correl.vector <- rep(1, Nh)
  else if (weighted.score.type == 1) 
  {
    tag.correl.vector <- correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } 
  else if (weighted.score.type == 2) 
  {
    tag.correl.vector <- 
      correl.vector[tag.loc.vector] * correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } 
  else 
  {
    tag.correl.vector <- correl.vector[tag.loc.vector] * weighted.score.type
    tag.correl.vector <- abs(tag.correl.vector)
  }
  
  norm.tag <- 1.0/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1.0/Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
  tag.diff.vector[2:Nh] <- 
    tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  
  return(list(ES = ES))
}



#convert cell-cell communication table to json
networkDataToJson<-function(network_nodes,network_edges){
  a<-toJSON(network_nodes)
  b<-toJSON(network_edges)
  d<-str_c('{"nodes":',a,',"links":',b,'}')
  final_data<-fromJSON(d,simplifyDataFrame = FALSE)
}

#convert bubble table to json
bubbleDataToJson<-function(bubbleData){
  toJSON(bubbleData)
}

#calculate cell distribution for every group
cellDistribution<-function(type_list,group_list){
  unique_type<-as.character(unique(type_list))
  num_type<-length(unique_type)
  if(is.null(group_list)){
    cell_count<-c()
    total<-length(type_list)
    for(i in 1:num_type){
      num_cell<-sum(type_list==unique_type[i])/total
      cell_count<-c(cell_count,num_cell)
    }
    cell_count_df<-data.frame(count=cell_count,type=unique_type)
    return(cell_count_df)
  }else{
    unique_group<-unique(group_list)
    num_group<-length(unique_group)
    cell_count_df<-c()
    for (i in 1:num_group){
      cell_count<-c()
      for(j in 1:num_type){
        total<-sum(group_list==unique_group[i])
        num_cell<-sum(type_list==unique_type[j]&group_list==unique_group[i])/total
        cell_count<-c(cell_count,num_cell)
      }
      if(i==1){
        cell_count_df<-cell_count
      }else{
        cell_count_df<-rbind(cell_count_df,cell_count)
      }
    }
    cell_count_df<-as.data.frame(cell_count_df)
    cell_count_df<-cbind(cell_count_df,unique_group)
    colnames(cell_count_df)<-c(as.character(unique_type),"group")
    melt_df<-melt(as.data.table(cell_count_df),id.vars="group")
  }
}


GONetworkGenerating<-function(GO_dn,GO_up,GO,GO_name){
  selected_GO_up<-GO_up[GO_up[,1]==GO,]
  selected_GO_dn<-GO_dn[GO_dn[,1]==GO,]
  selected_GO<-as.data.frame(rbind(selected_GO_up,selected_GO_dn))
  colnames(selected_GO)<-c("source","target")
  GO_nodes<-GO
  GO_type<-rep("GO",length(GO_nodes))
  gene_nodes<-as.character(unique(selected_GO[,2]))
  gene_names<-gene_nodes
  gene_type<-rep("gene",length(gene_nodes))
  nodes<-c(GO_nodes,gene_nodes)
  types<-c(GO_type,gene_type)
  names<-c(GO_name,gene_names)
  selected_GO_nodes<-data.frame(node=nodes,name=names,type=types)
  networkData<-networkDataToJson(selected_GO_nodes,selected_GO)
}


#if network have more than 15 pathways, filter some pathwayout
processActivateNetworkDisplay<-function(network_nodes,network_edges){
  important_pathway_list<-c(
    "Adipocytokine signaling pathway",
    "MAPK signaling pathway",
    "Fc epsilon RI signaling pathway",
    "TGF-beta signaling pathway",
    "T cell receptor signaling pathway",
    "ErbB signaling pathway",
    "Sphingolipid metabolism",
    "VEGF signaling pathway",
    "B cell receptor signaling pathway",
    "Oxytocin signaling pathway",
    "Ras signaling pathway",
    "FoxO signaling pathway",
    "Pancreatic secretion",
    "Calcium signaling pathway",
    "Cell cycle",
    "Chemokine signaling pathway",
    "Cytosolic DNA-sensing pathway" ,
    "ErbB signaling pathway",
    "Estrogen signaling pathway" ,
    "GnRH signaling pathway" ,
    "Hedgehog signaling pathway" ,
    "HIF-1 signaling pathway" ,
    "Hippo signaling pathway",
    "Insulin signaling pathway",
    "Jak-STAT signaling pathway",
    "mTOR signaling pathway",
    "Neurotrophin signaling pathway",
    "NF-kappa B signaling pathway" ,
    "NOD-like receptor signaling pathway" ,
    "Notch signaling pathway",
    "p53 signaling pathway" ,
    "PI3K-Akt signaling pathway" ,
    "Prolactin signaling pathway" ,
    "cGMP-PKG signaling pathway" ,
    "Rap1 signaling pathway",
    "TNF signaling pathway",
    "Wnt signaling pathway" ,
    "Toll-like receptor signaling pathway" ,
    "Thyroid hormone signaling pathway" ,
    "RIG-I-like receptor signaling pathway" ,
    "Retrograde endocannabinoid signaling"  
  )
  unique_pathway<-unique(as.character(network_nodes$parent))
  if(length(unique_pathway)>15){
    selected_pathway<-unique_pathway[unique_pathway%in%important_pathway_list]
    if(length(selected_pathway)>15){
      selected_pathway<-selected_pathway[c(1:15)]
    }
    new_network_nodes<-network_nodes[network_nodes$parent%in%selected_pathway,]
    new_network_edges<-network_edges[network_edges$PathwayName%in%selected_pathway,]
  }else{
    new_network_nodes<-network_nodes
    new_network_edges<-network_edges
  }
  return(list(new_network_nodes,new_network_edges))
}

processDrugClusteringDisplay<-function(drug_node,drug_network){
  num_cluster<-sum(drug_node$type=="exemplar")
  if(num_cluster>10){
    selected_clusters<-drug_node[which(drug_node$type=="exemplar")[c(1:10)],]$node
    new_drug_network<-drug_network[drug_network$source%in%selected_clusters|drug_network$target%in%selected_clusters,]
    selected_nodes<-unique(c(new_drug_network$source,new_drug_network$target))
    new_drug_node<-drug_node[drug_node$node%in%selected_nodes,]
  }else{
    new_drug_node<-drug_node
    new_drug_network<-drug_network
  }
  return(list(new_drug_node,new_drug_network))
}



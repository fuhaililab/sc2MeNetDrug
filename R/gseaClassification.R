#calculate gene set enrichment score for every cell type in every group 
#based on marker gene information user selected
#choose the cell type corresponding to the maximum es in each group as the 
#classification result
gseaClassification<-function(df,tsne_cluster,gene_symbol,marker_gene,progress,outputDir){
  rownames(marker_gene)<-marker_gene[,1]
  marker_gene<-marker_gene[,-1]
  cluster_num<-max(as.numeric(tsne_cluster))

  for(j in(1:cluster_num)){
    progress$set((j/cluster_num)*0.7+0.1, detail ="Compute GSEA score")
    if(is.null(dim(df[,tsne_cluster==j]))){
      cluster_mean<-df[,tsne_cluster==j]
    }else{
      cluster_mean<-rowMeans(df[,tsne_cluster==j])
    }
    other_mean<-rowMeans(df[,tsne_cluster!=j])
    fold_change<-cluster_mean-other_mean
    
   # cluster_normal_mean<-rowMeans(df_normal[,tsne_cluster==j])
    # pv_list<-rep(1.0,num_gene)
    # for(i in 1:num_gene){
    #   print(i)
    #   pv_list[i] <- wilcox.test(as.numeric(df[i,tsne_cluster==j]),as.numeric(df[i,tsne_cluster!=j]))$p.value
    # }
    # pv_list<-p.adjust(pv_list,method = "BH")
    # summary<-data.frame(fc=fold_change,pv=pv_list)
    # save(summary,file="summary.RData")
    
    es_list<-c()
    gene_list<-names(cluster_mean)
#    gene.list.exp<-order(cluster_normal_mean,decreasing = T)
    gene.list.fc<-order(fold_change,decreasing = T)
    for(i in 1:ncol(marker_gene)){
      marker_genes<-rownames(marker_gene)[marker_gene[,i]==1]
      tryCatch(
        {
          gene.set<-intersect(marker_genes,gene_list)
          if(length(gene.set)>0){
            gene.set<-which(names(cluster_mean)%in%gene.set)
         #   ES_exp<-unlist(GSEA.EnrichmentScore2(gene.list.exp,gene.set))
            ES_fc<-unlist(GSEA.EnrichmentScore2(gene.list.fc,gene.set))
            ES<-ES_fc
            if(is.na(ES)==T){
              ES<--Inf
            }
          }else{
            ES<--Inf
          }
          es_list<-c(es_list,ES)
          #result<-summary(gsea(fold_change,gsets=gene_list,minGenes = 1,maxGenes = 20,center=TRUE,B=500))
          #print(result)
          #pvalue<-result[4]
          #fdr<-result[6]
          #if(pvalue<=1&fdr<=1){
          #  es<-result[2]
          #  es_list<-c(es_list,es)
        #  }else{
         #   es_list<-c(es_list,-Inf)
        #  }
        },
        error=function(e){
          print(e)
        },
        finally={
          if(length(es_list)==i-1){
            es_list<-c(es_list,-Inf)
          }
        }
      )
    }
    if(j==1){
      classification_result<-es_list
    }
    else{
      classification_result<-rbind(classification_result,es_list)
    }
  }
  
  colnames(classification_result)<-colnames(marker_gene)
  rownames(classification_result)<-c(1:cluster_num)
  annotation_es<-classification_result
  save(annotation_es,file=paste0(outputDir,"/annotation_es.RData"))
  type_list<-c()
  cluster_type_list<-c()
  for(k in 1:cluster_num){
    max_es<-max(classification_result[k,])
    
    if(max_es>0){
      cell_type<-names(which(classification_result[k,]==max_es))
      cluster_type<-paste(k,cell_type,sep="-")
      cluster_type_list<-c(cluster_type_list,cluster_type)
      type_list<-c(type_list,cell_type)
    }else{
      cluster_type<-paste(k,"unknown",sep="-")
      cluster_type_list<-c(cluster_type_list,cluster_type)
      type_list<-c(type_list,"unknown")
    }
  }
  
  final_result<-cbind(c(1:cluster_num),type_list,cluster_type_list)
  outlier<-c(0,"outlier","0-outlier")
  final_result<-rbind(final_result,outlier)
  colnames(final_result)<-c("cluster","type","cluster_type")
  final_result<-as.data.frame(final_result)
  
}


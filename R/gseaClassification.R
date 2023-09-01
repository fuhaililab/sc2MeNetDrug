

#' Computes the weighted GSEA score of gene.set in gene.list. It is the same
#' calculation as in GSEA.EnrichmentScore but faster (x8) without producing the
#' RES, arg.RES and tag.indicator outputs.
#' This call is intended to be used to asses the enrichment of random
#' permutations rather than the observed one.
#' The weighted score type is the exponent of the correlation
#' weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted).
#' When the score type is 1 or 2 it is necessary to input the correlation vector
#' with the values in the same order as in the gene list.
#
#'
#' @param gene.list The ordered gene list (e.g. integers indicating the
#'                  original position in the input dataset)
#' @param gene.set A gene set (e.g. integers indicating the location of
#'                  those genes in the input dataset)
#' @param weighted.score.type Type of score: weight: 0
#'    (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)
#' @param correl.vector A vector with the coorelations (e.g. signal to noise scores)
#                   corresponding to the genes in the gene list
#'
#' @return ES: Enrichment score (real number between -1 and +1)
#' @export
#'
#' @examples
GSEA.EnrichmentScore2 <-
  function(gene.list,
           gene.set,
           weighted.score.type = 0,
           correl.vector = NULL) {
    N <- length(gene.list)
    Nh <- length(gene.set)
    Nm <-  N - Nh
    peak.res.vector <- valley.res.vector <-
      tag.diff.vector <- vector(length = Nh, mode = "numeric")
    
    loc.vector <- vector(length = N, mode = "numeric")
    loc.vector[gene.list] <- seq_len(N)
    tag.loc.vector <- loc.vector[gene.set]
    tag.loc.vector <- sort(tag.loc.vector, decreasing = FALSE)
    
    if (weighted.score.type == 0)
      tag.correl.vector <- rep(1, Nh)
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
      tag.correl.vector <-
        correl.vector[tag.loc.vector] * weighted.score.type
      tag.correl.vector <- abs(tag.correl.vector)
    }
    
    norm.tag <- 1.0 / sum(tag.correl.vector)
    tag.correl.vector <- tag.correl.vector * norm.tag
    norm.no.tag <- 1.0 / Nm
    tag.diff.vector[1] <- (tag.loc.vector[1] - 1)
    tag.diff.vector[2:Nh] <-
      tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
    tag.diff.vector <- tag.diff.vector * norm.no.tag
    peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
    valley.res.vector <- peak.res.vector - tag.correl.vector
    max.ES <- max(peak.res.vector)
    min.ES <- min(valley.res.vector)
    ES <-
      signif(ifelse(max.ES > -min.ES, max.ES, min.ES), digits = 5)
    
    return(list(ES = ES))
  }


#' calculate gene set enrichment score for every cell type in every group
#' based on marker gene information user selected
#' choose the cell type corresponding to the maximum es in each group as the
#' classification result
#'
#' @param df The scRNA expression matrix.
#' @param tsne_cluster Clustering result table.
#' @param gene_symbol Gene symbol list.
#' @param marker_gene Selected cell type and corresponding marker genes.
#' @param progress Progress bar for display.
#' @param outputDir Directory for saving results.
#'
#' @return
#' @export
#'
#' @examples
gseaClassification <-
  function(df,
           tsne_cluster,
           gene_symbol,
           marker_gene,
           progress,
           outputDir) {
    unique_cell_type <- unique(marker_gene$cell_type)
    cluster_num <- max(as.numeric(tsne_cluster))
    
    for (j in(1:cluster_num)) {
      progress$set((j / cluster_num) * 0.7 + 0.1, detail = "Compute GSEA score")
      if (is.null(dim(df[, tsne_cluster == j]))) {
        cluster_mean <- df[, tsne_cluster == j]
      } else{
        cluster_mean <- rowMeans(df[, tsne_cluster == j])
      }
      other_mean <- rowMeans(df[, tsne_cluster != j])
      fold_change <- cluster_mean - other_mean
      
      
      es_list <- c()
      gene_list <- names(cluster_mean)
      gene.list.fc <- order(fold_change, decreasing = T)
      for (i in 1:length(unique_cell_type)) {
        cell_type <- unique_cell_type[i]
        marker_genes <-
          marker_gene[marker_gene[, 1] == cell_type, 2]
        tryCatch({
          gene.set <- intersect(marker_genes, gene_list)
          if (length(gene.set) > 0) {
            gene.set <- which(names(cluster_mean) %in% gene.set)
            ES_fc <-
              unlist(GSEA.EnrichmentScore2(gene.list.fc, gene.set))
            ES <- ES_fc
            if (is.na(ES) == T) {
              ES <- -Inf
            }
          } else{
            ES <- -Inf
          }
          es_list <- c(es_list, ES)
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
        error = function(e) {
          print(e)
        },
        finally = {
          if (length(es_list) == i - 1) {
            es_list <- c(es_list, -Inf)
          }
        })
      }
      if (j == 1) {
        classification_result <- es_list
      }
      else{
        classification_result <- rbind(classification_result, es_list)
      }
    }
    
    colnames(classification_result) <- unique_cell_type
    rownames(classification_result) <- c(1:cluster_num)
    annotation_es <- classification_result
    save(annotation_es, file = paste0(outputDir, "/annotation_es.RData"))
    type_list <- c()
    cluster_type_list <- c()
    for (k in 1:cluster_num) {
      max_es <- max(classification_result[k, ])
      
      if (max_es > 0) {
        cell_type <- names(which(classification_result[k, ] == max_es)[1])
        cluster_type <- paste(k - 1, cell_type, sep = "-")
        cluster_type_list <- c(cluster_type_list, cluster_type)
        type_list <- c(type_list, cell_type)
      } else{
        cluster_type <- paste(k - 1, "unknown", sep = "-")
        cluster_type_list <- c(cluster_type_list, cluster_type)
        type_list <- c(type_list, "unknown")
      }
    }
    
    final_result <-
      cbind(c(0:(cluster_num - 1)), type_list, cluster_type_list)
    colnames(final_result) <- c("cluster", "type", "cluster_type")
    final_result <- as.data.frame(final_result)
    
  }

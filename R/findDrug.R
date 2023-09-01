#' target drug discovery function using drugbank.
#'
#' @param network_nodes The node data in the computed inter-cell signaling network.
#' @param network_edges The edge data in the computed inter-cell signaling network.
#' @param drugBankInteraction DrugBank drug-target information.
#' @param drugBankInformation DrugBank drug information.
#'
#' @return
#' @export
#'
#' @examples
findDrug <-
  function(network_nodes,
           network_edges,
           drugBankInteraction,
           drugBankInformation) {
    nodes_name <- as.character(network_nodes$name)
    drug_signal <- c()
    num_drug_list <- c()
    drug_gene_set <- c()
    drug_list <- c()
    drug_set_for_gene <- c()
    for (i in 1:length(nodes_name)) {
      node <- nodes_name[i]
      node_drug_list <-
        drugBankInteraction[which(drugBankInteraction$target == node), ]
      
      if (nrow(node_drug_list) != 0) {
        num_drug <- nrow(node_drug_list)
        num_drug_list <- c(num_drug_list, num_drug)
        drug_signal <- c(drug_signal, "haveDrug")
        drug_set_for_gene <-
          c(drug_set_for_gene,
            paste(node_drug_list$drug_name, collapse = ",  "))
        for (j in 1:num_drug) {
          drug <- as.character(node_drug_list[j, 1])
          drug_gene_list <-
            paste(drugBankInteraction[drugBankInteraction$drug_name == drug, 2], collapse =
                    ", ")
          drug_gene_set <- c(drug_gene_set, drug_gene_list)
          drug_list <- c(drug_list, drug)
        }
      } else{
        drug_signal <- c(drug_signal, "noDrug")
        drug_set_for_gene <- c(drug_set_for_gene, " ")
        num_drug_list <- c(num_drug_list, 0)
      }
    }
    drug_table <-
      data.frame(
        drug_name = drug_list,
        drug_gene_set = drug_gene_set,
        stringsAsFactors = F
      )
    drug_table <- unique(drug_table)
    drug_table <-
      merge(
        drug_table,
        drugBankInformation,
        by.x = "drug_name",
        by.y = "drug_name",
        all.x = T
      )
    new_network_nodes <-
      cbind(
        network_nodes,
        drugSignal = drug_signal,
        drugList = drug_set_for_gene,
        numDrug = num_drug_list
      )
    final_data <- networkDataToJson(new_network_nodes, network_edges)
    
    drug_smiles_list <- as.character(drug_table$drug_smile)
    names(drug_smiles_list) <- as.character(drug_table$drug_name)
    drugClusteringResult <-
      drugClusteringChemicalStructure(drug_smiles_list)
    
    list(
      drug_table,
      final_data,
      drugClusteringResult[[1]],
      drugClusteringResult[[2]],
      drugClusteringResult[[3]]
    )
  }

#' Drug clustering function
#'
#' @param drug_list Drug list for the clustering.
#'
#' @return
#' @export
#'
#' @examples
drugClusteringChemicalStructure <- function(drug_list) {
  num_drug <- length(drug_list)
  drug_name_list <- names(drug_list)
  mol_list <- parse.smiles(drug_list, kekulise = TRUE)
  fps <- lapply(mol_list, get.fingerprint, type = 'extended')
  fp.sim <- fp.sim.matrix(fps, method = 'tanimoto')
  rownames(fp.sim) <- drug_name_list
  colnames(fp.sim) <- drug_name_list
  drugCluster <- apcluster(fp.sim, seed = 123)
  exemplars <- labels(drugCluster, "exemplars")
  drug_index <- c(1:num_drug)
  drug_network <- data.frame(source = drug_index, target = exemplars)
  exemplars_names <- labels(drugCluster, "names")
  unique_exemplars <- unique(exemplars_names)
  drug_type <-
    sapply(drug_name_list, function(x) {
      ifelse(x %in% unique_exemplars, "exemplar", "normal")
    })
  drug_node <-
    data.frame(node = drug_index,
               type = drug_type,
               name = drug_name_list)
  return(list(drugCluster, drug_node, drug_network))
}



#' Signature based drug discovery function.
#'
#' @param network_nodes The node data in the computed inter-cell signaling network.
#' @param gSymZs Gene symbol list.
#' @param rankMatrix Drug rank matrix processed from the drug file.
#' @param drug_mapping Drug mapping data.
#' @param topDrug Number of top drug.
#' @param useFDAOnly If true, only include drug that been proved by the FDA.
#' @param drugBank DrugBank database.
#'
#' @return
#' @export
#'
#' @examples
findDrug2 <-
  function(network_nodes,
           gSymZs,
           rankMatrix,
           drug_mapping,
           topDrug,
           useFDAOnly,
           drugBank) {
    FDAdrugs <- unique(drugBank$drug_name)
    #fint top drug for two cell type
    n.probe.all <- dim(rankMatrix)[1]
    n.drg <- dim(rankMatrix)[2]
    
    #down_index1<-fc1<fc_downThres&pv1<=pv_thres
    geneSetUp <- unique(as.character(network_nodes$name))
    es_score <-
      getPatient2DrugReverseDis2(geneSetUp, NULL, gSymZs, rankMatrix, n.probe.all, n.drg)
    es_score_rank <-
      es_score[order(es_score, decreasing = F)] # up.ks - down.ks; then negative relationship
    name_list <- names(es_score_rank)
    name_list <-
      sapply(name_list, function(x) {
        substring(x, 1, nchar(x) - 3)
      })
    index <- c(1:length(es_score_rank))
    top <-
      data.frame(
        rankName = name_list,
        enrichment_score = as.vector(es_score_rank),
        drug_index = index
      )
    drug_table <-
      merge(top,
            drug_mapping,
            all.x = TRUE,
            by.x = "rankName",
            by.y = "rankName")
    drug_table <- drug_table[order(drug_table$drug_index), ]
    drug_table <- drug_table[, c(1, 2, 5, 10, 11)]
    rownames(drug_table) <- index
    if (useFDAOnly == 1) {
      drug_table <-
        drug_table[toupper(drug_table$pert_iname) %in% toupper(FDAdrugs), ]
      drug_table <- drug_table[c(1:topDrug), ]
      top_name <- drug_table$rankName
    } else{
      drug_table <- drug_table[c(1:topDrug), ]
      top_name <- drug_table$rankName
    }
    
    #drug clustering for two cell type
    drug_list <- colnames(rankMatrix)
    if (!is.null(drug_table)) {
      sub_drug_index <- which(drug_list %in% top_name)
      sub_rankMatrix <- rankMatrix[, sub_drug_index]
      drug_sim <-
        apply(sub_rankMatrix,
              2,
              cal_drug_es,
              gSymZs = gSymZs,
              rankMatrix = sub_rankMatrix)
      result <- apcluster(drug_sim, seed = 123)
      exemplars <- labels(result, "exemplars")
      drug_index <- c(1:topDrug)
      drug_network <- data.frame(source = drug_index, target = exemplars)
      exemplars_names <- labels(result)
      unique_exemplars <- unique(exemplars_names)
      drug_rank_names <- colnames(drug_sim)
      drug_names <- c()
      for (i in 1:length(drug_rank_names)) {
        infor <- drug_table[which(drug_table[, 1] == drug_rank_names[i]), ]
        drug_names <-
          c(drug_names, paste0(
            as.character(infor$pert_iname),
            ", ",
            as.character(round(infor[2], 3))
          ))
      }
      drug_type <-
        sapply(drug_rank_names, function(x) {
          ifelse(x %in% unique_exemplars, "exemplar", "normal")
        })
      drug_node <-
        data.frame(node = drug_index,
                   type = drug_type,
                   name = drug_names)
    } else{
      result <- NULL
      drug_node <- NULL
      drug_network <- NULL
    }
    
    
    list(drug_table, result, drug_node, drug_network)
  }

#' Compute signature drug score.
#'
#' @param x gene log-fold-change list.
#' @param gSymZs gene symbol list.
#' @param rankMatrix Drug rank matrix.
#'
#' @return
#' @export
#'
#' @examples
cal_drug_es <- function(x, gSymZs, rankMatrix) {
  x_order <- order(x)
  n.drg <- dim(rankMatrix)[2]
  geneSet_up <- gSymZs[x_order[1:150]]
  geneSet_down <- gSymZs[x_order[n.drg - 149:n.drg]]
  n.probe.all <- dim(rankMatrix)[1]
  
  n.drg <- dim(rankMatrix)[2]
  
  getPatient2DrugReverseDis2(geneSet_up,
                             geneSet_down,
                             gSymZs,
                             rankMatrix,
                             n.probe.all,
                             n.drg)
  
}



# -------------
getPatient2DrugReverseDis2 <-
  function(geneSetUp,
           geneSetDown,
           gSymZs,
           rankMatrix,
           n.probe.all,
           n.drg) {
    ###
    m.bigS <- matrix(-2.0, 1, n.drg)
    idxt <- 1:n.probe.all
    
    probe.set.up <- idxt[gSymZs %in% geneSetUp]
    
    if (length(probe.set.up) > 1 & sum(is.na(probe.set.up) == TRUE) < 1) {
      up.ks = apply(rankMatrix, 2, GSEA.EnrichmentScore2, gene.set = probe.set.up)
      up.ks <- unlist(up.ks) # newly added
    } else {
      up.ks = 0.0
      
    }
    probe.set.down <- idxt[gSymZs %in% geneSetDown]
    
    if (length(probe.set.down) > 1 &
        sum(is.na(probe.set.down) == TRUE) < 1) {
      down.ks = apply(rankMatrix, 2, GSEA.EnrichmentScore2, gene.set = probe.set.down)
      down.ks = unlist(down.ks)
      # newly added1
    } else {
      down.ks = 0.0
      
    }
    m.bigS = up.ks - down.ks
    return(m.bigS)
  }

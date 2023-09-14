#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Some Miscellaneous util functions.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert network data table to json.
#'
#' @param network_nodes Network node table.
#' @param network_edges Network edge table.
#'
#' @return
#' @export
#'
#' @examples
networkDataToJson <- function(network_nodes, network_edges) {
  a <- toJSON(network_nodes)
  b <- toJSON(network_edges)
  d <- str_c('{"nodes":', a, ',"links":', b, '}')
  final_data <- fromJSON(d, simplifyDataFrame = FALSE)
}

#' Convert marker gene bubble table to json.
#'
#' @param bubbleData data for the bubble plot.
#'
#' @return
#' @export
#'
#' @examples
bubbleDataToJson <- function(bubbleData) {
  toJSON(bubbleData)
}

#' Calculate cell distribution for every group.
#'
#' @param type_list Cell type list.
#' @param group_list Cell group list.
#'
#' @return
#' @export
#'
#' @examples
cellDistribution <- function(type_list, group_list) {
  unique_type <- as.character(unique(type_list))
  num_type <- length(unique_type)
  if (is.null(group_list)) {
    cell_count <- c()
    total <- length(type_list)
    for (i in 1:num_type) {
      num_cell <- sum(type_list == unique_type[i]) / total
      cell_count <- c(cell_count, num_cell)
    }
    cell_count_df <- data.frame(count = cell_count, type = unique_type)
    return(cell_count_df)
  } else{
    unique_group <- unique(group_list)
    num_group <- length(unique_group)
    cell_count_df <- c()
    for (i in 1:num_group) {
      cell_count <- c()
      for (j in 1:num_type) {
        total <- sum(group_list == unique_group[i])
        num_cell <-
          sum(type_list == unique_type[j] & group_list == unique_group[i]) / total
        cell_count <- c(cell_count, num_cell)
      }
      if (i == 1) {
        cell_count_df <- cell_count
      } else{
        cell_count_df <- rbind(cell_count_df, cell_count)
      }
    }
    cell_count_df <- as.data.frame(cell_count_df)
    cell_count_df <- cbind(cell_count_df, unique_group)
    colnames(cell_count_df) <- c(as.character(unique_type), "group")
    melt_df <- melt(as.data.table(cell_count_df), id.vars = "group")
  }
}


#' GO network generation function.
#'
#' @param GO_dn Down-regulated GOs table.
#' @param GO_up Up-regulated GOs table.
#' @param GO The GO terms.
#' @param GO_name All GO names.
#'
#' @return
#' @export
#'
#' @examples
GONetworkGenerating <- function(GO_dn, GO_up, GO, GO_name) {
  selected_GO_up <- GO_up[GO_up[, 1] == GO, ]
  selected_GO_dn <- GO_dn[GO_dn[, 1] == GO, ]
  selected_GO <- as.data.frame(rbind(selected_GO_up, selected_GO_dn))
  colnames(selected_GO) <- c("source", "target")
  GO_nodes <- GO
  GO_type <- rep("GO", length(GO_nodes))
  gene_nodes <- as.character(unique(selected_GO[, 2]))
  gene_names <- gene_nodes
  gene_type <- rep("gene", length(gene_nodes))
  nodes <- c(GO_nodes, gene_nodes)
  types <- c(GO_type, gene_type)
  names <- c(GO_name, gene_names)
  selected_GO_nodes <- data.frame(node = nodes,
                                  name = names,
                                  type = types)
  networkData <- networkDataToJson(selected_GO_nodes, selected_GO)
}


#' if generated signaling pathway network have more than 15 pathways,
#' filter some pathways out (Only display some important pathways).
#'
#' @param network_nodes Network node table.
#' @param network_edges Network edge table.
#'
#' @return
#' @export
#'
#' @examples
processActivateNetworkDisplay <-
  function(network_nodes, network_edges) {
    important_pathway_list <- c(
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
    unique_pathway <- unique(as.character(network_nodes$parent))
    if (length(unique_pathway) > 15) {
      selected_pathway <-
        unique_pathway[unique_pathway %in% important_pathway_list]
      if (length(selected_pathway) > 15) {
        selected_pathway <- selected_pathway[c(1:15)]
      }
      new_network_nodes <-
        network_nodes[network_nodes$parent %in% selected_pathway, ]
      new_network_edges <-
        network_edges[network_edges$PathwayName %in% selected_pathway, ]
    } else{
      new_network_nodes <- network_nodes
      new_network_edges <- network_edges
    }
    return(list(new_network_nodes, new_network_edges))
  }

#' Process drug clustering result for network generation (Only keep top 10 clusters).
#'
#' @param drug_node Drug network node table.
#' @param drug_network Drug network edge table.
#'
#' @return
#' @export
#'
#' @examples
processDrugClusteringDisplay <- function(drug_node, drug_network) {
  num_cluster <- sum(drug_node$type == "exemplar")
  if (num_cluster > 10) {
    selected_clusters <-
      drug_node[which(drug_node$type == "exemplar")[c(1:10)], ]$node
    new_drug_network <-
      drug_network[drug_network$source %in% selected_clusters |
                     drug_network$target %in% selected_clusters, ]
    selected_nodes <-
      unique(c(new_drug_network$source, new_drug_network$target))
    new_drug_node <- drug_node[drug_node$node %in% selected_nodes, ]
  } else{
    new_drug_node <- drug_node
    new_drug_network <- drug_network
  }
  return(list(new_drug_node, new_drug_network))
}


setRnaDataIdents <- function(rv) {
  data <- rv$rna_df
  if (!is.null(rv$rna_group_list)) {
    ident_list <-
      paste(rv$rna_group_list, rv$rna_type_list, sep = "_")

  } else{
    ident_list <- rv$rna_type_list
  }
  
  if (sum(Idents(data) == ident_list) == length(ident_list)) {
    rm(ident_list)
    return (0)
  }
  else {
    data <- SetIdent(data, value = ident_list)
    rv$rna_df <- data
    save(data, file = paste0(rv$outputDir, "/rna_df.RData"))
    rm(data, ident_list)
    
  }
  gc()
}
  

#' calculate foldchange and pvalue for all the cell type using Seurat DEG analysis method.
#'
#' @param data Seurat data.
#' @param normalGroup Group list used as the control.
#' @param testGroup Group list used as the text.
#' @param cell_type_list cell type list for cell in scRNA dataset.
#' @param cell_group_list cell group list for cell in scRNA dataset.
#' @param progress Progress bar for displaying.
#'
#' @return
#' @export
#'
#' @examples
upRegTest <-
  function(data,
           normalGroup = NULL,
           testGroup = NULL,
           cell_type_list,
           cell_group_list = NULL,
           progress) {
    if (!is.null(cell_group_list)) {
      normal_cell_type <-
        as.character(unique(cell_type_list[cell_group_list %in% normalGroup]))
      test_cell_type <-
        as.character(unique(cell_type_list[cell_group_list %in% testGroup]))
      unique_cell_type <-
        intersect(normal_cell_type, test_cell_type)
    } else{
      unique_cell_type <- as.character(unique(cell_type_list))
    }
    num_each_condition <- table(Idents(data))
    #exculde
    unique_cell_type <-
      unique_cell_type[!unique_cell_type %in% c("outlier", "unknown")]
    unique_filter_cell_type <- c()
    upRegGenes <- c()
    num_type <- length(unique_cell_type)
    genesInformation <- c()
    #use golbal mean
    for (i in 1:num_type) {
      progress$set((i / num_type) * 0.8, detail = "")
      cell_type <- unique_cell_type[i]
      if (is.null(cell_group_list)) {
        if (sum(num_each_condition[cell_type]) < 3) {
          next
        }
        unique_filter_cell_type <-
          c(unique_filter_cell_type, cell_type)
        wlicox_result <-
          FindMarkers(
            data,
            ident.1 = cell_type,
            ident.2 = NULL,
            logfc.threshold = 0,
            only.pos = T
          )
        wlicox_gene <- rownames(wlicox_result)
        
        bimod_result <-
          FindMarkers(
            data,
            ident.1 = cell_type,
            ident.2 = NULL,
            test.use = "bimod",
            logfc.threshold = 0,
            only.pos = T
          )
        bimod_gene <- rownames(bimod_result)
        exp_genes <- intersect(wlicox_gene, bimod_gene)
        exp_genes_table <-
          data.frame(
            logFc = wlicox_result[exp_genes, 2],
            pct1 = wlicox_result[exp_genes, 3],
            pct2 = wlicox_result[exp_genes, 4],
            wlicox_p_value = wlicox_result[exp_genes, 1],
            wlicox_p_adjust = wlicox_result[exp_genes, 5],
            bimod_p_value = bimod_result[exp_genes, 1],
            bimod_p_adjust = bimod_result[exp_genes, 5]
          )
        rownames(exp_genes_table) <- exp_genes
      } else{
        num_normal_group <- length(normalGroup)
        num_test_group <- length(testGroup)
        normal_group_type <-
          paste(normalGroup, rep(cell_type, num_normal_group), sep = "_")
        test_group_type <-
          paste(testGroup, rep(cell_type, num_test_group), sep = "_")
        each_condition <- unique(Idents(data))
        normal_group_type <-
          intersect(normal_group_type, each_condition)
        test_group_type <-
          intersect(test_group_type, each_condition)
        if (length(normal_group_type) == 0 ||
            length(test_group_type) == 0) {
          next
        }
        if (sum(num_each_condition[normal_group_type]) < 3 ||
            sum(num_each_condition[test_group_type]) < 3) {
          next
        }
        unique_filter_cell_type <-
          c(unique_filter_cell_type, cell_type)
        wlicox_result <-
          FindMarkers(
            data,
            ident.1 = test_group_type,
            ident.2 = normal_group_type,
            logfc.threshold = 0,
            only.pos = T
          )
        wlicox_gene <- rownames(wlicox_result)
        
        bimod_result <-
          FindMarkers(
            data,
            ident.1 = test_group_type,
            ident.2 = normal_group_type,
            test.use = "bimod",
            logfc.threshold = 0,
            only.pos = T
          )
        bimod_gene <- rownames(bimod_result)
        exp_genes <- intersect(wlicox_gene, bimod_gene)
        exp_genes_table <-
          data.frame(
            logFc = wlicox_result[exp_genes, 2],
            pct1 = wlicox_result[exp_genes, 3],
            pct2 = wlicox_result[exp_genes, 4],
            wlicox_p_value = wlicox_result[exp_genes, 1],
            wlicox_p_adjust = wlicox_result[exp_genes, 5],
            bimod_p_value = bimod_result[exp_genes, 1],
            bimod_p_adjust = bimod_result[exp_genes, 5]
          )
        rownames(exp_genes_table) <- exp_genes
        exp_genes_table <-
          exp_genes_table[order(exp_genes_table$logFc, decreasing = T),]
      }
      
      genesInformation <-
        c(genesInformation, list(exp_genes_table))
    }
    names(genesInformation) <- unique_filter_cell_type
    genesInformation
  }

#' find up-regluated ligands and receptors based on the threshold
#' user choose. then generate network json used for display
#'
#' @param genesInformation DEG results from upRegTest function.
#' @param gene_symbol Gene symbol list for all genes in scRNA dataset.
#' @param ligand_fc_thres Fold change threshold for ligands (log-scale).
#' @param rec_fc_thres Fold change threshold for receptors (log-scale).
#' @param pv_thres P-value threshold.
#' @param ligRecDatabase Ligand-receptor database.
#'
#' @return
#' @export
#'
#' @examples
up_regulated_network <-
  function(genesInformation,
           gene_symbol,
           ligand_fc_thres,
           rec_fc_thres,
           pv_thres,
           ligRecDatabase) {
    vLigand <- ligRecDatabase[, 1]
    Ligand <- unique(vLigand)
    Ligand <- toupper(Ligand)
    vReceptor <- ligRecDatabase[, 2]
    Receptor <- unique(vReceptor)
    Receptor <- toupper(Receptor)
    unique_cell_type <- names(genesInformation)
    
    cell_ligand_cell <- c()
    cell_ligand_ligand <- c()
    cell_ligand_node <- c()
    cell_ligand_type <- c()
    
    cell_receptor_cell <- c()
    cell_receptor_receptor <- c()
    cell_receptor_node <- c()
    cell_receptor_type <- c()
    for (i in 1:length(genesInformation)) {
      exp_genes_table <- genesInformation[[i]]
      cell_type <- names(genesInformation)[i]
      exp_gene_list <- toupper(rownames(exp_genes_table))
      ligand_up_genes <-
        exp_gene_list[exp_genes_table[, 1] >= ligand_fc_thres &
                        exp_genes_table[, 5] <= pv_thres &
                        exp_genes_table[, 7] <= pv_thres]
      #up-regulated ligand
      upRegLigand <- ligand_up_genes[ligand_up_genes %in% Ligand]
      num_ligand <- length(upRegLigand)
      if (num_ligand > 0) {
        cell_ligand_cell <-
          c(cell_ligand_cell, rep(cell_type, times = num_ligand))
        cell_ligand_ligand <- c(cell_ligand_ligand, upRegLigand)
        cell_ligand_node <- c(cell_ligand_node, upRegLigand)
        cell_ligand_type <-
          c(cell_ligand_type, rep("ligand", times = num_ligand))
        cell_ligand_node <- c(cell_ligand_node, cell_type)
        cell_ligand_type <- c(cell_ligand_type, "cell")
      }
      
      
      #up-regulated receptor
      rec_up_genes <-
        exp_gene_list[exp_genes_table[, 1] >= rec_fc_thres &
                        exp_genes_table[5] <= pv_thres &
                        exp_genes_table[, 7] <= pv_thres]
      #up-regulated ligand
      upRegReceptor <- rec_up_genes[rec_up_genes %in% Receptor]
      num_receptor <- length(upRegReceptor)
      if (num_receptor > 0) {
        cell_receptor_cell <-
          c(cell_receptor_cell, rep(cell_type, times = num_receptor))
        cell_receptor_receptor <-
          c(cell_receptor_receptor, upRegReceptor)
        cell_receptor_node <- c(cell_receptor_node, upRegReceptor)
        cell_receptor_type <-
          c(cell_receptor_type,
            rep("receptor", times = num_receptor))
        
        cell_receptor_node <- c(cell_receptor_node, cell_type)
        cell_receptor_type <- c(cell_receptor_type, "cell")
      }
    }
    cell_ligand_network <-
      data.frame(source = cell_ligand_cell,
                 target = cell_ligand_ligand,
                 stringsAsFactors = F)
    cell_ligand_nodes <-
      data.frame(id = cell_ligand_node,
                 parent = cell_ligand_type,
                 stringsAsFactors = F)
    if (dim(cell_ligand_nodes)[1] != 0 &&
        !is.null(dim(cell_ligand_nodes))) {
      cell_ligand_nodes <- unique(cell_ligand_nodes)
      cell_ligand_nodes <-
        cell_ligand_nodes[order(cell_ligand_nodes$parent),]
      ligand_network <-
        networkDataToJson(cell_ligand_nodes, cell_ligand_network)
    } else{
      ligand_network <- NULL
    }
    
    cell_receptor_network <-
      data.frame(source = cell_receptor_cell,
                 target = cell_receptor_receptor,
                 stringsAsFactors = F)
    cell_receptor_nodes <-
      data.frame(id = cell_receptor_node,
                 parent = cell_receptor_type,
                 stringsAsFactors = F)
    if (dim(cell_receptor_nodes)[1] != 0 &&
        !is.null(dim(cell_receptor_nodes))) {
      cell_receptor_nodes <- unique(cell_receptor_nodes)
      cell_receptor_nodes <-
        cell_receptor_nodes[order(cell_receptor_nodes$parent),]
      receptor_network <-
        networkDataToJson(cell_receptor_nodes, cell_receptor_network)
    } else{
      receptor_network <- NULL
    }
    
    
    cell_ligand_data <-
      list(ligand_network, cell_ligand_network, cell_ligand_nodes)
    cell_receptor_data <-
      list(receptor_network,
           cell_receptor_network,
           cell_receptor_nodes)
    final_data <- list(cell_ligand_data, cell_receptor_data)
  }


#' Generate ligand-receptor communication network given threshold by connecting all ligand-receptors pair
#' exists in the database and meet the threshold.
#'
#' @param genesInformation DEG results from upRegTest function.
#' @param gene_symbol Gene symbol list for all genes in scRNA dataset.
#' @param upReg_fc_thres Up-regulated fold-change threshold.
#' @param pv_thres P-value threshold.
#' @param ligRecDatabase ligand-receptor database.
#'
#' @return
#' @export
#'
#' @examples
upStreamNetwork <-
  function(genesInformation,
           gene_symbol,
           upReg_fc_thres,
           pv_thres,
           ligRecDatabase) {
    upLigand_expReceptor <-
      up_regulated_network(genesInformation,
                           gene_symbol,
                           upReg_fc_thres,
                           0,
                           pv_thres,
                           ligRecDatabase)
    expLigand_upReceptor <-
      up_regulated_network(genesInformation,
                           gene_symbol,
                           0,
                           upReg_fc_thres,
                           pv_thres,
                           ligRecDatabase)
    upLigand_upReceptor <-
      up_regulated_network(
        genesInformation,
        gene_symbol,
        upReg_fc_thres,
        upReg_fc_thres,
        pv_thres,
        ligRecDatabase
      )
    
    #up regulated ligand-expressed receptor
    ligand_links <- upLigand_expReceptor[[1]][[2]]
    ligand_nodes <- upLigand_expReceptor[[1]][[3]]
    ligand_type <-
      sapply(ligand_nodes[, 2], function(x) {
        ifelse(x == "cell", "ligand cell", x)
      })
    ligand_nodes[, 2] <- ligand_type
    ligand_gene_type <-
      apply(ligand_nodes, 1, function(x) {
        paste(x, collapse = "_")
      })
    ligand_nodes <-
      data.frame(id = ligand_nodes[, 1],
                 parent = ligand_type,
                 gene_type = ligand_gene_type)
    ligand_links <-
      merge(
        ligand_links,
        ligand_nodes,
        by.x = "source",
        by.y = "id",
        all.x = TRUE
      )
    ligand_links <-
      merge(
        ligand_links,
        ligand_nodes,
        by.x = "target",
        by.y = "id",
        all.x = TRUE
      )[, c(4, 6)]
    colnames(ligand_links) <- c("source", "target")

    receptor_links <- upLigand_expReceptor[[2]][[2]]
    receptor_nodes <- upLigand_expReceptor[[2]][[3]]
    receptor_type <-
      sapply(receptor_nodes[, 2], function(x) {
        ifelse(x == "cell", "receptor cell", x)
      })
    receptor_nodes[, 2] <- receptor_type
    receptor_gene_type <-
      apply(receptor_nodes, 1, function(x) {
        paste(x, collapse = "_")
      })
    receptor_nodes <-
      data.frame(id = receptor_nodes[, 1],
                 parent = receptor_type,
                 gene_type = receptor_gene_type)
    receptor_links <-
      merge(
        receptor_links,
        receptor_nodes,
        by.x = "source",
        by.y = "id",
        all.x = TRUE
      )
    receptor_links <-
      merge(
        receptor_links,
        receptor_nodes,
        by.x = "target",
        by.y = "id",
        all.x = TRUE
      )[, c(4, 6)]

    receptor_source <- receptor_links[, 2]
    receptor_links[, 2] <- receptor_links[, 1]
    receptor_links[, 1] <- receptor_source
    colnames(receptor_links) <- c("source", "target")
    
    ligand_list <- ligand_nodes[ligand_nodes[, 2] == "ligand", 1]
    receptor_list <-
      receptor_nodes[receptor_nodes[, 2] == "receptor", 1]
    
    upLigand_expReceptor_nodes <-
      unique(rbind(ligand_nodes, receptor_nodes))
    upLigand_expReceptor_links <-
      unique(rbind(ligand_links, receptor_links))
    for (i in 1:length(ligand_list)) {
      ligand = ligand_list[i]
      pair_list <- ligRecDatabase
      colnames(pair_list) <- c("source", "target")
      pair_list <- pair_list[pair_list[, 1] == ligand,]
      if (is.null(dim(pair_list))) {
        if (!pair_list[2] %in% receptor_list) {
          pair_list = NULL
        } else{
          pair_list <-
            data.frame(source = pair_list[1], target = pair_list[2])
        }
      } else{
        pair_list <- pair_list[pair_list[, 2] %in% receptor_list,]
      }
      if (!is.null(pair_list)) {
        if (!is.null(dim(pair_list))) {
          pair_list_ligand <-
            apply(pair_list, 1, function(x) {
              paste(x[1], "ligand", sep = "_")
            })
          pair_list_receptor <-
            apply(pair_list, 1, function(x) {
              paste(x[2], "receptor", sep = "_")
            })
          pair_list <-
            data.frame(source = pair_list_ligand, target = pair_list_receptor)
          upLigand_expReceptor_links <-
            rbind(upLigand_expReceptor_links, pair_list)
        } else{
          pair_list <-
            c(
              paste(pair_list[1], "ligand", sep = "_"),
              paste(pair_list[2], "receptor", sep = "_")
            )
          upLigand_expReceptor_links <-
            rbind(upLigand_expReceptor_links, pair_list)
        }
      }
    }
    upLigand_expReceptor_links <- unique(upLigand_expReceptor_links)
    
    #expressed ligand-upRegulated receptor
    ligand_links <- expLigand_upReceptor[[1]][[2]]
    ligand_nodes <- expLigand_upReceptor[[1]][[3]]
    ligand_type <-
      sapply(ligand_nodes[, 2], function(x) {
        ifelse(x == "cell", "ligand cell", x)
      })
    ligand_nodes[, 2] <- ligand_type
    ligand_gene_type <-
      apply(ligand_nodes, 1, function(x) {
        paste(x, collapse = "_")
      })
    ligand_nodes <-
      data.frame(id = ligand_nodes[, 1],
                 parent = ligand_type,
                 gene_type = ligand_gene_type)
    ligand_links <-
      merge(
        ligand_links,
        ligand_nodes,
        by.x = "source",
        by.y = "id",
        all.x = TRUE
      )
    ligand_links <-
      merge(
        ligand_links,
        ligand_nodes,
        by.x = "target",
        by.y = "id",
        all.x = TRUE
      )[, c(4, 6)]
    colnames(ligand_links) <- c("source", "target")
    receptor_links <- expLigand_upReceptor[[2]][[2]]
    receptor_nodes <- expLigand_upReceptor[[2]][[3]]
    receptor_type <-
      sapply(receptor_nodes[, 2], function(x) {
        ifelse(x == "cell", "receptor cell", x)
      })
    receptor_nodes[, 2] <- receptor_type
    receptor_gene_type <-
      apply(receptor_nodes, 1, function(x) {
        paste(x, collapse = "_")
      })
    receptor_nodes <-
      data.frame(id = receptor_nodes[, 1],
                 parent = receptor_type,
                 gene_type = receptor_gene_type)
    receptor_links <-
      merge(
        receptor_links,
        receptor_nodes,
        by.x = "source",
        by.y = "id",
        all.x = TRUE
      )
    receptor_links <-
      merge(
        receptor_links,
        receptor_nodes,
        by.x = "target",
        by.y = "id",
        all.x = TRUE
      )[, c(4, 6)]
    receptor_source <- receptor_links[, 2]
    receptor_links[, 2] <- receptor_links[, 1]
    receptor_links[, 1] <- receptor_source
    colnames(receptor_links) <- c("source", "target")
    
    ligand_list <- ligand_nodes[ligand_nodes[, 2] == "ligand", 1]
    receptor_list <-
      receptor_nodes[receptor_nodes[, 2] == "receptor", 1]
    expLigand_upReceptor_nodes <-
      unique(rbind(ligand_nodes, receptor_nodes))
    expLigand_upReceptor_links <-
      unique(rbind(ligand_links, receptor_links))
    
    for (i in 1:length(ligand_list)) {
      ligand = ligand_list[i]
      pair_list <- ligRecDatabase
      colnames(pair_list) <- c("source", "target")
      pair_list <- pair_list[pair_list[, 1] == ligand,]
      if (is.null(dim(pair_list))) {
        if (!pair_list[2] %in% receptor_list) {
          pair_list = NULL
        } else{
          pair_list <-
            data.frame(source = pair_list[1], target = pair_list[2])
        }
      } else{
        pair_list <- pair_list[pair_list[, 2] %in% receptor_list,]
      }
      if (!is.null(pair_list)) {
        if (!is.null(dim(pair_list))) {
          pair_list_ligand <-
            apply(pair_list, 1, function(x) {
              paste(x[1], "ligand", sep = "_")
            })
          pair_list_receptor <-
            apply(pair_list, 1, function(x) {
              paste(x[2], "receptor", sep = "_")
            })
          pair_list <-
            data.frame(source = pair_list_ligand, target = pair_list_receptor)
          expLigand_upReceptor_links <-
            rbind(expLigand_upReceptor_links, pair_list)
        } else{
          pair_list <-
            c(
              paste(pair_list[1], "ligand", sep = "_"),
              paste(pair_list[2], "receptor", sep = "_")
            )
          expLigand_upReceptor_links <-
            rbind(expLigand_upReceptor_links, pair_list)
        }
      }
    }
    expLigand_upReceptor_links <- unique(expLigand_upReceptor_links)
    
    #up regulated ligand-upRegulated receptor
    ligand_links <- upLigand_upReceptor[[1]][[2]]
    ligand_nodes <- upLigand_upReceptor[[1]][[3]]
    ligand_type <-
      sapply(ligand_nodes[, 2], function(x) {
        ifelse(x == "cell", "ligand cell", x)
      })
    ligand_nodes[, 2] <- ligand_type
    ligand_gene_type <-
      apply(ligand_nodes, 1, function(x) {
        paste(x, collapse = "_")
      })
    ligand_nodes <-
      data.frame(id = ligand_nodes[, 1],
                 parent = ligand_type,
                 gene_type = ligand_gene_type)
    ligand_links <-
      merge(
        ligand_links,
        ligand_nodes,
        by.x = "source",
        by.y = "id",
        all.x = TRUE
      )
    ligand_links <-
      merge(
        ligand_links,
        ligand_nodes,
        by.x = "target",
        by.y = "id",
        all.x = TRUE
      )[, c(4, 6)]
    colnames(ligand_links) <- c("source", "target")
    receptor_links <- upLigand_upReceptor[[2]][[2]]
    receptor_nodes <- upLigand_upReceptor[[2]][[3]]
    receptor_type <-
      sapply(receptor_nodes[, 2], function(x) {
        ifelse(x == "cell", "receptor cell", x)
      })
    receptor_nodes[, 2] <- receptor_type
    receptor_gene_type <-
      apply(receptor_nodes, 1, function(x) {
        paste(x, collapse = "_")
      })
    receptor_nodes <-
      data.frame(id = receptor_nodes[, 1],
                 parent = receptor_type,
                 gene_type = receptor_gene_type)
    receptor_links <-
      merge(
        receptor_links,
        receptor_nodes,
        by.x = "source",
        by.y = "id",
        all.x = TRUE
      )
    receptor_links <-
      merge(
        receptor_links,
        receptor_nodes,
        by.x = "target",
        by.y = "id",
        all.x = TRUE
      )[, c(4, 6)]
    receptor_source <- receptor_links[, 2]
    receptor_links[, 2] <- receptor_links[, 1]
    receptor_links[, 1] <- receptor_source
    colnames(receptor_links) <- c("source", "target")
    
    ligand_list <- ligand_nodes[ligand_nodes[, 2] == "ligand", 1]
    receptor_list <-
      receptor_nodes[receptor_nodes[, 2] == "receptor", 1]
    upLigand_upReceptor_nodes <-
      unique(rbind(ligand_nodes, receptor_nodes))
    upLigand_upReceptor_links <-
      unique(rbind(ligand_links, receptor_links))
    for (i in 1:length(ligand_list)) {
      ligand = ligand_list[i]
      pair_list <- ligRecDatabase
      colnames(pair_list) <- c("source", "target")
      pair_list <- pair_list[pair_list[, 1] == ligand,]
      if (is.null(dim(pair_list))) {
        if (!pair_list[2] %in% receptor_list) {
          pair_list = NULL
        } else{
          pair_list <-
            data.frame(source = pair_list[1], target = pair_list[2])
        }
      } else{
        pair_list <- pair_list[pair_list[, 2] %in% receptor_list,]
      }
      if (!is.null(pair_list)) {
        if (!is.null(dim(pair_list))) {
          pair_list_ligand <-
            apply(pair_list, 1, function(x) {
              paste(x[1], "ligand", sep = "_")
            })
          pair_list_receptor <-
            apply(pair_list, 1, function(x) {
              paste(x[2], "receptor", sep = "_")
            })
          pair_list <-
            data.frame(source = pair_list_ligand, target = pair_list_receptor)
          upLigand_upReceptor_links <-
            rbind(upLigand_upReceptor_links, pair_list)
        } else{
          pair_list <-
            c(
              paste(pair_list[1], "ligand", sep = "_"),
              paste(pair_list[2], "receptor", sep = "_")
            )
          upLigand_upReceptor_links <-
            rbind(upLigand_upReceptor_links, pair_list)
        }
      }
    }
    upLigand_upReceptor_links <- unique(upLigand_upReceptor_links)
    #combine 1 and 2
    combine_links <-
      unique(rbind(upLigand_expReceptor_links, expLigand_upReceptor_links))
    combine_nodes <-
      unique(rbind(upLigand_expReceptor_nodes, expLigand_upReceptor_nodes))
    
    if (nrow(upLigand_expReceptor_links) != 0) {
      upLigand_expReceptor_data <-
        upStreamNetworkGenerating(upLigand_expReceptor_nodes, upLigand_expReceptor_links)
    } else{
      upLigand_expReceptor_data <- NULL
    }
    if (nrow(expLigand_upReceptor_links) != 0) {
      expLigand_upReceptor_data <-
        upStreamNetworkGenerating(expLigand_upReceptor_nodes, expLigand_upReceptor_links)
      
    } else{
      expLigand_upReceptor_data <- NULL
    }
    if (nrow(combine_links) != 0) {
      combine_data <-
        upStreamNetworkGenerating(combine_nodes, combine_links)
    } else{
      combine_data <- NULL
    }
    if (nrow(upLigand_upReceptor_links) != 0) {
      upLigand_upReceptor_data <-
        upStreamNetworkGenerating(upLigand_upReceptor_nodes, upLigand_upReceptor_links)
      
    } else{
      upLigand_upReceptor_data <- NULL
    }
    
    
    list(
      upLigand_expReceptor_data,
      expLigand_upReceptor_data,
      combine_data,
      upLigand_upReceptor_data
    )
  }



#' upstream network generating function
#'
#' @param nodes Table for all nodes in the network.
#' @param links Table for all links in the network.
#'
#' @return
#' @export
#'
#' @examples
upStreamNetworkGenerating <- function(nodes, links) {
  network_nodes <-
    data.frame(
      id = c(1:nrow(nodes)),
      name = as.character(nodes$id),
      parent = as.character(nodes$parent),
      gene_type = nodes$gene_type
    )
  network_edges <- unique(links)
  network_edges <-
    merge(network_edges,
          network_nodes,
          by.x = "source",
          by.y = "gene_type")[,-c(4, 5)]
  network_edges <-
    merge(network_edges,
          network_nodes,
          by.x = "target",
          by.y = "gene_type")[,-c(5, 6)]
  network_edges <-
    data.frame(source = network_edges[, 3], target = network_edges[, 4])
  if (sum(is.na(network_edges$source) |
          is.na(network_edges$target)) > 0) {
    network_edges <-
      network_edges[-which(is.na(network_edges$source) |
                             is.na(network_edges$target)),]
  }
  network_nodes <- network_nodes[,-4]
  network_edges <- as.matrix(network_edges)
  network_graph <- graph.edgelist(network_edges, directed = T)
  #igraph only accept name attribute as what can return in edgelist. set id as name attribute
  vertex_attr(network_graph) <-
    list(name = network_nodes$id,
         symbol = network_nodes$name,
         parent = network_nodes$parent)
  node_outdegree <- degree(network_graph, mode = "out")
  node_indegree <- degree(network_graph, mode = "in")
  keep_nodes <-
    (network_nodes$parent == "ligand" &
       node_outdegree > 0) |
    (network_nodes$parent == "receptor" &
       node_indegree > 0) |
    network_nodes$parent %in% c("ligand cell", "receptor cell")
  network_subgraph <-
    induced.subgraph(network_graph, vids = which(keep_nodes))
  sub_degree <- degree(network_subgraph)
  network_subgraph <-
    induced.subgraph(network_subgraph, vids = which(sub_degree > 0))
  network_edges <- get.edgelist(network_subgraph, names = T)
  network_edges <- data.frame(network_edges)
  colnames(network_edges) <- c("source", "target")
  network_nodes <-
    data.frame(
      id = get.vertex.attribute(network_subgraph, "name"),
      name = get.vertex.attribute(network_subgraph, "symbol"),
      parent = get.vertex.attribute(network_subgraph, "parent")
    )
  final_data <- networkDataToJson(network_nodes, network_edges)
  list(final_data, network_edges, network_nodes)
}

#' Function for processing drug file for signature drug discovery.
#'
#' @param dataDir Directory for loading raw drug files.
#' @param drugBankInformation DrugBank database.
#' @param progress Progress bar for display
#'
#' @return
#' @export
#'
#' @examples
processDrugFile <- function(dataDir, drugBankInformation, progress) {
  progress$set(0.1, detail = "Process drug data")
  ioGctx <- 'GCTXio.R'
  utilsGctx <- 'GCTXutils.R'
  dataGctx <-  'GCTXdata.R'
  source(ioGctx)
  source(utilsGctx)
  source(dataGctx)
  gctx_file <-
    paste(dataDir,
          'GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx',
          sep = "/")
  sig_info_file <-
    paste(dataDir, 'GSE92742_Broad_LINCS_sig_info.txt', sep = '/')
  
  # experimental annotation file
  sig_info = read.delim(sig_info_file, sep = "\t")
  progress$set(0.2, detail = "Process drug data")
  pert_id = sig_info[, 2]
  cell_id0 = as.character(sig_info[, 5])
  pert_idose0 = round(as.numeric(as.character(sig_info[, 6])))
  
  pert_idose0[is.na(pert_idose0)] <- -10
  
  pert_itime0 = as.character(sig_info[, 9]) # hours: 24   6  96 144 120 168  48   1   3   2   4  72
  pert_dose_unit <-
    as.character(sig_info[, 7]) #  %     µM    -666  µL    ng/mL ng    ng/µL
  
  cell_IDs <-
    c("PC3",
      "VCAP",
      "A375",
      "A549",
      "HA1E",
      "HCC515",
      "HT29",
      "MCF7",
      "HEPG2")
  cidT <-
    which((pert_itime0 == "24") &
            (pert_idose0 == 10) &
            (pert_dose_unit == 'µM') & cell_id0 %in% cell_IDs
    )
  cidT <- intersect(cidT, grep("BRD-", pert_id)) # compounds
  rm(
    sig_info_file,
    sig_info,
    pert_id,
    cell_id0,
    pert_idose0,
    pert_itime0,
    pert_dose_unit,
    cell_IDs,
    ioGctx,
    utilsGctx,
    dataGctx
  )
  gc()
  progress$set(0.3, detail = "Process drug data")
  ds <- parse.gctx(gctx_file, cid = cidT) # 12328 x 38418
  ds <- ds@mat
  gc()
  gene_id1 <- row.names(ds)
  
  progress$set(0.5, detail = "Process drug data")
  num_col <- dim(ds)[2]
  iter_time = 8
  separate <- ceiling(num_col / iter_time)
  for (i in c(1:iter_time)) {
    if (i == 1) {
      rankMatrix <- apply(-1.0 * ds[, c(1:separate)], 2, order)
    } else if (i == iter_time) {
      rankMatrix <- cbind(rankMatrix, apply(-1.0 * ds, 2, order))
    }
    else{
      rankMatrix <- cbind(rankMatrix, apply(-1.0 * ds[, c(1:separate)], 2, order))
    }
    
    if (i != iter_time) {
      ds <- ds[, -c(1:separate)]
    }
    gc()
  }
  rm(ds, gctx_file, num_col, separate, iter_time)
  gc()
  # get the gene information of the each probe
  f_gene <-
    paste(dataDir, 'GSE92742_Broad_LINCS_gene_info.txt', sep = '/')
  gene_info <- read.delim(f_gene, sep = '\t')
  gene_id0 <- as.character(gene_info[, 1])
  nt <- length(gene_id1)
  idx_1a <- rep(0, nt)
  for (i in 1:nt) {
    st <- gene_id1[i]
    idxt <- which(gene_id0 %in% st)
    if (length(idxt) > 0) {
      idx_1a[i] <- idxt[1]
    }
  }
  gSymZs <- as.character(gene_info[idx_1a, 2])
  rm(f_gene, gene_info, gene_id0, nt, idx_1a, st, idxt)
  gc()
  
  progress$set(0.6, detail = "Process drug data")
  
  pert_iname <- colnames(rankMatrix)
  rankname_table <- data.frame(pert_iname)
  drug_name_split <- sapply(pert_iname, function(x) {
    strsplit(x, ":")
  })
  progress$set(0.7, detail = "Process drug data")
  drug_information_name <- c()
  for (i in 1:length(drug_name_split)) {
    a <- drug_name_split[[i]][2]
    drug_information_name <- c(drug_information_name, a)
  }
  
  for (i in 1:length(drug_information_name)) {
    if (nchar(drug_information_name[i]) == 22) {
      drug_information_name[i] <- substring(drug_information_name[i], 1, 13)
    }
  }
  drug_index <- c(1:length(drug_information_name))
  rankname_table <- cbind(rankname_table, drug_information_name)
  colnames(rankname_table) <- c("rankName", "pert_id")
  rankname_table <- cbind(rankname_table, drug_index)
  
  rm(pert_iname,
     drug_name_split,
     drug_information_name,
     drug_index)
  gc()
  
  progress$set(0.8, detail = "Process drug data")
  
  f_drug <-
    paste(dataDir, "GSE92742_Broad_LINCS_pert_info.txt", sep = "/")
  drug_information <- read.delim(f_drug)
  drug_mapping <-
    merge(rankname_table,
          drug_information,
          all.x = TRUE,
          by = "pert_id")
  rm(drug_information)
  gc()
  drug_mapping <- drug_mapping[order(drug_mapping$drug_index), ]
  pert_iname <- as.character(drug_mapping$pert_iname)
  
  
  drugBank_drug <-
    toupper(as.character(drugBankInformation$drug_name))
  drug_filter <- function(x, drugBank_drug) {
    x <- toupper(x)
    return(grepl("^BRD|TRCN|CMAP", x) || x %in% drugBank_drug)
    
  }
  keep_drug <-
    sapply(pert_iname, drug_filter, drugBank_drug = drugBank_drug)
  rm(drugBank_drug, pert_iname)
  gc()
  
  drug_mapping <- drug_mapping[keep_drug, ]
  drug_mapping <- drug_mapping[, -3]
  num_row <- nrow(drug_mapping)
  rownames(drug_mapping) <- c(1:num_row)
  rankMatrix <- rankMatrix[, keep_drug]
  progress$set(0.9, detail = "Save files to application")
  save(drug_mapping, file = "./data/drug_mapping.RData")
  save(rankMatrix, gSymZs, file = './data/rankMatrix92742a.RData')
  return(list(rankMatrix, drug_mapping, gSymZs))
}

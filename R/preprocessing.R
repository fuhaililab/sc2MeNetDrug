#' quality control for raw scRNA-seq data. This function will remove cells
#' in the data with extremely low/high read counts or mitochondrial contamination.
#'
#' @param seurat_data Seurat object with the raw read count data saved in the object.
#'
#' @return Seurat object with the filtered data.
#' @export
#'
#' @examples
quality_control <- function(seurat_data) {
  seurat_data[["percent.mt"]] <-
    PercentageFeatureSet(seurat_data, pattern = "^MT-")
  seurat_data <-
    subset(seurat_data, subset = nFeature_RNA > 200 &
             nFeature_RNA < 7500 & percent.mt < 10)
  seurat_data
}

#' read count normalization and scaling.
#'
#' @param seurat_data Seurat object with the raw read count data saved in the object.
#'
#' @return Seura object with the normalized data saved.
#' @export
#'
#' @examples
sctransform_seurat <- function(seurat_data) {
  seurat_data <- SCTransform(
    seurat_data,
    assay = seurat_data@active.assay,
    method = "glmGamPoi",
    vars.to.regress = "percent.mt",
    verbose = FALSE
  )
  seurat_data
}


#' Use the ALRA algorithm to do the imputation for scRNA-seq data.
#' Notes that this function will be used after the normalization and scaling.
#' @param seurat_data Seurat object with the raw read count data saved in the object.
#'
#' @return Seurat object with the imputed data saved.
#' @export
#'
#' @examples
imputation <- function(seurat_data) {
  seurat_data <- RunALRA(seurat_data, assay = seurat_data@active.assay)
  
  seurat_data <-
    FindVariableFeatures(seurat_data,
                         nfeatures = 3000,
                         assay = seurat_data@active.assay)
  seurat_data <-
    ScaleData(seurat_data, assay = seurat_data@active.assay)
  seurat_data
}


#' Convert the mouse gene symbol to the corresponding human gene symbol.
#'
#' @param data read count data with each row represent each gene.
#' @param gene_list The gene list that will be converted.
#'
#' @return read count data.
#' @export
#'
#' @examples
mouse_to_human_convert <- function(data, gene_list) {
  load("data/mouse_human_gene_pair.RData")
  gene_list <- data.frame(original_gene = gene_list)
  map_gene_result <-
    gene_list %>% inner_join(mouse_human_pair, by = c("original_gene" = "mouse_gene"))
  
  data <- data[map_gene_result$original_gene,]
  rownames(data) <- as.character(map_gene_result$human_gene)
  if (sum(duplicated(rownames(data))) > 0) {
    gene_variance <- apply(data, 1, var)
    data <- data[order(gene_variance, decreasing = T),]
    data <- data[!duplicated(rownames(data)),]
  }
  data
}

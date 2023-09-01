# Algorithms for dimension reduction

#' Function to do the PCA analysis of scRNA-seq data.
#'
#' @param seurat_data Seurat object
#'
#' @return Seurat object with the PCA results saved.
#' @export
#'
#' @examples
PCA <- function(seurat_data) {
  seurat_data <- RunPCA(seurat_data, verbose = FALSE)
  seurat_data
}

#' Function to do the UMAP analysis for scRNA-seq data.
#'
#' @param seurat_data Seurat object.
#' @param nPC The number of top dimensions used for UMAP algorithm.
#'
#' @return Seurat object with UMAP results saved.
#' @export
#'
#' @examples
UMAP <- function(seurat_data, nPC) {
  seurat_data <- RunUMAP(seurat_data, dims = 1:nPC, verbose = FALSE)
  seurat_data
}
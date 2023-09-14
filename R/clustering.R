# Function to find neighbor for each cell in scRNA-seq data.
#'
#' @param seurat_data Seurat object.
#' @param nPC The number of top dimensions used for the clustering algorithm.
#'
#' @return Seurat object with the neighbor results saved.
#' @export
#'
#' @examples
#'
find_neighbors <- function(seurat_data, nPC) {
  seurat_data <- FindNeighbors(seurat_data,
                               reduction="pca",
                               dims = 1:nPC,
                               verbose = FALSE)
  seurat_data
}


#' Function to do the clustering for scRNA-seq data.
#'
#' @param seurat_data Seurat object.
#' @param resolution Clustering resolution, a larger value with results in more clusters.
#'
#' @return Seurat object with clustering results saved.
#' @export
#'
#' @examples
clustering <- function(seurat_data, resolution) {
  seurat_data <- FindClusters(seurat_data,
                              resolution = resolution,
                              verbose = FALSE)
  seurat_data
}
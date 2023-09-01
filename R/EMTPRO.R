#' 0-x normalization.
#'
#' @param data data input for the normalization.
#' @param base The maximum number x.
#'
#' @return
#' @export
#'
#' @examples
normalization <- function(data, base = 2) {
  gene_max <- apply(data, base, max)
  center <- sweep(data, base, apply(data, base, min))
  R <- gene_max - apply(data, base, min)
  #prevent Nan
  R[R == 0] = 0.00001
  df_marker_normal <- sweep(center, base, R, "/")
  df_marker_normal
}



#' Calculate EMT-PRO score use corresponding marker gene.
#'
#' @param df scRNA expression matrix.
#' @param gene_list Gene symbol list of the scRNA dataset.
#' @param group_list cell group list of the scRNA dataset.
#' @param type_list cell type list of the scRNA dataset.
#' @param epGroup The group list included for EMT-PRO analysis
#' @param epType The cell type for computing EMT-PRO analysis.
#'
#' @return
#' @export
#'
#' @examples
emt_pro <-
  function(df,
           gene_list,
           group_list = NULL,
           type_list,
           epGroup,
           epType) {
    PRO <-
      as.vector(read.table("./data/HALLMARK/HALLMARK_E2F_TARGETS.txt")$V1)
    EMT <-
      as.vector(
        read.table(
          "./data/HALLMARK/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt"
        )$V1
      )
    intersect_gene <- intersect(gene_list, c(PRO, EMT))
    PRO_GENE_LIST <- intersect(gene_list, PRO)
    EMT_GENE_LIST <- intersect(gene_list, EMT)
    #prevent error if no gene in data
    if (length(EMT_GENE_LIST) == 0 || length(PRO_GENE_LIST) == 0) {
      stop(safeError("no enough gene"))
    }
    if (!is.null(group_list)) {
      selected_df <- df[, group_list %in% epGroup & type_list %in% epType]
    } else{
      selected_df <- df[, type_list %in% epType]
    }
    selected_df_gene <-
      selected_df[c(PRO_GENE_LIST, EMT_GENE_LIST),]
    normal_expression <- normalization(selected_df_gene, base = 1)
    normal_expression <- as.matrix(normal_expression)
    pro_score <- apply(normal_expression[PRO_GENE_LIST,], 2, mean)
    emt_score <- apply(normal_expression[EMT_GENE_LIST,], 2, mean)
    print(length(pro_score))
    print(length(emt_score))
    list(emt_score, pro_score)
  }

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# All reactive modules used in the app. Mainly contain functions that need reactive
# loading (Like loading table).
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Listen change in DEG result table data for up-regulated
#' ligand receptor analysis and dynamic show or hide table
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadTableReactive <- function(input, rv) {
  if (is.null(rv$cell_type1_ligands) ||
      is.null(rv$cell_type1_receptors) ||
      is.null(rv$cell_type2_ligands) ||
      is.null(rv$cell_type2_receptors)) {
    return (NULL)
  }
  if (input$networkTableSelection == "ligType1") {
    rv$df_table = rv$cell_type1_ligands
  } else if (input$networkTableSelection == "recType1") {
    rv$df_table = rv$cell_type1_receptors
  } else if (input$networkTableSelection == "ligType2") {
    rv$df_table = rv$cell_type2_ligands
  } else if (input$networkTableSelection == "recType2") {
    rv$df_table = rv$cell_type2_receptors
  }
  rv$networkTableDisplayUI <- {
    DTOutput("networkTable")
  }
  
}

#' Listen change in signature drug table data 1 for drug discovery
#'  and dynamic show or hide table.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadDrugTable1Reactive <- function(rv) {
  if (is.null(rv$drug_table1)) {
    return (NULL)
  }
  rv$drugTableTitle1UI <- {
    tags$h3(
      paste(
        "Signaling Signature Drug Discovering Result for Downstream Network from",
        paste(rv$cell_type1, collapse = "+"),
        "to"
        ,
        paste(rv$cell_type2, collapse = "+"),
        sep = " "
      )
    )
  }
  rv$drugTable1UI <- {
    DTOutput("drugTable1")
  }
}

#' Listen change in signature drug table data 2 for drug discovery
#'  and dynamic show or hide table.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadDrugTable2Reactive <- function(rv) {
  if (is.null(rv$drug_table2)) {
    return (NULL)
  }
  rv$drugTableTitle2UI <- {
    tags$h3(
      paste(
        "Signaling Signature Drug Discovering Result for Downstream Network from",
        paste(rv$cell_type2, collapse = "+"),
        "to"
        ,
        paste(rv$cell_type1, collapse = "+"),
        sep = " "
      )
    )
  }
  rv$drugTable2UI <- {
    DTOutput("drugTable2")
  }
}


#' Listen change in GO table data for GO analysis
#'  and dynamic show or hide table.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadGOTableReactive <- function(input, rv) {
  if (rv$GOData == 0) {
    return(NULL)
  }
  if (input$GOTableSelection == "Up-regulated GO Result") {
    rv$GO_table = rv$GO_up_table
  } else{
    rv$GO_table = rv$GO_dn_table
  }
  rv$GOTableDisplayUI <- {
    DTOutput("GOTable")
  }
}

#' Listen change in target drug table data 1 for drug discovery
#'  and dynamic show or hide table.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadDrugMappingTable1 <- function(input, rv) {
  if (is.null(rv$drug_mapping_table1)) {
    return(NULL)
  }
  rv$targetDrugTableTitle1UI <- {
    tags$h3(
      paste(
        "Targets Drug Discovering Result for Downstream Network from",
        paste(rv$cell_type1, collapse = "+"),
        "to"
        ,
        paste(rv$cell_type2, collapse = "+"),
        sep = " "
      )
    )
  }
  DTOutput("drug_mapping_table1")
}


#' Listen change in target drug table data 2 for drug discovery
#'  and dynamic show or hide table.
#'
#' @param input Input variables of the Shiny app.
#' @param rv Variables list in the Shiny app.
#'
#' @return
#' @export
#'
#' @examples
loadDrugMappingTable2 <- function(input, rv) {
  if (is.null(rv$drug_mapping_table2)) {
    return(NULL)
  }
  rv$targetDrugTableTitle2UI <- {
    tags$h3(
      paste(
        "Targets Drug Discovering Result for Downstream Network from",
        paste(rv$cell_type2, collapse = "+"),
        "to"
        ,
        paste(rv$cell_type1, collapse = "+"),
        sep = " "
      )
    )
  }
  DTOutput("drug_mapping_table2")
}
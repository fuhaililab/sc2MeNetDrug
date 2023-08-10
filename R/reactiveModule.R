#listen change in table data and 
#dynamic show or hide table
loadTableRactive<-function(input,rv){
  if(is.null(rv$cell_type1_ligands)||is.null(rv$cell_type1_receptors)||is.null(rv$cell_type2_ligands)||is.null(rv$cell_type2_receptors)){
    return (NULL)
  }
  if(input$networkTableSelection=="ligType1"){
    rv$df_table=rv$cell_type1_ligands
  }else if(input$networkTableSelection=="recType1"){
    rv$df_table=rv$cell_type1_receptors
  }else if(input$networkTableSelection=="ligType2"){
    rv$df_table=rv$cell_type2_ligands
  }else if(input$networkTableSelection=="recType2"){
    rv$df_table=rv$cell_type2_receptors
  }
    rv$networkTableDisplayUI<-{
      DTOutput("networkTable")}

}

#listen change in table data and 
#dynamic show or hide table
loadDrugTable1Ractive<-function(rv){
  if(is.null(rv$drug_table1)){
    return (NULL)
  }
  rv$drugTableTitle1UI<-{
    tags$h3(paste("Signaling Signature Drug Discovering Result for Downstream Network from",paste(rv$cell_type1,collapse = "+"),"to" 
                  ,paste(rv$cell_type2,collapse = "+"), sep=" "))
  }
  rv$drugTable1UI<-{
    DTOutput("drugTable1")}
}
loadDrugTable2Ractive<-function(rv){

  if(is.null(rv$drug_table2)){
    return (NULL)
  }
  rv$drugTableTitle2UI<-{
    tags$h3(paste("Signaling Signature Drug Discovering Result for Downstream Network from",paste(rv$cell_type2,collapse = "+"),"to" 
                  ,paste(rv$cell_type1,collapse = "+"), sep=" "))
  }
  rv$drugTable2UI<-{
    DTOutput("drugTable2")}
}

loadGOTableRactive<-function(input,rv){
  if(rv$GOData==0){
    return(NULL)
  }
  if(input$GOTableSelection=="Up-regulated GO Result"){
    rv$GO_table=rv$GO_up_table
  }else{
    rv$GO_table=rv$GO_dn_table
  }
  rv$GOTableDisplayUI<-{DTOutput("GOTable")}
}
loadDrugMappingTable1<-function(input,rv){
  if(is.null(rv$drug_mapping_table1)){
    return(NULL)
  }
  rv$targetDrugTableTitle1UI<-{
    tags$h3(paste("Targets Drug Discovering Result for Downstream Network from",paste(rv$cell_type1,collapse = "+"),"to" 
                  ,paste(rv$cell_type2,collapse = "+"), sep=" "))
  }
  DTOutput("drug_mapping_table1")
}
loadDrugMappingTable2<-function(input,rv){
  if(is.null(rv$drug_mapping_table2)){
    return(NULL)
  }
  rv$targetDrugTableTitle2UI<-{
    tags$h3(paste("Targets Drug Discovering Result for Downstream Network from",paste(rv$cell_type2,collapse = "+"),"to" 
                  ,paste(rv$cell_type1,collapse = "+"), sep=" "))
  }
  DTOutput("drug_mapping_table2")
}
#calculate EMT-PRO score use corresponding marker gene
emt_pro<-function(df,gene_list,group_list=NULL,type_list,epGroup,epType){
  PRO<-as.vector(read.table("./data/HALLMARK/HALLMARK_E2F_TARGETS.txt")$V1)
  EMT<-as.vector(read.table("./data/HALLMARK/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt")$V1)
  intersect_gene<-intersect(gene_list,c(PRO,EMT))
  PRO_GENE_LIST<-intersect(gene_list,PRO)
  EMT_GENE_LIST<-intersect(gene_list,EMT)
  #prevent error if no gene in data
  if(length(EMT_GENE_LIST)==0||length(PRO_GENE_LIST)==0){
    stop(safeError("no enough gene"))
  }
  if(!is.null(group_list)){
    selected_df<-df[,group_list%in%epGroup&type_list%in%epType]
  }else{
    selected_df<-df[,type_list%in%epType]
  }
  selected_df_gene<-selected_df[c(PRO_GENE_LIST,EMT_GENE_LIST),]
  normal_expression<-normalization(selected_df_gene,base=1)
  normal_expression<-as.matrix(normal_expression)
  pro_score<-apply(normal_expression[PRO_GENE_LIST,],2,mean)
  emt_score<-apply(normal_expression[EMT_GENE_LIST,],2,mean)
  print(length(pro_score))
  print(length(emt_score))
  list(emt_score,pro_score)
}



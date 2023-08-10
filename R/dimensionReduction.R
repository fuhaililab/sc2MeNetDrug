#mutli-core tsne function
runTsne<-function(df,max_iter,perplexity){
  df_tsne=Rtsne(df,initial_config = NULL, dims = 2,  
               max_iter = max_iter,perplexity =perplexity, pca=FALSE,num_threads=2,verbose=TRUE,theta=0.3,check_duplicates=FALSE)
  
  return(as.data.frame(df_tsne$Y))
}


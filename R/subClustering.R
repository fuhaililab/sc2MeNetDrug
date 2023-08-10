#sub clustering, loop for every main cluster,
#if the shape of cluster is not good, use gmm to discover 
#sub cluster in the main cluster group
sub_cluster<-function(df_cluster,progress){
  pre_cluster<-df_cluster[,3]
  #exclude outlier 
  num_clu<-length(unique(pre_cluster))-1
  lx<-max(df_cluster[,1])-min(df_cluster[,1])
  ly<-max(df_cluster[,2])-min(df_cluster[,2])
  n<-nrow(df_cluster)
  
  for(i in c(1:num_clu)){
    progress$set((i/num_clu)*0.95,detail="")
    sub_df<-df_cluster[pre_cluster==i,]
    sub_n<-nrow(sub_df)
    sub_lx<-max(sub_df[,1])-min(sub_df[,1])
    sub_ly<-max(sub_df[,2])-min(sub_df[,2])
    if((sub_n>0.2*n)|(sub_lx>0.3*lx)|(sub_ly>0.3*ly)){
      num_sub_clu<-max(round((sub_lx/lx)*7),round((sub_ly/ly)*7))
      print(num_sub_clu)
      if(num_sub_clu<2|sub_n<=2|sub_n<=num_sub_clu){
        next
      }
      else{
        sub_clu<-do_gmm(sub_df[,c(1,2)],num_sub_clu)
        sub_clu<-sapply(sub_clu,function(x){paste(i,x,sep="-")})
        df_cluster[pre_cluster==i,3]<-sub_clu
      }
    }
  }
  
  new_clu<-unique(df_cluster[,3])
  new_clu<-new_clu[which(new_clu!="0")]
  num_new_clu<-length(new_clu)
  new_clu_df<-data.frame(index=c(1:num_new_clu),cluster=new_clu,stringsAsFactors = FALSE)
  new_clu_df<-rbind(c(0,"0"),new_clu_df)
  colnames(df_cluster)[3]<-"cluster"
  #prevent sort in merge
  df_cluster$id<-c(1:nrow(df_cluster))
  sub_df_cluster<-merge(df_cluster,new_clu_df,by="cluster",all.x=TRUE,sort=FALSE)
  sub_df_cluster<-sub_df_cluster[order(sub_df_cluster$id),]
  sub_df_cluster<-sub_df_cluster[,-c(1,4)]
  colnames(sub_df_cluster)[which(colnames(sub_df_cluster)=="index")]="cluster"
  sub_df_cluster
}

do_gmm<-function(df,G){
  if(G<=0){
    Mclust(df, G=1, model="EII")
  }
  gmm_result<-Mclust(df, G=G, model="EII")
  if(is.null(gmm_result)){
    do_gmm(df,G-1)
  }else{
    gmm_result$classification
  }
}


set.seed(123)
#optics clustering method
findOrder<-function(df,eps=NULL,minPts=5){
  result<-optics(df,eps=NULL,minPts=5)
  result
}

extract<-function(result,xi){
  cl<-extractDBSCAN(result,xi)
}


# 
# result<-findOrder(df_tsne,eps=1,minPts = 50)
# dfClust<-extractDBSCAN(result,0.8)
# unique(dfClust$cluster)
# df_cluster<-cbind(df_tsne,dfClust$cluster)
# 
# colnames(df_cluster)[3]<-"cluster"
# 
# plot_ly(df_cluster,x=~V1,y=~V2,type="scatter",mode="markers"
#         ,color=~factor(cluster),colors="Accent",width=700,height=600,
#         text=~factor(cluster),hovertemplate = paste('<b>Cluster</b>: %{text}'))%>%
#   layout(legend=c(itemdoubleclick="toggleothers"))

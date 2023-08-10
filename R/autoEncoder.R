#auto encoder, the input layer should be 2048
autoEncoder<-function(df_normal,input_dim,progress){
  #AutoEncoder
  input<-layer_input(shape=c(input_dim))
  #encoder
  encoder<-input%>%
    layer_dense(units=1024,activation="relu")%>%
    layer_batch_normalization()%>%
    layer_activation("relu")%>%
    
    
    layer_dense(units=512,activation="relu")%>%
    layer_batch_normalization()%>%
    layer_activation("relu")%>%
    layer_dropout(0.2)%>%
    
    layer_dense(units=128,activation="relu")%>%
    layer_batch_normalization()%>%
    layer_activation("relu")%>%
    layer_dropout(0.3)%>%
    
    layer_dense(units=64,activation="relu")
  
  encoder_model<-keras_model(input=input,output=encoder)
  
  decoder<-encoder%>%
    layer_dense(units=128,activation="relu")%>%
    layer_batch_normalization()%>%
    layer_activation("relu")%>%
    #layer_dropout(0.3)%>%
    
    layer_dense(units=512,activation="relu")%>%
    layer_batch_normalization()%>%
    layer_activation("relu")%>%
    #layer_dropout(0.2)%>%
    
    layer_dense(units=1024,activation="relu")%>%
    layer_batch_normalization()%>%
    layer_activation("relu")%>%
    
    layer_dense(units=2048,activation="relu")
  
  AE_model<-keras_model(input=input,output=decoder)
  AE_model%>%compile(
    loss = 'mse',
    optimizer = optimizer_adam()
  )

  #listen the training progress and update the 
  #progress in application UI
  epochMonitor <- R6::R6Class("epochMonitor",
                             inherit = KerasCallback,
                             public = list(
                               on_epoch_begin = function(epoch,logs=list()) {
                                 progress$inc(0.04,detail=paste0("epoch: ",epoch+1))
                               }
                             ))
  
  
  epochMonitor<-epochMonitor$new()
  result<-AE_model%>%fit(x=df_normal,y=df_normal,batch_size=128,epochs=15,callbacks=list(epochMonitor),verbose=0)
  
  encodered_df<-encoder_model%>%predict(df_normal,batch_size=128)
}



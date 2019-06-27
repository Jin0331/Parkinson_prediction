source("keras_function_multiple_NN.R")
# source("keras_function.R")
library_load()

setname <- geneset_load()[[2]]

# callback list
callback_list <- list(
  callback_early_stopping(monitor = "val_acc", patience = 20),
  callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.15, patience = 5)
)

build_model <- function(){
  model <- keras_model_sequential() %>%
    layer_dense(units = 1000, activation = "relu",input_shape = DNN_input$common$dim,
                kernel_regularizer = regularizer_l2(0.0009))  %>%
    layer_dropout(0.5) %>%
    layer_dense(units = 1000, activation = "relu",
                kernel_regularizer = regularizer_l2(0.0009)) %>%
    layer_dropout(0.5) %>%
    layer_dense(units = 1, activation = "sigmoid")
  
  model %>% compile(
    optimizer = "rmsprop",#optimizer_rmsprop(lr = 1e-5)
    loss = "binary_crossentropy",
    metrics = "accuracy"
  )
}

for(gene_name in c(2,3,14:15, 17, 22, 24)){
  DNN_input <- gene_row_extract(set_name = setname[gene_name])
  model <- build_model()
  
  # cross validation
  k <- 4
  indices <- sample(1:nrow(DNN_input$common$x_train))
  folds <- cut(1:length(indices), breaks = k, labels = F)
  
  num_epoch <- 1000
  all_score <- tibble(.rows = 0, loss = 0, acc = 0)
  
  for(index in 1:k){
    cat("processing fold #", index, "\n")
    
    val_indices <- which(folds == index, arr.ind = T)
    val_data <- DNN_input$common$x_train[val_indices,]
    val_targets <- DNN_input$common$y_train[val_indices]
    
    partial_train_data <- DNN_input$common$x_train[-val_indices,]
    partial_train_targets <- DNN_input$common$y_train[-val_indices]
    
    model <- build_model()
    model %>% fit(partial_train_data, partial_train_targets,
                  epochs = num_epoch, batch_size = 30, verbose = 0)
    
    results <- model %>% evaluate(val_data, val_targets, verbose = 0)
    all_score <- all_score %>% bind_rows(., tibble(loss = results$loss, acc = results$acc))
    
  }
  
  write_delim(x = all_score, path = paste0(setname[gene_name], "_val_score_common_snp.txt"), delim = "\t")
}


# history <-  model %>% fit(DNN_input$common$x_train, 
#                           DNN_input$common$y_train, 
#                           epochs = 100, batch_size = 10, shuffle = T, validation_split = 0.2, callbacks = callback_list) 

source("C://Users/JINOO/Desktop/Parkinson_prediction/ver_1/keras_function_multiple_NN.R")
# source("keras_function.R")
library_load()

# save image 생성 ====
for(test_type in c("fisher", "logit")){
  for(p_value_type in c("UNADJ", "BONF")){
    for(dim_type in c("raw", "dim")){
      train_data <- varinat_selection(test_type = test_type, p_value_type = p_value_type, dim_type = dim_type)
      save.image(file = paste("common",test_type,p_value_type, dim_type, ".RData",sep = "_"))
    }
  }
}

# build_model ====
callback_list <- list(
  callback_early_stopping(monitor = "val_loss", patience = 4)
  # ,callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.01, patience = 5, verbose = 1)
)

build_model_FNN <- function(){
  model <- keras_model_sequential() %>%
    layer_dense(units = 500, activation = "relu", input_shape = train_data$dim) %>%
    layer_dropout(0.3) %>% 
    layer_dense(units = 64, activation = "relu") %>%
    layer_dropout(0.3) %>% 
    
    layer_dense(units = 1, activation = "sigmoid")
  
  model %>% compile(
    optimizer = optimizer_adam(),
    loss = "binary_crossentropy",
    metrics = "accuracy"
  )
}

# build_model_CNN <- function(){
#   model <- keras_model_sequential() %>% 
#     
# }



# fit ====
model_path <- "C://Users/JINOO/Desktop/Parkinson_prediction/R_data/model_save/"
model <- build_model_FNN()
history <-  model %>% fit(train_data$x_train,
                          train_data$y_train,
                          epochs = 20, batch_size = 20, shuffle = T, validation_split = 0.2
                          ,callbacks = callback_list
)

plot(history)

name <- "fisher_FNN(3)_0807"
save_model_hdf5(object = model, filepath = paste0(model_path, name, ".hdf5"))



# k-fold ======
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

ggsave()


# PCA
install.packages("DAAG")
library(DAAG)



minmax <- function(x) (x - min(x) / (max(x) - min(x)))
pca_data <- bind_cols(tibble(pheno = train_data$y_train), as_tibble(train_data$x_train)) 
pca_data_scale <- apply(pca_data[, 2:ncol(pca_data)], 2, minmax)

pca <- prcomp(x = pca_data_scale)

qplot(x = 1:4222, y = cumsum(pca$sdev)/sum(pca$sdev), geom = "line")
ggplot(as_tibble(pca$x), aes(x = PC1231, y = PC122, col = factor(pca_data$pheno))) + geom_point()

trans
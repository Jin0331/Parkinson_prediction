# library function
library_load <- function(){
  library(data.table)
  library(tidyverse)
  library(keras)
  library(SKAT);library(tidyverse);library(data.table)
  library(parallel);library(doParallel);library(progress);library(glue)
}

# pre-processing function
geneset_load <- function(){
  print("Geneset load_0428!")
  return_list <- list()
  geneset_merge <- fread(file = "NeuroX/1 gene sets 190428_exclude_ARX_for SKAT_O.txt", header = T, 
                         sep = "\t", stringsAsFactors = F, data.table = F) %>% as_tibble()
  
  one_hot <- lapply(X=select(geneset_merge, -HGNC, -M2, -'M2-GBA-LRRK2', -'M2_Brain_gene', -'M2_NOT_brain'), 
                    FUN = function(x){ifelse(nchar(x) > 0, 1, 0)}) %>% bind_cols()
  geneset_onehot <- cbind(GENE=geneset_merge$HGNC, one_hot) %>% mutate(GENE = as.character(GENE)) %>% 
    cbind(M2 = as.numeric(geneset_merge$M2)) %>% cbind('M2-GBA-LRRK2'= as.numeric(geneset_merge$`M2-GBA-LRRK2`)) %>%
    cbind('M2_Brain_gene' = as.numeric(geneset_merge$`M2_Brain_gene`)) %>% 
    cbind('M2_NOT_brain'= as.numeric(geneset_merge$`M2_NOT_brain`))
  
  geneset_onehot <- geneset_onehot %>% mutate(., 'M2_Brain_gene-LRRK2-GBA' = M2_Brain_gene)
  geneset_onehot$`M2_Brain_gene-LRRK2-GBA`[[554]] <- 0
  geneset_onehot$`M2_Brain_gene-LRRK2-GBA`[[755]] <- 0
  
  for(geneset in 2:ncol(geneset_onehot)){
    for(k in 1:nrow(geneset_onehot)){
      if(geneset_onehot[k, geneset] >= 1){
        geneset_onehot[k,geneset] <- geneset_onehot[k,1]
      } else{
        geneset_onehot[k,geneset] <- NA
      }
    }
  }
  
  # CHD, DEG adding
  
  CHD <- fread(file = "NeuroX/CHD_0424.txt", header = T, 
               sep = "\t", stringsAsFactors = F, data.table = F) %>% as_tibble() %>% 
    bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.)), CHD = NA))
  geneset_onehot <- bind_cols(geneset_onehot, CHD)
  
  DEG <- fread(file = "NeuroX/DEG_geneset.txt", header = T, sep = "\t", stringsAsFactors = F, data.table = F)
  
  DEG1 <- DEG$DEG1 %>% as_tibble() %>% bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(DEG1) <- "DEG1"
  DEG2 <- DEG$DEG1YJK[DEG$DEG1YJK != ""] %>% as_tibble() %>% bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(DEG2) <- "DEG1YJK"
  DEG3 <- DEG$DEG1Mito[DEG$DEG1Mito != ""] %>% as_tibble() %>% bind_rows(., tibble(.rows = (nrow(geneset_merge) - nrow(.))));colnames(DEG3) <- "DEG1Mito"
  
  geneset_onehot <- geneset_onehot %>% bind_cols(., DEG1) %>% bind_cols(., DEG2) %>% bind_cols(., DEG3)
  
  
  return_list[[1]] <- geneset_onehot %>% as.list()
  return_list[[2]]<- colnames(geneset_onehot)
  
  return(return_list)
} ## 1 = geneset_list, 2 = col_name
fix_load <- function(data_name){
  print(paste0(data_name, " fix load!"))
  test_fix <- fread(file = "NeuroX/NeuroX_fix.txt", sep = "\t", header = T, stringsAsFactors = F, data.table = F) %>% 
    as_tibble()
  test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(X)", replacement = "23")
  test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(Y)", replacement = "24")
  test_fix <- test_fix %>% mutate(ID2 = paste(CHROM, POS, sep = ":"))
  
  test_id <- test_fix$ID;test_ch <- test_fix$CHROM;test_pos <- test_fix$POS
  for(change in 1:length(test_id)){
    if(test_id[change] == ""){
      test_id[change] <- paste0(test_ch[change], ":", test_pos[change])
    }
  }
  test_fix$ID <- test_id
  
  return(test_fix)
} # data_name, "IPDGC", "NeuroX"

gene_row_extract <- function(set_name, maf=1){

  load("NeuroX/data_test_fix_0329.RData")
  fix <- test_fix[[2]]
  geneset <- geneset_load()
  geneset<- unlist(geneset[[1]][names(geneset[[1]]) == set_name]) %>% .[!is.na(.)]
  
  if(maf <= 0.03){
    snp_ <- lapply(X = geneset, FUN = function(gene){
      temp <- filter(fix, Gene.knownGene  == gene & ExAC_AF <= maf & CADD13_PHRED > 12.37)
      return(temp)
    }) %>% bind_rows() ###  rare-variant
    
  }else{
    snp_ <- lapply(X = geneset, FUN = function(gene){
      temp <- filter(fix, Gene.knownGene  == gene & ExAC_AF >= 0.05)
      return(temp)
    }) %>% bind_rows() ###  common-variant
  }
  
  
  fwrite(x = tibble(snp_$ID), file = paste0(set_name,"_snp_extract.txt"), sep = "\t", col.names = F)
  system(glue("D://tool/plink.exe --file NeuroX/NeuroX_row --extract {name}_snp_extract.txt --recode A --out {name}_NeuroX_dosage", 
              name = set_name), show.output.on.console = F)
  system(glue("D://tool/plink.exe --file NeuroX/NeuroX_row --logistic beta --out {name}_NeuroX_dosage_logit", 
              name = set_name), show.output.on.console = F)
  
  temp <- fread(file = "NeuroX/NeuroX.cov", header = T, data.table = F) %>% select(., -FID, -IID) # phenotype
  another_p1 <- as.data.frame(scale(temp)) # phenotype
  
  dosage_raw <- fread(file = paste0(set_name,"_NeuroX_dosage.raw"), header = T, nThread = 20) %>% 
    as_tibble() %>% select(.,-SEX, -FID, -IID, -PAT, -MAT)
  
  phenotype <- dosage_raw %>% select(., PHENOTYPE) %>% bind_cols(., as.data.frame(another_p1))
  dosage_raw <- dosage_raw %>% select(., -PHENOTYPE)
  
  ## logit(pij)
  dosage_beta <- fread(file = paste0(set_name, "_NeuroX_dosage_logit.assoc.logistic"), header = T) %>%
    filter(., TEST == "ADD")
  
  for(col_len in 1:ncol(dosage_raw)){
    index <- which(str_detect(dosage_beta$SNP, pattern = paste0(str_sub(colnames(dosage_raw)[col_len], end = -3), "$")))
    beta <- dosage_beta[index, ]$BETA
    dosage_raw[col_len] <- dosage_raw[col_len] * beta
  }
  
  dosage_raw[is.na(dosage_raw)] <- 0.000
  dosage_raw_scale <- as.data.frame(scale(dosage_raw))
  dosage_raw_scale[is.nan(dosage_raw_scale)] <- 0
  
  result_raw <- bind_cols(phenotype, dosage_raw)
  
  
  return(train_test_split(result_raw))
}

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))


train_test_split <- function(dataDF){
  # dataDF <- fread(file = file, header = T, stringsAsFactors = F)
  control <- dataDF %>% filter(., PHENOTYPE == 1)
  case <- dataDF %>% filter(., PHENOTYPE == 2)
  
  # case control split
  control_sample <- sample.int(n = nrow(control), size = floor(.7 * nrow(control)), replace = F)
  control_train <- control[control_sample, ] %>% mutate(., PHENOTYPE = 0)
  control_test <- control[-control_sample, ] %>% mutate(., PHENOTYPE = 0)
  
  case_sample <- sample.int(n = nrow(case), size = floor(.7 * nrow(case)), replace = F)
  case_train <- case[case_sample, ] %>% mutate(., PHENOTYPE = 1)
  case_test <- case[-case_sample, ] %>% mutate(., PHENOTYPE = 1)
  
  train <- rbind(control_train, case_train)
  train <- train[sample(1:nrow(train)), ] # shuffle
  
  test <- rbind(control_test, case_test)
  test <- test[sample(1:nrow(test)), ] # shuffle
  
  
  # x, y split
  x_train <- train %>% select(., -PHENOTYPE)
  y_train <- train %>% select(., PHENOTYPE)
  
  x_test <- test %>% select(., -PHENOTYPE)
  y_test <- test %>% select(., PHENOTYPE)
  
  
  return(list(x_train = as.matrix(x_train), y_train = as.integer(y_train$PHENOTYPE), 
              x_test = as.matrix(x_test), y_test = as.integer(y_test$PHENOTYPE),
              dim = ncol(x_train)))
}


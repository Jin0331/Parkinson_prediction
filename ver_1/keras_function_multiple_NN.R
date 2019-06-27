# library function
library_load <- function(){
  library(data.table)
  library(tidyverse)
  library(keras)
  library(SKAT);library(tidyverse);library(data.table)
  library(parallel);library(doParallel);library(progress);library(glue)
}

# association test
assoc_test <- function(){
  system("D://tool/plink.exe --file NeuroX/NeuroX_row --logistic beta --out NeuroX/NeuroX_dosage_logit", show.output.on.console = F)
}

# min max scale
min_max_scale <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}

# pre-processing function
geneset_load <- function(){
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
is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))
}
fix_load <- function(data_name){
  # print(paste0(data_name, " fix load!"))
  # test_fix <- fread(file = "NeuroX/NeuroX_fix.txt", sep = "\t", header = T, stringsAsFactors = F, data.table = F) %>%
  #   as_tibble()
  # test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(X)", replacement = "23")
  # test_fix$CHROM <- gsub(x = test_fix$CHROM, pattern = "(Y)", replacement = "24")
  # test_fix <- test_fix %>% mutate(ID2 = paste(CHROM, POS, sep = ":"))
  # 
  # test_id <- test_fix$ID;test_ch <- test_fix$CHROM;test_pos <- test_fix$POS
  # for(change in 1:length(test_id)){
  #   if(test_id[change] == ""){
  #     test_id[change] <- paste0(test_ch[change], ":", test_pos[change])
  #   }
  # }
  # test_fix$ID <- test_id
  # 
  # return(test_fix)
} # data_name, "IPDGC", "NeuroX"

gene_row_extract <- function(set_name, rare_maf = 0.03, common_maf = 0.05, seed = sample(x = 1:1000, size = 1)){
  
  set.seed(seed) # set seed
  
  load("NeuroX/data_test_fix_0329.RData")
  fix <- test_fix[[2]]
  geneset <- geneset_load()
  geneset<- unlist(geneset[[1]][names(geneset[[1]]) == set_name]) %>% .[!is.na(.)]
  
  variant_rare <- lapply(X = geneset, FUN = function(gene){
      temp <- filter(fix, Gene.knownGene  == gene & plink_maf <= rare_maf)
      return(temp)
  }) %>% bind_rows() ###  rare-variant
    
  variant_common <- lapply(X = geneset, FUN = function(gene){
      temp <- filter(fix, Gene.knownGene  == gene & plink_maf >= common_maf)
      return(temp)
  }) %>% bind_rows() ###  common-variant

  # case control split for variant(rare, common, age_sex, mds_fmiss)
  rare_snp <- dosage_calc(snp_ = variant_rare, state = "rare", set_name = set_name)
  common_snp <- dosage_calc(snp_ = variant_common, state = "common", set_name = set_name)
  
  rare_control <- rare_snp[[1]] %>% filter(., PHENOTYPE == 1)
  rare_case <- rare_snp[[1]] %>% filter(., PHENOTYPE == 2)
  
  common_control <- common_snp[[1]] %>% filter(., PHENOTYPE == 1)
  common_case <- common_snp[[1]] %>% filter(., PHENOTYPE == 2)

  age_mds_fmiss_control <- rare_snp[[2]] %>% filter(., PHENOTYPE == 1)
  age_mds_fmiss_case <- rare_snp[[2]] %>% filter(., PHENOTYPE == 2)
  
  ### 
  control_sample <- sample.int(n = nrow(rare_control), size = floor(.7 * nrow(rare_control)), replace = F)
  case_sample <- sample.int(n = nrow(rare_case), size = floor(.7 * nrow(rare_case)), replace = F)
  ### 
  
  rare_result <- train_test_split(control_df = rare_control, case_df = rare_case, 
                                  control_sp_num = control_sample, case_sp_num = case_sample)
  
  common_result <- train_test_split(control_df = common_control, case_df = common_case, 
                                    control_sp_num = control_sample, case_sp_num = case_sample)
  
  age_mds_fmiss_result <- train_test_split(control_df = age_mds_fmiss_control, case_df = age_mds_fmiss_case, 
                                     control_sp_num = control_sample, case_sp_num = case_sample)
  
  return(list(rare = rare_result, common = common_result, age_mds_fmiss = age_mds_fmiss_result))
}
dosage_calc <- function(snp_, state, set_name){
  fwrite(x = tibble(snp_$ID), file = paste0(set_name,"_",state,"_snp_extract.txt"), sep = "\t", col.names = F)
  system(glue("D://tool/plink.exe --file NeuroX/NeuroX_row --extract {name}_{state}_snp_extract.txt --recode A --out {name}_{state}_NeuroX_dosage", 
              name = set_name, state = state), show.output.on.console = F)
  
  dosage_raw <- fread(file = paste0(set_name,"_", state,"_NeuroX_dosage.raw"), header = T, nThread = 20) %>% 
    as_tibble() %>% select(., -FID, -IID, -PAT, -MAT) 
  
  COVARIATE <- fread(file = "NeuroX/NeuroX.cov", header = T, data.table = F) %>% select(., -FID, -IID) # AGE, 4 MDS
  AGE_MDS_F_MISS <- COVARIATE %>% scale() %>% as.data.frame() %>% cbind(PHENOTYPE = dosage_raw$PHENOTYPE, .) # 4 MDS & F_MISS
  
  
  PHENOTYPE<- dosage_raw %>% select(., PHENOTYPE)
  dosage_raw <- dosage_raw %>% select(., -PHENOTYPE, -SEX)
  
  ## logit(pij)
  dosage_beta <- fread(file = "NeuroX/NeuroX_dosage_logit.assoc.logistic", header = T) %>%
    filter(., TEST == "ADD")
  
  # dosage_beta$BETA <- min_max_scale(dosage_beta$BETA)
  
  for(col_len in 1:ncol(dosage_raw)){
    index <- which(str_detect(dosage_beta$SNP, pattern = paste0(str_sub(colnames(dosage_raw)[col_len], end = -3), "$")))
    beta <- dosage_beta[index, ]$BETA
    dosage_raw[col_len] <- dosage_raw[col_len] * beta
  }
  
  delete_count <- c()
  for(na_col_len in 1:ncol(dosage_raw)){
    temp <- dosage_raw[na_col_len]
    temp[is.na(temp)] <- -999
    colnames(temp) <- "temp"
    
    if(nrow(filter(temp, temp == -999)) >= 10000){
      delete_count <- c(delete_count, na_col_len)
    }
  }
  
  if(length(delete_count) == 0){
    print("no snp")
  }else{
    dosage_raw <- dosage_raw[,-delete_count]  
  }
  
  dosage_raw[is.na(dosage_raw)] <- 0.000
  # dosage_raw_scale <- as.data.frame(scale(dosage_raw))
  # dosage_raw_scale[is.nan(dosage_raw_scale)] <- 0.000
  
  result_raw <- list()
  result_raw[[1]] <- bind_cols(PHENOTYPE, dosage_raw)
  result_raw[[2]] <- AGE_MDS_F_MISS
  
  return(result_raw)
   
}
train_test_split <- function(control_df, case_df, control_sp_num, case_sp_num){

  control_train <- control_df[control_sp_num, ] %>% mutate(., PHENOTYPE = 0)
  control_test <- control_df[-control_sp_num, ] %>% mutate(., PHENOTYPE = 0)
  case_train <- case_df[case_sp_num, ] %>% mutate(., PHENOTYPE = 1)
  case_test <- case_df[-case_sp_num, ] %>% mutate(., PHENOTYPE = 1)
  
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
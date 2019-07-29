library(tidyverse)
library(data.table)



# 0725 variant 조사 및 variant 셀렉션
{
  # (1) variant 수 확인
  data_path <- "C://Users/JINOO/Desktop/Parkinson_prediction/data/NeuroX_QC/"
  NeuroX_freq <- fread(paste(data_path,"NeuroX_freq.frq", sep = "\\")) %>% as_tibble()
  NeuroX_freq %>% filter(MAF >= 0.05) %>% nrow() # 35,642 variants
  NeuroX_freq %>% filter(MAF >= 0.005 & MAF < 0.05) %>% nrow() # 15,635 variants
  NeuroX_freq %>% filter(MAF < 0.005) %>% nrow() # 127,110 variants
  
  # (2) Fisher's exact test
  NeuroX_fisher <- fread(paste(data_path, "NeuroX_fisher.assoc.fisher", sep = "\\")) %>% as_tibble()
  View(NeuroX_fisher) 
  
  # (3) join___(1),(2)
  NeuroX_maf_fisher <- load_freq_maf()
  
  # (4) 
}


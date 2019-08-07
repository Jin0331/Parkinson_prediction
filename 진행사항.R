# 0725 variant 조사 및 variant 셀렉션 ======

  # (1) variant 수 확인
  data_path <- "C://Users/JINOO/Desktop/Parkinson_prediction/data/NeuroX_QC/"
  NeuroX_freq <- fread(paste(data_path,"NeuroX_freq.frq", sep = "\\")) %>% as_tibble()
  NeuroX_freq %>% filter(MAF >= 0.05) %>% nrow() # 35,642 variants
  NeuroX_freq %>% filter(MAF >= 0.005 & MAF < 0.05) %>% nrow() # 15,635 variants
  NeuroX_freq %>% filter(MAF < 0.005) %>% nrow() # 127,110 variants
  
  # (2) Fisher's exact test
  NeuroX_fisher <- fread(paste(data_path, "NeuroX_fisher.assoc.fisher", sep = "\\")) %>% as_tibble()
  NeuroX_fisher <- fread(paste(data_path, "NeuroX_fisher_adjust.assoc.fisher.adjusted", sep = "\\")) %>% as_tibble()
  View(NeuroX_fisher)
  
  # (3) join___(1),(2)
  NeuroX_maf_fisher <- load_freq_maf()
  
  # (4) 


# 0729 assoc과 logistic??? =====
  # 다수의 논문의 경우 assoc이 아닌 logistic regression 으로 하는 경우가 많음. 
  # unadjust p-value가 아닌 adjust p-value로, bonferroni 방법을 적용하여 QC를 진행함.
  # inner join을 통해 각각의 파일에서 필요한 unadjust p-value와 bonferroni p-value, OR 구하기(완)
  
  load_freq_maf(type = "fisher") %>% View()
  


# 0730 dosage 어떻게 할 것인가 =====
  # dosage_clac 수정


# 0801 ===== 
  # case 5384, control 5689


# 0802 ======
  # train 2D -> 3D
  # wild type, Heterozygous mutant, Homozygous Mutant


# 0807  =====
  # additive type (0, 1, 2)와 wild type, Hetero, Homo -> one-hot encoding 2가지 버젼추가
  
  
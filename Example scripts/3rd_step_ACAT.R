AD_DPR <- $path_to_TWAS_DPR
AD_EN <- $path_to_TWAS_EN
AD_FUSION <- $path_to_TWAS_FUSION

AD_DPR <- AD_DPR[,c(1:5,10)]
AD_EN <- AD_EN[,c(1:5,10)]
AD_FUSION <- AD_FUSION[,c(1:5,10)]

colnames(AD_DPR)[6] <- "p_DPR"
colnames(AD_EN)[6] <- "p_EN"
colnames(AD_FUSION)[6] <- "p_FUSION"

library(dplyr)
result <- AD_DPR %>%
  full_join(AD_EN, by = c("CHROM", "GeneStart", "GeneEnd", "TargetID", "GeneName")) %>%
  full_join(AD_FUSION, by = c("CHROM", "GeneStart", "GeneEnd", "TargetID", "GeneName"))

result <- result %>%
  mutate(across(6:8, ~ifelse(. == 1, 2/3, .)))

devtools::install_github("yaowuliu/ACAT")

library(ACAT)
ACAT_withNA = function(p_vec){
  p_vec_noNA = p_vec[is.na(p_vec) == F]
  ACAT(p_vec_noNA)
}

result$p_ACAT <- NA
for (i in 1:nrow(result)){
ps <- as.vector(t(result[i,6:8]))
p_new <- ACAT_withNA(ps)
result[i,9] <- p_new
}

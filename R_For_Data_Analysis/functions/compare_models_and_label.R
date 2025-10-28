#-------------------------------------------------------------------------------
# compare models, add significant and foldchange labels
#-------------------------------------------------------------------------------
compare_models_and_label <- function(res1, res2, res3,
                                     key_col = "Genus",
                                     fc_cutoff = 0.01, 
                                     fdr_cutoff = 0.05) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(rlang)
  })
  
  key_sym <- sym(key_col)
  
  # combine 3 models
  res_compare <- res1 %>%
    select(!!key_sym, Estimate1 = Estimate, FDR1 = FDR) %>%
    left_join(
      res2 %>% select(!!key_sym, Estimate2 = Estimate, FDR2 = FDR),
      by = key_col
    ) %>%
    left_join(
      res3 %>% select(!!key_sym, Estimate3 = Estimate, FDR3 = FDR),
      by = key_col
    )
  
  # add labels
  res_compare <- res_compare %>%
    mutate(
      Sig1 = FDR1 < fdr_cutoff,
      Sig2 = FDR2 < fdr_cutoff,
      Sig3 = FDR3 < fdr_cutoff,
      Category = case_when(
        Sig1 & Sig2 & Sig3 ~ "All_Significant",
        Sig1 & !Sig2 & !Sig3 ~ "Only_res1",
        !Sig1 & Sig2 & !Sig3 ~ "Only_res2",
        !Sig1 & !Sig2 & Sig3 ~ "Only_res3",
        Sig1 & Sig2 & !Sig3 ~ "Res1_Res2",
        Sig1 & !Sig2 & Sig3 ~ "Res1_Res3",
        !Sig1 & Sig2 & Sig3 ~ "Res2_Res3",
        TRUE ~ "None"
      ),
      FC_label = ifelse(Estimate1 > fc_cutoff | Estimate1 < -fc_cutoff, 
                        "FoldChange > 2", 
                        "Stable")
    )
  
  # results
  sig_fc <- subset(res_compare, 
                   FC_label == "FoldChange > 2" & Category == "All_Significant")
  
  return(list(full = res_compare, filtered = sig_fc))
}


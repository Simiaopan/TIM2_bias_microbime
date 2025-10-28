
prm.ls <- list()

prm.ls$general <- list(seed_all = 20021, 
                    lib = c("tidyverse"),
                    data_path = "data/",
                    meta_path = "data/Cleaned_Meta.xlsx",
                    report_path = "report/denoising_stats.tsv"
                    )

prm.ls$sample_QC <- list(
                     functions = c("functions/veiw_report.R"),
                     final_reads_threshold = 1000,
                     final_percent_threshold = 30
                     )

prm.ls$data_pre <- list(tax_lvls = c("ASV", "Genus", "Family", "Phylum"),
                       min_otu = 10,
                       min_sample =  20
                        )


prm.ls$pcoa <- list(norm_m = "Rare",
                    distance = c("bray", "jaccard", "unifrac", "wunifrac")
                             )
prm.ls$RPCA <- list(dataset_names = c("_all","dtt1","t1","dtt1ts", 
                                      "t1_T0","t1_T24","t1_T48","t1_T72"),
                    input_dir = "data/rpca_input",
                    output_dir = "data/rpca_output_90")
 
prm.ls$LMM <- list(mean_ab = "asv_mean_abund_raw")

prm.ls$prediction <- list(ref_path = "tools/Tax4Fun2_ReferenceData_v2",
                          TaxFun2_dir = "D:/Internship/TIM2_bias_microbime/R_For_Data_Analysis/data/Tax4Fun2_tmp")

save(file = "rdata/PRM.Rdata", list = c("prm.ls"))


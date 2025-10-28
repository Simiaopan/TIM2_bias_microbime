library(dplyr)
library(stringr)
library(readxl)
library(writexl)

##------------------ set seed -----------------------------
set.seed(20021)

##------------------- Clean and generate columns like fastq_code ---------------------

# Read the raw Excel data
df_raw <- read_excel("D:/internship/TIM2_bias_microbime/R_For_Data_Analysis/Raw_Meta.xlsx", sheet = "All_Meta", col_types = "text")

# 1️⃣ Remove rows where Seq_Nr is empty and only keep numbers
df_raw <- df_raw %>% filter(!is.na(Seq_Nr) & Seq_Nr != "")

df_raw <- df_raw %>%
  mutate(Seq_Nr_clean = gsub("[^0-9]", "", Seq_Nr))


# ##----- Detailed replacement for SIEM ----------------------------
# # Replace "control" etc. in df_raw$Testproduct with "SIEM"
# df_raw$Testproduct <- df_raw$Testproduct %>%
#   str_replace_all(pattern = "^control (no treatment)$|^Control$|^Cont1$|^SEIM$", replacement = "SIEM")
# 
# # For specific df_raw$Project, replace "SIEM" in df_raw$Testproduct with "SIEMNC"
# df_raw$Testproduct <- ifelse(
#   df_raw$Project %in% c("P0076-24", "P0042") & df_raw$Testproduct == "SIEM",
#   "SIEM-NC",
#   df_raw$Testproduct
# )

##------------------ Detailed replacement for P0087 --------------
# Only operate on rows where Project equals P0087
df_raw$Donor <- ifelse(
  df_raw$Project == "P0087", 
  paste(df_raw$Donor, gsub("[^0-9]", "", df_raw$Sample), sep = "-"), 
  df_raw$Donor
)


##-- Cleaning and modification Process --------------------------------------
df_cleaned <- df_raw %>%
  # Step 1: Run_Project
  mutate(Run_Project = paste(RUN, Project, sep = "_")) %>%
  
  # Step 2-3: Extract Time and Unit numbers
  mutate(
    Time_clean = str_extract(Time, "\\d+"),
    Time_clean = ifelse(!is.na(Time_clean), paste0("T", Time_clean), NA),
    Unit_clean = str_extract(Unit, "\\d+"),
    Unit_clean = ifelse(!is.na(Unit_clean), paste0("U", Unit_clean), NA)
  ) %>%
  
  # Create IDs
  mutate(Sample_ID = ifelse(is.na(Testproduct) & is.na(Time_clean),
                            NA,
                            paste(Project, Donor, Testproduct, Time_clean, sep = "_"))) %>%
  
  mutate(Single_ID = ifelse(is.na(Sample_ID),
                            NA,
                            paste(Project, Donor, Testproduct, Time_clean, batch, sep = "_"))) %>%
  
  mutate(Batch_ID = ifelse(!is.na(Single_ID), paste(Project, Donor, Testproduct, batch, sep="_"),
                            NA)) %>%
  
  mutate(Treat_ID = ifelse(is.na(Single_ID), NA, 
                           paste(Project, Donor, Testproduct, sep = "_"))) %>%
  
  # Step 5: Custom logic for generating fastq_code
  rowwise() %>%
  mutate(
    fastq_code = {
      run_clean <- as.character(RUN)
      snr <- as.character(Seq_Nr_clean)
      
      if (nchar(snr) < 4) {
        snr
      } else {
        stripped <- str_remove(snr, paste0("^", run_clean))
        stripped <- str_replace(stripped, "^0+", "")  # Remove leading zeros
        ifelse(stripped == "", snr, stripped)
      }
    }
  ) %>%
  ungroup() %>%
  
  # Pad RUN to three digits
  mutate(RUN = str_pad(RUN, width = 3, side = "left", pad = "0")) %>%
  
  # Initialize the redo replacement source column
  mutate(Replaced_From = NA_character_) %>%
  
  # Select the output columns
  select(
    RUN,
    Run_Project,
    Project,
    Testproduct,
    Batch = batch,
    Duplication,
    TimePattern,
    Batch_ID,
    Treat_ID,
    Sample_ID,
    Single_ID,
    fastq_code,
    Time = Time_clean,
    Unit = Unit_clean,
    Donor,
    Date,
    Seq_Nr = Seq_Nr_clean,
    Redo,
    Replaced_From
  )

# Optional output intermediate results
#write.csv(df_cleaned, "D:/internship/Meta_data/Meta_Data_Clean.csv", row.names = FALSE)

##--------------------- Handle redo replacement -------------------------------------

##=================== Chain Replacement Process ===================##

# Ensure Seq_Nr and Redo are characters
df_cleaned <- df_cleaned %>%
  mutate(
    Seq_Nr = as.character(Seq_Nr),
    Redo   = as.character(Redo),
    Redo   = ifelse(Redo == "" | is.na(Redo), NA, Redo)
  )

# 1️⃣ Build redo_map: old_sample → new_sample
redo_links <- df_cleaned %>% filter(!is.na(Redo) & Redo != "")
redo_map <- setNames(redo_links$Seq_Nr, redo_links$Redo)  # Redo is the old sample, who replaced it

# 2️⃣ Find chain starts (old samples that have never been replaced)
chain_starts <- setdiff(redo_links$Redo, redo_links$Seq_Nr)
cat("✅ Found redo chain starts: ", paste(chain_starts, collapse = ", "), "\n\n")

# 3️⃣ Build replacement chains (old → new)
replacement_chains <- list()
for (start in chain_starts) {
  chain <- start
  current <- start
  while (current %in% names(redo_map)) {
    current <- redo_map[[current]]
    chain <- c(chain, current)
  }
  if (length(chain) > 1) {
    replacement_chains[[start]] <- chain
  }
}

cat("✅ Number of replacement chains: ", length(replacement_chains), "\n")
for (start in names(replacement_chains)) {
  cat("🔁 Replacement chain: ", paste(replacement_chains[[start]], collapse = " -> "), "\n")
}
cat("\n==== Replacement chain execution begins ====\n\n")

# Initialize
rows_to_delete <- c()
df_final <- df_cleaned

for (chain in replacement_chains) {
  final_seq <- tail(chain, 1)
  rep_seqs <- head(chain, -1)
  
  cat("🔁 Current chain: ", paste(chain, collapse = " -> "), "\n")
  cat("  ✅ Final to keep: ", final_seq, "\n")
  
  tgt_idx <- which(df_final$Seq_Nr == final_seq)
  if (length(tgt_idx) != 1) {
    cat("  ⚠️ Final sample ", final_seq, " not found or not unique in df_final, skipping\n\n")
    next
  }
  
  # Get the oldest and newest samples
  oldest_row <- df_cleaned %>% filter(Seq_Nr == chain[1])
  newest_row <- df_cleaned %>% filter(Seq_Nr == final_seq)
  
  if (nrow(oldest_row) == 1 && nrow(newest_row) == 1) {
    oldest_row <- oldest_row[1, ]
    newest_row <- newest_row[1, ]
    
    # ✅ Keep the information of the oldest sample
    cols_keep_old <- c("Time","Project", "Unit", "Testproduct", "Date", "Donor","Batch","Duplication","TimePattern","Batch_ID", "Treat_ID","Sample_ID","Single_ID")
    for (col in cols_keep_old) {
      df_final[tgt_idx, col] <- oldest_row[[col]]
    }
    
    # ✅ Use the information of the newest replacer
    cols_use_new <- c("Seq_Nr", "fastq_code", "Run_Project", "RUN")
    for (col in cols_use_new) {
      df_final[tgt_idx, col] <- newest_row[[col]]
    }
    
    # ✅ Mark source (the second-to-last)
    if (length(chain) >= 2) {
      df_final[tgt_idx, "Replaced_From"] <- chain[length(chain) - 1]
    }
    
    cat("    ✅ Final sample updated: ", final_seq, "\n")
  }
  
  # Delete intermediate samples
  for (rep_seq in rep_seqs) {
    del_idx <- which(df_final$Seq_Nr == rep_seq)
    if (length(del_idx) == 1) {
      rows_to_delete <- c(rows_to_delete, del_idx)
      cat("    ❌ Marked for deletion: ", rep_seq, "\n")
    }
  }
  
  cat("\n")
}

# Delete the intermediate replaced samples
rows_to_delete <- unique(rows_to_delete[rows_to_delete > 0 & rows_to_delete <= nrow(df_final)])
if (length(rows_to_delete) > 0) {
  df_final <- df_final[-rows_to_delete, ]
  message("✅ redo replacement successful, ", length(rows_to_delete), " old samples deleted.")
} else {
  message("⚠️ No valid redo rows to delete, rows_to_delete is empty or invalid.")
}

# Sort
df_final <- df_final %>%
  mutate(RUN_numeric = as.numeric(RUN)) %>%
  arrange(RUN_numeric) %>%
  select(-RUN_numeric)

# Tracking table
replacement_df <- lapply(replacement_chains, function(chain) {
  tibble(
    Final_Seq_Nr = chain[1],
    Replacement_Chain = paste(rev(chain), collapse = " → ")
  )
}) %>% bind_rows()

# Log
redo_log <- df_cleaned %>%
  filter(!is.na(Redo)) %>%
  select(
    Seq_Nr, Redo, Run_Project, Batch_ID, Treat_ID, Sample_ID, Single_ID, fastq_code,
    Unit, Time, Donor, Testproduct, Date
  ) %>%
  arrange(Redo)

#filter in the end
df_final <- df_final %>%
    filter(!is.na(Time) & !is.na(Testproduct) & !is.na(Project))

#Add global ID
n <- nrow(df_final)
df_final$Global_ID <- sprintf("%04d", sample(1:n,n))
# Export
write_xlsx(
  list(
    "Meta_Data_Final"    = df_final,
    "Replacement_Chain"  = replacement_df,
    "Redo_Log_Raw"        = redo_log
  ),
  path = "D:/internship/TIM2_bias_microbime/R_For_Data_Analysis/data/Cleaned_Meta.xlsx"
)

# # Check Single_ID
# Single_Check <- as.data.frame(table(df_final$Single_ID))

# # Check Run_ID
# Batch_Check <- as.data.frame(table(df_final$Batch_ID))

# #Check Sample_ID
# Sample_Check <- as.data.frame(table(df_final$Sample_ID))

# #Check trend_ID
# Treat_Check <- as.data.frame(table(df_final$Treat_ID))

# #Testproject_sum
# Testproduct_sum <- as.data.frame(table(df_final$Testproduct))

# #Time_sum
# Time_sum <- as.data.frame(table(df_final$Time))

#export report
write_xlsx(list(
  "Batch_ID" = Batch_Check,
  "Sample_ID" = Sample_Check,
  "Single_ID" = Single_Check,
  "Treat_ID" = Treat_Check,
  "Testproduct" = Testproduct_sum,
  "Time" = Time_sum
  ),
  path = "data/Meta_Report.xlsx"
)


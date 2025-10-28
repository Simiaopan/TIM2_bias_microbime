
load("rdata/prm.Rdata")

invisible(lapply(prm.ls$general$lib, require, character.only = TRUE))
library(readxl)

Meta_path <- prm.ls$general$meta_path
Report_path <- prm.ls$general$report_path
functions <- prm.ls$sample_QC$functions
lapply(functions,source)


##----view the dada report--------------------------
meta_for_QC <- read_xlsx(Meta_path,sheet=1)
Report <- read.delim(Report_path, comment.char = "#", check.names = FALSE)

view_report(Report_path)
print("Please change the filter parameters based on the report")

low_quality_samples <- Report[Report$`non-chimeric` < 5000, 1]
low_quality_no_prefix <- sub("^np", "", low_quality_samples)

# find Batch_ID
bad_batches <- meta_for_QC %>%
  filter(Global_ID %in% low_quality_no_prefix) %>%
  pull(Batch_ID) %>%
  unique()

drop_ids <- meta_for_QC %>%
  filter(Batch_ID %in% bad_batches) %>%
  pull(Global_ID)

drop_ids_with_prefix <- paste0("np", drop_ids)

# change Treat_ID to False
bad_trends <- meta_for_QC %>%
  filter(Global_ID %in% drop_ids) %>%
  pull(Treat_ID) %>%
  unique()

meta_for_QC <- meta_for_QC %>%
  mutate(Duplication = ifelse(Treat_ID %in% bad_trends, FALSE, Duplication))

# remove drop_ids
message("Removed ", length(drop_ids), " samples in ", length(bad_batches), " batches.")

meta <- meta_for_QC %>%
  filter(!(Global_ID %in% drop_ids))


##------------save the variables---------------------
save(meta, file = "rdata/filtered_meta.Rdata")

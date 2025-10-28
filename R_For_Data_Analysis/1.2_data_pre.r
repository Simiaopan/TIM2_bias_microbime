#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("rdata/prm.Rdata")
load("rdata/filtered_meta.Rdata")

seed = prm.ls$general$seed_all
set.seed(seed)

suppressPackageStartupMessages({
  library(phyloseq)
  library(qiime2R)
  library(dplyr)
  library(colorspace)
  library(Biostrings)
})

#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("functions/phy_shorten_tax_names.R")
source("functions/phy_css_norm.R")
source("functions/phy_clr_norm.R")
source("functions/make_splits.R")

#-------------------------------------------------------------------------------
# Variables 
#-------------------------------------------------------------------------------
data.path  <- prm.ls$general$data_path
meta.path <- prm.ls$general$meta_path

gr.col <- "Testproduct"
time.col <- "Time"
time.p  <- "TimePoint"
tax.lvls <- prm.ls$data_pre$tax_lvls
min.otu   <- prm.ls$data_pre$min_otu
min.s <- prm.ls$data_pre$min_sample

data.ls <- list()

#-------------------------------------------------------------------------------
# Read data/ Create phyloseq
#-------------------------------------------------------------------------------
ps <- qza_to_phyloseq(
  features = paste0(data.path, "asv_table.qza"), 
  tree     = paste0(data.path, "trees/rooted-tree.qza"), 
  taxonomy = paste0(data.path, "taxonomy_07.qza")
)

seqs <- readDNAStringSet(paste0(data.path, "dna-sequences.fasta"))
names(seqs) <- gsub(" .*", "", names(seqs))
seqs <- seqs[taxa_names(ps)]
ps <- merge_phyloseq(ps, refseq(seqs))

#-------------------------------------------------------------------------------
# Filter out taxa
#-------------------------------------------------------------------------------
otu <- as(otu_table(ps), "matrix")
if(!taxa_are_rows(ps)) otu <- t(otu)
keep_taxa <- rowSums(otu > min.otu) >= min.s
ps <- prune_taxa(keep_taxa, ps)

ps <- prune_taxa(!tax_table(ps)[, "Genus"] %in% "Mitochondria", ps)
ps <- prune_taxa(!tax_table(ps)[, "Genus"] %in% "Chloroplast", ps)
ps <- prune_taxa(!is.na(tax_table(ps)[, "Phylum"])[,1], ps)
ps <- prune_taxa(tax_table(ps)[, "Kingdom"] %in% c("d__Bacteria","d__Archaea"), ps)

#-------------------------------------------------------------------------------
# Add metadata
#-------------------------------------------------------------------------------
meta$Global_ID <- paste0("np",meta$Global_ID)
meta <- as.data.frame(meta)
rownames(meta) <- meta[["Global_ID"]]

over.samp <- intersect(sample_names(ps), rownames(meta))

# based on RPCA result, RUN 060 need to be filtered out!
if ("RUN" %in% colnames(meta)) {
  rm_samples <- rownames(meta[meta$RUN == "060", ])
  over.samp <- setdiff(over.samp, rm_samples)
}

ps   <- prune_samples(over.samp, ps)
ps   <- prune_taxa(taxa_sums(ps) > 0, ps)



meta <- meta[sample_names(ps), ] %>% 
  mutate(
    Time = as.numeric(gsub("^T", "", TimePoint)),
    
    TimePoint   = factor(TimePoint, levels = paste0("T", sort(unique(as.numeric(gsub("^T", "", TimePoint)))))),
    Batch_ID    = factor(Batch_ID),
    Treat_ID    = factor(Treat_ID),
    Sample_ID   = factor(Sample_ID),
    Single_ID   = factor(Single_ID),
    Run_Project = factor(Run_Project),
    Project     = factor(Project),
    
    Testproduct = factor(
      Testproduct,
      levels = c("SIEM", sort(unique(Testproduct[Testproduct != "SIEM"])))
    )
  )

data.ls[["meta"]]  <- meta
sample_data(ps) <- meta


################################################################################
# Prepare list with data 
################################################################################

for(i.lvl in tax.lvls)  {
  if(i.lvl == "ASV") { 
    ps.inst <- ps
  } else { 
    ps.inst <- tax_glom(ps, i.lvl)
  }
  taxa_names(ps.inst) <- phy_shorten_tax_names(ps.inst) %>% make.unique()

  data.ls[["PS"]][[i.lvl]] <- list(
    "Raw"  = ps.inst, 
    "Rare" = rarefy_even_depth(ps.inst, rngseed = 347), 
    "CSS"  = phy_css_norm(ps.inst),
    "CLR"  = phy_clr_norm(ps.inst)
  )

  data.ls[[i.lvl]] <- list()

  for (nm in c("Raw","Rare","CSS","CLR")) {
    phy.v <- data.ls[["PS"]][[i.lvl]][[nm]]
    splits <- make_splits(phy.v)
    data.ls[["Ps"]][[i.lvl]][[toupper(nm)]] <- c(
      list(`_all` = phy.v),
      splits
    )
  }
}

#-------------------------------------------------------------------------------
# Aesthetics 
#-------------------------------------------------------------------------------
aes.ls <- list(col = list(), shape = list(), col_alt = list())

col.vec <- qualitative_hcl(100, palette="Dark 3") 

shape.vec <- c(16, 17, 18, 15, 8, 1:7)

aes.ls[["col"]][[gr.col]] <- col.vec[1:length(levels(meta[[gr.col]]))] %>% 
  setNames(levels(meta[[gr.col]]))

aes.ls[["col"]][[time.p]] <- col.vec[sample(length(col.vec), 
                                            length(levels(meta[[time.p]])))] %>% 
  setNames(levels(meta[[time.p]]))

aes.ls[["shape"]][[time.p]] <- shape.vec[1:length(levels(meta[[time.p]]))] %>% 
  setNames(levels(meta[[time.p]]))
# same project has similar color
aes.ls[["col"]][["Treat_ID"]] <- col.vec[1:length(levels(meta$Treat_ID))] %>% 
  setNames(levels(meta$Treat_ID))
# same project has diferent color
nTreat <- length(levels(meta$Treat_ID))
col.vec <- qualitative_hcl(max(200, nTreat), palette="Dark 3")
idx <- round(seq(1, length(col.vec), length.out = nTreat))
aes.ls[["col_alt"]][["Treat_ID"]] <- col.vec[idx] %>%
  setNames(levels(meta$Treat_ID))
# Project color
# nProj <- length(levels(meta$Project))
# col.vec <- qualitative_hcl(nProj, palette = "Dark 3")
# aes.ls[["col_alt"]][["Project"]] <- col.vec %>%
#   setNames(levels(meta$Project))
col.vec <- c(
  "#264653", "#2A9D8F", "#8AB17D", "#E9C46A",
  "#F4A261", "#E76F51", "#6D597A", "#B56576",
  "#E56B6F", "#EAAC8B", "#669BBC", "#9B2226"
)
aes.ls[["col_alt"]][["Project"]] <- setNames(col.vec, levels(meta$Project))

#-------------------------------------------------------------------------------
save(list = c("data.ls", "aes.ls", "meta"), 
     file = "rdata/data.Rdata")
save(list = c("aes.ls"), file = "rdata/aes.Rdata")

rm(list = ls())
gc()



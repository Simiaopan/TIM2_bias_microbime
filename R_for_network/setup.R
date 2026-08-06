install.packages("renv")
renv::init()
renv::status()
install.packages(c(
  "remotes",
  "BiocManager",
  "igraph",
  "ggplot2",
  "dplyr",
  "tibble",
  "readr"
))
BiocManager::install(c(
  "phyloseq",
  "microbiome",
  "WGCNA",
  "DESeq2",
  "multtest",
  "preprocessCore"
), update = FALSE, ask = FALSE)

install.packages("pkgbuild")
Sys.setenv(PATH = paste(
  "C:/RBuildTools/4.4/usr/bin",
  "C:/RBuildTools/4.4/x86_64-w64-mingw32.static.posix/bin",
  Sys.getenv("PATH"),
  sep = ";"
))
pkgbuild::check_build_tools(debug = TRUE)

install.packages("remotes")
remotes::install_github("zdk123/SpiecEasi", update = "never")
library(SpiecEasi)
remotes::install_github("GraceYoon/SPRING", update = "never")
library(SPRING)
remotes::install_github("stefpeschel/NetCoMi", dependencies = TRUE, update = "never")
library(NetCoMi)

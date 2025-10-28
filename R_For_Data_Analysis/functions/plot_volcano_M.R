#-------------------------------------------------------------------------------
# Volcano plot for MaAsLin2 and quality results
#-------------------------------------------------------------------------------
plot_volcano_M <- function(df, n_split = 6) {
  library(ggplot2)
  library(dplyr)
  library(scales)
  
  #--------------------------------
  # Remove rows with any NA and report
  #--------------------------------
  n_before <- nrow(df)
  df <- df %>% na.omit()
  n_after <- nrow(df)
  message(sprintf("Removed %d rows with NA values; %d rows remain.", 
                  n_before - n_after, n_after))
  
  #--------------------------------
  # Data preparation
  #--------------------------------
  df <- df %>%
    mutate(
      log10FDR = -log10(FDR_lmm),
      Remarkable = as.logical(Remarkable),
      Both_Remarkable = as.logical(Both_Remarkable)
    )
  
  # Split points for mean relative abundance (used for point shape)
  abund_breaks <- quantile(df$RelAbundance_mean,
                           probs = seq(0, 1, length.out = n_split + 1),
                           na.rm = TRUE)
  abund_labels <- paste0(
    signif(abund_breaks[-(n_split + 1)], 2), "–", signif(abund_breaks[-1], 2)
  )
  df$AbundanceGroup <- cut(
    df$RelAbundance_mean,
    breaks = abund_breaks,
    labels = abund_labels,
    include.lowest = TRUE
  )
  shape_vals <- c(0, 1, 2, 15, 16, 17)
  
  #--------------------------------
  # Volcano plot 1 – significance & remarkable logic
  #--------------------------------
  df <- df %>%
    mutate(
      ColorGroup1 = case_when(
        Both_Remarkable ~ "Both Remarkable",
        Remarkable & !Both_Remarkable ~ "Remarkable only",
        FDR_lmm > 0.05 ~ "FDR > 0.05",
        TRUE ~ "Significant"
      )
    )
  
  color_map1 <- c(
    "Both Remarkable" = "#FF8C00",
    "Remarkable only" = "deepskyblue",
    "Significant" = "forestgreen",
    "FDR > 0.05" = "gray70"
  )
  
  p1 <- ggplot(df, aes(x = log2FC72, y = log10FDR)) +
    geom_point(aes(color = ColorGroup1, shape = AbundanceGroup),
               size = 2, alpha = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
    scale_color_manual(values = color_map1, name = "Category") +
    scale_shape_manual(values = shape_vals, name = "Mean relative abundance") +
    theme_bw(base_size = 13) +
    labs(
      title = "Volcano plot 1: significance & remarkables",
      x = "log2FC72",
      y = "-log10(FDR_lmm)"
    ) +
    theme(legend.position = "right")
  
  print(p1)
  
  #--------------------------------
  # Volcano plot 2 – color only by model_quality
  #--------------------------------
  df <- df %>%
    mutate(
      ColorGroup2 = case_when(
        model_quality == "Suspect" ~ "Suspect",
        model_quality == "Reliable" ~ "Reliable",
        TRUE ~ "Other"
      )
    )
  
  color_map2 <- c(
    "Suspect" = "red",
    "Reliable" = "gray70",
    "Other" = "forestgreen"
  )
  
  p2 <- ggplot(df, aes(x = log2FC72, y = log10FDR)) +
    geom_point(aes(color = ColorGroup2, shape = AbundanceGroup),
               size = 2, alpha = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
    scale_color_manual(values = color_map2, name = "Model quality") +
    scale_shape_manual(values = shape_vals, name = "Mean relative abundance") +
    theme_bw(base_size = 13) +
    labs(
      title = "Volcano plot 2: model quality only",
      x = "log2FC72",
      y = "-log10(FDR_lmm)"
    ) +
    theme(legend.position = "right")
  
  print(p2)
}

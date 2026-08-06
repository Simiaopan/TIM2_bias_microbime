analyze_candidate_node_metrics <- function(node_table_list,
                                           metrics = c(
                                             "degree",
                                             "betweenness",
                                             "closeness",
                                             "eigenvector"
                                           ),
                                           group_col = "candidate",
                                           method = c("wilcox", "welch"),
                                           p_adjust_method = "BH") {
  
  method <- match.arg(method)
  
  res <- lapply(names(node_table_list), function(subset) {
    
    df <- node_table_list[[subset]]
    
    # Check group column
    if (!group_col %in% colnames(df)) {
      stop("group_col not found in subset: ", subset)
    }
    
    # Check metric columns
    missing_metrics <- setdiff(metrics, colnames(df))
    if (length(missing_metrics) > 0) {
      stop(
        "These metrics are missing in subset ",
        subset,
        ": ",
        paste(missing_metrics, collapse = ", ")
      )
    }
    
    out <- lapply(metrics, function(metric) {
      
      tmp <- df |>
        dplyr::select(
          dplyr::all_of(c(group_col, metric))
        ) |>
        dplyr::filter(
          !is.na(.data[[group_col]]),
          !is.na(.data[[metric]])
        )
      
      # Ensure group variable contains exactly two levels (TRUE/FALSE)
      groups <- unique(tmp[[group_col]])
      
      if (length(groups) != 2) {
        return(data.frame(
          metric = metric,
          method = method,
          n_candidate = sum(tmp[[group_col]] == TRUE, na.rm = TRUE),
          n_background = sum(tmp[[group_col]] == FALSE, na.rm = TRUE),
          mean_candidate = NA_real_,
          mean_background = NA_real_,
          median_candidate = NA_real_,
          median_background = NA_real_,
          statistic = NA_real_,
          p_value = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      
      x_candidate <- tmp[[metric]][tmp[[group_col]] == TRUE]
      x_background <- tmp[[metric]][tmp[[group_col]] == FALSE]
      
      if (length(x_candidate) < 2 || length(x_background) < 2) {
        return(data.frame(
          metric = metric,
          method = method,
          n_candidate = length(x_candidate),
          n_background = length(x_background),
          mean_candidate = mean(x_candidate, na.rm = TRUE),
          mean_background = mean(x_background, na.rm = TRUE),
          median_candidate = median(x_candidate, na.rm = TRUE),
          median_background = median(x_background, na.rm = TRUE),
          statistic = NA_real_,
          p_value = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      
      test_res <- tryCatch({
        if (method == "wilcox") {
          stats::wilcox.test(
            x_candidate,
            x_background,
            exact = FALSE
          )
        } else {
          stats::t.test(
            x_candidate,
            x_background,
            var.equal = FALSE
          )
        }
      }, error = function(e) NULL)
      
      if (is.null(test_res)) {
        statistic <- NA_real_
        p_value <- NA_real_
      } else {
        statistic <- unname(test_res$statistic)
        p_value <- test_res$p.value
      }
      
      data.frame(
        metric = metric,
        method = method,
        n_candidate = length(x_candidate),
        n_background = length(x_background),
        mean_candidate = mean(x_candidate, na.rm = TRUE),
        mean_background = mean(x_background, na.rm = TRUE),
        median_candidate = median(x_candidate, na.rm = TRUE),
        median_background = median(x_background, na.rm = TRUE),
        statistic = statistic,
        p_value = p_value,
        stringsAsFactors = FALSE
      )
    })
    
    out <- dplyr::bind_rows(out)
    
    out |>
      dplyr::mutate(
        p_adj = p.adjust(p_value, method = p_adjust_method),
        signif = dplyr::case_when(
          p_adj < 0.001 ~ "***",
          p_adj < 0.01  ~ "**",
          p_adj < 0.05  ~ "*",
          p_adj < 0.1   ~ ".",
          TRUE          ~ "ns"
        )
      )
  })
  
  names(res) <- names(node_table_list)
  
  res
}
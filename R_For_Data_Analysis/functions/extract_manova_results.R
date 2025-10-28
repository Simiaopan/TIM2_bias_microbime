#-------------------------------------------------------------------------------
# extract_manova_results(): extract all multivariate test results
#   - Extracts results from summary.Anova.mlm$multivariate.tests
#   - Supports linearHypothesis.mlm printed structures
#-------------------------------------------------------------------------------
extract_manova_results <- function(manova_res) {
  suppressPackageStartupMessages({
    library(stringr)
    library(dplyr)
    library(purrr)
  })
  
  if (!"multivariate.tests" %in% names(manova_res)) {
    stop("Input must be a summary.Anova.mlm object with $multivariate.tests.")
  }
  
  tests <- manova_res$multivariate.tests
  
  extract_term <- function(term, obj) {
    lines <- capture.output(print(obj))
    start <- grep("^Multivariate Tests:", lines)
    if (length(start) == 0) return(NULL)
    
    tbl_lines <- lines[(start + 1):length(lines)]
    tbl_lines <- tbl_lines[grepl("Pillai|Wilks|Hotelling|Roy", tbl_lines)]
    if (length(tbl_lines) == 0) return(NULL)
    
    df_term <- purrr::map_dfr(tbl_lines, function(line) {
      parts <- strsplit(line, "\\s+")[[1]]
      parts <- parts[parts != ""]
      if (length(parts) < 8) return(NULL)
      data.frame(
        Term = term,
        Test = parts[1],
        Df = as.numeric(parts[2]),
        Stat = as.numeric(parts[3]),
        F_value = as.numeric(parts[4]),
        num_Df = as.numeric(parts[5]),
        den_Df = as.numeric(parts[6]),
        p_value = suppressWarnings(as.numeric(parts[7])),
        sig = parts[8],
        stringsAsFactors = FALSE
      )
    })
    df_term
  }
  
  df <- purrr::map_dfr(names(tests), function(term) extract_term(term, tests[[term]]))
  
  if (nrow(df) == 0) stop("No valid multivariate test results parsed.")
  
  # Clean term names
  df$Term <- gsub("biplot_df\\[\\[time_var\\]\\]", "Time", df$Term)
  df$Term <- gsub("biplot_df\\[\\[project_var\\]\\]", "Project", df$Term)
  df$Term <- gsub("biplot_df\\[\\[time_var\\]\\]:biplot_df\\[\\[project_var\\]\\]", 
                  "Time:Project", df$Term)
  
  # Create Pillai-specific summary table (explicit dplyr calls)
  pillai_df <- df %>%
    dplyr::filter(Test == "Pillai") %>%
    dplyr::mutate(
      partial_eta2 = Stat / (Stat + 1),
      contribution = round(partial_eta2 / sum(partial_eta2), 3)
    ) %>%
    dplyr::select(Term, Pillai = Stat, partial_eta2, contribution, F_value, p_value, sig)
  
  return(list(
    manova_df = df,
    pillai_df = pillai_df
  ))
}

build_netcomi_networks <- function(ps_list,
                                   candidate_features = NULL,
                                   n_background = 400,
                                   measure = "spearman",
                                   normMethod = "clr",
                                   zeroMethod = "pseudo",
                                   sparsMethod = "threshold",
                                   thresh = 0.3,
                                   dataType = "counts",
                                   weighted = TRUE,
                                   cores = 1L,
                                   seed = 123,
                                   graphlet = FALSE,
                                   gcmHeat = FALSE,
                                   gcmHeatLCC = TRUE,
                                   centrLCC = TRUE,
                                   weightDeg = FALSE,
                                   verbose = 2) {
  
  if (!requireNamespace("NetCoMi", quietly = TRUE)) stop("Package 'NetCoMi' is required.")
  if (!requireNamespace("phyloseq", quietly = TRUE)) stop("Package 'phyloseq' is required.")
  if (!requireNamespace("glue", quietly = TRUE)) stop("Package 'glue' is required.")
  
  flatten_ps_list <- function(node, current_name = "") {
    items <- list()
    
    if (inherits(node, "phyloseq")) {
      items[[current_name]] <- node
      return(items)
    }
    
    if (is.list(node)) {
      nm <- names(node)
      if (is.null(nm)) nm <- paste0("Item", seq_along(node))
      
      for (i in seq_along(node)) {
        sub_name <- if (current_name == "") nm[i] else paste0(current_name, "__", nm[i])
        items <- c(items, flatten_ps_list(node[[i]], sub_name))
      }
    }
    
    items
  }
  
  get_candidates_for_task <- function(candidate_features, task_name) {
    if (is.null(candidate_features)) return(character(0))
    
    if (is.list(candidate_features)) {
      if (!task_name %in% names(candidate_features)) {
        warning(glue::glue(
          "No candidate_features found for subset [{task_name}]. Using empty candidate set."
        ))
        return(character(0))
      }
      return(as.character(candidate_features[[task_name]]))
    }
    
    as.character(candidate_features)
  }
  
  select_taxa_keep <- function(ps, candidate_vec, n_background) {
    otu <- as(phyloseq::otu_table(ps), "matrix")
    
    if (!phyloseq::taxa_are_rows(ps)) {
      otu <- t(otu)
    }
    
    sample_sums <- colSums(otu)
    otu_rel <- sweep(otu, 2, sample_sums, "/")
    otu_rel[is.na(otu_rel)] <- 0
    
    taxa_var <- apply(otu_rel, 1, var, na.rm = TRUE)
    var_asvs <- names(sort(taxa_var, decreasing = TRUE))
    var_asvs <- head(var_asvs, min(n_background, length(var_asvs)))
    
    candidate_present <- intersect(candidate_vec, phyloseq::taxa_names(ps))
    candidate_missing <- setdiff(candidate_vec, phyloseq::taxa_names(ps))
    
    keep_taxa <- union(var_asvs, candidate_present)
    ps_keep <- phyloseq::prune_taxa(keep_taxa, ps)
    
    list(
      ps = ps_keep,
      keep_taxa = keep_taxa,
      background_taxa = var_asvs,
      candidate_present = candidate_present,
      candidate_missing = candidate_missing,
      n_input_taxa = phyloseq::ntaxa(ps),
      n_background = length(var_asvs),
      n_candidate_input = length(candidate_vec),
      n_candidate_present = length(candidate_present),
      n_candidate_missing = length(candidate_missing),
      n_final_taxa = phyloseq::ntaxa(ps_keep)
    )
  }
  
  flat_tasks <- flatten_ps_list(ps_list)
  
  if (length(flat_tasks) == 0) {
    stop("No valid phyloseq objects found in ps_list.")
  }
  
  construct_required <- c(
    "data", "dataType", "measure", "zeroMethod", "normMethod",
    "filtTax", "filtTaxPar", "sparsMethod", "thresh",
    "weighted", "cores", "seed", "verbose"
  )
  
  missing_construct <- setdiff(
    construct_required,
    names(formals(NetCoMi::netConstruct))
  )
  
  if (length(missing_construct) > 0) {
    stop(glue::glue(
      "Current NetCoMi::netConstruct() lacks required arguments: {paste(missing_construct, collapse = ', ')}"
    ))
  }
  
  analyze_required <- c(
    "net", "centrLCC", "weightDeg", "normDeg", "normBetw",
    "normClose", "graphlet", "gcmHeat", "gcmHeatLCC",
    "connectivity", "verbose"
  )
  
  missing_analyze <- setdiff(
    analyze_required,
    names(formals(NetCoMi::netAnalyze))
  )
  
  if (length(missing_analyze) > 0) {
    stop(glue::glue(
      "Current NetCoMi::netAnalyze() lacks required arguments: {paste(missing_analyze, collapse = ', ')}"
    ))
  }
  
  if (is.list(candidate_features)) {
    missing_candidate_sets <- setdiff(names(flat_tasks), names(candidate_features))
    if (length(missing_candidate_sets) > 0) {
      warning(glue::glue(
        "These ps subsets have no matching candidate list: {paste(missing_candidate_sets, collapse = ', ')}"
      ))
    }
  }
  
  constructed <- list()
  analyzed <- list()
  filter_reports <- list()
  status <- list()
  
  cat(glue::glue(
    "\nNetCoMi version: {as.character(packageVersion('NetCoMi'))}\n",
    "Network Construction Pipeline Initialized ({length(flat_tasks)} subsets found)\n",
    "Settings: {measure} | dataType: {dataType} | norm: {normMethod} | zero: {zeroMethod} | ",
    "background: {n_background} | spars: {sparsMethod} ({thresh}) | ",
    "graphlet: {graphlet} | gcmHeat: {gcmHeat}\n\n"
  ))
  
  for (task_name in names(flat_tasks)) {
    cat(glue::glue("[{task_name}] Building network...\n"))
    start_time <- Sys.time()
    
    current_candidates <- get_candidates_for_task(candidate_features, task_name)
    
    filter_obj <- select_taxa_keep(
      ps = flat_tasks[[task_name]],
      candidate_vec = current_candidates,
      n_background = n_background
    )
    
    current_ps <- filter_obj$ps
    filter_reports[[task_name]] <- filter_obj[names(filter_obj) != "ps"]
    
    if (filter_obj$n_candidate_present > 0) {
      lost_after_prune <- setdiff(
        filter_obj$candidate_present,
        phyloseq::taxa_names(current_ps)
      )
      
      if (length(lost_after_prune) > 0) {
        stop(glue::glue(
          "[{task_name}] Candidate ASVs were lost after pruning: {paste(lost_after_prune, collapse = ', ')}"
        ))
      }
    }
    
    net_out <- tryCatch({
      NetCoMi::netConstruct(
        data = current_ps,
        dataType = dataType,
        measure = measure,
        zeroMethod = zeroMethod,
        normMethod = normMethod,
        filtTax = "none",
        filtTaxPar = NULL,
        sparsMethod = sparsMethod,
        thresh = thresh,
        weighted = weighted,
        cores = cores,
        seed = seed,
        verbose = verbose
      )
    }, error = function(e) {
      status[[task_name]] <<- paste("construct_failed:", e$message)
      cat(glue::glue("FAILED at netConstruct: {e$message}\n\n"))
      return(NULL)
    })
    
    if (is.null(net_out)) next
    
    net_anal <- tryCatch({
      NetCoMi::netAnalyze(
        net = net_out,
        centrLCC = centrLCC,
        weightDeg = weightDeg,
        normDeg = TRUE,
        normBetw = TRUE,
        normClose = TRUE,
        graphlet = graphlet,
        gcmHeat = gcmHeat,
        gcmHeatLCC = gcmHeatLCC,
        connectivity = TRUE,
        verbose = verbose
      )
    }, error = function(e) {
      status[[task_name]] <<- paste("analyze_failed:", e$message)
      cat(glue::glue("FAILED at netAnalyze: {e$message}\n\n"))
      return(NULL)
    })
    
    if (is.null(net_anal)) next
    
    constructed[[task_name]] <- net_out
    analyzed[[task_name]] <- net_anal
    status[[task_name]] <- "success"
    
    elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
    
    cat(glue::glue(
      "SUCCESS | Input taxa: {filter_obj$n_input_taxa} | Final taxa: {filter_obj$n_final_taxa} | ",
      "Candidate kept: {filter_obj$n_candidate_present}/{filter_obj$n_candidate_input} | ",
      "{elapsed} min\n\n"
    ))
  }
  
  n_success <- sum(unlist(status) == "success")
  n_total <- length(flat_tasks)
  
  if (n_success == n_total) {
    cat(glue::glue("Completed successfully: {n_success}/{n_total} subsets.\n"))
  } else {
    cat(glue::glue("Finished with incomplete results: {n_success}/{n_total} subsets succeeded.\n"))
    print(status)
  }
  
  list(
    constructed = constructed,
    analyzed = analyzed,
    filter_reports = filter_reports,
    status = status,
    settings = list(
      NetCoMi_version = as.character(packageVersion("NetCoMi")),
      SpiecEasi_version = if (requireNamespace("SpiecEasi", quietly = TRUE)) {
        as.character(packageVersion("SpiecEasi"))
      } else {
        NA_character_
      },
      SPRING_version = if (requireNamespace("SPRING", quietly = TRUE)) {
        as.character(packageVersion("SPRING"))
      } else {
        NA_character_
      },
      measure = measure,
      dataType = dataType,
      normMethod = normMethod,
      zeroMethod = zeroMethod,
      n_background = n_background,
      sparsMethod = sparsMethod,
      thresh = thresh,
      weighted = weighted,
      seed = seed,
      graphlet = graphlet,
      gcmHeat = gcmHeat,
      gcmHeatLCC = gcmHeatLCC,
      centrLCC = centrLCC,
      weightDeg = weightDeg
    )
  )
}
compute_SBC_single <- function(vars_and_generated, backend, cores,
                               keep_fit, thin_ranks,
                               ensure_num_ranks_divisor,
                               dquants,
                               var_attributes,
                               save_samples_vars = NULL,
                               n_samples_keep = 100,
                               progressor = NULL) {
  variables <- vars_and_generated$variables
  generated <- vars_and_generated$generated
  
  result_with_output <- SBC:::capture_all_outputs({
    res <- tryCatch({
      fit <- SBC_fit(backend, generated, cores = cores)
      c(list(fit = fit, error = NULL))
    }, error = function(e) { list(fit = NULL, error = e) })
  })
  
  res <- result_with_output$result
  res$output <- result_with_output$output
  res$messages <- result_with_output$messages
  res$warnings <- result_with_output$warnings
  
  if(is.null(res$error)) {
    error_stats <- SBC:::capture_all_outputs({
      tryCatch( {
        res$stats <- SBC::SBC_statistics_from_single_fit(
          res$fit, variables = variables, thin_ranks = thin_ranks,
          ensure_num_ranks_divisor = ensure_num_ranks_divisor,
          generated = generated, dquants = dquants,
          var_attributes = var_attributes,
          backend = backend)
        res$backend_diagnostics <- SBC::SBC_fit_to_diagnostics(
          res$fit, res$output, res$messages, res$warnings)
        NULL
      }, error = identity)
    })
    
    if(!is.null(error_stats$result)) {
      res$error <- error_stats$result
    }
    if(!is.null(error_stats$output) && length(error_stats$output) > 0) {
      res$output <- c(res$output, "\n== Output from computing statistics ==\n", error_stats$output)
    }
    if(!is.null(error_stats$messages) && length(error_stats$messages) > 0) {
      res$messages <- c(res$messages, "== Messages from computing statistics ==", error_stats$messages)
    }
    if(!is.null(error_stats$warnings) && length(error_stats$warnings) > 0) {
      res$warnings <- c(res$warnings, "== Warnings from computing statistics ==", error_stats$warnings)
    }
    
    # NEW: Extract posterior samples BEFORE removing fit object
    if(!is.null(save_samples_vars) && !is.null(res$fit)) {
      res$posterior_samples <- list()
      
      for(var in save_samples_vars) {
        tryCatch({
          # Extract draws for the variable
          draws <- as_draws_matrix(res$fit$draws(var)) %>% as.data.frame()
          
          # Use seq with length.out for even spacing across posterior samples
          if(nrow(draws) > 0) {
            # Create indices with even spacing to get desired number of samples
            idx <- seq(1, nrow(draws), length.out = min(n_samples_keep, nrow(draws))) %>%
              as.integer()
            
            # Extract subset - save as data.frame for easier handling later
            res$posterior_samples[[var]] <- draws[idx, , drop = FALSE]
          }
        }, error = function(e) {
          warning(paste("Failed to extract samples for variable", var, ":", e$message))
        })
      }
    }
    
  } else {
    res$stats <- NULL
    res$backend_diagnostics <- NULL
    res$posterior_samples <- NULL
  }
  
  # CRITICAL: Remove fit object HERE inside compute_SBC_single to save memory
  if(!keep_fit) {
    res$fit <- NULL
  }
  
  if(!is.null(progressor)) {
    progressor()
  }
  
  res
}


# Complete compute_SBC_with_sample_saving function
compute_SBC <- function(datasets, backend,
                        cores_per_fit = default_cores_per_fit(length(datasets)),
                        keep_fits = TRUE,
                        thin_ranks = SBC_backend_default_thin_ranks(backend),
                        ensure_num_ranks_divisor = 2,
                        chunk_size = default_chunk_size(length(datasets)),
                        dquants = NULL,
                        cache_mode = "none",
                        cache_location = NULL,
                        globals = list(),
                        gen_quants = NULL,
                        n_samples_keep = 100,
                        save_samples_vars = NULL
                        ) {
  stopifnot(length(datasets) > 0)
  
  if(!is.null(gen_quants)) {
    warning("gen_quants argument is deprecated, use dquants")
    if(is.null(dquants)) {
      dquants <- gen_quants
    }
  }
  
  datasets <- validate_SBC_datasets(datasets)
  if(!is.null(dquants)) {
    dquants <- validate_derived_quantities(dquants)
  }
  
  ## Handle caching
  if(cache_mode == "results") {
    if(is.null(cache_location) || !dir.exists(dirname(cache_location))) {
      stop(SBC_error("SBC_invalid_argument_error",
                     "When using cache_mode == 'results', the cache_location argument must provide a filename in an existing directory"))
    }
    cache_basename <- basename(cache_location)
    if(!endsWith(cache_basename, ".rds")) {
      cache_location <- file.path(dirname(cache_location), paste0(cache_basename, ".rds"))
    }
    
    backend_hash <- SBC_backend_hash_for_cache(backend)
    
    # Standardize datasets for hash
    if(is.null(datasets$var_attributes)) {
      datasets$var_attributes <- NULL
    }
    data_hash <- rlang::hash(datasets)
    
    # Ensure backwards compatibility of cache
    datasets_old <- datasets
    names(datasets_old)[names(datasets) == "variables"] <- "parameters"
    data_hash_old <- rlang::hash(datasets_old)
    
    if(file.exists(cache_location)) {
      results_from_cache <- readRDS(cache_location)
      # Ensure backwards compatibility of cache
      if(!("dquants" %in% names(results_from_cache)) && ("gen_quants" %in% names(results_from_cache))) {
        # This type of assignment necessary to preserve NULL values
        results_from_cache["dquants"] <- list(results_from_cache$gen_quants)
      }
      if(!is.list(results_from_cache) ||
         !all(
           c("result", "backend_hash", "data_hash", "thin_ranks", "dquants","keep_fits")
           %in% names(results_from_cache))) {
        warning("Cache file exists but is in invalid format. Will recompute.")
      } else if(results_from_cache$backend_hash != backend_hash) {
        message("Cache file exists but the backend hash differs. Will recompute.")
      } else if(results_from_cache$data_hash != data_hash && results_from_cache$data_hash != data_hash_old) {
        message("Cache file exists but the datasets hash differs. Will recompute.")
      } else {
        if(is.null(results_from_cache$ensure_num_ranks_divisor)) {
          results_from_cache$ensure_num_ranks_divisor <- 1
        }
        
        result <- tryCatch(validate_SBC_results(results_from_cache$result),
                           error = function(e) { NULL })
        
        error_dquants <- "error dquants"
        if(!is.null(results_from_cache$dquants)) {
          results_from_cache$dquants <-
            tryCatch(validate_derived_quantities(results_from_cache$dquants),
                     error = function(e) { error_dquants })
          
        }
        if(is.null(result)) {
          warning("Cache file contains invalid SBC_results object. Will recompute.")
        } else if(results_from_cache$thin_ranks != thin_ranks ||
                  !identical(results_from_cache$dquants, dquants) ||
                  results_from_cache$ensure_num_ranks_divisor != ensure_num_ranks_divisor)  {
          if(identical(results_from_cache$dquants, error_dquants)) {
            warning("dquants loaded from cache are invalid")
          }
          if(!results_from_cache$keep_fits) {
            message("Cache file exists, but was computed with different thin_ranks/dquants/ensure_num_ranks_divisor and keep_fits == FALSE. Will recompute.")
          } else {
            message(paste0("Results loaded from cache file '", cache_basename,
                           "' but it was computed with different thin_ranks/dquants/ensure_num_ranks_divisor.\n",
                           "Calling recompute_SBC_statistics."))
            return(recompute_SBC_statistics(old_results = result, datasets = datasets,
                                            thin_ranks = thin_ranks,
                                            ensure_num_ranks_divisor = ensure_num_ranks_divisor,
                                            dquants = dquants,
                                            backend = backend))
          }
        } else {
          message(paste0("Results loaded from cache file '", cache_basename, "'"))
          check_all_SBC_diagnostics(result)
          
          return(result)
        }
      }
    }
  } else if(cache_mode == "none") {
    if(!is.null(cache_location)) {
      warning("cache_location is provided, but cache_mode is set to 'none' - no caching will take place.")
    }
  } else {
    stop(SBC_error("SBC_invalid_argument_error", "Unrecognized cache mode"))
  }
  ## End of caching
  
  
  # Create combined data for computation
  vars_and_generated_list <- list()
  for(i in 1:length(datasets)) {
    vars_and_generated_list[[i]] <- list(
      variables = datasets$variables[i,],
      generated = datasets$generated[[i]]
    )
  }
  if(is.null(dquants)) {
    future.globals <- globals
  } else {
    dq_globals <- attr(dquants, "globals")
    future.globals <- bind_globals(globals, dq_globals)
  }
  
  if(requireNamespace("progressr", quietly = TRUE)) {
    progressor <- progressr::progressor(along = vars_and_generated_list)
  } else {
    progressor <- NULL
  }
  
  results_raw <- future.apply::future_lapply(
    vars_and_generated_list, compute_SBC_single,
    backend = backend, cores = cores_per_fit,
    keep_fit = keep_fits, thin_ranks = thin_ranks,
    ensure_num_ranks_divisor = ensure_num_ranks_divisor,
    dquants = dquants,
    var_attributes = datasets$var_attributes,
    future.seed = TRUE,
    future.globals = future.globals,
    future.chunk.size = chunk_size,
    progressor = progressor, 
    save_samples_vars = save_samples_vars,
    n_samples_keep = n_samples_keep
    )
  
  # Combine, check and summarise
  fits <- rep(list(NULL), length(datasets))
  outputs <- rep(list(NULL), length(datasets))
  messages <- rep(list(NULL), length(datasets))
  warnings <- rep(list(NULL), length(datasets))
  errors <- rep(list(NULL), length(datasets))
  stats_list <- list()
  backend_diagnostics_list <- list()
  n_errors <- 0
  max_errors_to_show <- 5
  for(i in 1:length(datasets)) {
    if(!is.null(results_raw[[i]]$fit)) {
      fits[[i]] <- results_raw[[i]]$fit
    }
    if(is.null(results_raw[[i]]$error)) {
      stats_list[[i]] <- results_raw[[i]]$stats
      stats_list[[i]]$sim_id <- i
      stats_list[[i]] <- dplyr::select(stats_list[[i]], sim_id, tidyselect::everything())
      backend_diagnostics_list[[i]] <- results_raw[[i]]$backend_diagnostics
      if(!is.null(results_raw[[i]]$backend_diagnostics)){
        backend_diagnostics_list[[i]]$sim_id <- i
        backend_diagnostics_list[[i]] <- dplyr::select(backend_diagnostics_list[[i]], sim_id, tidyselect::everything())
      }
    }
    else {
      if(n_errors < max_errors_to_show) {
        if(is.null(results_raw[[i]]$fit)) {
          message("Simulation ", i, " resulted in error when fitting.\n")
          message(results_raw[[i]]$error, "\n")
          if(!is.null(results_raw[[i]]$warnings)) {
            message(" --- Warnings for sim ", i, " ----")
            message(paste0(results_raw[[i]]$warnings, collapse = "\n"))
          }
          if(!is.null(results_raw[[i]]$messages)) {
            message(" --- Messages for sim ", i, " ----")
            message(paste0(results_raw[[i]]$messages, collapse = "\n"))
          }
          if(is.null(results_raw[[i]]$output)) {
            message(" --- Nothing in stdout ---")
          } else {
            message(" ---- Model output ----")
            cat(paste0(results_raw[[i]]$output, collapse = "\n"))
          }
          message("\n ---- End of output for simulation ", i, " -----")
        } else {
          message("Simulation ", i, " resulted in error when post-processing the fit.\n",
                  "Calling `recompute_SBC_statistics` after you've found and fixed the problem could ",
                  "let you move further without refitting")
          message(results_raw[[i]]$error, "\n")
        }
        
      } else if(n_errors == max_errors_to_show) {
        message("Too many simulations produced errors. Further error messages not shown.\n")
      }
      n_errors <- n_errors + 1
      errors[[i]] <- results_raw[[i]]$error
    }
    if(!is.null(results_raw[[i]]$output)) {
      outputs[[i]] <- results_raw[[i]]$output
    }
    if(!is.null(results_raw[[i]]$messages)) {
      messages[[i]] <- results_raw[[i]]$messages
    }
    if(!is.null(results_raw[[i]]$warnings)) {
      warnings[[i]] <- results_raw[[i]]$warnings
    }
  }
  
  if(n_errors == length(datasets)) {
    warning("All simulations produced error when fitting")
  } else if(n_errors > 0) {
    warning("Total of ", n_errors, " simulations produced errors.")
  }
  
  stats <- do.call(rbind, stats_list)
  backend_diagnostics <- do.call(rbind, backend_diagnostics_list)
  
  if(!is.null(stats)) {
    check_stats(stats, datasets, thin_ranks = thin_ranks,
                ensure_num_ranks_divisor = ensure_num_ranks_divisor,
                iid_draws = SBC_backend_iid_draws(backend))
  } else {
    # Return dummy stats that let the rest of the code work.
    stats <- data.frame(sim_id = integer(0), rhat = numeric(0), ess_bulk = numeric(0),
                        ess_tail = numeric(0), has_na = logical(0),
                        rank = integer(0), simulated_value = numeric(0), max_rank = integer(0),
                        attributes = character(0))
  }
  
  default_diagnostics <-  tryCatch(
    { compute_default_diagnostics(stats) },
    error = function(e) { warning("Error when computing default per-variable diagnostics. ", e); NULL })
  
  
  
  # Collect posterior samples if they were saved
  posterior_samples_list <- list()
  if(!is.null(save_samples_vars)) {
    for(i in 1:length(datasets)) {
      if(!is.null(results_raw[[i]]$posterior_samples)) {
        posterior_samples_list[[i]] <- results_raw[[i]]$posterior_samples
      }
    }
    # Only keep non-NULL entries
    posterior_samples_list <- posterior_samples_list[!sapply(posterior_samples_list, is.null)]
  }
  
  
  
  res <- SBC_results(stats = stats, fits = fits, outputs = outputs,
                     messages = messages,
                     warnings = warnings,
                     backend_diagnostics = backend_diagnostics,
                     default_diagnostics = default_diagnostics,
                     errors = errors)
  
  if(cache_mode == "results") {
    results_for_cache <- list(result = res, backend_hash = backend_hash,
                              data_hash = data_hash, thin_ranks = thin_ranks,
                              ensure_num_ranks_divisor = ensure_num_ranks_divisor,
                              dquants = dquants, keep_fits = keep_fits)
    tryCatch(saveRDS(results_for_cache, file = cache_location),
             error = function(e) { warning("Error when saving cache file: ", e) })
  }
  
  check_all_SBC_diagnostics(res)
  
  res$saved_samples <- posterior_samples_list
  
  res
}


# Helper function to create plots from saved samples
plot_variable_density_from_saved <- function(varname, results) {
  post_df <- lapply(seq_along(results$saved_samples), function(i) {
    if(!is.null(results$saved_samples[[i]][[varname]])) {
      return(results$saved_samples[[i]][[varname]])
    } else {
      return(NULL)
    }
  }) %>%
    purrr::discard(is.null) %>%
    bind_rows()
  
  true_df <- results$stats %>% 
    filter(variable == varname) %>% 
    select(sim_id, simulated_value)
  
  ggplot() +
    geom_density(data = post_df, aes(.data[[varname]]),  # Use .data for safe evaluation
                 fill = "skyblue", alpha = 0.4, colour = "black") +
    geom_density(data = true_df, aes(x = simulated_value),
                 fill = "red", alpha = 0.3, colour = "red") +
    labs(x = varname, y = "Density") +
    theme_minimal()
}

# Function to create all plots at once
create_all_density_plots <- function(variables, results) {
  if(is.null(results$saved_samples)) {
    stop("No saved samples found in results. Set save_samples_vars when running compute_SBC_with_sample_saving.")
  }
  
  plots <- lapply(variables, function(var) {
    tryCatch({
      plot_variable_density_from_saved(var, results)
    }, error = function(e) {
      warning(paste("Failed to create plot for variable", var, ":", e$message))
      return(NULL)
    })
  })
  
  plots <- purrr::discard(plots, is.null)
  
  if(length(plots) > 0) {
    return(wrap_plots(plots))
  } else {
    stop("No plots could be created. Check variable names and saved samples.")
  }
}




check_stats <- function(stats, datasets, thin_ranks,
                        ensure_num_ranks_divisor, iid_draws) {
  unique_max_ranks <- unique(stats$max_rank)
  if(length(unique_max_ranks) != 1) {
    warning("Differening max_rank across fits")
  }
  
  if(min(unique_max_ranks) < 49) {
    if(iid_draws) {
      message_end = " (the backend produces i.i.d. samples so thin_ranks = 1 is the most sensible)."
    } else {
      message_end = "."
    }
    warning("Ranks were computed from fewer than 50 samples, the SBC checks will have low ",
            "precision.\nYou may need to increase the number of samples from the backend and make sure that ",
            "the combination of thinning in the backend, `thin_ranks` and `ensure_num_ranks_divisor` is sensible.\n",
            "Currently thin_ranks = ", thin_ranks, ", ensure_num_ranks_divisor = ",
            ensure_num_ranks_divisor,
            message_end)
    
  }
  
  all_vars <- dplyr::summarise(
    dplyr::group_by(stats, sim_id),
    all_vars = paste0(variable, collapse = ","), .groups = "drop")
  if(length(unique(all_vars$all_vars)) != 1) {
    warning("Not all fits share the same variables")
  }
  
  missing_vars <- setdiff(posterior::variables(datasets$variables), stats$variable)
  if(length(missing_vars) > 0) {
    warning("Some variables missing in fits: ", paste0(missing_vars, collapse = ", "))
    
  }
}



compute_default_diagnostics <- function(stats) {
  if(is.null(stats$attributes)) {
    stats$attributes <- ""
  }
  eligible_for_check <- function(value, attributes) {
    value[!is.na(value)
          | !attribute_present_stats(possibly_constant_var_attribute(), attributes)
    ]
  }
  
  should_check_na <- function(attributes) {
    !attribute_present_stats(na_valid_var_attribute(), attributes)
  }
  
  val <- dplyr::summarise(dplyr::group_by(stats, sim_id),
                          n_vars = dplyr::n(),
                          n_has_na = sum(has_na & should_check_na(attributes)),
                          n_na_rhat = sum(is.na(eligible_for_check(rhat, attributes))),
                          n_na_ess_bulk = sum(is.na(eligible_for_check(ess_bulk, attributes))),
                          n_na_ess_tail = sum(is.na(eligible_for_check(ess_tail, attributes))),
                          max_rhat = max(c(-Inf, eligible_for_check(rhat, attributes)), na.rm = TRUE),
                          min_ess_bulk = min(c(Inf, eligible_for_check(ess_bulk, attributes))),
                          min_ess_tail = min(c(Inf, eligible_for_check(ess_tail, attributes))),
                          min_ess_to_rank = min(c(Inf, eligible_for_check(ess_tail / max_rank, attributes))),
                          .groups = "drop"
  )
  
  return(val)
}

validate_SBC_results <- function(x) {
  stopifnot(is.list(x))
  stopifnot(inherits(x, "SBC_results"))
  if(!is.data.frame(x$stats)) {
    stop("SBC_results object has to have a 'stats' field of type data.frame")
  }
  
  # Ensure backwards compatibility of results (important to keep caches
  # valid and intact)
  # From 0.3
  if("dataset_id" %in% names(x$stats)) {
    x$stats <- dplyr::rename(x$stats, sim_id = dataset_id)
  }
  
  if("parameter" %in% names(x$stats)) {
    x$stats <- dplyr::rename(x$stats, variable = parameter)
  }
  
  # From 0.4
  if(!("has_na" %in% names(x$stats))) {
    x$stats$has_na <- FALSE
  }
  
  if(!("attributes" %in% names(x$stats))) {
    x$stats$attributes <- NA_character_
  }
  
  if("mad" %in% names(x$stats)) {
    x$stats <- dplyr::select(x$stats, -mad)
  }
  
  
  # Check validity
  if(!is.list(x$fits)) {
    stop("SBC_results object has to have a 'fits' field of type list")
  }
  
  if(!is.null(x$backend_diagnostics) && !is.data.frame(x$backend_diagnostics)) {
    stop("If the SBC_results object has a 'backend_diagnostics' field, it has to inherit from data.frame")
  }
  
  if(!is.data.frame(x$default_diagnostics)) {
    stop("The SBC_results needs a 'default_diagnostics' field, and it has to inherit from data.frame")
  }
  
  # Ensure backwards compatibility
  # From 0.3
  if("parameter" %in% names(x$default_diagnostics)) {
    x$stats <- dplyr::rename(x$stats, variable = parameter)
  }
  
  # From 0.4
  new_default_diagnostics <- c("n_has_na", "n_na_rhat", "n_na_ess_bulk", "n_na_ess_tail")
  for(nd in new_default_diagnostics) {
    if(!(nd %in% names(x$default_diagnostics))) {
      x$default_diagnostics[[nd]] <- 0
    }
  }
  
  if(!is.list(x$errors)) {
    stop("SBC_results object has to have an 'errors' field of type list")
  }
  
  if(nrow(x$stats) > 0) {
    if(!is.numeric(x$stats$sim_id)) {
      stop("The sim_id column of stats needs to be a number.")
    }
    
    
    if(min(x$stats$sim_id) < 1 || max(x$stats$sim_id) > length(x$fits)) {
      stop("stats$sim_id values must be between 1 and number of fits")
    }
  }
  
  if(!is.null(x$outputs)) {
    if(!is.list(x$outputs) || length(x$outputs) != length(x$fits)) {
      stop("outputs can only be a list of the same length as fits")
    }
  }
  
  if(!is.null(x$messages)) {
    if(!is.list(x$messages) || length(x$messages) != length(x$fits)) {
      stop("messages can only be a list of the same length as fits")
    }
  }
  
  if(!is.null(x$warnings)) {
    if(!is.list(x$warnings) || length(x$warnings) != length(x$fits)) {
      stop("warnings can only be a list of the same length as fits")
    }
  }
  
  if(!is.null(x$backend_diagnostics) && nrow(x$backend_diagnostics) > 0) {
    
    # Ensure backwards compatibility
    if("dataset_id" %in% names(x$backend_diagnostics)) {
      x$backend_diagnostics <- dplyr::rename(x$backend_diagnostics, sim_id = dataset_id)
    }
    
    
    if(!is.numeric(x$backend_diagnostics$sim_id)) {
      stop("The sim_id column of 'backend_diagnostics' needs to be a number.")
    }
    
    
    if(min(x$backend_diagnostics$sim_id) < 1 || max(x$backend_diagnostics$sim_id > length(x$fits))) {
      stop("backend_diagnostics$sim_id values must be between 1 and number of fits")
    }
  }
  
  if(nrow(x$default_diagnostics) > 0) {
    # Ensure backwards compatibility
    if("dataset_id" %in% names(x$default_diagnostics)) {
      x$default_diagnostics <- dplyr::rename(x$default_diagnostics, sim_id = dataset_id)
    }
    
    if(!is.numeric(x$default_diagnostics$sim_id)) {
      stop("The sim_id column of 'default_diagnostics' needs to be a number.")
    }
    
    
    if(min(x$default_diagnostics$sim_id) < 1 || max(x$default_diagnostics$sim_id > length(x$fits))) {
      stop("default_diagnostics$sim_id values must be between 1 and number of fits")
    }
  }
  
  
  if(length(x$fits) != length(x$errors)) {
    stop("Needs equal no. of fits and errors")
  }
  
  #TODO check identical var names
  x
}



cv_prior_rate <- function(cv_lim, p) -log(1 - p) / cv_lim




# Helper function to create pairs plot comparing posterior and true values
plot_variable_pairs_from_saved <- function(variables, results, n_samples_display = 1000) {
  if(is.null(results$saved_samples)) {
    stop("No saved samples found in results.")
  }
  
  # 데이터 수집
  post_data <- list()
  for(var in variables) {
    var_samples <- lapply(seq_along(results$saved_samples), function(i) {
      if(!is.null(results$saved_samples[[i]][[var]])) {
        samples <- results$saved_samples[[i]][[var]]
        if(nrow(samples) > n_samples_display) {
          idx <- sample(nrow(samples), n_samples_display)
          samples <- samples[idx, , drop = FALSE]
        }
        return(samples[[var]])
      } else {
        return(NULL)
      }
    }) %>% purrr::discard(is.null) %>% unlist()
    post_data[[var]] <- var_samples
  }
  
  true_data <- list()
  for(var in variables) {
    true_vals <- results$stats %>% filter(variable == var) %>% pull(simulated_value)
    if(length(true_vals) > 0) {
      n_post <- length(post_data[[var]])
      n_sims <- length(true_vals)
      n_per_sim <- n_post / n_sims
      true_data[[var]] <- rep(true_vals, each = ceiling(n_per_sim))[1:n_post]
    }
  }
  
  post_df <- data.frame(post_data)
  true_df <- data.frame(true_data)
  
  # 명확하게 factor level 설정
  post_df$type <- factor("Posterior", levels = c("Posterior", "True"))
  true_df$type <- factor("True", levels = c("Posterior", "True"))
  
  combined_df <- bind_rows(post_df, true_df)
  # 다시 한번 factor 확실히 설정
  combined_df$type <- factor(as.character(combined_df$type), levels = c("Posterior", "True"))
  
  # alpha를 aes에서 제거하고 geom에서 설정
  GGally::ggpairs(
    combined_df, 
    columns = variables,
    mapping = aes(color = type),  # alpha 제거
    upper = list(continuous = GGally::wrap("cor", alpha = 0.6)),
    lower = list(continuous = GGally::wrap("points", alpha = 0.6)),
    diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.6))
  ) +
    scale_color_manual(
      name = "Type",
      values = c("Posterior" = "skyblue", "True" = "red"),
      labels = c("Posterior", "True"),
      drop = FALSE
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 8),
      axis.text = element_text(size = 6)
    )
}


# Function to create correlation comparison matrix
create_correlation_comparison <- function(variables, results) {
  if(is.null(results$saved_samples)) {
    stop("No saved samples found in results.")
  }
  
  # Collect all posterior samples
  post_data <- list()
  for(var in variables) {
    var_samples <- lapply(seq_along(results$saved_samples), function(i) {
      if(!is.null(results$saved_samples[[i]][[var]])) {
        return(results$saved_samples[[i]][[var]][[var]])
      } else {
        return(NULL)
      }
    }) %>%
      purrr::discard(is.null) %>%
      unlist()
    
    post_data[[var]] <- var_samples
  }
  
  # Collect true values
  true_data <- list()
  for(var in variables) {
    true_vals <- results$stats %>% 
      filter(variable == var) %>% 
      pull(simulated_value)
    
    if(length(true_vals) > 0) {
      n_post <- length(post_data[[var]])
      n_sims <- length(true_vals)
      n_per_sim <- n_post / n_sims
      true_data[[var]] <- rep(true_vals, each = ceiling(n_per_sim))[1:n_post]
    }
  }
  
  # Calculate correlation matrices
  post_cor <- cor(data.frame(post_data), use = "complete.obs")
  true_cor <- cor(data.frame(true_data), use = "complete.obs")
  
  # Create comparison data frame
  cor_comparison <- data.frame(
    Variable_1 = character(0),
    Variable_2 = character(0),
    Posterior_Correlation = numeric(0),
    True_Correlation = numeric(0),
    Difference = numeric(0)
  )
  
  for(i in 1:(length(variables)-1)) {
    for(j in (i+1):length(variables)) {
      var1 <- variables[i]
      var2 <- variables[j]
      post_val <- post_cor[i, j]
      true_val <- true_cor[i, j]
      
      cor_comparison <- rbind(cor_comparison, data.frame(
        Variable_1 = var1,
        Variable_2 = var2,
        Posterior_Correlation = post_val,
        True_Correlation = true_val,
        Difference = post_val - true_val
      ))
    }
  }
  
  return(cor_comparison)
}

# Function to plot correlation comparison
plot_correlation_comparison <- function(variables, results) {
  cor_comp <- create_correlation_comparison(variables, results)
  
  if(nrow(cor_comp) == 0) {
    stop("No correlation pairs found.")
  }
  
  cor_comp$Pair <- paste(cor_comp$Variable_1, "vs", cor_comp$Variable_2)
  
  # Reshape for plotting
  cor_long <- cor_comp %>%
    select(Pair, Posterior_Correlation, True_Correlation) %>%
    pivot_longer(cols = c(Posterior_Correlation, True_Correlation),
                 names_to = "Type", values_to = "Correlation") %>%
    mutate(Type = str_remove(Type, "_Correlation"))
  
  ggplot(cor_long, aes(x = Pair, y = Correlation, fill = Type)) +
    geom_col(position = "dodge", alpha = 0.7) +
    scale_fill_manual(values = c("Posterior" = "skyblue", "True" = "red")) +
    labs(
      x = "Variable Pairs",
      y = "Correlation Coefficient",
      title = "Posterior vs True Correlations"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Master function to create all correlation plots
create_all_correlation_plots <- function(variables, results) {
  if(length(variables) < 2) {
    stop("Need at least 2 variables for correlation analysis.")
  }
  
  if(is.null(results$saved_samples)) {
    stop("No saved samples found in results. Set save_samples_vars when running compute_SBC.")
  }
  
  plots <- list()
  
  # 1. Pairs plot
  tryCatch({
    plots$pairs <- plot_variable_pairs_from_saved(variables, results)
  }, error = function(e) {
    warning(paste("Failed to create pairs plot:", e$message))
  })
  
  # 2. Correlation comparison plot
  tryCatch({
    plots$correlation_comparison <- plot_correlation_comparison(variables, results)
  }, error = function(e) {
    warning(paste("Failed to create correlation comparison plot:", e$message))
  })
  
  # 3. Correlation comparison table
  tryCatch({
    plots$correlation_table <- create_correlation_comparison(variables, results)
  }, error = function(e) {
    warning(paste("Failed to create correlation table:", e$message))
  })
  
  return(plots)
}




# This file defines re-themed versions of the bayesplot::ppc_pit_ecdf_grouped
# and the SBC::plot_sim_estimated functions.

modified_ppc_pit_ecdf_grouped <- function(
    y, yrep, group, ..., K = NULL, pit = NULL, prob = 0.99,
    plot_diff = FALSE, interpolate_adj = NULL, facet_labeller = NULL, nrow = 2) {
  bayesplot:::check_ignored_arguments(..., ok_args = c(
    "K", "pit", "prob",
    "plot_diff", "interpolate_adj"
  ))
  require("dplyr")
  if (is.null(pit)) {
    pit <- bayesplot::ppc_data(y, yrep, group) %>%
      group_by(.data$y_id) %>%
      group_map(~ mean(.x$value[.x$is_y] > .x$value[!.x$is_y]) +
                  runif(1, max = mean(.x$value[.x$is_y] == .x$value[!.x$is_y]))) %>%
      unlist()
    if (is.null(K)) {
      K <- min(nrow(yrep) + 1, 1000)
    }
  } else {
    rlang::inform("'pit' specified so ignoring 'y' and 'yrep' if specified.")
    pit <- bayesplot:::validate_pit(pit)
  }
  N <- length(pit)
  gammas <- lapply(unique(group), function(g) {
    N_g <- sum(group == g)
    bayesplot:::adjust_gamma(
      N = N_g, K = ifelse(is.null(K), N_g, K),
      prob = prob, interpolate_adj = interpolate_adj
    )
  })
  names(gammas) <- unique(group)
  data <- data.frame(pit = pit, group = group) %>%
    group_by(group) %>%
    group_map(~ data.frame(
      ecdf_value = ecdf(.x$pit)(seq(0,
                                    1,
                                    length.out = ifelse(is.null(K), nrow(.x), K)
      )),
      group = .y[1], lims_upper = bayesplot:::ecdf_intervals(
        gamma = gammas[[unlist(.y[1])]],
        N = nrow(.x), K = ifelse(is.null(K), nrow(.x),
                                 K
        )
      )$upper[-1] / nrow(.x), lims_lower = bayesplot:::ecdf_intervals(
        gamma = gammas[[unlist(.y[1])]],
        N = nrow(.x), K = ifelse(is.null(K), nrow(.x),
                                 K
        )
      )$lower[-1] / nrow(.x), x = seq(0, 1, length.out = ifelse(is.null(K),
                                                                nrow(.x), K
      ))
    )) %>%
    bind_rows()
  ggplot(data) +
    aes(x = .data$x, y = .data$ecdf_value - (plot_diff ==
                                               TRUE) * .data$x, group = .data$group, color = "y") +
    geom_step(show.legend = FALSE) +
    geom_step(
      aes(y = .data$lims_upper -
            (plot_diff == TRUE) * .data$x, color = "yrep"),
      linetype = 2,
      show.legend = FALSE
    ) +
    geom_step(
      aes(y = .data$lims_lower -
            (plot_diff == TRUE) * .data$x, color = "yrep"),
      linetype = 2,
      show.legend = FALSE
    ) +
    labs(y = ifelse(plot_diff, "ECDF difference",
                    "ECDF"
    ), x = "PIT") +
    bayesplot:::yaxis_ticks(FALSE) +
    bayesplot::bayesplot_theme_get() +
    facet_wrap("group", labeller = facet_labeller, nrow = nrow) +
    bayesplot:::scale_color_ppc() +
    bayesplot:::force_axes_in_facets()
}

modified_plot_sim_estimated <- function(
    x,
    variables = NULL,
    estimate = "mean",
    uncertainty = c("q5", "q95"),
    alpha = NULL,
    parameters = NULL,
    facet_labeller = NULL,
    nrow = 2) {
  require(ggplot2)
  if (!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if (is.null(variables)) {
      variables <- parameters
    }
  }
  if ("parameter" %in% names(x)) {
    if (!("variable" %in% names(x))) {
      warning("The x parameter contains a `parameter` column, which is deprecated, use `variable` instead.")
      x$variable <- x$parameter
    }
  }
  required_columns <- c("variable", estimate, uncertainty)
  if (!all(required_columns %in% names(x))) {
    stop(
      "The data.frame needs to have the following columns: ",
      paste0("'", required_columns, "'", collapse = ", ")
    )
  }
  if (!is.null(variables)) {
    x <- dplyr::filter(x, variable %in% variables)
  }
  if (is.null(alpha)) {
    n_points <- dplyr::summarise(dplyr::group_by(x, variable),
                                 count = dplyr::n()
    )
    max_points <- max(n_points$count)
    alpha_guess <- 1 / ((max_points * 0.06) + 1)
    alpha <- max(0.05, alpha_guess)
  }
  x$estimate__ <- x[[estimate]]
  if (length(uncertainty) != 2) {
    stop("'uncertainty' has to be null or a character vector of length 2")
  }
  x$low__ <- x[[uncertainty[1]]]
  x$high__ <- x[[uncertainty[2]]]
  all_aes <- aes(
    x = simulated_value, y = estimate__, ymin = low__,
    ymax = high__
  )
  y_label <- paste0(
    "posterior ", estimate, " (", uncertainty[1],
    " - ", uncertainty[2], ")"
  )
  if (nrow(x) == 0) {
    stop("No data to plot.")
  }
  ggplot(x, all_aes) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed"
    ) +
    geom_linerange(
      alpha = alpha,
      linewidth = .6,
      color = "#AAAAAA"
    ) +
    geom_point(
      stroke = 0,
      shape = 21,
      fill = "#386cb0",
      color = "#386cb0",
      size = 1.2
    ) +
    labs(y = y_label) +
    facet_wrap(
      ~variable,
      scales = "free",
      labeller = facet_labeller,
      nrow = nrow
    )
}

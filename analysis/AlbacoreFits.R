#------------------------------------------------------------------------------
#------------------------------ Preamble --------------------------------------
#------------------------------------------------------------------------------

library(cmdstanr)
library(here)
library(purrr)
library(posterior)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggridges)
source("utils/pipeline_utils.R")

if (!exists("ntimes", inherits = FALSE) || !exists("It", inherits = FALSE) || !exists("Ct", inherits = FALSE)) {
  source("data/Albacore_Data.R")
}

if (!exists("nchains", inherits = FALSE)) nchains <- env_get_int("ALBACOREFITS_NCHAINS", 2L)
if (!exists("nwarmup", inherits = FALSE)) nwarmup <- env_get_int("ALBACOREFITS_NWARMUP", 200L)
if (!exists("nsample", inherits = FALSE)) nsample <- env_get_int("ALBACOREFITS_NSAMPLE", 400L)
if (!exists("adaptDelta", inherits = FALSE)) adaptDelta <- env_get_num("ALBACOREFITS_ADAPT_DELTA", 0.9)
if (!exists("thin", inherits = FALSE)) thin <- env_get_int("ALBACOREFITS_THIN", 1L)
if (!exists("maxTree", inherits = FALSE)) maxTree <- env_get_int("ALBACOREFITS_MAX_TREE", 10L)
if (!exists("seed", inherits = FALSE)) seed <- env_get_int("ALBACOREFITS_SEED", 123L)

albacore_selected_models <- env_get_vec("ALBACOREFITS_MODELS", character())
albacore_save_plot_data <- tolower(env_get_chr("ALBACOREFITS_SAVE_PLOT_DATA", "true")) %in% c("1","true","yes","y")
albacore_inline_plots <- tolower(env_get_chr("ALBACOREFITS_INLINE_PLOTS", "false")) %in% c("1","true","yes","y")

normalize_model_key <- function(x) {
  x <- trimws(x)
  x <- basename(x)
  x <- sub("\\.stan$", "", x, ignore.case = TRUE)
  tolower(x)
}

#------------------------------------------------------------------------------
#------------------------------ Compile models --------------------------------
#------------------------------------------------------------------------------

stan_dir <- here("models")
stan_files <- list.files(stan_dir, pattern = "\\.stan$", full.names = FALSE)
if (length(stan_files) == 0) stop("No .stan files found in models/")

models <- setNames(
  lapply(stan_files, function(f) {
    cmdstanr::cmdstan_model(
      stan_file = file.path(stan_dir, f),
      pedantic = FALSE,
      force_recompile = FALSE
    )
  }),
  tools::file_path_sans_ext(stan_files)
)

if (length(albacore_selected_models) > 0) {
  selected_keys <- unique(normalize_model_key(albacore_selected_models))
  model_keys <- normalize_model_key(names(models))
  keep <- model_keys %in% selected_keys
  models <- models[keep]
  if (length(models) == 0) {
    stop(
      "No models left after ALBACOREFITS_MODELS filter. ",
      "Selected: ", paste(albacore_selected_models, collapse = ", "),
      " | Available: ", paste(stan_files, collapse = ", ")
    )
  }
}

#------------------------------------------------------------------------------
#----------------------------- Step 1: Shared data ----------------------------
#------------------------------------------------------------------------------

## Priors
logrPrior   <- c(-1.38, 0.2600780)
logKPrior   <- c(5.042905, 0.2659315)
sigma2Prior <- c(3.785518, 0.010223)
logqPrior   <- c(log(1e-6), log(1))

## tau2Prior: third element is a flag (0 prior-only, 1 fitted)
tau2Prior_base <- c(1.708603, 0.008613854, 0)

InputData_base <- list(
  ntimes      = ntimes,
  It          = It,
  Ct          = Ct,
  logKPrior   = logKPrior,
  logrPrior   = logrPrior,
  sigma2Prior = sigma2Prior,
  tau2Prior   = tau2Prior_base,
  logqPrior   = logqPrior,
  cdiff_scale = 0.01
)

#------------------------------------------------------------------------------
#----------------------------- Step 2: Unified random init --------------------
#------------------------------------------------------------------------------

make_init_list <- function(ntimes, nchains, seed) {
  set.seed(seed)
  lapply(seq_len(nchains), function(i) {
    
    ut0 <- rnorm(ntimes + 1, mean = 0, sd = 0.05)
    Pt0 <- rlnorm(ntimes + 1, meanlog = log(0.7), sdlog = 0.4) # >0 always
    
    list(
      logK   = log(300) + rnorm(1, 0, 0.2),
      logr   = log(0.3) + rnorm(1, 0, 0.2),
      logq   = log(1e-3) + rnorm(1, 0, 0.2),
      sigma2 = abs(rnorm(1, 0.05, 0.02)),
      tau2   = abs(rnorm(1, 0.05, 0.02)),
      ut     = ut0,
      Pt     = Pt0
    )
  })
}

init_list <- make_init_list(
  ntimes  = InputData_base$ntimes,
  nchains = nchains,
  seed    = seed
)

#------------------------------------------------------------------------------
#----------------------------- Step 3: Run both scenarios ---------------------
#------------------------------------------------------------------------------

run_models_tauflag <- function(tau_flag, InputData_base, init_list, out_root = "fitcsv") {
  
  InputData <- InputData_base
  InputData$tau2Prior[3] <- tau_flag
  
  fits <- lapply(names(models), function(mname) {
    
    message("Sampling (tau2Prior[3]=", tau_flag, "): ", mname)
    
    out_dir <- file.path(out_root, paste0("tau", tau_flag), mname)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    models[[mname]]$sample(
      data = InputData,
      init = init_list,
      seed = seed,
      output_dir = out_dir,
      chains = nchains,
      parallel_chains = nchains,
      iter_warmup = nwarmup,
      iter_sampling = as.integer(nsample / nchains),
      adapt_delta = adaptDelta,
      thin = thin,
      max_treedepth = maxTree,
      refresh = 100
    )
  })
  
  names(fits) <- names(models)
  fits
}

fits_tau0 <- run_models_tauflag(tau_flag = 0, InputData_base = InputData_base, init_list = init_list, out_root = "fitcsv")
fits_tau1 <- run_models_tauflag(tau_flag = 1, InputData_base = InputData_base, init_list = init_list, out_root = "fitcsv")

# Quick divergence check
divs_tau0 <- sapply(fits_tau0, function(f) f$diagnostic_summary()$num_divergent)
divs_tau1 <- sapply(fits_tau1, function(f) f$diagnostic_summary()$num_divergent)
print(divs_tau0)
print(divs_tau1)

#------------------------------------------------------------------------------
#----------------------------- Step 4: Build plotting data --------------------
#------------------------------------------------------------------------------

pick_biomass_cols <- function(nm) {
  pats <- c(
    "^Bt\\[\\d+\\]$",
    "^B\\[\\d+\\]$",
    "^BStatus\\[\\d+\\]$",
    "^BtoverBmsy\\[\\d+\\]$",
    "^B_t\\[\\d+\\]$"
  )
  for (p in pats) {
    hit <- grep(p, nm, value = TRUE)
    if (length(hit) > 0) return(hit)
  }
  character(0)
}

extract_bt_long <- function(fits_obj, scenario_label) {
  bind_rows(lapply(names(fits_obj), function(mname) {
    
    dr <- fits_obj[[mname]]$draws()
    if (posterior::ndraws(dr) == 0) return(NULL)
    
    df_all <- posterior::as_draws_df(dr)
    nm <- names(df_all)
    
    bio_cols <- pick_biomass_cols(nm)
    if (length(bio_cols) == 0) return(NULL)
    
    df_all %>%
      select(all_of(bio_cols)) %>%
      pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
      mutate(
        scenario = scenario_label,
        model = mname,
        time  = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])")),
        var   = str_replace(param, "\\[\\d+\\]$", "")
      ) %>%
      select(scenario, model, var, time, value)
  }))
}

extract_param_long <- function(fits_obj, scenario_label) {
  bind_rows(lapply(names(fits_obj), function(mname) {
    
    cand <- c("r","K","q","logr","logK","logq","sigma2","tau2")
    dr <- fits_obj[[mname]]$draws(variables = cand)
    if (posterior::ndraws(dr) == 0) return(NULL)
    
    df <- posterior::as_draws_df(dr)
    
    out <- list()
    if ("r" %in% names(df))      out$r      <- df$r      else if ("logr" %in% names(df)) out$r <- exp(df$logr)
    if ("K" %in% names(df))      out$K      <- df$K      else if ("logK" %in% names(df)) out$K <- exp(df$logK)
    if ("q" %in% names(df))      out$q      <- df$q      else if ("logq" %in% names(df)) out$q <- exp(df$logq)
    if ("sigma2" %in% names(df)) out$sigma2 <- df$sigma2
    if ("tau2"   %in% names(df)) out$tau2   <- df$tau2
    if (length(out) == 0) return(NULL)
    
    as.data.frame(out) %>%
      pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
      mutate(scenario = scenario_label, model = mname) %>%
      select(scenario, model, param, value)
  })) %>%
    filter(is.finite(value))
}

bt_long_all <- bind_rows(
  extract_bt_long(fits_tau0, "tau2Prior[3]=0 (prior-only)"),
  extract_bt_long(fits_tau1, "tau2Prior[3]=1 (fitted)")
) %>%
  filter(time <= ntimes) %>%
  mutate(year = factor(years[time], levels = years))

param_long_all <- bind_rows(
  extract_param_long(fits_tau0, "tau2Prior[3]=0 (prior-only)"),
  extract_param_long(fits_tau1, "tau2Prior[3]=1 (fitted)")
)

if (albacore_save_plot_data) {
  out_dir <- ensure_dir(file.path("results", "albacore_fits"))
  model_tag <- paste(unique(names(models)), collapse = "-")
  model_tag <- clean_name_for_path(model_tag)
  if (nchar(model_tag) > 80) {
    model_tag <- paste0(substr(model_tag, 1, 80), "_plus", length(unique(names(models))), "models")
  }
  out_file <- file.path(
    out_dir,
    paste0(
      "albacore_fits_plot_data_",
      model_tag, "_",
      "nsample", nsample, "_",
      "ad", clean_name_for_path(as.character(adaptDelta)), "_",
      format(Sys.time(), "%Y%m%d-%H%M%S"),
      ".rds"
    )
  )
  saveRDS(
    list(
      metadata = list(
        created_at = as.character(Sys.time()),
        models = names(models),
        settings = list(
          nchains = nchains,
          nwarmup = nwarmup,
          nsample = nsample,
          adaptDelta = adaptDelta,
          thin = thin,
          maxTree = maxTree,
          seed = seed
        )
      ),
      diagnostics = list(
        divs_tau0 = divs_tau0,
        divs_tau1 = divs_tau1
      ),
      bt_long_all = bt_long_all,
      param_long_all = param_long_all
    ),
    out_file
  )
  cat("Saved AlbacoreFits plot data:", out_file, "\n")
}

if (!albacore_inline_plots) {
  cat("Skipping inline AlbacoreFits plots (set ALBACOREFITS_INLINE_PLOTS=true to show).\n")
} else {
  #------------------------------------------------------------------------------
  #----------------------------- Step 5: Plots ----------------------------------
  #------------------------------------------------------------------------------

  ## Bt boxplot: x=time, boxplots by model, facet by scenario (and by var if multiple)
  p_bt_box_time_both <- ggplot(bt_long_all, aes(x = year, y = value, fill = model)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.size = 0.35) +
    facet_grid(scenario ~ var, scales = "free_y") +
    labs(x = "Year", y = "Biomass state") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  p_bt_box_time_both

  p_params_density_both_line <- ggplot(
    param_long_all,
    aes(x = value, colour = model, group = model)
  ) +
    geom_density(adjust = 1.2, linewidth = 0.7) +
    facet_grid(param ~ scenario, scales = "free") +
    labs(x = "Parameter value", y = "Density") +
    theme_bw() +
    theme(legend.position = "bottom")

  p_params_density_both_line
}

## ============================================================
## Prior-only model comparison:
## - Fit the selected Stan models
## - Simulate prior-predictive "truth" via SimData()
## - Save posterior draws and truth draws for the quantities used
##   in the study under results/prior_only/
##
## This reproducibility repository does not rebuild manuscript figures.
## ============================================================

source("utils/simulation.R")
source("data/Albacore_Data.R")
source("utils/pipeline_utils.R")

library(here)
library(cmdstanr)
library(posterior)

env_int <- function(name, default) {
  x <- Sys.getenv(name, "")
  if (x == "") return(default)
  as.integer(x)
}

env_num <- function(name, default) {
  x <- Sys.getenv(name, "")
  if (x == "") return(default)
  as.numeric(x)
}

selected_models <- env_get_vec("PRIORONLY_MODELS", character())

normalize_model_key <- function(x) {
  x <- trimws(x)
  x <- basename(x)
  x <- sub("\\.stan$", "", x, ignore.case = TRUE)
  tolower(x)
}

## -----------------------------
## 1) Model list
## -----------------------------
stan_dir <- if (dir.exists(here("stan", "priorSBC"))) {
  here("stan", "priorSBC")
} else {
  here("models")
}

stan_files <- list.files(stan_dir, pattern = "\\.stan$", full.names = FALSE)
if (length(stan_files) == 0) stop("No .stan files found in: ", stan_dir)

SSPM <- lapply(stan_files, function(m) {
  cmdstan_model(
    stan_file = file.path(stan_dir, m),
    pedantic = FALSE,
    force_recompile = FALSE
  )
})
names(SSPM) <- tools::file_path_sans_ext(stan_files)
models <- names(SSPM)

if (length(selected_models) > 0) {
  selected_keys <- unique(normalize_model_key(selected_models))
  model_keys <- normalize_model_key(names(SSPM))
  keep <- model_keys %in% selected_keys
  SSPM <- SSPM[keep]
  models <- names(SSPM)
  if (length(models) == 0) {
    stop(
      "No models left after PRIORONLY_MODELS filter. ",
      "Selected: ", paste(selected_models, collapse = ", "),
      " | Available: ", paste(stan_files, collapse = ", ")
    )
  }
}

cat("Models to run:\n")
print(models)

## -----------------------------
## 2) Settings
## -----------------------------
nchains    <- env_int("PRIORONLY_NCHAINS", 10)
nwarmup    <- env_int("PRIORONLY_NWARMUP", 10000)
nsample    <- env_int("PRIORONLY_NSAMPLE", 60000)
adaptDelta <- env_num("PRIORONLY_ADAPT_DELTA", 0.99)
thin       <- env_int("PRIORONLY_THIN", 1)
maxTree    <- env_int("PRIORONLY_MAX_TREE", 10)
seed       <- env_int("PRIORONLY_SEED", 123)

n_draws <- env_int("PRIORONLY_N_DRAWS",2000000)
n_sims  <- env_int("PRIORONLY_N_SIMS", 60000)
save_plot_data <- tolower(env_get_chr("PRIORONLY_SAVE_PLOT_DATA", "true")) %in% c("1", "true", "yes", "y")

if (!is.na(n_draws) && !is.na(n_sims) && n_draws < n_sims) {
  warning(
    "PRIORONLY_N_DRAWS (", n_draws, ") is smaller than PRIORONLY_N_SIMS (", n_sims,
    "). Increasing n_draws to ", n_sims, " so SimData() can run."
  )
  n_draws <- n_sims
}

## -----------------------------
## 3) Priors and data
## -----------------------------
tau2Prior   <- c(10, 10, 0)
logrPrior   <- c(-1.38, 0.2600780)
logKPrior   <- c(5.042905, 0.2659315)
sigma2Prior <- c(3.785518, 0.010223)
logqPrior   <- c(log(1e-6), log(1))

make_input_data <- function() {
  list(
    ntimes      = ntimes,
    It          = It,
    Ct          = Ct,
    logKPrior   = logKPrior,
    logrPrior   = logrPrior,
    sigma2Prior = sigma2Prior,
    tau2Prior   = tau2Prior,
    logqPrior   = logqPrior
  )
}

## Generate per-chain initial values (matched to model parameterisation)
make_chain_inits <- function(model_name, nchains, ntimes, seed = 1) {
  set.seed(seed)
  is_nct <- grepl("^nct", model_name)
  
  inits <- vector("list", nchains)
  for (ch in seq_len(nchains)) {
    init_ch <- list(
      logK   = log(300) + rnorm(1, 0, 0.05),
      logr   = log(0.3) + rnorm(1, 0, 0.05),
      logq   = log(0.001) + rnorm(1, 0, 0.10),
      sigma2 = exp(log(0.05) + rnorm(1, 0, 0.02)),
      tau2   = exp(log(0.05) + rnorm(1, 0, 0.02))
    )
    
    if (is_nct) {
      init_ch$ut <- rep(0.01, ntimes + 1) + rnorm(ntimes + 1, 0, 0.05)
    } else {
      init_ch$Pt <- pmin(0.999, pmax(1e-6, rep(0.9, ntimes + 1) + rnorm(ntimes + 1, 0, 0.05)))
    }
    
    inits[[ch]] <- init_ch
  }
  inits
}

run_one_model <- function(model, model_name) {
  input_data <- make_input_data()
  inits <- make_chain_inits(model_name, nchains = nchains, ntimes = ntimes, seed = seed)
  
  out_dir <- file.path("fitcsv", model_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  model$sample(
    data = input_data,
    init = inits,
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
}

## -----------------------------
## 4) Fit all models
## -----------------------------
fits <- mapply(
  FUN = run_one_model,
  model = SSPM,
  model_name = names(SSPM),
  SIMPLIFY = FALSE
)

## -----------------------------
## 5) Simulate prior predictive "truth"
## -----------------------------
input_data <- make_input_data()

RandomDraws <- list(
  K      = as.matrix(exp(rnorm(n_draws, logKPrior[1], logKPrior[2]))),
  r      = as.matrix(exp(rnorm(n_draws, logrPrior[1], logrPrior[2]))),
  q      = as.matrix(exp(runif(n_draws, logqPrior[1], logqPrior[2]))),
  et     = matrix(rnorm(n_draws * ntimes, 0, 1), nrow = n_draws, ncol = ntimes),
  ut     = matrix(rnorm(n_draws * (ntimes + 1), 0, 1), nrow = n_draws, ncol = ntimes + 1),
  sigma2 = as.matrix(1 / rgamma(n_draws, sigma2Prior[1], sigma2Prior[2])),
  tau2   = as.matrix(1 / rgamma(n_draws, tau2Prior[1], tau2Prior[2]))
)

save_samples_vars <- c(
  "K", "r", "q", "sigma2", "tau2", "MSY",
  "Pt[23]", "BtoverBmsy[23]", "HtoverHmsy[23]",
  "loglik"
)
prior_pt_targets <- paste0("Pt[", 21:23, "]")

sspm_dataset <- SimData(
  draws = RandomDraws,
  SBCtype = "priorSBC",
  ntimes = ntimes,
  Ct = Ct,
  ngen = n_sims,
  Input = input_data,
  # Ensure selected Pt indices are available for prior-only plot-data output even when
  # the saved SBC variable set is reduced.
  save_samples_vars = unique(c(save_samples_vars, prior_pt_targets))
)

truth_df <- as.data.frame(sspm_dataset$variables)

## -----------------------------
## 6) Extract posterior draws for selected saved variables
## -----------------------------
plot_param_priority <- c("K", "r", "q", "sigma2", "tau2", prior_pt_targets)
prior_plot_vars <- unique(c(save_samples_vars, prior_pt_targets))
plot_params <- intersect(plot_param_priority, prior_plot_vars)

post_parts <- lapply(models, function(mname) {
  dr <- fits[[mname]]$draws(variables = plot_params)
  if (posterior::ndraws(dr) == 0) return(NULL)
  d <- as_draws_df(dr)
  
  do.call(rbind, lapply(plot_params, function(p) {
    if (!p %in% names(d)) return(NULL)
    data.frame(model = mname, param = p, value = d[[p]])
  }))
})
post_parts <- Filter(Negate(is.null), post_parts)
if (length(post_parts) == 0) stop("No posterior draws available for selected plot parameters.")
post_long <- do.call(rbind, post_parts)

## Truth draws (replicated across models so each facet has the same truth reference)
available_truth_params <- intersect(plot_params, names(truth_df))
if (length(available_truth_params) == 0) {
  stop("Truth data missing all selected plot parameters.")
}
sim_long_one <- do.call(rbind, lapply(available_truth_params, function(p) {
  if (!p %in% names(truth_df)) stop("Truth data missing column: ", p)
  data.frame(param = p, value = truth_df[[p]])
}))
sim_long <- do.call(rbind, lapply(models, function(mname) {
  transform(sim_long_one, model = mname)
}))

if (save_plot_data) {
  out_dir <- ensure_dir(file.path("results", "prior_only"))
  model_tag <- paste(unique(models), collapse = "-")
  model_tag <- clean_name_for_path(model_tag)
  if (nchar(model_tag) > 80) {
    model_tag <- paste0(substr(model_tag, 1, 80), "_plus", length(unique(models)), "models")
  }
  out_file <- file.path(
    out_dir,
    paste0(
      "prior_only_plot_data_",
      model_tag, "_",
      "nsims", n_sims, "_",
      "ndraws", n_draws, "_",
      format(Sys.time(), "%Y%m%d-%H%M%S"),
      ".rds"
    )
  )
  saveRDS(
    list(
      metadata = list(
        created_at = as.character(Sys.time()),
        models = models,
        settings = list(
          nchains = nchains,
          nwarmup = nwarmup,
          nsample = nsample,
          adaptDelta = adaptDelta,
          thin = thin,
          maxTree = maxTree,
          seed = seed,
          n_draws = n_draws,
          n_sims = n_sims
        )
      ),
      post_long = post_long,
      sim_long = sim_long
    ),
    out_file
  )
  cat("Saved prior-only plot data:", out_file, "\n")
}

quit(save = "no", status = 0)

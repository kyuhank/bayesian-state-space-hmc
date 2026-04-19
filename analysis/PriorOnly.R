## ============================================================
## Prior-only model comparison:
## - Fit all Stan models in the model directory
## - Simulate prior-predictive "truth" via SimData()
## - Ridge-style densities for Pt[21:23] by model
## - Posterior density: filled (model colour)
## - Truth density: red outline
## - Highlight ONLY the area where posterior and truth densities differ (shaded)
##
## Notes:
## - Difference shading uses densities computed on a common grid.
## - Each (model, param) density is normalised to max=1 so ridge heights are comparable.
## - Lines are omitted (no posterior line, no baseline line); only fills are shown.
## ============================================================

source("utils/simulation.R")
source("data/Albacore_Data.R")
source("utils/pipeline_utils.R")

library(here)
library(cmdstanr)
library(posterior)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(ggridges)

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
skip_plots <- tolower(env_get_chr("PRIORONLY_SKIP_PLOTS", "false")) %in% c("1", "true", "yes", "y")

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

## -----------------------------
## 7) Ordering for facets and ridge stacking
## -----------------------------
param_levels <- intersect(plot_param_priority, unique(as.character(post_long$param)))
post_long$param <- factor(post_long$param, levels = rev(param_levels))
sim_long$param  <- factor(sim_long$param,  levels = rev(param_levels))

ct_models   <- models[grepl("^ct",  models)]
nct_models  <- models[grepl("^nct", models)]
model_order <- c(ct_models, nct_models)

post_long$model <- factor(post_long$model, levels = model_order)
sim_long$model  <- factor(sim_long$model,  levels = model_order)

ncol_facets <- max(length(ct_models), length(nct_models))

plot_dir <- file.path("plot", "PriorOnly")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (skip_plots) {
  cat("Skipping plot generation because PRIORONLY_SKIP_PLOTS is enabled.\n")
  quit(save = "no", status = 0)
}

## -----------------------------
## 8) Build densities on a common grid (join-by-idx) + normalise (max=1)
## -----------------------------
xgrid <- seq(0, 1, length.out = 512)

dens_long <- function(df, out_col) {
  do.call(rbind, lapply(split(df, list(df$model, df$param), drop = TRUE), function(dd) {
    if (nrow(dd) == 0) return(NULL)
    d <- density(dd$value, from = 0, to = 1, n = length(xgrid), na.rm = TRUE)
    
    out <- data.frame(
      model = as.character(dd$model[1]),
      param = as.character(dd$param[1]),
      idx   = seq_along(d$x),
      x     = d$x,
      dens  = d$y
    )
    names(out)[names(out) == "dens"] <- out_col
    out
  }))
}

post_den  <- dens_long(post_long, "dens_post")
truth_den <- dens_long(sim_long,  "dens_truth")

den_df <- left_join(
  post_den,
  truth_den %>% select(model, param, idx, dens_truth),
  by = c("model", "param", "idx")
)

den_df$dens_post  <- ifelse(is.na(den_df$dens_post),  0, den_df$dens_post)
den_df$dens_truth <- ifelse(is.na(den_df$dens_truth), 0, den_df$dens_truth)

den_df <- den_df %>%
  mutate(
    param = factor(param, levels = levels(post_long$param)),
    model = factor(model, levels = model_order),
    y0    = as.numeric(param)
  ) %>%
  group_by(model, param) %>%
  mutate(
    dens_post_s  = dens_post  / max(dens_post,  na.rm = TRUE),
    dens_truth_s = dens_truth / max(dens_truth, na.rm = TRUE)
  ) %>%
  ungroup()

## -----------------------------
## 9) Convert densities to ridge coordinates and define the difference band
## -----------------------------
scale_ridge <- 0.85

den_df <- den_df %>%
  mutate(
    y_post  = y0 + scale_ridge * dens_post_s,
    y_truth = y0 + scale_ridge * dens_truth_s,
    ymin = pmin(y_post, y_truth),
    ymax = pmax(y_post, y_truth)
  )

## -----------------------------
## 10) Plot: posterior fill + truth outline + shaded difference only
## -----------------------------
p_highlight <- ggplot() +
  ## Posterior density (filled)
  geom_ribbon(
    data = den_df,
    aes(x = x, ymin = y0, ymax = y_post, fill = model, group = interaction(model, param)),
    alpha = 0.45,
    colour = NA
  ) +
  ## Difference band between posterior and truth densities (shaded)
  geom_ribbon(
    data = den_df,
    aes(x = x, ymin = ymin, ymax = ymax, group = interaction(model, param)),
    fill = "grey30",
    alpha = 0.90
  ) +
  ## Truth density (red outline)
  geom_line(
    data = den_df,
    aes(x = x, y = y_truth, group = interaction(model, param)),
    colour = "red",
    linewidth = 0.5
  ) +
  facet_wrap(~ model, nrow = 2, ncol = ncol_facets) +
  scale_y_continuous(
    breaks = seq_along(levels(post_long$param)),
    labels = levels(post_long$param),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x     = "Value",
    y     = NULL,
    fill  = "Model",
    title = "Prior-only: posterior (filled) vs truth (red), shaded = difference (normalised densities)"
  ) +
  theme_bw()

print(p_highlight)

ggsave(
  filename = file.path(plot_dir, "ridge_pt21_23_shaded_difference_norm_nolines.png"),
  plot     = p_highlight,
  width    = 12,
  height   = 7,
  dpi      = 300
)

cat("\nSaved shaded ridge plot to: ",
    file.path(plot_dir, "ridge_pt21_23_shaded_difference_norm_nolines.png"),
    "\n", sep = "")

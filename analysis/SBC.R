
#------------------------------------------------------------------------------
#------------------------------ Preamble --------------------------------------
#------------------------------------------------------------------------------

## compile stan models
if (!exists("resolve_sbc_model_path")) {
  source("utils/pipeline_utils.R")
}

# Internal implementation mode is fixed for this pipeline.
sbc_model_subdir <- "priorSBC"
stan_model_file <- resolve_sbc_model_path(sbc_type = sbc_model_subdir, sbc_model = SBCmodel)
message("Using Stan model file: ", stan_model_file)

SSPM <- cmdstanr::cmdstan_model(stan_model_file, 
                                pedantic=F,  
                                force_recompile=F)

#------------------------------------------------------------------------------
#----------------------------- Step 1: Fit to data ----------------------------
#------------------------------------------------------------------------------


## Prior
logrPrior=c(-1.38, 0.2600780)
logKPrior=c(5.042905, 0.2659315)
sigma2Prior=c(3.785518, 0.010223)
logqPrior   <- c(log(1e-6), log(1))

## Input data 
InputData=list("ntimes"=ntimes,
               "It"=It,
               "Ct"=Ct,
               "logKPrior"=logKPrior,
               "logrPrior"=logrPrior,
               "sigma2Prior"=sigma2Prior,
               "tau2Prior"=tau2Prior,
               "logqPrior"=logqPrior)


init_bounds <- list(
  logK = c(log(250), log(350)),
  logr = c(log(0.25), log(0.35)),
  logq = c(log(1e-5), log(1e-3)),
  sigma2 = c(1e-3, 1e-1),
  tau2 = c(1e-3, 1e-1),
  ut = c(-0.2, 0.2),
  Pt = c(0.7, 0.9)
)

draw_bounded <- function(n, bounds) {
  runif(n, min = bounds[1], max = bounds[2])
}

make_fixed_sbc_init_list <- function(ntimes, nchains, init_seed, init_bounds) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(init_seed)

  lapply(seq_len(nchains), function(i) {
    list(
      logK   = draw_bounded(1, init_bounds$logK),
      logr   = draw_bounded(1, init_bounds$logr),
      logq   = draw_bounded(1, init_bounds$logq),
      sigma2 = draw_bounded(1, init_bounds$sigma2),
      tau2   = draw_bounded(1, init_bounds$tau2),
      ut     = draw_bounded(ntimes + 1, init_bounds$ut),
      Pt     = draw_bounded(ntimes + 1, init_bounds$Pt)
    )
  })
}

fixed_init_seed <- seed + 999
fixed_init_list <- make_fixed_sbc_init_list(
  ntimes = ntimes,
  nchains = nchains,
  init_seed = fixed_init_seed,
  init_bounds = init_bounds
)


#------------------------------------------------------------------------------
#----------------------------- Step 2: Perform priorSBC -----------------------
#------------------------------------------------------------------------------

n_draws <- 1000000

RandomDraws <- list(
  K   = as.matrix(exp(rnorm(n_draws, logKPrior[1], logKPrior[2]))),
  r   = as.matrix(exp(rnorm(n_draws, logrPrior[1], logrPrior[2]))),
  q   = as.matrix(exp(runif(n_draws, logqPrior[1], logqPrior[2]))),
 et = matrix(rnorm(n_draws * (ntimes), 0, 1 ), nrow = n_draws, ncol = ntimes),
  ut = matrix(rnorm(n_draws * (ntimes + 1), 0, 1 ), nrow = n_draws, ncol = ntimes+1),
 sigma2 = as.matrix(1 / rgamma(n_draws, sigma2Prior[1], sigma2Prior[2])),
 tau2   = as.matrix(1 / rgamma(n_draws, tau2Prior[1], tau2Prior[2]))
)



sspm_dataset <- SimDataPriorSBC(
  draws = RandomDraws,
  ntimes = ntimes,
  Ct = Ct,
  ngen = n_sims,
  Input = InputData,
  save_samples_vars = save_samples_vars
)

# SBC backend
sspm_backend <- SBC_backend_cmdstan_sample(SSPM,
                                           iter_warmup     = nwarmup,
                                           iter_sampling   = as.integer(nsample/nchains), 
                                           adapt_delta     = adaptDelta,
                                           chains          = nchains,
                                           init_factory = \(stan_data) {
                                             fixed_init_list
                                           }
                                           )


results<-compute_SBC(sspm_dataset,
                     sspm_backend,
                     thin_ranks = thin,
                     save_samples_vars = save_samples_vars,
                     n_samples_keep = 10,
                     keep_fits = FALSE)

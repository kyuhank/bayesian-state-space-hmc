library(posterior)
library(bayesplot)
library(cmdstanr)
library(SBC)
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(patchwork)
library(tidyr)
library(GGally)


setwd(here::here())
progressr::handlers(global = TRUE)

source("utils/pipeline_utils.R")

n_sims <- as.integer(Sys.getenv("n_sims", 100))
thin <- as.integer(Sys.getenv("thin", 10))
adaptDelta <- as.numeric(Sys.getenv("adaptDelta", 0.99))
nchains <- as.integer(Sys.getenv("nchains", 10))
nwarmup <- as.integer(Sys.getenv("nwarmup", 5000))
nsample <- as.integer(Sys.getenv("nsample", 10000))
seed <- as.integer(Sys.getenv("seed", 12345))
maxTree <- as.integer(Sys.getenv("maxTree", 10))
SBCmodel <- Sys.getenv("SBCmodel", "ct_soft2.stan")

# tau2Prior<-c("no"="13.11 0.1272 0",
#              "0.1"="13.11 0.1272 1",
#              "0.3"="13.11 1.1013 1",
#              "0.5"="13.11 2.8516 1")


tau2Prior <- Sys.getenv("tau2Prior", "13.11 2.8516 1")
tau2Prior <- as.numeric(strsplit(tau2Prior, " ")[[1]])
tau2Prior_name <- derive_tau2_prior_name(tau2Prior)



cat("SBC model:", SBCmodel, "\n")
cat("Total sims:", n_sims, "\n")
cat("adaptDelta:", adaptDelta, "\n")
cat("nchains:", nchains, "\n")
cat("nwarmup:", nwarmup, "\n")
cat("nsample:", nsample, "\n")
cat("thin:", thin, "\n")
cat("maxTree:", maxTree, "\n")
cat("seed:", seed, "\n")
cat("tau2Prior:", tau2Prior, "\n")


set.seed(seed)

if (!interactive()) {
  env_cmdstan <- Sys.getenv("CMDSTAN_PATH", "")
  if (nzchar(env_cmdstan) && dir.exists(env_cmdstan)) {
    set_cmdstan_path(env_cmdstan)
  } else if (dir.exists("/usr/local/cmdstan/cmdstan-2.37.0")) {
    set_cmdstan_path("/usr/local/cmdstan/cmdstan-2.37.0")
  }
}

source('data/Albacore_Data.R')
source('utils/simulation.R')
source('utils/sbc_utils.R')
source('utils/helpers.R')

###################
## Run the model ##
###################

dir.create("fitcsv", showWarnings = FALSE)


save_samples_vars <- c("K", 
                       "r", 
                       "q", 
                       "sigma2", 
                       "tau2", 
                       "MSY",
                       "Pt[23]",
                       "BtoverBmsy[23]",
                       "HtoverHmsy[23]",
                       "loglik")


source("analysis/SBC.R")

cat("Fixed init seed:", fixed_init_seed, "\n")
cat("Fixed init chains:", length(fixed_init_list), "\n")


dir.create("results", showWarnings = FALSE)

result_stem <- make_sbc_result_stem(
  SBCmodel = SBCmodel,
  adaptDelta = adaptDelta,
  tau2Prior_name = tau2Prior_name,
  n_sims = n_sims
)

result_file <- file.path("results", paste0(result_stem, ".rds"))

saveRDS(
  list(
    run_type = "sbc",
    results = results,
    metadata = list(
      SBCmodel = SBCmodel,
      n_sims = n_sims,
      thin = thin,
      adaptDelta = adaptDelta,
      nchains = nchains,
      nwarmup = nwarmup,
      nsample = nsample,
      seed = seed,
      maxTree = maxTree,
      tau2Prior = tau2Prior,
      tau2Prior_name = tau2Prior_name,
      fixed_init_seed = fixed_init_seed,
      init_bounds = init_bounds,
      fixed_init_list = fixed_init_list
    )
  ),
  file = result_file
)

cat("Saved SBC result bundle:", result_file, "\n")

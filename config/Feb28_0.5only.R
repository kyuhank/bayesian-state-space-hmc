# Active Condor profile.
# Run with: CONDOR_CONFIG_NAME=prod make condor-submit
# This is to test the full set of models with a quick sampler configuration.
# PriorOnly / AlbacoreFits use the same profile but different Make targets.

model_pair <- c("ct_fmax.stan",
                "ct_if.stan",
                "ct_soft.stan",
                "ct_soft2.stan",
                
                "nct_fmax.stan",
                "nct_if.stan",
                "nct_soft.stan",
                "nct_soft2.stan")

condor_user_defaults <- list(
  repo = list(
    # Used in remote_dir path: <github_repo>/<folder_name>/...
    branch = "paper_prelim2",
    folder_name = "Feb28_0.5only_quick"
  ),
  condor = list(
    # Default only when calling `Rscript launch_condor.R` directly.
    # Makefile targets override this with CONDOR_MAKE_TARGET.
    make_target = "runSBC",
    cpus = 12L,
    memory = "16GB",
    disk = "20GB"
  ),
  sampler = list(
    # Increase these for production Condor runs
    n_sims = "200",
    nchains = "10",
    nwarmup = "5000",
    nsample = "20000",
    thin = "10",
    maxTree = "10",
    seed = "123"
  ),
  grid = list(
    # Keep this aligned with local smoke-test first, then scale sampler only.
    sbc_models = model_pair,
    adapt_deltas = c("0.99"),
    tau2_keys = c("0.1","0.5")
  ),
  prior_only = list(
    # Models are submitted one model per Condor job
    models = model_pair,
    n_sims = "60000",
    n_draws = "200000",
    nchains = "10",
    nwarmup = "5000",
    nsample = "60000",
    adapt_delta = "0.99"
  ),
  albacore_fits = list(
    # Models are submitted one model per Condor job
    models = model_pair,
    nchains = "10",
    nwarmup = "5000",
    nsample = "20000",
    adapt_delta = "0.99"
  )
)

# Copy this file to `config/<profile>.R` (e.g., `config/prod.R`) and edit values below.
# Run with: CONDOR_CONFIG_NAME=<profile> make condor-submit
# For PriorOnly / AlbacoreFits use:
#   CONDOR_CONFIG_NAME=<profile> make condor-submit-prioronly
#   CONDOR_CONFIG_NAME=<profile> make condor-submit-albacorefits

model_pair <- c("ct_soft.stan", "nct_soft.stan")

condor_user_defaults <- list(
  repo = list(
    # Used in remote_dir path: <github_repo>/<folder_name>/...
    branch = "nofix",
    folder_name = "23_Feb_2026"
  ),
  condor = list(
    # Default only when calling `Rscript launch_condor.R` directly.
    # Makefile targets override this with CONDOR_MAKE_TARGET.
    make_target = "runSBC",
    # Public GHCR image, so login is usually not needed.
    ghcr_login = FALSE,
    cpus = 12L,
    memory = "20GB",
    disk = "40GB"
  ),
  sampler = list(
    # Increase these for production Condor runs
    n_sims = "2000",
    nchains = "10",
    nwarmup = "5000",
    nsample = "10000",
    thin = "10",
    maxTree = "10",
    seed = "123"
  ),
  grid = list(
    # Keep this aligned with local smoke-test first, then scale sampler only.
    sbc_models = model_pair,
    adapt_deltas = c("0.8"),
    tau2_keys = c("no", "0.1", "0.3")
  ),
  prior_only = list(
    # Models are submitted one model per Condor job
    models = model_pair,
    n_sims = "60000",
    n_draws = "200000",
    nchains = "10",
    nwarmup = "10000",
    nsample = "60000",
    adapt_delta = "0.99"
  ),
  albacore_fits = list(
    # Models are submitted one model per Condor job
    models = model_pair,
    nchains = "10",
    nwarmup = "5000",
    nsample = "10000",
    adapt_delta = "0.99"
  )
)

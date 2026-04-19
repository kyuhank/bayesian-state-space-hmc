#!/usr/bin/env bash
set -euo pipefail

main_models=(
  ct_fmax.stan
  ct_if.stan
  ct_soft.stan
  ct_soft2.stan
  nct_fmax.stan
  nct_if.stan
  nct_soft.stan
  nct_soft2.stan
)

: "${PRIORONLY_MODELS:=$(IFS=,; echo "${main_models[*]}")}"
: "${PRIORONLY_NCHAINS:=10}"
: "${PRIORONLY_NWARMUP:=10000}"
: "${PRIORONLY_NSAMPLE:=60000}"
: "${PRIORONLY_ADAPT_DELTA:=0.99}"
: "${PRIORONLY_THIN:=1}"
: "${PRIORONLY_MAX_TREE:=10}"
: "${PRIORONLY_SEED:=123}"
: "${PRIORONLY_N_DRAWS:=1500000}"
: "${PRIORONLY_N_SIMS:=50000}"
: "${PRIORONLY_SAVE_PLOT_DATA:=true}"
: "${PRIORONLY_SKIP_PLOTS:=false}"

export PRIORONLY_MODELS PRIORONLY_NCHAINS PRIORONLY_NWARMUP PRIORONLY_NSAMPLE
export PRIORONLY_ADAPT_DELTA PRIORONLY_THIN PRIORONLY_MAX_TREE PRIORONLY_SEED
export PRIORONLY_N_DRAWS PRIORONLY_N_SIMS PRIORONLY_SAVE_PLOT_DATA PRIORONLY_SKIP_PLOTS

echo "[prior-only] Running the main eight-model prior-only comparison"
echo "[prior-only] models: ${PRIORONLY_MODELS}"
echo "[prior-only] n_sims=${PRIORONLY_N_SIMS}, n_draws=${PRIORONLY_N_DRAWS}, nchains=${PRIORONLY_NCHAINS}, nwarmup=${PRIORONLY_NWARMUP}, nsample=${PRIORONLY_NSAMPLE}"

Rscript analysis/PriorOnly.R

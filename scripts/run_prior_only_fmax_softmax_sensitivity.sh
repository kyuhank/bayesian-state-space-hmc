#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

export PRIORONLY_MODELS="${PRIORONLY_MODELS:-ct_fmax.stan,nct_fmax.stan,ct_softmax.stan,nct_softmax.stan}"
export PRIORONLY_NCHAINS="${PRIORONLY_NCHAINS:-10}"
export PRIORONLY_NWARMUP="${PRIORONLY_NWARMUP:-10000}"
export PRIORONLY_NSAMPLE="${PRIORONLY_NSAMPLE:-60000}"
export PRIORONLY_ADAPT_DELTA="${PRIORONLY_ADAPT_DELTA:-0.99}"
export PRIORONLY_MAX_TREE="${PRIORONLY_MAX_TREE:-10}"
export PRIORONLY_N_DRAWS="${PRIORONLY_N_DRAWS:-2000000}"
export PRIORONLY_N_SIMS="${PRIORONLY_N_SIMS:-60000}"
export PRIORONLY_SAVE_PLOT_DATA="${PRIORONLY_SAVE_PLOT_DATA:-true}"
export PRIORONLY_SKIP_PLOTS="${PRIORONLY_SKIP_PLOTS:-false}"

echo "[prior-only] running fmax vs softmax sensitivity comparison"
echo "[prior-only] models: ${PRIORONLY_MODELS}"
echo "[prior-only] n_sims=${PRIORONLY_N_SIMS}, n_draws=${PRIORONLY_N_DRAWS}, nchains=${PRIORONLY_NCHAINS}, nwarmup=${PRIORONLY_NWARMUP}, nsample=${PRIORONLY_NSAMPLE}"

Rscript analysis/PriorOnly.R

echo "[prior-only] outputs written under results/prior_only/ and plot/PriorOnly/"

#!/usr/bin/env bash
set -euo pipefail

# Optional shared env file (same style used by launch_condor.R).
if [[ -z "${CONDOR_ENV_FILE:-}" ]]; then
  for candidate in condor_submit.env condor_env.txt config/condor_submit.env; do
    if [[ -f "${candidate}" ]]; then
      CONDOR_ENV_FILE="${candidate}"
      break
    fi
  done
fi
if [[ -n "${CONDOR_ENV_FILE:-}" && -f "${CONDOR_ENV_FILE}" ]]; then
  # shellcheck disable=SC1090
  set -a
  source "${CONDOR_ENV_FILE}"
  set +a
  echo "[smoke] Loaded env file: ${CONDOR_ENV_FILE}"
fi

# Small local defaults for fast pipeline validation.
: "${SMOKE_MODELS:=${CONDOR_MODELS:-ct_soft.stan,nct_soft.stan}}"
: "${SMOKE_MODEL:=ct_soft.stan}"
: "${SMOKE_RUN_SBC:=true}"
: "${SMOKE_SBC_CV_KEYS:=${CONDOR_TAU_KEYS:-no,0.1,0.3,0.5}}"
: "${SMOKE_SBC_ADAPT_DELTAS:=${CONDOR_ADAPT_DELTAS:-0.8}}"

IFS=',' read -r -a smoke_models_array <<< "${SMOKE_MODELS}"
if [[ ${#smoke_models_array[@]} -eq 0 || -z "${smoke_models_array[0]}" ]]; then
  smoke_models_array=("ct_soft.stan" "nct_soft.stan")
fi

# Keep compatibility with SMOKE_MODEL while preferring SMOKE_MODELS.
if [[ -n "${SMOKE_MODEL}" && "${SMOKE_MODELS}" == "ct_soft.stan,nct_soft.stan" ]]; then
  if [[ "${SMOKE_MODEL}" != "ct_soft.stan" ]]; then
    smoke_models_array=("${SMOKE_MODEL}")
  fi
fi

# Minimal SBC setup.R smoke run (produces one .RData for plots.rmd SBC sections)
: "${SBCmodel:=ct_soft.stan}"
: "${n_sims:=${CONDOR_NSIMS:-2}}"
: "${nchains:=${CONDOR_NCHAINS:-1}}"
: "${nwarmup:=${CONDOR_NWARMUP:-50}}"
: "${nsample:=${CONDOR_NSAMPLE:-100}}"
: "${adaptDelta:=${CONDOR_ADAPT_DELTAS:-0.8}}"
: "${thin:=${CONDOR_THIN:-1}}"
: "${maxTree:=${CONDOR_MAX_TREE:-8}}"
: "${seed:=${CONDOR_SEED:-123}}"
: "${tau2Prior:=13.11 0.1272 0}"

IFS=',' read -r -a smoke_sbc_cv_keys_array <<< "${SMOKE_SBC_CV_KEYS}"
if [[ ${#smoke_sbc_cv_keys_array[@]} -eq 0 || -z "${smoke_sbc_cv_keys_array[0]}" ]]; then
  smoke_sbc_cv_keys_array=("no" "0.1" "0.3" "0.5")
fi

IFS=',' read -r -a smoke_sbc_adapt_deltas_array <<< "${SMOKE_SBC_ADAPT_DELTAS}"
if [[ ${#smoke_sbc_adapt_deltas_array[@]} -eq 0 || -z "${smoke_sbc_adapt_deltas_array[0]}" ]]; then
  smoke_sbc_adapt_deltas_array=("0.8")
fi

tau2_prior_from_key() {
  case "$1" in
    no)  echo "13.11 0.1272 0" ;;
    0.1) echo "13.11 0.1272 1" ;;
    0.3) echo "13.11 1.1013 1" ;;
    0.5) echo "13.11 2.8516 1" ;;
    *)
      echo "Unsupported SMOKE_SBC_CV_KEYS value: $1" >&2
      return 1
      ;;
  esac
}

: "${PRIORONLY_MODELS:=$(IFS=,; echo "${smoke_models_array[*]}")}"
: "${PRIORONLY_NCHAINS:=1}"
: "${PRIORONLY_NWARMUP:=50}"
: "${PRIORONLY_NSAMPLE:=100}"
: "${PRIORONLY_ADAPT_DELTA:=0.8}"
: "${PRIORONLY_THIN:=1}"
: "${PRIORONLY_MAX_TREE:=8}"
: "${PRIORONLY_SEED:=123}"
: "${PRIORONLY_N_DRAWS:=2000}"
: "${PRIORONLY_N_SIMS:=50}"
: "${PRIORONLY_SAVE_PLOT_DATA:=true}"
: "${PRIORONLY_SKIP_PLOTS:=true}"

: "${ALBACOREFITS_MODELS:=$(IFS=,; echo "${smoke_models_array[*]}")}"
: "${ALBACOREFITS_NCHAINS:=1}"
: "${ALBACOREFITS_NWARMUP:=50}"
: "${ALBACOREFITS_NSAMPLE:=100}"
: "${ALBACOREFITS_ADAPT_DELTA:=0.8}"
: "${ALBACOREFITS_THIN:=1}"
: "${ALBACOREFITS_MAX_TREE:=8}"
: "${ALBACOREFITS_SEED:=123}"
: "${ALBACOREFITS_SAVE_PLOT_DATA:=true}"
: "${ALBACOREFITS_INLINE_PLOTS:=false}"

export PRIORONLY_MODELS PRIORONLY_NCHAINS PRIORONLY_NWARMUP PRIORONLY_NSAMPLE
export PRIORONLY_ADAPT_DELTA PRIORONLY_THIN PRIORONLY_MAX_TREE PRIORONLY_SEED
export PRIORONLY_N_DRAWS PRIORONLY_N_SIMS PRIORONLY_SAVE_PLOT_DATA PRIORONLY_SKIP_PLOTS

export ALBACOREFITS_MODELS ALBACOREFITS_NCHAINS ALBACOREFITS_NWARMUP ALBACOREFITS_NSAMPLE
export ALBACOREFITS_ADAPT_DELTA ALBACOREFITS_THIN ALBACOREFITS_MAX_TREE ALBACOREFITS_SEED
export ALBACOREFITS_SAVE_PLOT_DATA ALBACOREFITS_INLINE_PLOTS

export SBCmodel n_sims nchains nwarmup nsample adaptDelta thin maxTree seed tau2Prior

if [[ "${SMOKE_RUN_SBC}" == "true" ]]; then
  for sbc_model_item in "${smoke_models_array[@]}"; do
    for sbc_adapt_delta in "${smoke_sbc_adapt_deltas_array[@]}"; do
      for sbc_cv_key in "${smoke_sbc_cv_keys_array[@]}"; do
        export SBCmodel="${sbc_model_item}"
        export adaptDelta="${sbc_adapt_delta}"
        export tau2Prior="$(tau2_prior_from_key "${sbc_cv_key}")"
        echo "[smoke] Running setup.R (SBC) with model: ${SBCmodel}, adaptDelta: ${adaptDelta}, cv-key: ${sbc_cv_key}"
        Rscript setup.R
      done
    done
  done
fi

echo "[smoke] Running PriorOnly with models: ${PRIORONLY_MODELS}"
Rscript analysis/PriorOnly.R

echo "[smoke] Running AlbacoreFits with models: ${ALBACOREFITS_MODELS}"
Rscript -e "source('analysis/AlbacoreFits.R')"

echo "[smoke] Building manifests"
Rscript scripts/build_results_manifest.R
Rscript scripts/build_fitcsv_manifest.R

echo "[smoke] Rendering plots/plots.rmd"
Rscript -e "rmarkdown::render('plot/plots.rmd')"

echo "[smoke] Done"

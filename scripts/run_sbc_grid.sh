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

tau2_prior_from_key() {
  case "$1" in
    no)  echo "13.11 0.1272 0" ;;
    0.1) echo "13.11 0.1272 1" ;;
    0.3) echo "13.11 1.1013 1" ;;
    0.5) echo "13.11 2.8516 1" ;;
    *)
      echo "Unsupported SBC_CV_KEYS value: $1" >&2
      return 1
      ;;
  esac
}

: "${SBC_MODELS:=$(IFS=,; echo "${main_models[*]}")}"
: "${SBC_CV_KEYS:=no,0.1,0.3,0.5}"
: "${SBC_ADAPT_DELTAS:=0.99}"

: "${n_sims:=2000}"
: "${nchains:=10}"
: "${nwarmup:=5000}"
: "${nsample:=10000}"
: "${thin:=10}"
: "${maxTree:=10}"
: "${seed:=12345}"

IFS=',' read -r -a sbc_models_array <<< "${SBC_MODELS}"
IFS=',' read -r -a sbc_cv_keys_array <<< "${SBC_CV_KEYS}"
IFS=',' read -r -a sbc_adapt_deltas_array <<< "${SBC_ADAPT_DELTAS}"

echo "[sbc] Running the paper SBC grid"
echo "[sbc] models: ${SBC_MODELS}"
echo "[sbc] settings: ${SBC_CV_KEYS}"
echo "[sbc] adapt_delta: ${SBC_ADAPT_DELTAS}"
echo "[sbc] n_sims=${n_sims}, nchains=${nchains}, nwarmup=${nwarmup}, nsample=${nsample}, thin=${thin}, maxTree=${maxTree}"

for sbc_model_item in "${sbc_models_array[@]}"; do
  for sbc_adapt_delta in "${sbc_adapt_deltas_array[@]}"; do
    for sbc_cv_key in "${sbc_cv_keys_array[@]}"; do
      export SBCmodel="${sbc_model_item}"
      export adaptDelta="${sbc_adapt_delta}"
      export tau2Prior="$(tau2_prior_from_key "${sbc_cv_key}")"
      echo "[sbc] setup.R with model=${SBCmodel}, adaptDelta=${adaptDelta}, setting=${sbc_cv_key}"
      Rscript setup.R
    done
  done
done

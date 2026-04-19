# Stan model variants

This directory contains the eight Stan models used in the realised-prior, observed-data, and SBC comparisons, plus a small number of sensitivity-check variants.

## Naming convention

- `ct_*`: centred parameterisation
- `nct_*`: non-centred parameterisation
- `if`: hard, non-differentiable branch rule that clips catch
- `fmax`: hard, non-differentiable lower floor on the biomass update
- `soft`: smooth catch cap with a soft biomass floor
- `soft2`: soft biomass floor with a smooth harvest-rate penalty
- `softmax`: shorthand for ``soft `fmax`'', meaning a smooth softplus replacement for the hard `fmax` floor while retaining observed catch directly

## Model list

- `ct_if.stan`
- `nct_if.stan`
- `ct_fmax.stan`
- `nct_fmax.stan`
- `ct_soft.stan`
- `nct_soft.stan`
- `ct_soft2.stan`
- `nct_soft2.stan`

## Sensitivity-check variants

These variants are not part of the main eight-model comparison. They were added to isolate how the hard `fmax` floor behaves under centred and non-centred constructions.

- `ct_softmax.stan`: centred version with observed catch retained and the hard `fmax` floor replaced by a smooth softplus floor
- `nct_softmax.stan`: non-centred analogue of `ct_softmax.stan`

The manuscript uses these models only in a supplementary prior-only sensitivity check comparing:

- `ct_fmax.stan`
- `nct_fmax.stan`
- `ct_softmax.stan`
- `nct_softmax.stan`

## Shared conventions

- Biomass state is represented as relative biomass, `Pt = Bt / K`.
- `tau2Prior[3] = 1` includes the abundance-index likelihood; other values switch it off for prior-only runs.
- `soft` models define an effective catch and include a penalty term linking model-effective catch to input catch.
- `fmax` and `soft2` retain observed catch directly and modify the state update through a floor or penalty.

## Parameterisation

- Centred models sample latent biomass states directly.
- Non-centred models sample standard-normal process innovations and reconstruct the states deterministically.

## Reproducing the softmax sensitivity check

For the four-model prior-only sensitivity comparison used in the supplement, run:

- `bash scripts/run_prior_only_fmax_softmax_sensitivity.sh`

By default, this script uses the same prior-only settings described in the manuscript:

- `10` chains
- `10,000` warmup iterations
- `60,000` post-warmup draws in total
- up to `2,000,000` candidate truth draws
- `60,000` retained feasible simulations

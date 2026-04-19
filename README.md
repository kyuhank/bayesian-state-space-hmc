# Surplus-SBC-repro

Reproducibility repository for the study of parameterisation and feasibility treatment in a Bayesian surplus production model fitted with Hamiltonian Monte Carlo.

This repository contains:

- analysis scripts
- Stan model files
- saved results used in the paper
- figure-generation scripts
- helper utilities

This repository does **not** include the manuscript source files.

## Repository layout

- `analysis/`: scripts for SBC, prior-only runs, and observed-data fits
- `config/`: local and Condor-style run profiles
- `data/`: input data objects
- `models/`: Stan model files
- `plot/`: figure-generation notebooks
- `results/`: saved outputs used to reproduce the paper figures
- `scripts/`: helper scripts for manifests and sensitivity-figure export
- `utils/`: shared helper functions
- `figures_generated/`: output directory for regenerated figures

## Quick start

The saved results are included, so the paper figures can be rebuilt without rerunning the models.

1. Build the results manifest:
   - `make collect-results`
2. Rebuild the paper figures:
   - `make paper-figures`

Figures will be written to:

- `figures_generated/`

## Full reruns

To rerun the analyses from code:

- SBC grid:
  - `make run-sbc`
- Prior-only analysis:
  - `make run-prior-only`
- Observed-data fits:
  - `make run-observed-fits`

After rerunning analyses, rebuild the manifest and figures:

- `make collect-results`
- `make paper-figures`

## Notes

- `fitcsv/` is excluded by default because it is large and can be regenerated locally.
- The workflow was developed in the same software environment used for the paper repository.
- The original project used the Docker image `ghcr.io/pacificcommunity/bayes:v1.5` for containerised runs.

## Suggested archive workflow

For publication:

1. Push this repository to GitHub.
2. Create a tagged release.
3. Archive that release in Zenodo.


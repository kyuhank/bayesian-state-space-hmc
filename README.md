# Surplus-SBC-repro

This repository is provided to support transparency of the methods used in the study and to give other analysts a reproducible code base that can be inspected, rerun, and adapted. It accompanies the article *Parameterisation and feasibility treatment in a Bayesian stock assessment model fitted with Hamiltonian Monte Carlo*.

This repository contains:

- analysis scripts
- Stan model files
- saved results used in the study
- helper utilities

This repository does **not** include the manuscript source files.

## Repository layout

- `analysis/`: scripts for SBC and prior-only runs
- `data/`: input data objects
- `models/`: Stan model files
- `results/`: saved outputs used in the study
- `scripts/`: helper scripts for manifests and sensitivity exports
- `utils/`: shared helper functions

## Quick start

The saved results are included, and the main analyses can also be rerun from code.

1. Show the available commands:
   - `make help`
2. Rebuild the results manifest:
   - `make collect-results`

## Main make targets

- `make help`
  - Show the available commands.
- `make collect-results`
  - Rebuild `results/results_manifest.rds` from the saved SBC outputs.
- `make sbc`
  - Run the SBC analysis grid via `setup.R`.
- `make prior-only`
  - Run the prior-only analysis.
- `make sensitivity-softmax`
  - Run the four-model `fmax` versus `softmax` sensitivity analysis.
- `make sensitivity-figures`
  - Export the manuscript-style sensitivity figures from the saved four-model bundle.
- `make sensitivity-trajectory`
  - Export the sensitivity trajectory boxplot after `make sensitivity-softmax`.
- `make smoke-local`
  - Run a reduced local smoke test.
- `make rerun-all`
  - Run SBC and prior-only analyses, then rebuild the results manifest.

## Docker

For a containerised run that matches the software environment used in the original project, use:

- `ghcr.io/pacificcommunity/bayes:v1.5`

Pull the image:

```bash
docker pull ghcr.io/pacificcommunity/bayes:v1.5
```

Open an interactive shell in that environment:

```bash
make docker-shell
```

Run the main analyses in Docker:

```bash
make docker-sbc
make docker-prior-only
make docker-rerun-all
```

Equivalent direct `docker run` form:

```bash
docker run --rm -v "$(pwd):/work" -w /work ghcr.io/pacificcommunity/bayes:v1.5 make rerun-all
```

This gives a more reproducible execution environment than a local run, provided the saved results, input files, and container tag are kept fixed.

## Notes

- `fitcsv/` is excluded by default because it is large and can be regenerated locally.

## Suggested archive workflow

For publication:

1. Push this repository to GitHub.
2. Create a tagged release.
3. Archive that release in Zenodo.

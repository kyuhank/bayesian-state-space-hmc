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
- `scripts/`: helper scripts for manifests and sensitivity analyses
- `utils/`: shared helper functions

## Quick start

The saved results used in the study are already included under `results/`. In practice, most users will want to inspect those saved outputs rather than rerun everything locally.

1. Show the available commands without starting a rerun:
   - `make help`
2. Rebuild the results manifest from the saved outputs:
   - `make collect-results`
3. Run the full paper workflow locally, one step at a time:
   - `make`

The default `make` target runs the full paper workflow sequentially. This can take a very long time on a local machine.

## Main make targets

- `make`
  - Run the full paper workflow locally and sequentially.
- `make help`
  - Show the available commands.
- `make collect-results`
  - Rebuild `results/results_manifest.rds` from the saved SBC outputs.
- `make sbc`
  - Run the full SBC grid used in the study:
    eight main models under `no`, `0.1`, `0.3`, and `0.5`.
- `make prior-only`
  - Run the independent prior-only comparison for the main eight models.
- `make sensitivity-softmax`
  - Run the four-model `fmax` versus `softmax` sensitivity analysis.
- `make rerun-all`
  - Run all paper analyses, then rebuild the results manifest.

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

- The original study used `HTCondor` to run analyses in parallel. A local `make` run executes them sequentially and can take a very long time.
- `fitcsv/` is excluded by default because it is large and can be regenerated locally.
- Saved outputs used in the study are included under `results/`.

## Suggested archive workflow

For publication:

1. Push this repository to GitHub.
2. Create a tagged release.
3. Archive that release in Zenodo.

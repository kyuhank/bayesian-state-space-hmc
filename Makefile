.DEFAULT_GOAL := help

help:
	@printf "\nSurplus-SBC-repro targets\n\n"
	@printf "  make collect-results        Build results/results_manifest.rds from saved SBC outputs\n"
	@printf "  make paper-figures          Rebuild the main paper figures into figures_generated/\n"
	@printf "  make all-figures            Rebuild the full figure set into figures_generated/\n"
	@printf "  make sbc                    Run the SBC analysis grid (setup.R)\n"
	@printf "  make prior-only             Run the prior-only analysis\n"
	@printf "  make observed-fits          Run the observed-data albacore fits\n"
	@printf "  make sensitivity-softmax    Run the four-model fmax/softmax sensitivity analysis\n"
	@printf "  make sensitivity-figures    Export the manuscript-style sensitivity figures\n"
	@printf "  make sensitivity-trajectory Export the sensitivity trajectory boxplot\n"
	@printf "  make smoke-local            Run a reduced local smoke test\n"
	@printf "  make rerun-all              Run SBC, prior-only, observed fits, then rebuild figures\n\n"

collect-results:
	Rscript scripts/build_results_manifest.R

paper-figures:
	mkdir -p figures_generated
	Rscript -e "rmarkdown::render('plot/paper_plots.rmd')"

all-figures:
	mkdir -p figures_generated
	Rscript -e "rmarkdown::render('plot/plots.rmd')"

sbc: run-sbc

run-sbc:
	Rscript setup.R

prior-only: run-prior-only

run-prior-only:
	Rscript analysis/PriorOnly.R

observed-fits: run-observed-fits

run-observed-fits:
	Rscript -e "source('analysis/AlbacoreFits.R')"

sensitivity-softmax:
	bash scripts/run_prior_only_fmax_softmax_sensitivity.sh

sensitivity-figures:
	Rscript scripts/export_prior_only_sensitivity_figures.R

sensitivity-trajectory:
	Rscript scripts/export_prior_only_sensitivity_trajectory_boxplot.R

smoke-local:
	bash scripts/smoke_local.sh

rerun-all:
	$(MAKE) run-sbc
	$(MAKE) run-prior-only
	$(MAKE) run-observed-fits
	$(MAKE) collect-results
	$(MAKE) paper-figures

.PHONY: help collect-results paper-figures all-figures sbc run-sbc prior-only run-prior-only observed-fits run-observed-fits sensitivity-softmax sensitivity-figures sensitivity-trajectory smoke-local rerun-all

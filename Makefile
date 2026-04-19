.DEFAULT_GOAL := help

DOCKER_IMAGE ?= ghcr.io/pacificcommunity/bayes:v1.5
DOCKER_WORKDIR ?= /work

help:
	@printf "\nSurplus-SBC-repro targets\n\n"
	@printf "  make collect-results        Build results/results_manifest.rds from saved SBC outputs\n"
	@printf "  make sbc                    Run the SBC analysis grid (setup.R)\n"
	@printf "  make prior-only             Run the prior-only analysis\n"
	@printf "  make observed-fits          Run the observed-data albacore fits\n"
	@printf "  make sensitivity-softmax    Run the four-model fmax/softmax sensitivity analysis\n"
	@printf "  make sensitivity-figures    Export the manuscript-style sensitivity figures\n"
	@printf "  make sensitivity-trajectory Export the sensitivity trajectory boxplot\n"
	@printf "  make smoke-local            Run a reduced local smoke test\n"
	@printf "  make rerun-all              Run SBC, prior-only, observed fits, then rebuild the manifest\n\n"
	@printf "  make docker-pull            Pull the Docker image used for reproducible runs\n"
	@printf "  make docker-shell           Open an interactive shell in the Docker image\n"
	@printf "  make docker-sbc             Run the SBC analysis grid in Docker\n"
	@printf "  make docker-prior-only      Run the prior-only analysis in Docker\n"
	@printf "  make docker-observed-fits   Run the observed-data fits in Docker\n"
	@printf "  make docker-rerun-all       Run the full analysis sequence in Docker\n\n"

collect-results:
	Rscript scripts/build_results_manifest.R

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

docker-pull:
	docker pull $(DOCKER_IMAGE)

docker-shell:
	docker run --rm -it -v "$(CURDIR):$(DOCKER_WORKDIR)" -w $(DOCKER_WORKDIR) $(DOCKER_IMAGE) bash

docker-sbc:
	docker run --rm -v "$(CURDIR):$(DOCKER_WORKDIR)" -w $(DOCKER_WORKDIR) $(DOCKER_IMAGE) make sbc

docker-prior-only:
	docker run --rm -v "$(CURDIR):$(DOCKER_WORKDIR)" -w $(DOCKER_WORKDIR) $(DOCKER_IMAGE) make prior-only

docker-observed-fits:
	docker run --rm -v "$(CURDIR):$(DOCKER_WORKDIR)" -w $(DOCKER_WORKDIR) $(DOCKER_IMAGE) make observed-fits

docker-rerun-all:
	docker run --rm -v "$(CURDIR):$(DOCKER_WORKDIR)" -w $(DOCKER_WORKDIR) $(DOCKER_IMAGE) make rerun-all

.PHONY: help collect-results sbc run-sbc prior-only run-prior-only observed-fits run-observed-fits sensitivity-softmax sensitivity-figures sensitivity-trajectory smoke-local rerun-all docker-pull docker-shell docker-sbc docker-prior-only docker-observed-fits docker-rerun-all

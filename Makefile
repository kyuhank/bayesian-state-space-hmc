collect-results:
	Rscript scripts/build_results_manifest.R

paper-figures:
	mkdir -p figures_generated
	Rscript -e "rmarkdown::render('plot/paper_plots.rmd')"

all-figures:
	mkdir -p figures_generated
	Rscript -e "rmarkdown::render('plot/plots.rmd')"

run-sbc:
	Rscript setup.R

run-prior-only:
	Rscript analysis/PriorOnly.R

run-observed-fits:
	Rscript -e "source('analysis/AlbacoreFits.R')"

.PHONY: collect-results paper-figures all-figures run-sbc run-prior-only run-observed-fits

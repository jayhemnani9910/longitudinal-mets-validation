.PHONY: all data scores analyses dca subgroups manuscript clean test

# End-to-end
all: data scores analyses dca subgroups manuscript

# Phase 1: Data + scores
data:
	Rscript scripts/01_download_nhanes.R
	Rscript scripts/02_link_mortality.R
	Rscript scripts/03_apply_inclusion.R

scores:
	Rscript scripts/04_compute_scores.R

# Phase 2: Survival analyses
analyses:
	Rscript scripts/05_survival_allcause.R
	Rscript scripts/06_survival_cv.R
	Rscript scripts/07_survival_diabetes.R
	Rscript scripts/07b_calibration.R

# Phase 3: DCA + comparisons + ML
dca:
	Rscript scripts/08_dca.R
	Rscript scripts/09_pairwise_comparisons.R
	uv run python scripts/10_xgboost_sensitivity.py

# Phase 6: Subgroups
subgroups:
	Rscript scripts/11_subgroup_sex.R
	Rscript scripts/12_subgroup_race.R
	Rscript scripts/13_subgroup_age.R

# Phase 7: Manuscript
manuscript:
	cd manuscript && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex

# Run all R tests
test:
	Rscript -e 'testthat::test_dir("tests")'

clean:
	rm -rf results/cache/
	rm -f manuscript/*.aux manuscript/*.log manuscript/*.out manuscript/*.bbl manuscript/*.blg

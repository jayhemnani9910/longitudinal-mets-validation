.PHONY: all data scores analyses dca subgroups bootstrap manuscript clean test

# Use --no-init-file to skip .Rprofile / renv auto-activation so that
# the explicit .libPaths() call in each script takes effect.
RSCRIPT := Rscript --no-init-file

# End-to-end
all: data scores analyses dca subgroups manuscript

# Phase 1: Data + scores
data:
	$(RSCRIPT) scripts/01_download_nhanes.R
	$(RSCRIPT) scripts/02_link_mortality.R
	$(RSCRIPT) scripts/03_apply_inclusion.R

scores:
	$(RSCRIPT) scripts/04_compute_scores.R

# Phase 2: Survival analyses
analyses:
	$(RSCRIPT) scripts/05_survival_allcause.R
	$(RSCRIPT) scripts/06_survival_cv.R
	$(RSCRIPT) scripts/07_survival_diabetes.R
	$(RSCRIPT) scripts/07b_calibration.R

# Phase 3: DCA + comparisons + ML
dca:
	$(RSCRIPT) scripts/08_dca.R
	$(RSCRIPT) scripts/09_pairwise_comparisons.R
	uv run python scripts/10_xgboost_sensitivity.py

# Phase 6: Subgroups
subgroups:
	$(RSCRIPT) scripts/11_subgroup_sex.R
	$(RSCRIPT) scripts/12_subgroup_race.R
	$(RSCRIPT) scripts/13_subgroup_age.R

# Phase 8: Bootstrap sensitivity (500 PSU-cluster resamples, ~4h wall)
# Toggle rep count via N_REPS env var. Default in script is 10 for smoke tests.
bootstrap:
	mkdir -p results/cache
	N_REPS=500 $(RSCRIPT) scripts/14_bootstrap_cis.R 2>&1 | tee results/cache/bootstrap.log

bootstrap-smoke:
	mkdir -p results/cache
	N_REPS=10 $(RSCRIPT) scripts/14_bootstrap_cis.R 2>&1 | tee results/cache/bootstrap_smoke.log

# Phase 7: Manuscript
manuscript:
	cd manuscript && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex

# Run all R tests
test:
	$(RSCRIPT) -e 'testthat::test_dir("tests")'

clean:
	rm -rf results/cache/
	rm -f manuscript/*.aux manuscript/*.log manuscript/*.out manuscript/*.bbl manuscript/*.blg

# Contributing

Thanks for looking. This is a research project, so the most useful help is checking the score code against the papers, spotting errors in the analysis, and small fixes to the scripts and docs.

I try to reply to issues within a few days.

## Setup

You need R. The pipeline was run on R 4.3.3 with the packages pinned in `renv.lock`.

```bash
git clone https://github.com/jayhemnani9910/longitudinal-mets-validation.git
cd longitudinal-mets-validation
make test
```

Install `testthat` first, for example with `install.packages("testthat")` in R. `make test` runs the unit tests in `tests/scores/`. They cover the five score implementations and need only `testthat`. CI runs the same tests on every pull request.

The full pipeline also needs the R packages in `renv.lock` and `uv` for the XGBoost step. See Reproducing in the README for the `make` targets. The data download needs internet access, and the bootstrap step (500 resamples) is slow.

## Data

The NHANES tables and the public-use Linked Mortality File are downloaded by `scripts/01` and `scripts/02` into `data/`, which is gitignored. Never commit raw data or any file with one row per participant. Only aggregate tables and plots go in `results/`.

## Changes

- Branch from `master` and open a pull request against it. Keep one change per PR.
- If you change a file in `R/scores/`, add or update a test in `tests/scores/`, and say which paper equation or table it follows.
- If your change moves a number in `results/`, say which scripts you re-ran and how the number changed. The README, the manuscript and `docs/index.html` quote those numbers.
- Match the style of the file you are editing: two-space indents, `<-` for assignment, short comments that explain why.
- Commit messages follow the history: `fix: ...`, `feat: ...`, `docs: ...`, `test: ...`, `chore: ...`.

## Bugs and questions

Open an issue with the bug report template. Include the script or function, the command you ran, what you expected and what happened, plus your R version and OS.

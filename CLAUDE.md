# longitudinal-mets-validation — Claude Code context

> Project-scoped instructions. Read this file before doing anything. It captures rules, decisions, and gotchas so you don't re-derive them every session.

## Project at a glance

External longitudinal validation of two MetS risk scoring methods (RMRS and the B9 decision tree, both from Shin/Shim/Oh) on NHANES + NHANES Linked Mortality File. Benchmarked against ACC/AHA Pooled Cohort Equations (PCE), Framingham 2008, and FINDRISC. Target venue: JAMA Network Open or npj Digital Medicine.

**Spec** (the locked decisions): `/home/po/projects/work/research-on-simonshim/specs/2026-05-17-longitudinal-mets-validation-design.md`
**Plan** (the task breakdown): `/home/po/projects/work/research-on-simonshim/plans/2026-05-17-longitudinal-mets-validation-plan.md`
**Source papers**: `/home/po/projects/work/research-on-simonshim/pdfs/` (B8 PeerJ 2024, B9 PLoS One 2023, B12 ICCT 2022 are the relevant Shin/Shim/Oh papers)

## Who is Jay (the user)

Former student of Prof. Simon Shim (San Jose State University). No longer formally his student. The goal of this project is to (a) extend Shim's metabolic syndrome research line with rigorous longitudinal validation, (b) publish in a JAMA Net Open / npj Digital Medicine tier journal, and (c) email Shin/Shim/Oh with preliminary results in week 7-8 of the implementation timeline to offer co-authorship.

**Do NOT** assume Jay is a current Shim co-author. He is implementing this independently. He is NOT "DK Gajulamandyam" — that is a different person on several Shim papers.

## Locked decisions (do not re-litigate)

- **Outcomes (3)**: cardiovascular mortality, diabetes-related mortality (UCOD113 = 7 proxy for T2D), all-cause mortality. Capped at 10-year follow-up.
- **Dataset (1)**: NHANES 1999-2018 + NHANES Linked Mortality File 2019 release (both public-use, no IRB needed).
- **Risk scores under test (5)**: RMRS (B8), B9 decision tree, PCE (primary CVD baseline), Framingham 2008 (secondary CVD baseline), FINDRISC (T2D screening).
- **Methodology**: NHANES survey-weighted Fine-Gray competing-risks for CV/diabetes outcomes, Cox for all-cause. IPCW time-dependent AUC at 5y and 10y. TRIPOD+AI + STROBE reporting. OSF pre-registration before submission.
- **Engagement**: silent build for ~6-8 weeks, then email Shin + Shim + Oh with preliminary results and offer of co-authorship.
- **Out of scope**: UK Biobank, KoGES Korean cohort, the B12 two-path model (subsumed by RMRS/B9), SHAP, real-time decision support apps. Don't pull these in without explicit user approval.

## Repo layout

```
longitudinal-mets-validation/
├── R/scores/        # 5 risk score implementations + b9 features
├── R/utils/         # NHANES loader, survey design constructor
├── scripts/         # Numbered pipeline: 01_download → 10_xgboost
├── tests/scores/    # testthat tests (all pass as of 2026-05-17)
├── data/raw/        # NHANES + LMF (gitignored)
├── data/processed/  # cohort + scores (gitignored, regeneratable)
├── results/         # Output tables, plots, cached models
├── prereg/          # OSF pre-reg draft, literature scan, pause notes
├── renv.lock        # R package lockfile (rms pinned to 6.7-1)
├── pyproject.toml   # Python env (xgboost ML sensitivity only)
└── Makefile         # make data | scores | analyses | dca
```

## How to resume execution

```bash
cd /home/po/projects/work/longitudinal-mets-validation
make data         # NHANES + LMF download + cohort assembly  (~10-15 min, idempotent)
make scores       # compute 5 risk scores per subject         (~1 min)
make analyses     # all-cause Cox + CV/diabetes Fine-Gray     (~10 min)
make dca          # DCA + pairwise + XGBoost                  (~5 min)
```

Then update `prereg/osf-preregistration-draft.md` with actual N (already 17,031 from first run), submit to OSF, push the OSF URL to README.

## Tooling rules (don't relearn these)

1. **`make` uses `Rscript --no-init-file`** to bypass renv auto-activation. The hardcoded `.libPaths()` at the top of each script depends on this.
2. **`rms` is pinned to 6.7-1** (CRAN archive). rms 6.8+ requires R 4.4. We have R 4.3.3. Don't `renv::install("rms")` blindly — use the archive URL.
3. **Skip the `arrow` R package.** 30-45 min Apache Arrow C++ compile blocks everything. Use `pyreadstat` Python-side and `.rds` R-side for interop.
4. **NHANES `haven_labelled` columns** must be stripped before `bind_rows`. `R/utils/load_nhanes.R` handles this.
5. **HDL column name varies by cycle**: `LBDHDL` (1999-2002) / `LBXHDD` (2003-04) / `LBDHDD` (2005+). Coalesce all three. `scripts/03_apply_inclusion.R` does this.
6. **Fasting subsample restricts cohort to ~half** of adult NHANES. Plan for 17-30K, not 40-50K.
7. **timeROC `iid=TRUE` OOMs at N=17k.** Use `iid = FALSE` and compute bootstrap CIs separately. Already applied in scripts/05-07.
8. **10-year AUC at t=10 cap is degenerate** (IPCW weights blow up at exactly the cap). Use t=9.5 as the 10-year horizon. Not yet fixed as of 2026-05-18.

System dependencies (already installed via apt): `libxml2-dev libcurl4-openssl-dev libssl-dev libgsl-dev libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev libpng-dev libtiff5-dev libjpeg-dev`.

## Open questions (resolve before final analysis)

1. **Diabetes mortality only has 6 events at 10y** — too few for Fine-Gray. Options: broaden to include diabetes-contributing causes (UCOD113 flagged DIABETES=1), drop diabetes outcome entirely, or extend follow-up. **Decide before OSF submission.**
2. **10y AUC NA across the board** — use t=9.5 horizon. Update `scripts/05_survival_allcause.R`, `06_survival_cv.R`, `07_survival_diabetes.R`.
3. **Initial results pattern**: Framingham (0.792 5y all-cause AUC) >> RMRS (0.582) >> B9 (0.525). Clinical scores dominate. Reframe paper accordingly — the story shifts to "MetS scores don't replace clinical risk equations for mortality prediction but offer at-home screening value."
4. **FINDRISC proxies are crude.** Refine before final analysis by mapping NHANES PAQ + DR1TOT + MCQ300C to FINDRISC items properly.
5. **B9 left-subtree approximated.** Add the US-refit sensitivity per spec Task 5.2b.

Full open-questions list lives in user memory at `~/.claude/projects/-home-po-projects-work-research-on-simonshim/memory/project_open_questions.md`.

## Writing-style ban list (inherited from Jay's global CLAUDE.md, applied everywhere)

Any prose written for the user (commit messages, README copy, manuscript drafts, OSF text, outreach emails, plot captions) must NOT contain:

- **Em-dashes** (`—`, `---`, `--`). Use periods, commas, colons, parens, or restructure.
- **AI ban-words**: `thrilled, seamlessly, leveraged, leverage (verb), utilize, cutting-edge, state-of-the-art, world-class, best-in-class, robust, passionate, deeply, truly, absolutely, incredibly, furthermore, moreover, additionally, in addition, delighted, excited, game-changing, transformative, next-generation, innovative, quietly`.
- **LLM cadence patterns**: 3-word punchy closer sentences, "X, not Y" punchline endings, "There is..." openers, triadic parallels, "X sits in Y", "real X" as ersatz precision, "in under N minutes" tidy time quantifications.
- **Never mention SJSU / San Jose State University** under Jay's name. Jay attended PDEU, not SJSU. Shim's affiliation can be mentioned factually inside research notes but not in outward-facing prose authored by Jay.

Code comments and internal task tracking are exempt.

## Engagement strategy reminder

Build silently weeks 1-7. Email Shin / Shim / Oh in week 7-8 with preliminary results and offer of co-authorship. Email template lives in `specs/2026-05-17-longitudinal-mets-validation-design.md` section 8.3. Do NOT email earlier without explicit user approval.

## Interaction style

- Ask clarifying questions only when the next step is a genuine fork. Don't stack premature multi-choice questions.
- For ambiguous decisions, ask via `AskUserQuestion` with 2-4 concrete options.
- Break technical concepts into plain language first, then add depth.
- Mark every commit with a clear subject line under 72 chars.

## Quick "don't-need-to-ask" facts

- Author full name: **Jay Hemnani**
- GitHub username: **jayhemnani9910**
- Email on this machine: jayhemnani992000@gmail.com
- Repo is public on GitHub
- License: MIT for code, CC-BY for manuscript
- Working dir: `/home/po/projects/work/longitudinal-mets-validation`
- Master branch (not main) — `git push origin master`
- R 4.3.3 (Ubuntu 24.04 default). Don't upgrade unless explicitly needed.
- Python 3.11.14 via uv

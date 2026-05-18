# longitudinal-mets-validation — Claude Code context

> Project-scoped instructions. Read this file before doing anything. It captures stable rules, locked decisions, and gotchas so you don't re-derive them every session. Current results, numbers, and open questions live in user memory and in `results/` — see the "Current state" section below.

## Project at a glance

External longitudinal validation of two MetS risk scoring methods (RMRS and the B9 decision tree, both from Shin/Shim/Oh) on NHANES + NHANES Linked Mortality File. Benchmarked against ACC/AHA Pooled Cohort Equations (PCE), Framingham 2008, and FINDRISC. Target venue: JAMA Network Open or npj Digital Medicine.

**Spec** (the locked design decisions): `/home/po/projects/work/research-on-simonshim/specs/2026-05-17-longitudinal-mets-validation-design.md`
**Plan** (task breakdown): `/home/po/projects/work/research-on-simonshim/plans/2026-05-17-longitudinal-mets-validation-plan.md`
**Source papers**: `/home/po/projects/work/research-on-simonshim/pdfs/` (B8 PeerJ 2024, B9 PLoS One 2023, B12 ICCT 2022 are the relevant Shin/Shim/Oh papers)

## Who is Jay (the user)

Former student of Prof. Simon Shim. No longer formally his student. The goal of this project is to (a) extend Shim's metabolic syndrome research line with rigorous longitudinal validation, (b) publish in a JAMA Net Open / npj Digital Medicine tier journal, and (c) email Shin/Shim/Oh with preliminary results in week 7-8 of the implementation timeline to offer co-authorship.

**Do NOT** assume Jay is a current Shim co-author. He is implementing this independently. He is NOT "DK Gajulamandyam" — that name appears as a co-author on several Shim papers and refers to a different person.

## Locked decisions (do not re-litigate)

- **Outcomes (3)**: cardiovascular mortality, diabetes-related mortality, all-cause mortality.
  - Diabetes-related mortality uses the **broadened** definition: UCOD_LEADING == 7 OR UCOD_LEADING == 9 (nephritis) OR LMF DIABETES contributing-cause flag == 1.
  - Follow-up cap: 10 years for CV + all-cause; 15 years for diabetes-related (broader definition needs the longer window to accumulate events).
- **Dataset (1)**: NHANES 1999-2018 + NHANES Linked Mortality File 2019 release (both public-use, no IRB needed).
- **Risk scores under test (5)**: RMRS (B8), B9 decision tree, PCE (primary CVD baseline), Framingham 2008 (secondary CVD baseline), FINDRISC (T2D screening).
- **Methodology**:
  - NHANES survey-weighted Fine-Gray competing-risks for CV and diabetes outcomes; Cox for all-cause.
  - IPCW time-dependent AUC at 5y and **9.5y** for CV/all-cause; at 5y, 10y, and **14.5y** for diabetes. (The off-cap horizons avoid IPCW degeneracy at exactly the follow-up cap.)
  - Delta-AUC inference uses `pROC::roc.test` DeLong on horizon-cap binary outcome, because timeROC `iid=TRUE` OOMs at N>=17k (see tooling rule 7). Bootstrap CIs at 500 PSU-cluster resamples are the pre-registered definitive sensitivity.
  - DCA recalibrates each score to a horizon-specific predicted probability via per-score Cox, so all 5 scores enter DCA on the same probability scale.
  - TRIPOD+AI + STROBE reporting. OSF pre-registration before manuscript submission.
- **Manuscript narrative**: outcome-specific value. MetS scores do not replace clinical risk equations for CV or all-cause mortality; RMRS has a defined diabetes-mortality lane. The full reasoning lives in user memory (`project_direction.md`) and in the OSF draft at `prereg/osf-preregistration-draft.md`. Don't restate the narrative in code or commit messages without checking the current memory entry first.
- **Engagement**: silent build for ~6-8 weeks, then email Shin + Shim + Oh with preliminary results and offer of co-authorship.
- **Out of scope**: UK Biobank, KoGES Korean cohort, the B12 two-path model (subsumed by RMRS/B9), SHAP, real-time decision support apps. Don't pull these in without explicit user approval.

## Repo layout

```
longitudinal-mets-validation/
├── R/scores/        # 5 risk score implementations + b9 features
├── R/utils/         # NHANES loader, survey design constructor
├── scripts/         # Numbered pipeline: 01_download to 10_xgboost
├── tests/scores/    # testthat tests
├── data/raw/        # NHANES + LMF (gitignored)
├── data/processed/  # cohort + scores (gitignored, regeneratable)
├── results/         # Output tables, plots, cached models
├── prereg/          # OSF pre-reg draft, literature scan, session pause notes
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

Before re-running, **read user memory** to see what's been done and what the current open questions are. Don't repeat work that already landed.

## Current state (read these instead of hard-coding numbers here)

Numbers, event counts, AUC tables, and the current open-questions list move every session. The structured place for them:

- **User memory**: `~/.claude/projects/-home-po-projects-work-research-on-simonshim/memory/` — start with `MEMORY.md`, then `project_direction.md` for the latest narrative + AUC tables, `project_open_questions.md` for what's next.
- **Result tables**: `results/*_summary.csv`, `results/pairwise_comparisons.csv`, `results/incremental_dm.csv`, `results/xgboost_summary.csv`, `results/dca_*_netbenefit.csv`.
- **OSF draft**: `prereg/osf-preregistration-draft.md` (the registered analysis plan + the registered findings).
- **Commit history**: `git log --oneline` will show the most recent phase completions.

Do not re-derive cohort N, event counts, or AUCs by reading and reasoning about scripts. Read the result CSVs.

## Tooling rules (don't relearn these)

1. **`make` uses `Rscript --no-init-file`** to bypass renv auto-activation. The hardcoded `.libPaths()` at the top of each script depends on this.
2. **`rms` is pinned to 6.7-1** (CRAN archive). rms 6.8+ requires R 4.4. We have R 4.3.3. Don't `renv::install("rms")` blindly. Use the archive URL.
3. **Skip the `arrow` R package.** 30-45 min Apache Arrow C++ compile blocks everything. Use `pyreadstat` Python-side and `.rds` R-side. For R-to-Python flow the pattern is: R writes `.rds`, a small R helper writes `.csv` from the `.rds`, Python reads the `.csv`.
4. **NHANES `haven_labelled` columns** must be stripped before `bind_rows`. `R/utils/load_nhanes.R` handles this.
5. **HDL column name varies by cycle**: `LBDHDL` (1999-2002) / `LBXHDD` (2003-04) / `LBDHDD` (2005+). Coalesce all three. `scripts/03_apply_inclusion.R` does this.
6. **Fasting subsample restricts cohort to ~half** of adult NHANES. Plan for cohort N in the 17-30K range, not 40-50K.
7. **timeROC `iid=TRUE` OOMs at N>=17k** (verified 14 GB RSS OOM-kill). Use `iid = FALSE`. For delta-AUC inference fall back to `pROC::roc.test` on the horizon-cap binary outcome (already wired in `scripts/09_pairwise_comparisons.R`). Bootstrap 500 PSU-cluster resamples for the definitive CIs.
8. **AUC horizons at the follow-up cap are degenerate.** Use t=9.5 (for 10y cap) and t=14.5 (for 15y cap) — never the exact cap value.
9. **XGBoost grouping**: use `sdmvstra x sdmvpsu` (~300 clusters), not `sdmvpsu` alone (only 3 levels, breaks GroupKFold).
10. **DCA cross-score recalibration**: FINDRISC is integer 0-26; the others are probability-style. `scripts/08_dca.R` Cox-recalibrates each score to a horizon-specific predicted probability before DCA. Document this in manuscript methods.

System dependencies (already installed via apt): `libxml2-dev libcurl4-openssl-dev libssl-dev libgsl-dev libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev libpng-dev libtiff5-dev libjpeg-dev`.

Full tooling notes with context and links live in user memory: `~/.claude/projects/-home-po-projects-work-research-on-simonshim/memory/feedback_r_tooling.md`.

## Open questions

Maintained in user memory at `~/.claude/projects/-home-po-projects-work-research-on-simonshim/memory/project_open_questions.md`. Read that file at the start of any session that touches this project. Do not duplicate the list here; it moves too fast.

## Writing-style ban list (inherited from Jay's global CLAUDE.md, applied everywhere)

Any prose written for the user (commit messages, README copy, manuscript drafts, OSF text, outreach emails, plot captions) must NOT contain:

- **Em-dashes** (`—`, `---`, `--`). Use periods, commas, colons, parens, or restructure.
- **AI ban-words**: `thrilled, seamlessly, leveraged, leverage (verb), utilize, cutting-edge, state-of-the-art, world-class, best-in-class, robust, passionate, deeply, truly, absolutely, incredibly, furthermore, moreover, additionally, in addition, delighted, excited, game-changing, transformative, next-generation, innovative, quietly`.
  - Exception: "Robust" inside the proper noun "Robust Metabolic Syndrome Risk Score" (the published name of the B8 score) is allowed.
- **LLM cadence patterns**: 3-word punchy closer sentences, "X, not Y" punchline endings, "There is..." openers, triadic parallels for rhythmic effect, "X sits in Y", "real X" as ersatz precision, "in under N minutes" tidy time quantifications.
- **Never mention SJSU / San Jose State University** under Jay's name. Jay attended PDEU, not SJSU. Shim's institutional affiliation can appear in third-person factual citations where required, but Jay must not be presented as an SJSU author.

Code comments and internal task tracking are exempt.

## Engagement strategy reminder

Build silently weeks 1-7. Email Shin / Shim / Oh in week 7-8 with preliminary results and offer of co-authorship. Email template lives in `specs/2026-05-17-longitudinal-mets-validation-design.md` section 8.3. Do NOT email earlier without explicit user approval.

## Interaction style

- Ask clarifying questions only when the next step is a genuine fork. Don't stack premature multi-choice questions.
- For ambiguous decisions, ask via `AskUserQuestion` with 2-4 concrete options.
- Break technical concepts into plain language first, then add depth. Gloss every acronym the first time it appears in user-facing prose.
- Mark every commit with a clear subject line under 72 chars.

## Subagent dispatch pattern

Long-running, well-defined engineering work (script edits + rerun + commit + push) should be dispatched to an `implementer` subagent with a self-contained prompt that includes the relevant locked decisions and the writing-style ban list. Background mode is fine for runs > ~5 minutes.

## Quick "don't-need-to-ask" facts

- Author full name: **Jay Hemnani**
- GitHub username: **jayhemnani9910**
- Email on this machine: jayhemnani992000@gmail.com
- Repo is public on GitHub at `https://github.com/jayhemnani9910/longitudinal-mets-validation`
- License: MIT for code, CC-BY for manuscript
- Working dir: `/home/po/projects/work/longitudinal-mets-validation`
- Master branch (not main). `git push origin master`. The Claude auto-mode classifier may block direct pushes to master and need either inline user authorization ("push to master") or the user pushing manually.
- R 4.3.3 (Ubuntu 24.04 default). Don't upgrade unless explicitly needed.
- Python 3.11.14 via uv.

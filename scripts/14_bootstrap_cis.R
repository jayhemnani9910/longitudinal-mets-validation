#!/usr/bin/env Rscript
# scripts/14_bootstrap_cis.R
#
# Bootstrap 95% percentile CIs at PSU-cluster resamples for the registered
# primary AUC and delta-AUC quantities. This is the registered definitive
# interval estimator for the H3 non-inferiority test (RMRS vs FINDRISC on
# diabetes-related mortality at t = 14.5, registered margin -0.05 AUC).
#
# Resampling unit: NHANES PSU cluster defined as the cross of sdmvstra and
# sdmvpsu (about 301 clusters). Within each stratum, sample its clusters with
# replacement, then take all subjects in those clusters. Do not sample
# individuals.
#
# Per resample, recompute:
#   1. AUC at primary horizons for each registered score-outcome pair
#      via timeROC::timeROC(iid = FALSE). iid = TRUE OOMs at N >= 17k.
#   2. Delta-AUCs for the registered pairwise contrasts: primary is the paired
#      IPCW time-dependent delta (timeROC), secondary is pROC DeLong on the
#      horizon-cap binary outcome (matches script 09).
#   3. DCA net benefit at one threshold per outcome (5% all-cause, 2% CV,
#      1% diabetes), per score, with Cox-recalibrated horizon-specific
#      predicted risk (matches script 08).
#
# Toggle the rep count via env var N_REPS (default 10 for smoke test).
# Seed = 20260519 for reproducibility.
#
# Inputs:  data/processed/cohort_with_scores.rds
# Outputs:
#   results/bootstrap_auc_distributions.rds
#   results/bootstrap_deltaauc_distributions.rds
#   results/bootstrap_dca_netbenefit_distributions.rds
#   results/bootstrap_summary.csv
#   results/bootstrap_dca_netbenefit.csv
#   results/cache/bootstrap.log (via tee in the calling Makefile target)

.libPaths("renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(timeROC)
  library(survival)
  library(riskRegression)
  library(prodlim)
  library(pROC)
  library(survIDINRI)
  library(dplyr)
})

source("R/utils/dca_competing.R")

# Optional parallel backend; only enabled if N_WORKERS env var is set > 1.
HAS_FUTURE <- requireNamespace("future.apply", quietly = TRUE) &&
              requireNamespace("future", quietly = TRUE)

# Tunables via env vars
N_REPS    <- as.integer(Sys.getenv("N_REPS", "10"))
N_WORKERS <- as.integer(Sys.getenv("N_WORKERS", "1"))
SEED      <- 20260519

if (is.na(N_REPS) || N_REPS < 1) N_REPS <- 10
if (is.na(N_WORKERS) || N_WORKERS < 1) N_WORKERS <- 1

message(sprintf("Bootstrap: N_REPS=%d, N_WORKERS=%d, SEED=%d",
                N_REPS, N_WORKERS, SEED))

#-----------------------------------------------------------------------
# Load cohort and define resampling unit
#-----------------------------------------------------------------------

df <- readRDS("data/processed/cohort_with_scores.rds")
df$cluster_id <- paste(df$sdmvstra, df$sdmvpsu, sep = "_")
clusters_all <- sort(unique(df$cluster_id))
n_clusters <- length(clusters_all)

# Pre-split row indices by cluster so resampling is index-only.
cluster_rows <- split(seq_len(nrow(df)), df$cluster_id)

# Map each design stratum (sdmvstra) to its PSU clusters, so the bootstrap can
# resample PSUs WITHIN strata (the NHANES-appropriate scheme) rather than from
# the pooled cluster set.
uc <- unique(df[, c("sdmvstra", "cluster_id")])
strata_to_clusters <- split(uc$cluster_id, uc$sdmvstra)
n_strata <- length(strata_to_clusters)

message(sprintf("Cohort N = %d, PSU clusters = %d, design strata = %d",
                nrow(df), n_clusters, n_strata))

#-----------------------------------------------------------------------
# Registered AUC quantities (score x outcome x horizon)
#-----------------------------------------------------------------------

# Outcome metadata mirrors scripts 05/06/07.
# competing/cause_coded mirror scripts 06-08: cause-specific outcomes use the
# Fine-Gray CIF + competing-risks net benefit on the 1999-2014 cause-coded
# cohort; all-cause uses every cycle. DCA thresholds sit on each outcome's
# achievable predicted-risk range (CV mortality never reaches the old 7.5%).
outcomes <- list(
  allcause = list(time_col   = "followup_years",
                  event_col  = "event_allcause",
                  delta_col  = "event_allcause",
                  status_col = "event_allcause",
                  competing  = FALSE, cause_coded = FALSE,
                  cause      = 1,
                  horizons   = 9.5,
                  scores     = c("rmrs_score", "b9_score", "pce_score",
                                 "framingham_score", "findrisc_score"),
                  dca_thresh = 0.05),
  cv       = list(time_col   = "followup_years",
                  event_col  = "event_cv",
                  delta_col  = "competing_cv",
                  status_col = "competing_cv",
                  competing  = TRUE, cause_coded = TRUE,
                  cause      = 1,
                  horizons   = 9.5,
                  scores     = c("rmrs_score", "b9_score", "pce_score",
                                 "framingham_score"),
                  dca_thresh = 0.02),
  dm       = list(time_col   = "followup_years_dm",
                  event_col  = "event_dm",
                  delta_col  = "competing_dm",
                  status_col = "competing_dm",
                  competing  = TRUE, cause_coded = TRUE,
                  cause      = 1,
                  horizons   = 14.5,
                  scores     = c("rmrs_score", "b9_score", "findrisc_score"),
                  dca_thresh = 0.01)
)

# Registered pairwise contrasts. The H3 row is flagged below.
pairs <- list(
  list(label = "RMRS vs FINDRISC on diabetes mortality (H3)",
       score1 = "rmrs_score", score2 = "findrisc_score",
       outcome = "dm", registered_h3 = TRUE),
  list(label = "RMRS vs B9 on all-cause mortality",
       score1 = "rmrs_score", score2 = "b9_score",
       outcome = "allcause", registered_h3 = FALSE),
  list(label = "RMRS vs B9 on CV mortality",
       score1 = "rmrs_score", score2 = "b9_score",
       outcome = "cv", registered_h3 = FALSE),
  list(label = "RMRS vs B9 on diabetes mortality",
       score1 = "rmrs_score", score2 = "b9_score",
       outcome = "dm", registered_h3 = FALSE),
  list(label = "PCE vs Framingham on CV mortality",
       score1 = "pce_score", score2 = "framingham_score",
       outcome = "cv", registered_h3 = FALSE)
)

#-----------------------------------------------------------------------
# Helper functions
#-----------------------------------------------------------------------

# IPCW time-dependent AUC at one horizon, last element (matches script 09).
timeroc_point_auc <- function(time_v, delta_v, marker, cause, horizon) {
  ok <- !is.na(time_v) & !is.na(delta_v) & !is.na(marker)
  if (sum(ok) < 50) return(NA_real_)
  roc <- tryCatch(
    timeROC(T = time_v[ok], delta = delta_v[ok], marker = marker[ok],
            cause = cause, times = horizon, iid = FALSE),
    error = function(e) NULL
  )
  if (is.null(roc)) return(NA_real_)
  v <- if (!is.null(roc$AUC_1)) roc$AUC_1
       else if (!is.null(roc$AUC)) roc$AUC
       else return(NA_real_)
  as.numeric(v[length(v)])
}

# Paired IPCW time-dependent delta-AUC (timeROC) on the common non-missing
# subsample: delta = AUC(marker1) - AUC(marker2) at the horizon, both fit on the
# same rows. This is the registered PRIMARY discrimination metric and the basis
# of the H3 non-inferiority test; delong_delta() below is the SECONDARY
# binary-outcome check, kept for comparison only.
timeroc_delta <- function(df_p, marker1, marker2, time_col, delta_col,
                          cause, horizon) {
  time_v  <- df_p[[time_col]]
  delta_v <- df_p[[delta_col]]
  keep <- !is.na(time_v) & !is.na(delta_v) &
          !is.na(df_p[[marker1]]) & !is.na(df_p[[marker2]])
  if (sum(keep) < 50) return(NA_real_)
  a1 <- timeroc_point_auc(time_v[keep], delta_v[keep], df_p[[marker1]][keep],
                          cause, horizon)
  a2 <- timeroc_point_auc(time_v[keep], delta_v[keep], df_p[[marker2]][keep],
                          cause, horizon)
  if (is.na(a1) || is.na(a2)) return(NA_real_)
  a1 - a2
}

# Delta-AUC point estimate via DeLong on horizon-cap binary outcome, paired.
# Returns numeric scalar (delta = AUC1 - AUC2) or NA.
delong_delta <- function(df_p, marker1, marker2, time_col, delta_col,
                         cause, horizon) {
  time_v  <- df_p[[time_col]]
  delta_v <- df_p[[delta_col]]
  early_censor <- (time_v < horizon) & (delta_v == 0)
  keep <- !early_censor & !is.na(time_v) & !is.na(delta_v) &
          !is.na(df_p[[marker1]]) & !is.na(df_p[[marker2]])
  y    <- as.integer((delta_v == cause) & (time_v <= horizon))
  y_k  <- y[keep]
  m1_k <- df_p[[marker1]][keep]
  m2_k <- df_p[[marker2]][keep]
  if (sum(y_k) < 5 || sum(y_k) == length(y_k)) return(NA_real_)
  r1 <- tryCatch(pROC::roc(y_k, m1_k, quiet = TRUE, direction = "<"),
                 error = function(e) NULL)
  r2 <- tryCatch(pROC::roc(y_k, m2_k, quiet = TRUE, direction = "<"),
                 error = function(e) NULL)
  if (is.null(r1) || is.null(r2)) return(NA_real_)
  as.numeric(pROC::auc(r1) - pROC::auc(r2))
}

# Tie-aware survey-weighted AUC (weighted c-statistic) on a binary label. The
# primary IPCW timeROC AUC is unweighted; this gives the survey-weighted
# discrimination as a sensitivity, to show the design effect on a rank metric.
weighted_auc <- function(marker, y, w) {
  ok <- !is.na(marker) & !is.na(y) & !is.na(w) & w > 0
  marker <- marker[ok]; y <- y[ok]; w <- w[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  Wp <- sum(w[y == 1]); Wn <- sum(w[y == 0])
  if (Wp == 0 || Wn == 0) return(NA_real_)
  ord <- order(marker)
  mv <- marker[ord]; yv <- y[ord]; wv <- w[ord]
  grp <- cumsum(!duplicated(mv))            # ascending unique-value groups
  posW <- tapply(wv * (yv == 1), grp, sum)
  negW <- tapply(wv * (yv == 0), grp, sum)
  posW[is.na(posW)] <- 0; negW[is.na(negW)] <- 0
  cum_neg_below <- cumsum(c(0, negW))[seq_along(negW)]  # neg weight strictly below
  sum(posW * (cum_neg_below + 0.5 * negW)) / (Wp * Wn)
}

# Weighted binary AUC for one score on the horizon-capped outcome (same recoding
# as delong_delta, but survey-weighted by wt_fast).
weighted_binary_auc <- function(df_p, marker, time_col, delta_col, cause,
                                horizon, weight_col = "wt_fast") {
  time_v  <- df_p[[time_col]]; delta_v <- df_p[[delta_col]]
  early_censor <- (time_v < horizon) & (delta_v == 0)
  keep <- !early_censor
  y <- as.integer((delta_v == cause) & (time_v <= horizon))
  weighted_auc(df_p[[marker]][keep], y[keep], df_p[[weight_col]][keep])
}

#-----------------------------------------------------------------------
# One bootstrap replicate
#-----------------------------------------------------------------------

# Returns a list with auc (named numeric), delta_auc (named numeric),
# dca_nb (named numeric). DCA recalibration and net benefit reuse the shared
# competing-risks helpers (recalibrated_risk + cr_net_benefit) so the bootstrap
# matches the point analysis in script 08.
one_rep <- function(rep_idx, boot_rows) {
  df_b <- df[boot_rows, , drop = FALSE]
  frame_for <- function(o) if (isTRUE(o$cause_coded))
    df_b[df_b$cause_coded, , drop = FALSE] else df_b

  auc_out <- list()
  for (out_name in names(outcomes)) {
    o <- outcomes[[out_name]]
    df_o <- frame_for(o)
    for (s in o$scores) {
      for (h in o$horizons) {
        key <- sprintf("%s__%s__t%.1f", s, out_name, h)
        auc_out[[key]] <- timeroc_point_auc(
          df_o[[o$time_col]], df_o[[o$delta_col]],
          df_o[[s]], o$cause, h
        )
      }
    }
  }

  # Survey-weighted (wt_fast) binary AUC, same score x outcome grid (sensitivity).
  wauc_out <- list()
  for (out_name in names(outcomes)) {
    o <- outcomes[[out_name]]
    df_o <- frame_for(o)
    for (s in o$scores) {
      for (h in o$horizons) {
        key <- sprintf("%s__%s__t%.1f", s, out_name, h)
        wauc_out[[key]] <- weighted_binary_auc(
          df_o, s, o$time_col, o$delta_col, o$cause, h
        )
      }
    }
  }

  # Primary delta-AUC = paired IPCW timeROC delta; secondary = binary DeLong.
  delta_out <- list()
  delta_bin_out <- list()
  for (p in pairs) {
    o <- outcomes[[p$outcome]]
    df_o <- frame_for(o)
    for (h in o$horizons) {
      key <- sprintf("%s_vs_%s__%s__t%.1f",
                     p$score1, p$score2, p$outcome, h)
      delta_out[[key]] <- timeroc_delta(
        df_o, p$score1, p$score2,
        o$time_col, o$delta_col, o$cause, h
      )
      delta_bin_out[[key]] <- delong_delta(
        df_o, p$score1, p$score2,
        o$time_col, o$delta_col, o$cause, h
      )
    }
  }

  # Competing-risks DCA net benefit at one threshold per outcome, per score.
  # Recalibrate via the Fine-Gray CIF (Cox for all-cause) on the resampled
  # cohort, then score net benefit with the Aalen-Johansen incidence.
  nb_out <- list()
  for (out_name in names(outcomes)) {
    o <- outcomes[[out_name]]
    df_o <- frame_for(o)
    h <- o$horizons[1]
    for (s in o$scores) {
      key <- sprintf("%s__%s__t%.1f", s, out_name, h)
      risk <- recalibrated_risk(df_o, s, o$time_col, o$status_col,
                                o$cause, h, o$competing)
      nb_out[[key]] <- cr_net_benefit(risk, df_o[[o$time_col]],
                                      df_o[[o$status_col]], o$cause, h,
                                      o$dca_thresh)
    }
  }

  # H4 incremental value: IDI + continuous NRI for adding RMRS to FINDRISC on
  # diabetes mortality (t=14.5). Computed inside the cluster bootstrap so the CI
  # respects the survey design; survIDINRI's own perturbation (used for the
  # point estimate in script 09) treats rows as i.i.d. and is anticonservative.
  # npert=0 returns the point estimates only (the bootstrap supplies the CI).
  idinri_out <- c(idi = NA_real_, nri = NA_real_)
  do <- frame_for(outcomes$dm)
  ok <- !is.na(do$rmrs_score) & !is.na(do$findrisc_score) &
        !is.na(do$followup_years_dm) & !is.na(do$competing_dm)
  do <- do[ok, , drop = FALSE]
  if (nrow(do) > 50 && sum(do$competing_dm == 1) >= 10) {
    indata <- as.matrix(data.frame(
      time   = do$followup_years_dm,
      status = as.integer(do$competing_dm == 1)))   # cause-specific: competing censored
    covs0 <- as.matrix(do[, "findrisc_score", drop = FALSE])
    covs1 <- as.matrix(do[, c("findrisc_score", "rmrs_score")])
    res <- tryCatch(
      survIDINRI::IDI.INF(indata, covs0, covs1, t0 = 14.5, npert = 0),
      error = function(e) NULL)
    if (!is.null(res)) {
      idinri_out["idi"] <- as.numeric(res$m1[1])
      idinri_out["nri"] <- as.numeric(res$m2[1])
    }
  }

  list(auc = unlist(auc_out, use.names = TRUE),
       wauc = unlist(wauc_out, use.names = TRUE),
       delta_auc = unlist(delta_out, use.names = TRUE),
       delta_auc_binary = unlist(delta_bin_out, use.names = TRUE),
       dca_nb = unlist(nb_out, use.names = TRUE),
       idinri = idinri_out)
}

#-----------------------------------------------------------------------
# Point estimates on the original (unresampled) cohort
#-----------------------------------------------------------------------

message("\nComputing point estimates on the original cohort ...")
t_point_start <- Sys.time()
point_est <- one_rep(0L, seq_len(nrow(df)))
t_point_end <- Sys.time()
message(sprintf("  point estimates done in %.1f s",
                as.numeric(difftime(t_point_end, t_point_start, units = "secs"))))

message("\nPoint AUCs (sample):")
print(round(head(point_est$auc, 20), 3))
message("Point delta-AUCs:")
print(round(point_est$delta_auc, 3))
message("Point DCA net benefit (sample):")
print(round(head(point_est$dca_nb, 20), 4))

#-----------------------------------------------------------------------
# Resample plan
#-----------------------------------------------------------------------

set.seed(SEED)
# Stratified PSU bootstrap: within each design stratum, resample its PSU clusters
# with replacement (the textbook NHANES scheme). Singleton strata are kept as-is
# (no within-stratum variance to resample). Each element is one resample's vector
# of cluster ids.
resample_clusters <- function() {
  unlist(lapply(strata_to_clusters, function(cl) {
    if (length(cl) <= 1) cl else sample(cl, length(cl), replace = TRUE)
  }), use.names = FALSE)
}
boot_cluster_lists <- replicate(N_REPS, resample_clusters(), simplify = FALSE)

# Materialize row indices for a resample's cluster ids (a cluster drawn k times
# contributes its rows k times).
build_rows <- function(cluster_vec) {
  unlist(cluster_rows[cluster_vec], use.names = FALSE)
}

#-----------------------------------------------------------------------
# Run bootstrap
#-----------------------------------------------------------------------

run_one <- function(i) {
  t0 <- Sys.time()
  rows <- build_rows(boot_cluster_lists[[i]])
  res <- tryCatch(one_rep(i, rows), error = function(e) {
    message(sprintf("  rep %d failed: %s", i, e$message))
    NULL
  })
  t1 <- Sys.time()
  list(rep = i,
       wall_s = as.numeric(difftime(t1, t0, units = "secs")),
       result = res)
}

message(sprintf("\nRunning %d bootstrap reps ...", N_REPS))
t_boot_start <- Sys.time()

if (HAS_FUTURE && N_WORKERS > 1) {
  message(sprintf("  parallel backend: future.apply multisession workers=%d",
                  N_WORKERS))
  future::plan(future::multisession, workers = N_WORKERS)
  # future_lapply distributes the RNG safely via L'Ecuyer streams seeded
  # from the global seed set above.
  reps_out <- future.apply::future_lapply(
    seq_len(N_REPS), run_one,
    future.seed = TRUE,
    future.globals = c("df", "boot_cluster_lists", "clusters_all",
                       "cluster_rows", "strata_to_clusters", "outcomes", "pairs",
                       "build_rows", "one_rep",
                       "timeroc_point_auc", "timeroc_delta", "delong_delta",
                       "weighted_auc", "weighted_binary_auc",
                       "recalibrated_risk", "cr_net_benefit",
                       "cr_net_benefit_all", ".cif_at"),
    future.packages = c("timeROC", "survival", "riskRegression",
                        "prodlim", "pROC", "survIDINRI")
  )
  future::plan(future::sequential)
} else {
  reps_out <- vector("list", N_REPS)
  for (i in seq_len(N_REPS)) {
    reps_out[[i]] <- run_one(i)
    if (i %% max(1, N_REPS %/% 20) == 0 || i == 1 || i == N_REPS) {
      elapsed <- as.numeric(difftime(Sys.time(), t_boot_start, units = "secs"))
      eta <- elapsed * (N_REPS - i) / i
      message(sprintf("  rep %3d/%d  wall_s=%.1f  elapsed=%.0fs  ETA=%.0fs",
                      i, N_REPS, reps_out[[i]]$wall_s, elapsed, eta))
    }
  }
}

t_boot_end <- Sys.time()
total_wall <- as.numeric(difftime(t_boot_end, t_boot_start, units = "secs"))
message(sprintf("\nBootstrap wall time: %.1f s (%.2f min)",
                total_wall, total_wall / 60))

#-----------------------------------------------------------------------
# Assemble distributions
#-----------------------------------------------------------------------

# Names from point estimate define the column order.
auc_names <- names(point_est$auc)
wauc_names <- names(point_est$wauc)
delta_names <- names(point_est$delta_auc)
delta_bin_names <- names(point_est$delta_auc_binary)
nb_names <- names(point_est$dca_nb)
idinri_names <- names(point_est$idinri)

stack_field <- function(field, expected_names) {
  mat <- matrix(NA_real_, nrow = N_REPS, ncol = length(expected_names),
                dimnames = list(NULL, expected_names))
  for (i in seq_along(reps_out)) {
    r <- reps_out[[i]]$result
    if (is.null(r)) next
    v <- r[[field]]
    if (is.null(v)) next
    common <- intersect(expected_names, names(v))
    mat[i, common] <- v[common]
  }
  mat
}

auc_mat       <- stack_field("auc",              auc_names)
wauc_mat      <- stack_field("wauc",             wauc_names)
delta_mat     <- stack_field("delta_auc",        delta_names)
delta_bin_mat <- stack_field("delta_auc_binary", delta_bin_names)
nb_mat        <- stack_field("dca_nb",           nb_names)
idinri_mat    <- stack_field("idinri",           idinri_names)

#-----------------------------------------------------------------------
# Save full distributions
#-----------------------------------------------------------------------

dir.create("results", showWarnings = FALSE)

auc_dist <- as.list(as.data.frame(auc_mat))
wauc_dist <- as.list(as.data.frame(wauc_mat))
delta_dist <- as.list(as.data.frame(delta_mat))
delta_bin_dist <- as.list(as.data.frame(delta_bin_mat))
nb_dist <- as.list(as.data.frame(nb_mat))
idinri_dist <- as.list(as.data.frame(idinri_mat))

saveRDS(auc_dist,       "results/bootstrap_auc_distributions.rds")
saveRDS(wauc_dist,      "results/bootstrap_wauc_distributions.rds")
saveRDS(delta_dist,     "results/bootstrap_deltaauc_distributions.rds")
saveRDS(delta_bin_dist, "results/bootstrap_deltaauc_binary_distributions.rds")
saveRDS(nb_dist,        "results/bootstrap_dca_netbenefit_distributions.rds")
saveRDS(idinri_dist,    "results/bootstrap_idinri_distributions.rds")

#-----------------------------------------------------------------------
# Build summary CSV
#-----------------------------------------------------------------------

pct <- function(x, q) as.numeric(quantile(x, probs = q, na.rm = TRUE))

summarize_block <- function(mat, point_vec, kind, h3_keys = character(0)) {
  do.call(rbind, lapply(colnames(mat), function(nm) {
    x <- mat[, nm]
    data.frame(
      kind          = kind,
      quantity      = nm,
      point         = as.numeric(point_vec[nm]),
      pct_2_5       = pct(x, 0.025),
      pct_50        = pct(x, 0.50),
      pct_97_5      = pct(x, 0.975),
      n_valid       = sum(!is.na(x)),
      registered_h3 = nm %in% h3_keys,
      stringsAsFactors = FALSE
    )
  }))
}

h3_key <- "rmrs_score_vs_findrisc_score__dm__t14.5"

auc_summary       <- summarize_block(auc_mat,    point_est$auc,       "auc")
wauc_summary      <- summarize_block(wauc_mat,   point_est$wauc,      "weighted_auc")
delta_summary     <- summarize_block(delta_mat,  point_est$delta_auc, "delta_auc",
                                     h3_keys = h3_key)
delta_bin_summary <- summarize_block(delta_bin_mat, point_est$delta_auc_binary,
                                     "delta_auc_binary")
nb_summary        <- summarize_block(nb_mat,     point_est$dca_nb,    "dca_nb")
idinri_summary    <- summarize_block(idinri_mat, point_est$idinri,    "idi_nri")

summary_tbl <- rbind(auc_summary, wauc_summary, delta_summary, delta_bin_summary,
                     nb_summary, idinri_summary)
write.csv(summary_tbl, "results/bootstrap_summary.csv", row.names = FALSE)

# Separate DCA CSV as a convenience for manuscript tables.
nb_summary_clean <- nb_summary
nb_summary_clean$registered_h3 <- NULL
write.csv(nb_summary_clean, "results/bootstrap_dca_netbenefit.csv",
          row.names = FALSE)

#-----------------------------------------------------------------------
# Report the H3 result explicitly
#-----------------------------------------------------------------------

# Primary inference is the IPCW timeROC delta (kind == "delta_auc"); the binary
# DeLong delta (kind == "delta_auc_binary") is reported only as a secondary check.
h3_row <- summary_tbl[summary_tbl$quantity == h3_key &
                        summary_tbl$kind == "delta_auc", , drop = FALSE]
h3_bin <- summary_tbl[summary_tbl$quantity == h3_key &
                        summary_tbl$kind == "delta_auc_binary", , drop = FALSE]
if (nrow(h3_row) == 1) {
  margin <- -0.05
  crosses <- h3_row$pct_2_5 < margin
  message(sprintf(
    "\n=== H3 bootstrap CI (RMRS vs FINDRISC, diabetes mortality, t=14.5) ===\n  PRIMARY (IPCW timeROC delta-AUC):\n    point  = %+0.4f\n    2.5%%   = %+0.4f\n    median = %+0.4f\n    97.5%%  = %+0.4f\n  registered non-inferiority margin: delta-AUC > %+0.2f\n  95%% CI lower bound %s the margin (%s)",
    h3_row$point, h3_row$pct_2_5, h3_row$pct_50, h3_row$pct_97_5,
    margin,
    ifelse(crosses, "CROSSES", "remains above"),
    ifelse(crosses, "non-inferiority NOT established",
                    "non-inferiority established")
  ))
  if (nrow(h3_bin) == 1) {
    message(sprintf(
      "  SECONDARY (binary DeLong delta-AUC, for comparison only):\n    point = %+0.4f  95%% CI (%+0.4f, %+0.4f)",
      h3_bin$point, h3_bin$pct_2_5, h3_bin$pct_97_5))
  }
}

message(sprintf("\nResults written:\n  results/bootstrap_summary.csv (%d rows)\n  results/bootstrap_auc_distributions.rds\n  results/bootstrap_deltaauc_distributions.rds\n  results/bootstrap_deltaauc_binary_distributions.rds\n  results/bootstrap_dca_netbenefit_distributions.rds\n  results/bootstrap_dca_netbenefit.csv",
                nrow(summary_tbl)))

message("\nDone.")

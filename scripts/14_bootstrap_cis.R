#!/usr/bin/env Rscript
# scripts/14_bootstrap_cis.R
#
# Bootstrap 95% percentile CIs at PSU-cluster resamples for the registered
# primary AUC and delta-AUC quantities. This is the registered definitive
# interval estimator for the H3 non-inferiority test (RMRS vs FINDRISC on
# diabetes-related mortality at t = 14.5, registered margin -0.05 AUC).
#
# Resampling unit: NHANES PSU cluster defined as the cross of sdmvstra and
# sdmvpsu (about 301 clusters). Sample whole clusters with replacement, then
# take all subjects in those clusters. Do not sample individuals.
#
# Per resample, recompute:
#   1. AUC at primary horizons for each registered score-outcome pair
#      via timeROC::timeROC(iid = FALSE). iid = TRUE OOMs at N >= 17k.
#   2. Delta-AUCs for the registered pairwise contrasts on the horizon-cap
#      binary outcome via pROC::roc.test DeLong (matches script 09).
#   3. DCA net benefit at one threshold per outcome (7.5% CV, 5% all-cause,
#      2% diabetes), per score, with Cox-recalibrated horizon-specific
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

.libPaths("/home/po/projects/work/longitudinal-mets-validation/renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(timeROC)
  library(survival)
  library(riskRegression)
  library(prodlim)
  library(pROC)
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

message(sprintf("Cohort N = %d, PSU clusters = %d", nrow(df), n_clusters))

# Pre-split row indices by cluster so resampling is index-only.
cluster_rows <- split(seq_len(nrow(df)), df$cluster_id)

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

  delta_out <- list()
  for (p in pairs) {
    o <- outcomes[[p$outcome]]
    df_o <- frame_for(o)
    for (h in o$horizons) {
      key <- sprintf("%s_vs_%s__%s__t%.1f",
                     p$score1, p$score2, p$outcome, h)
      delta_out[[key]] <- delong_delta(
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

  list(auc = unlist(auc_out, use.names = TRUE),
       delta_auc = unlist(delta_out, use.names = TRUE),
       dca_nb = unlist(nb_out, use.names = TRUE))
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
# Each row of boot_cluster_idx is one resample: indices into clusters_all,
# sampled with replacement, length = n_clusters.
boot_cluster_idx <- replicate(N_REPS,
                              sample.int(n_clusters, n_clusters, replace = TRUE),
                              simplify = FALSE)

# Materialize row indices per rep up front.
build_rows <- function(idx_vec) {
  unlist(cluster_rows[clusters_all[idx_vec]], use.names = FALSE)
}

#-----------------------------------------------------------------------
# Run bootstrap
#-----------------------------------------------------------------------

run_one <- function(i) {
  t0 <- Sys.time()
  rows <- build_rows(boot_cluster_idx[[i]])
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
    future.globals = c("df", "boot_cluster_idx", "clusters_all",
                       "cluster_rows", "outcomes", "pairs",
                       "build_rows", "one_rep",
                       "timeroc_point_auc", "delong_delta",
                       "recalibrated_risk", "cr_net_benefit",
                       "cr_net_benefit_all", ".cif_at"),
    future.packages = c("timeROC", "survival", "riskRegression",
                        "prodlim", "pROC")
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
delta_names <- names(point_est$delta_auc)
nb_names <- names(point_est$dca_nb)

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

auc_mat   <- stack_field("auc",       auc_names)
delta_mat <- stack_field("delta_auc", delta_names)
nb_mat    <- stack_field("dca_nb",    nb_names)

#-----------------------------------------------------------------------
# Save full distributions
#-----------------------------------------------------------------------

dir.create("results", showWarnings = FALSE)

auc_dist <- as.list(as.data.frame(auc_mat))
delta_dist <- as.list(as.data.frame(delta_mat))
nb_dist <- as.list(as.data.frame(nb_mat))

saveRDS(auc_dist,   "results/bootstrap_auc_distributions.rds")
saveRDS(delta_dist, "results/bootstrap_deltaauc_distributions.rds")
saveRDS(nb_dist,    "results/bootstrap_dca_netbenefit_distributions.rds")

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

auc_summary   <- summarize_block(auc_mat,   point_est$auc,       "auc")
delta_summary <- summarize_block(delta_mat, point_est$delta_auc, "delta_auc",
                                 h3_keys = h3_key)
nb_summary    <- summarize_block(nb_mat,    point_est$dca_nb,    "dca_nb")

summary_tbl <- rbind(auc_summary, delta_summary, nb_summary)
write.csv(summary_tbl, "results/bootstrap_summary.csv", row.names = FALSE)

# Separate DCA CSV as a convenience for manuscript tables.
nb_summary_clean <- nb_summary
nb_summary_clean$registered_h3 <- NULL
write.csv(nb_summary_clean, "results/bootstrap_dca_netbenefit.csv",
          row.names = FALSE)

#-----------------------------------------------------------------------
# Report the H3 result explicitly
#-----------------------------------------------------------------------

h3_row <- summary_tbl[summary_tbl$quantity == h3_key, , drop = FALSE]
if (nrow(h3_row) == 1) {
  margin <- -0.05
  crosses <- h3_row$pct_2_5 < margin
  message(sprintf(
    "\n=== H3 bootstrap CI (RMRS vs FINDRISC, diabetes mortality, t=14.5) ===\n  point     = %+0.4f\n  2.5%%      = %+0.4f\n  median    = %+0.4f\n  97.5%%     = %+0.4f\n  registered non-inferiority margin: delta-AUC > %+0.2f\n  95%% CI lower bound %s the margin (%s)",
    h3_row$point, h3_row$pct_2_5, h3_row$pct_50, h3_row$pct_97_5,
    margin,
    ifelse(crosses, "CROSSES", "remains above"),
    ifelse(crosses, "non-inferiority NOT established",
                    "non-inferiority established")
  ))
}

message(sprintf("\nResults written:\n  results/bootstrap_summary.csv (%d rows)\n  results/bootstrap_auc_distributions.rds\n  results/bootstrap_deltaauc_distributions.rds\n  results/bootstrap_dca_netbenefit_distributions.rds\n  results/bootstrap_dca_netbenefit.csv",
                nrow(summary_tbl)))

message("\nDone.")

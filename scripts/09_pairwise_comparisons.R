#!/usr/bin/env Rscript
# scripts/09_pairwise_comparisons.R
#
# Pairwise discrimination + reclassification comparisons between risk scores
# within each outcome.
#
# Pairs covered:
#   1. RMRS vs PCE on CV mortality (9.5y)
#   2. RMRS vs Framingham on CV mortality (9.5y)
#   3. PCE vs Framingham on CV mortality (9.5y) clinical head-to-head
#   4. RMRS vs B9 on CV mortality (9.5y)
#   5. RMRS vs B9 on all-cause mortality (9.5y)
#   6. RMRS vs B9 on diabetes mortality (14.5y) broadened tree test
#   7. RMRS vs FINDRISC on diabetes mortality (14.5y) head-to-head
#   8. B9 vs FINDRISC on diabetes mortality (14.5y)
#
# Plus an incremental value test: RMRS added on top of FINDRISC for diabetes
# mortality. Cox linear predictors are compared via pROC at the 14.5y horizon
# and survIDINRI::IDI.INF supplies category-free NRI and IDI.
#
# Methodology notes:
#   * Primary AUC point estimates come from timeROC with iid=FALSE at the
#     same horizons used in the Phase 2 survival scripts (9.5y for allcause
#     / CV, 14.5y for diabetes). IPCW handles right-censoring properly.
#   * Delta-AUC inference (95% CI + DeLong p-value) is computed on the
#     horizon-cap binary outcome via pROC::roc.test. We tried IPCW-based
#     delta-AUC inference with timeROC iid=TRUE but it OOMs at N=17k (the
#     CLAUDE.md tooling rule documents this). The horizon-cap pROC version
#     loses some censoring information but is the standard fallback used
#     in the validation literature.
#   * Horizon-cap encoding for the binary AUC: event = (delta == cause) &
#     (followup <= horizon). Subjects censored before the horizon are
#     excluded to avoid biased AUC estimates.
#
# Inputs:  data/processed/cohort_with_scores.rds
# Outputs: results/pairwise_comparisons.{rds,csv}, results/incremental_dm.{rds,csv}

.libPaths("/home/po/projects/work/longitudinal-mets-validation/renv/library/R-4.3/x86_64-pc-linux-gnu")

suppressMessages({
  library(timeROC)
  library(survival)
  library(survIDINRI)
  library(pROC)
})

df <- readRDS("data/processed/cohort_with_scores.rds")

outcomes <- list(
  allcause = list(time_col = "followup_years",
                  delta_col = "event_allcause",
                  cause = 1,
                  horizon = 9.5,
                  is_competing = FALSE),
  cv       = list(time_col = "followup_years",
                  delta_col = "competing_cv",
                  cause = 1,
                  horizon = 9.5,
                  is_competing = TRUE),
  dm       = list(time_col = "followup_years_dm",
                  delta_col = "competing_dm",
                  cause = 1,
                  horizon = 14.5,
                  is_competing = TRUE)
)

pairs <- list(
  list(label = "RMRS vs PCE on CV mortality",
       score1 = "rmrs_score", score2 = "pce_score", outcome = "cv"),
  list(label = "RMRS vs Framingham on CV mortality",
       score1 = "rmrs_score", score2 = "framingham_score", outcome = "cv"),
  list(label = "PCE vs Framingham on CV mortality",
       score1 = "pce_score", score2 = "framingham_score", outcome = "cv"),
  list(label = "RMRS vs B9 on CV mortality",
       score1 = "rmrs_score", score2 = "b9_score", outcome = "cv"),
  list(label = "RMRS vs B9 on all-cause mortality",
       score1 = "rmrs_score", score2 = "b9_score", outcome = "allcause"),
  list(label = "RMRS vs B9 on diabetes mortality",
       score1 = "rmrs_score", score2 = "b9_score", outcome = "dm"),
  list(label = "RMRS vs FINDRISC on diabetes mortality",
       score1 = "rmrs_score", score2 = "findrisc_score", outcome = "dm"),
  list(label = "B9 vs FINDRISC on diabetes mortality",
       score1 = "b9_score", score2 = "findrisc_score", outcome = "dm")
)

# Extract AUC at the single fitted horizon. timeROC pads the AUC vector with
# a leading t=0 NA so we always take the last (largest-time) element.
roc_auc_last <- function(roc) {
  if (is.null(roc)) return(NA_real_)
  v <- if (!is.null(roc$AUC_1)) roc$AUC_1
       else if (!is.null(roc$AUC)) roc$AUC
       else return(NA_real_)
  as.numeric(v[length(v)])
}

# IPCW time-dependent AUC at the horizon. iid=FALSE keeps memory bounded.
timeroc_point_auc <- function(time_v, delta_v, marker, cause, horizon) {
  roc <- tryCatch(
    timeROC(T = time_v, delta = delta_v, marker = marker,
            cause = cause, times = horizon, iid = FALSE),
    error = function(e) { message("    timeROC failed: ", e$message); NULL }
  )
  roc_auc_last(roc)
}

# DeLong test for paired AUCs on the horizon-cap binary outcome. Returns
# delta = AUC1 - AUC2, 95% CI, and p-value. Subjects with follow-up < horizon
# but no event are dropped (administrative censors only; this avoids the
# IPCW degeneracy that would otherwise need iid=TRUE inference).
delong_paired <- function(df_p, marker1, marker2, time_col, delta_col,
                          cause, horizon) {
  time_v <- df_p[[time_col]]
  delta_v <- df_p[[delta_col]]
  # In-horizon event vs in-horizon non-event. Drop early-censored.
  early_censor <- (time_v < horizon) & (delta_v == 0)
  keep <- !early_censor & !is.na(time_v) & !is.na(delta_v)
  y <- as.integer((delta_v == cause) & (time_v <= horizon))
  y_k <- y[keep]
  m1_k <- df_p[[marker1]][keep]
  m2_k <- df_p[[marker2]][keep]
  if (sum(y_k) < 5 || sum(y_k) == length(y_k)) {
    return(list(delta = NA_real_, lo = NA_real_, hi = NA_real_,
                p = NA_real_, n = length(y_k), events = sum(y_k),
                auc1 = NA_real_, auc2 = NA_real_))
  }
  r1 <- pROC::roc(y_k, m1_k, quiet = TRUE, direction = "<")
  r2 <- pROC::roc(y_k, m2_k, quiet = TRUE, direction = "<")
  cmp <- tryCatch(pROC::roc.test(r1, r2, method = "delong"),
                  error = function(e) NULL)
  if (is.null(cmp)) {
    return(list(delta = as.numeric(pROC::auc(r1) - pROC::auc(r2)),
                lo = NA_real_, hi = NA_real_, p = NA_real_,
                n = length(y_k), events = sum(y_k),
                auc1 = as.numeric(pROC::auc(r1)),
                auc2 = as.numeric(pROC::auc(r2))))
  }
  delta_v <- as.numeric(cmp$estimate[1] - cmp$estimate[2])
  # roc.test returns CI for the difference if computed; otherwise build it.
  # pROC::roc.test stores the difference SE in cmp$statistic / cmp$z; reverse.
  z <- as.numeric(cmp$statistic)
  se <- if (!is.na(z) && z != 0) abs(delta_v) / abs(z) else NA_real_
  list(delta = delta_v,
       lo = delta_v - 1.96 * se,
       hi = delta_v + 1.96 * se,
       p  = as.numeric(cmp$p.value),
       n  = length(y_k),
       events = sum(y_k),
       auc1 = as.numeric(cmp$estimate[1]),
       auc2 = as.numeric(cmp$estimate[2]))
}

results <- list()

for (p in pairs) {
  message(sprintf("\n--- %s ---", p$label))
  meta <- outcomes[[p$outcome]]

  rows <- !is.na(df[[p$score1]]) & !is.na(df[[p$score2]]) &
          !is.na(df[[meta$time_col]]) & !is.na(df[[meta$delta_col]])
  df_p <- df[rows, ]

  # IPCW time-dependent AUC (point estimate) for each marker
  auc1_ipcw <- timeroc_point_auc(
    df_p[[meta$time_col]], df_p[[meta$delta_col]],
    df_p[[p$score1]], meta$cause, meta$horizon
  )
  auc2_ipcw <- timeroc_point_auc(
    df_p[[meta$time_col]], df_p[[meta$delta_col]],
    df_p[[p$score2]], meta$cause, meta$horizon
  )

  # DeLong test on horizon-cap binary AUC (paired)
  dl <- delong_paired(df_p, p$score1, p$score2,
                      meta$time_col, meta$delta_col,
                      meta$cause, meta$horizon)

  results[[p$label]] <- list(
    pair = p, horizon = meta$horizon,
    n_total = nrow(df_p),
    n_eval = dl$n,
    events_in_horizon = dl$events,
    auc_score1_ipcw = auc1_ipcw,
    auc_score2_ipcw = auc2_ipcw,
    auc_score1_binary = dl$auc1,
    auc_score2_binary = dl$auc2,
    delta_auc_binary = dl$delta,
    delta_ci_lo = dl$lo,
    delta_ci_hi = dl$hi,
    delong_p = dl$p
  )

  message(sprintf("  n_total=%d  n_eval=%d  events_in_horizon=%d",
                  nrow(df_p), dl$n, dl$events))
  message(sprintf("  IPCW    AUC(%s)=%.3f  AUC(%s)=%.3f",
                  p$score1, auc1_ipcw, p$score2, auc2_ipcw))
  message(sprintf("  Binary  AUC(%s)=%.3f  AUC(%s)=%.3f",
                  p$score1, dl$auc1, p$score2, dl$auc2))
  message(sprintf("  delta-AUC (binary) = %.3f (95%% CI %.3f, %.3f)  DeLong p=%.4g",
                  dl$delta, dl$lo, dl$hi, dl$p))
}

dir.create("results", showWarnings = FALSE)
saveRDS(results, "results/pairwise_comparisons.rds")

summary_tbl <- do.call(rbind, lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(
    comparison        = nm,
    horizon_y         = r$horizon,
    n_total           = r$n_total,
    n_eval            = r$n_eval,
    events_in_horizon = r$events_in_horizon,
    auc_score1_ipcw   = r$auc_score1_ipcw,
    auc_score2_ipcw   = r$auc_score2_ipcw,
    auc_score1_binary = r$auc_score1_binary,
    auc_score2_binary = r$auc_score2_binary,
    delta_auc_binary  = r$delta_auc_binary,
    ci_lo             = r$delta_ci_lo,
    ci_hi             = r$delta_ci_hi,
    delong_p          = r$delong_p
  )
}))
write.csv(summary_tbl, "results/pairwise_comparisons.csv", row.names = FALSE)

message("\n=== Pairwise comparisons ===")
print(summary_tbl)

#-----------------------------------------------------------------------
# Incremental value: RMRS added on top of FINDRISC for diabetes mortality.
#
# Approach: fit two Cox models on the cause-1 diabetes-related death outcome.
# Compare time-dependent AUC of the linear predictors at the 14.5y horizon
# (IPCW point estimate plus DeLong CI / p-value on the horizon-cap binary
# outcome). Run survIDINRI::IDI.INF for category-free IDI and NRI at the
# same horizon.
#-----------------------------------------------------------------------

message("\n--- Incremental value: RMRS added to FINDRISC on diabetes mortality ---")

meta_dm <- outcomes$dm
inc_rows <- !is.na(df$findrisc_score) & !is.na(df$rmrs_score) &
            !is.na(df[[meta_dm$time_col]]) & !is.na(df[[meta_dm$delta_col]])
df_i <- df[inc_rows, ]
df_i$event_cause1 <- as.integer(df_i[[meta_dm$delta_col]] == 1)

cox_base <- coxph(Surv(followup_years_dm, event_cause1) ~ findrisc_score, data = df_i)
cox_full <- coxph(Surv(followup_years_dm, event_cause1) ~ findrisc_score + rmrs_score,
                  data = df_i)
df_i$lp_base <- predict(cox_base, type = "lp")
df_i$lp_full <- predict(cox_full, type = "lp")

auc_base_ipcw <- timeroc_point_auc(
  df_i$followup_years_dm, df_i$event_cause1, df_i$lp_base,
  cause = 1, horizon = meta_dm$horizon
)
auc_full_ipcw <- timeroc_point_auc(
  df_i$followup_years_dm, df_i$event_cause1, df_i$lp_full,
  cause = 1, horizon = meta_dm$horizon
)

# Reuse the binary DeLong helper on df_i with event_cause1 as the cause
dl_inc <- delong_paired(df_i, "lp_full", "lp_base",
                        "followup_years_dm", "event_cause1",
                        cause = 1, horizon = meta_dm$horizon)

# IDI / NRI via survIDINRI
indicator <- as.matrix(df_i[, c("followup_years_dm", "event_cause1")])
covs_base <- as.matrix(df_i[, "findrisc_score", drop = FALSE])
covs_full <- as.matrix(df_i[, c("findrisc_score", "rmrs_score")])

set.seed(42)
message("  running IDI.INF (200 perturbation reps) ...")
idi_nri <- tryCatch(
  IDI.INF(indicator, covs_base, covs_full,
          t0 = meta_dm$horizon, npert = 200),
  error = function(e) { message("  IDI.INF failed: ", e$message); NULL }
)

inc_result <- list(
  n          = nrow(df_i),
  horizon    = meta_dm$horizon,
  auc_base_ipcw   = auc_base_ipcw,
  auc_full_ipcw   = auc_full_ipcw,
  auc_base_binary = dl_inc$auc2,
  auc_full_binary = dl_inc$auc1,
  delta_auc_binary = dl_inc$delta,
  ci_lo       = dl_inc$lo,
  ci_hi       = dl_inc$hi,
  delong_p    = dl_inc$p,
  idi_nri     = idi_nri
)

saveRDS(inc_result, "results/incremental_dm.rds")

inc_row <- data.frame(
  comparison        = "FINDRISC + RMRS vs FINDRISC on diabetes mortality",
  horizon_y         = meta_dm$horizon,
  n                 = inc_result$n,
  auc_base_ipcw     = inc_result$auc_base_ipcw,
  auc_full_ipcw     = inc_result$auc_full_ipcw,
  auc_base_binary   = inc_result$auc_base_binary,
  auc_full_binary   = inc_result$auc_full_binary,
  delta_auc_binary  = inc_result$delta_auc_binary,
  ci_lo             = inc_result$ci_lo,
  ci_hi             = inc_result$ci_hi,
  delong_p          = inc_result$delong_p,
  idi_est = if (!is.null(idi_nri)) idi_nri$m1[1] else NA_real_,
  idi_lo  = if (!is.null(idi_nri)) idi_nri$m1[2] else NA_real_,
  idi_hi  = if (!is.null(idi_nri)) idi_nri$m1[3] else NA_real_,
  nri_est = if (!is.null(idi_nri)) idi_nri$m2[1] else NA_real_,
  nri_lo  = if (!is.null(idi_nri)) idi_nri$m2[2] else NA_real_,
  nri_hi  = if (!is.null(idi_nri)) idi_nri$m2[3] else NA_real_
)
write.csv(inc_row, "results/incremental_dm.csv", row.names = FALSE)

message(sprintf("\n  n=%d", inc_result$n))
message(sprintf("  IPCW    AUC(FINDRISC)=%.3f  AUC(FINDRISC + RMRS)=%.3f",
                auc_base_ipcw, auc_full_ipcw))
message(sprintf("  Binary  AUC(FINDRISC)=%.3f  AUC(FINDRISC + RMRS)=%.3f",
                dl_inc$auc2, dl_inc$auc1))
message(sprintf("  delta-AUC (binary) = %.3f (95%% CI %.3f, %.3f)  DeLong p=%.4g",
                dl_inc$delta, dl_inc$lo, dl_inc$hi, dl_inc$p))
if (!is.null(idi_nri)) {
  message(sprintf("  IDI = %.4f (95%% CI %.4f, %.4f)",
                  idi_nri$m1[1], idi_nri$m1[2], idi_nri$m1[3]))
  message(sprintf("  Continuous NRI = %.4f (95%% CI %.4f, %.4f)",
                  idi_nri$m2[1], idi_nri$m2[2], idi_nri$m2[3]))
}

message("\nDone.")

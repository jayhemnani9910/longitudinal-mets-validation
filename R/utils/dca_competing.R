# Competing-risks-aware decision-curve helpers, shared by scripts/08_dca.R and
# scripts/14_bootstrap_cis.R so the main analysis and its bootstrap CIs use one
# definition of net benefit.
#
# Two pieces:
#   recalibrated_risk() turns a raw score into an absolute predicted risk of the
#     cause of interest by the horizon. Competing outcomes (CV, diabetes) use the
#     Fine-Gray subdistribution CIF; all-cause uses 1 - Cox survival. Recalibrating
#     a competing-risks outcome with an ordinary Cox model (treating competing
#     deaths as censoring) overstates absolute risk and distorts thresholding.
#   cr_net_benefit() computes net benefit at one decision threshold, estimating
#     the event probability among the test-positive group with the Aalen-Johansen
#     cumulative incidence so that right-censoring and competing deaths are handled
#     rather than ignored.

suppressMessages({
  library(survival)
})

# Aalen-Johansen cumulative incidence of `cause` at `horizon` for a sample.
# status integer: 0 = censored, 1 = cause of interest, 2+ = competing event.
# A two-state (0/1) status reduces to 1 - Kaplan-Meier.
.cif_at <- function(time, status, cause, horizon) {
  ok <- !is.na(time) & !is.na(status)
  time <- time[ok]; status <- status[ok]
  if (length(time) < 1 || all(status == 0)) return(0)
  n_states <- length(unique(status[status != 0]))
  if (n_states <= 1) {
    km <- survival::survfit(survival::Surv(time, as.integer(status == cause)) ~ 1)
    idx <- findInterval(horizon, km$time)
    if (idx < 1) return(0)
    return(1 - km$surv[idx])
  }
  st <- factor(status, levels = sort(unique(c(0, status))))
  fit <- tryCatch(survival::survfit(survival::Surv(time, st) ~ 1),
                  error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  cause_state <- as.character(cause)
  if (is.null(fit$states) || !(cause_state %in% fit$states)) return(NA_real_)
  col <- which(fit$states == cause_state)
  idx <- findInterval(horizon, fit$time)
  if (idx < 1) return(0)
  as.numeric(fit$pstate[idx, col])
}

# Net benefit at one decision threshold, competing-risks aware.
#   NB = TP/N - (FP/N) * (threshold / (1 - threshold))
# with TP/N = P(high) * CIF_high and FP/N = P(high) * (1 - CIF_high), where
# CIF_high is the Aalen-Johansen incidence of the cause within the test-positive
# group at the horizon.
cr_net_benefit <- function(risk, time, status, cause, horizon, threshold) {
  ok <- !is.na(risk) & !is.na(time) & !is.na(status)
  risk <- risk[ok]; time <- time[ok]; status <- status[ok]
  if (length(risk) < 50) return(NA_real_)
  high <- risk >= threshold
  p_high <- mean(high)
  if (p_high == 0) return(0)
  cif_high <- .cif_at(time[high], status[high], cause, horizon)
  if (is.na(cif_high)) return(NA_real_)
  tp <- p_high * cif_high
  fp <- p_high * (1 - cif_high)
  tp - fp * (threshold / (1 - threshold))
}

# Net benefit of the treat-all strategy at a threshold, using the overall CIF.
cr_net_benefit_all <- function(time, status, cause, horizon, threshold) {
  ok <- !is.na(time) & !is.na(status)
  if (sum(ok) < 50) return(NA_real_)
  cif <- .cif_at(time[ok], status[ok], cause, horizon)
  if (is.na(cif)) return(NA_real_)
  cif - (1 - cif) * (threshold / (1 - threshold))
}

# Absolute predicted risk of `cause` by `horizon` from a one-covariate model.
# competing = TRUE uses the Fine-Gray CIF (riskRegression::predictRisk on FGR);
# competing = FALSE uses 1 - Cox survival. status_col is the competing-risk
# status (0/1/2) when competing, or the 0/1 event indicator otherwise.
recalibrated_risk <- function(df_in, score, time_col, status_col, cause,
                              horizon, competing) {
  ok <- !is.na(df_in[[score]]) & !is.na(df_in[[time_col]]) &
        !is.na(df_in[[status_col]])
  out <- rep(NA_real_, nrow(df_in))
  if (sum(ok) < 50) return(out)
  df_s <- df_in[ok, ]
  if (competing) {
    f <- as.formula(sprintf("prodlim::Hist(%s, %s) ~ %s",
                            time_col, status_col, score))
    fit <- tryCatch(riskRegression::FGR(f, data = df_s, cause = cause),
                    error = function(e) NULL)
    if (is.null(fit)) return(out)
    r <- tryCatch(riskRegression::predictRisk(fit, newdata = df_s, times = horizon),
                  error = function(e) NULL)
    if (is.null(r)) return(out)
    out[ok] <- pmin(pmax(as.numeric(r), 0), 1)
  } else {
    f <- as.formula(sprintf("survival::Surv(%s, %s) ~ %s",
                            time_col, status_col, score))
    fit <- tryCatch(survival::coxph(f, data = df_s), error = function(e) NULL)
    if (is.null(fit)) return(out)
    sf <- tryCatch(survival::survfit(fit, newdata = df_s), error = function(e) NULL)
    if (is.null(sf)) return(out)
    idx <- findInterval(horizon, sf$time)
    if (idx < 1) idx <- 1
    out[ok] <- pmin(pmax(1 - sf$surv[idx, ], 0), 1)
  }
  out
}

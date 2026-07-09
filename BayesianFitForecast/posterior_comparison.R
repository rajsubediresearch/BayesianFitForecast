# ============================================================================
# posterior_comparison.R   (for BayesianFitForecast toolbox)
# Raj Subedi
# ============================================================================
# Additive diagnostic script. Reads the fit .Rdata that run_MCMC(_Parallel).R
# saves and produces, in a dedicated output subfolder:
#
#   1. Posterior histograms for every estimated parameter, with the PRIOR
#      overlaid where drawable (identifiability read: where posterior tracks
#      the prior, the data isn't informing that parameter).
#   2. Composite / ratio histograms (e.g. R0 = beta/gamma) with 95% CrI.
#   3. Pairwise scatter for estimated-parameter pairs, annotated with the
#      correlation -- this reveals identifiability RIDGES (|rho| ~ 1 means
#      only the combination, e.g. R0, is identified, not the parts).
#
# Drop this file in the toolbox folder (next to run_analyzeResults.R) and:
#   1. Set OPTIONS_FILE below to the options_*.R for your run.
#   2. source("posterior_comparison.R")
#
# It reads the SAME saved objects the analysis script uses (param_samples,
# composite_samples, pars, composite_expressions), so no other changes needed.
# ============================================================================

rm(list = ls())
suppressMessages({ library(ggplot2) })
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ---- USER CHOICE -----------------------------------------------------------
OPTIONS_FILE <- "options_SISMID_P2_Switzerland.R"   # <- your options file
# ----------------------------------------------------------------------------

source(OPTIONS_FILE)
errorstructure <- c("negativebinomial", "normal", "poisson")
err_name <- errorstructure[errstrc]

# figure out which parameters are ESTIMATED (paramsfix == 0) and their priors
est_idx   <- which(paramsfix == 0)
est_names <- params[est_idx]

# prior string for an estimated parameter, e.g. params1_prior
get_prior <- function(k) {
  nm <- paste0("params", k, "_prior")
  if (exists(nm)) get(nm) else NULL
}
prior_density_fun <- function(prior_str) {
  if (is.null(prior_str) || is.numeric(prior_str)) return(NULL)
  s <- gsub("\\s+", "", prior_str); s <- sub("T\\[[^]]*\\]", "", s)
  dist <- sub("\\(.*", "", s)
  args <- suppressWarnings(as.numeric(strsplit(sub(".*\\((.*)\\)", "\\1", s), ",")[[1]]))
  switch(dist,
    uniform     = function(x) dunif(x, args[1], args[2]),
    normal      = function(x) dnorm(x, args[1], args[2]),
    lognormal   = function(x) dlnorm(x, args[1], args[2]),
    exponential = function(x) dexp(x, args[1]),
    gamma       = function(x) dgamma(x, args[1], args[2]),
    beta        = function(x) dbeta(x, args[1], args[2]),
    cauchy      = function(x) dcauchy(x, args[1], args[2]),
    NULL)
}
est_prior <- setNames(lapply(est_idx, get_prior), est_names)

dir.create("output", showWarnings = FALSE)

for (calibrationperiod in calibrationperiods) {
  # load the fit saved by run_MCMC (same naming the analysis script uses)
  fit_rdata <- paste(model_name, "cal", calibrationperiod, "fcst",
                     forecastinghorizon, err_name, caddisease,
                     "fit.Rdata", sep = "-")
  if (!file.exists(fit_rdata)) {
    message("Fit file not found, skipping: ", fit_rdata); next
  }
  load(fit_rdata)                      # param_samples, composite_samples, pars, ...
  pars <- unlist(pars)

  out_dir <- file.path("output",
                       paste("posterior_comparison", model_name, caddisease,
                             "cal", calibrationperiod, sep = "-"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  tag <- paste0("cal", calibrationperiod)
  cat("Calibration", calibrationperiod, "-> writing to", out_dir, "\n")

  # ---- 1. posterior histograms w/ prior overlay ----
  for (p in pars) {
    draws <- param_samples[[p]]
    if (is.null(draws) || length(draws) == 0) next
    df <- data.frame(value = draws)
    g <- ggplot(df, aes(value)) +
      geom_histogram(aes(y = after_stat(density)), bins = 40,
                     fill = "#2a78d6", alpha = 0.55, color = NA) +
      labs(title = sprintf("Posterior: %s  (%s, %s)", p, model_name, tag),
           x = p, y = "density") +
      theme_minimal(base_size = 12)
    pf <- if (!is.null(est_prior[[p]])) prior_density_fun(est_prior[[p]]) else NULL
    if (!is.null(pf)) {
      xr <- range(draws); xs <- seq(xr[1], xr[2], length.out = 200)
      pri <- data.frame(x = xs, d = sapply(xs, pf)); pri <- pri[is.finite(pri$d), ]
      g <- g + geom_line(data = pri, aes(x, d), color = "#eb6834",
                         linewidth = 1, inherit.aes = FALSE) +
        labs(subtitle = "orange = prior; where posterior tracks prior, data is not informing that parameter")
    }
    ggsave(file.path(out_dir, sprintf("posterior_%s_%s.pdf", p, tag)),
           g, width = 7, height = 4.5)
  }

  # ---- 2. composite / ratio histograms ----
  if (exists("composite_samples")) {
    for (nm in names(composite_samples)) {
      cd <- composite_samples[[nm]]; cd <- cd[is.finite(cd)]
      if (!length(cd)) next
      q <- quantile(cd, c(0.025, 0.5, 0.975))
      g <- ggplot(data.frame(value = cd), aes(value)) +
        geom_histogram(aes(y = after_stat(density)), bins = 40,
                       fill = "#1baf7a", alpha = 0.6, color = NA) +
        geom_vline(xintercept = q, linetype = c("dashed","solid","dashed"),
                   color = "#0f6e56") +
        labs(title = sprintf("Composite: %s  (%s)", nm, tag),
             subtitle = sprintf("median = %.3f   95%% CrI [%.3f, %.3f]", q[2], q[1], q[3]),
             x = nm, y = "density") +
        theme_minimal(base_size = 12)
      ggsave(file.path(out_dir, sprintf("composite_%s_%s.pdf", nm, tag)),
             g, width = 7, height = 4.5)
    }
  }

  # ---- 3. pairwise ridge scatter ----
  present <- intersect(est_names, names(param_samples))
  if (length(present) >= 2) {
    for (cc in combn(present, 2, simplify = FALSE)) {
      x <- param_samples[[cc[1]]]; y <- param_samples[[cc[2]]]
      n <- min(length(x), length(y))
      df <- data.frame(x = x[seq_len(n)], y = y[seq_len(n)])
      rho <- suppressWarnings(cor(df$x, df$y))
      g <- ggplot(df, aes(x, y)) +
        geom_point(alpha = 0.25, size = 0.8, color = "#4a3aa7") +
        labs(title = sprintf("Joint posterior: %s vs %s  (%s)", cc[1], cc[2], tag),
             subtitle = sprintf("correlation = %.3f  %s", rho,
               if (abs(rho) > 0.9) "-> strong ridge: only jointly identified"
               else "-> reasonably separated"),
             x = cc[1], y = cc[2]) +
        theme_minimal(base_size = 12)
      ggsave(file.path(out_dir, sprintf("joint_%s_vs_%s_%s.pdf", cc[1], cc[2], tag)),
             g, width = 6.5, height = 5.5)
    }
  }
}

cat("\nDone. See output/posterior_comparison-* folders.\n")

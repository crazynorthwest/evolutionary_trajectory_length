# Set manually if you want to force a specific xmin.
# Otherwise leave as NA, and xmin will be estimated by KS minimization.
xmin_user <- NA

# Set manually if there is a known physical/experimental upper bound.
# Otherwise leave as NA, and xmax will be the largest observed value.
xmax_user <- generation_total

# Minimum number of observations retained in the tail for candidate xmin.
min_tail_n <- 50

# Candidate xmin values above this quantile are ignored to avoid very tiny tails.
candidate_max_quantile <- 0.90

# Parametric bootstrap repetitions for KS goodness-of-fit p values.
# Larger values are more stable but slower. Use 0 to skip p-value calculation.
ks_bootstrap_B <- 1000
ks_bootstrap_seed <- 123

# Numerical safety constant.
PENALTY <- 1e100
EPS <- .Machine$double.xmin

############################
# 1. Utility functions
############################
log_sum_exp <- function(logv) {
  logv <- logv[is.finite(logv)]
  if (length(logv) == 0) return(-Inf)
  m <- max(logv)
  m + log(sum(exp(logv - m)))
}

safe_value <- function(v) {
  if (length(v) != 1 || is.na(v) || !is.finite(v)) return(PENALTY)
  if (v < 0) return(PENALTY)
  min(as.numeric(v), PENALTY)
}

make_safe_fn <- function(fn) {
  force(fn)
  function(par) {
    v <- tryCatch(fn(par), error = function(e) PENALTY)
    safe_value(v)
  }
}

safe_optim_nm <- function(starts, fn, maxit = 5000) {
  safe_fn <- make_safe_fn(fn)
  best <- NULL
  best_val <- Inf
  for (i in seq_along(starts)) {
    st <- starts[[i]]
    fit <- tryCatch(
      optim(st, safe_fn, method = "Nelder-Mead", control = list(maxit = maxit)),
      error = function(e) NULL
    )
    if (!is.null(fit) && is.finite(fit$value) && fit$value < best_val) {
      best <- fit
      best_val <- fit$value
    }
  }
  if (is.null(best)) stop("Optimization failed for all starting values.")
  best
}

safe_optimize_1d <- function(fn, interval) {
  safe_fn <- function(z) safe_value(tryCatch(fn(z), error = function(e) PENALTY))
  optimize(safe_fn, interval = interval)
}

read_integer_data <- function(file_path) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  dat <- read.csv(file_path, check.names = FALSE)
  if (ncol(dat) < 1) stop("The CSV file has no columns.")

  numeric_cols <- names(dat)[sapply(dat, is.numeric)]
  if (length(numeric_cols) == 0) {
    # Try coercing the first column to numeric.
    x <- suppressWarnings(as.numeric(dat[[1]]))
    col_used <- names(dat)[1]
  } else {
    x <- dat[[numeric_cols[1]]]
    col_used <- numeric_cols[1]
  }

  x <- x[!is.na(x) & is.finite(x)]
  if (length(x) == 0) stop("No valid numeric values were found.")

  if (any(abs(x - round(x)) > 1e-8)) {
    warning("Some values are not integers. They were rounded to the nearest integer.")
  }

  x <- as.integer(round(x))
  x <- x[x > 0]
  if (length(x) == 0) stop("No positive integer values remain after cleaning.")

  list(x = x, column = col_used)
}

# Riemann/Hurwitz-like tail sum for pure discrete Pareto.
# Computes sum_{j=start}^{infinity} j^(-mu).
power_tail_sum <- function(mu, start) {
  if (!is.finite(mu) || mu <= 1 || start < 1) return(Inf)

  # Prefer R's zeta() if available.
  if (exists("zeta", mode = "function")) {
    z <- tryCatch(zeta(mu), error = function(e) NA_real_)
    if (is.finite(z)) {
      before <- if (start > 1) sum((1:(start - 1))^(-mu)) else 0
      val <- z - before
      if (is.finite(val) && val > 0) return(val)
    }
  }

  # Fallback: finite sum plus integral approximation.
  n_end <- max(start + 200000, ceiling(start * 100))
  s <- sum((start:n_end)^(-mu))
  tail <- (n_end + 0.5)^(1 - mu) / (mu - 1)
  s + tail
}

power_range_sum <- function(mu, start, end) {
  if (!is.finite(mu) || start > end) return(0)
  sum((start:end)^(-mu))
}

# Discrete probabilities for continuous distributions by unit-width bins.
# P(X=k) = F(k+0.5) - F(k-0.5)
disc_probs <- function(k, dist, pars) {
  lower <- pmax(k - 0.5, 0)
  upper <- k + 0.5

  if (dist == "lognormal") {
    p <- plnorm(upper, meanlog = pars[1], sdlog = pars[2]) -
      plnorm(lower, meanlog = pars[1], sdlog = pars[2])
  } else if (dist == "weibull") {
    p <- pweibull(upper, shape = pars[1], scale = pars[2]) -
      pweibull(lower, shape = pars[1], scale = pars[2])
  } else if (dist == "gamma") {
    p <- pgamma(upper, shape = pars[1], scale = pars[2]) -
      pgamma(lower, shape = pars[1], scale = pars[2])
  } else {
    stop("Unknown continuous distribution: ", dist)
  }

  pmax(p, EPS)
}

############################
# 2. Model probability functions
############################

finite_model_probs <- function(model, params, xmin, xmax) {
  support <- xmin:xmax

  if (model == "Poisson") {
    logp <- dpois(support, lambda = params$lambda, log = TRUE)
  } else if (model == "Negative binomial") {
    logp <- dnbinom(support, size = params$size, mu = params$mu, log = TRUE)
  } else if (model == "Lognormal") {
    logp <- log(disc_probs(support, "lognormal", c(params$meanlog, params$sdlog)))
  } else if (model == "Weibull") {
    logp <- log(disc_probs(support, "weibull", c(params$shape, params$scale)))
  } else if (model == "Gamma") {
    logp <- log(disc_probs(support, "gamma", c(params$shape, params$scale)))
  } else if (model == "Truncated Pareto") {
    logp <- -params$mu * log(support)
  } else {
    stop("finite_model_probs() does not handle model: ", model)
  }

  logZ <- log_sum_exp(logp)
  probs <- exp(logp - logZ)
  probs / sum(probs)
}

model_cdf_at <- function(fit, q, xmin, xmax) {
  q <- floor(q)
  q[q < xmin] <- xmin - 1

  if (fit$model == "Pure Pareto") {
    mu <- fit$params$mu
    z <- power_tail_sum(mu, xmin)
    out <- numeric(length(q))
    for (i in seq_along(q)) {
      if (q[i] < xmin) {
        out[i] <- 0
      } else {
        out[i] <- power_range_sum(mu, xmin, q[i]) / z
      }
    }
    return(pmin(pmax(out, 0), 1))
  }

  support <- xmin:xmax
  probs <- finite_model_probs(fit$model, fit$params, xmin, xmax)
  cdf <- cumsum(probs)

  idx <- pmin(pmax(q - xmin + 1, 0), length(support))
  out <- ifelse(idx <= 0, 0, cdf[idx])
  pmin(pmax(out, 0), 1)
}

model_ccdf_on_support <- function(fit, support, xmin, xmax) {
  if (fit$model == "Pure Pareto") {
    mu <- fit$params$mu
    z <- power_tail_sum(mu, xmin)
    return(sapply(support, function(k) power_tail_sum(mu, k) / z))
  }

  probs <- finite_model_probs(fit$model, fit$params, xmin, xmax)
  ccdf <- rev(cumsum(rev(probs)))
  ccdf
}

ks_statistic <- function(x, fit, xmin, xmax) {
  points <- sort(unique(x))
  emp_cdf <- sapply(points, function(z) mean(x <= z))
  mod_cdf <- model_cdf_at(fit, points, xmin, xmax)
  max(abs(emp_cdf - mod_cdf))
}

# Simulate data from a fitted model.
# For finite-support models, samples are drawn from [xmin, xmax].
# For pure Pareto, the infinite support is approximated by a large finite support.
simulate_from_fit <- function(fit, n, xmin, xmax) {
  if (fit$model == "Pure Pareto") {
    mu <- fit$params$mu
    # Large approximation support for the discrete pure Pareto.
    # Increase this if your fitted mu is very close to 1 and extremely long simulated tails are expected.
    max_support <- max(xmax * 20, xmin + 10000, 50000)
    support <- xmin:max_support
    logp <- -mu * log(support)
    logZ <- log_sum_exp(logp)
    probs <- exp(logp - logZ)
    return(sample(support, size = n, replace = TRUE, prob = probs))
  }

  support <- xmin:xmax
  probs <- finite_model_probs(fit$model, fit$params, xmin, xmax)
  sample(support, size = n, replace = TRUE, prob = probs)
}

# Refit the same model to a bootstrap sample.
refit_same_model <- function(model_name, x_boot, xmin, xmax) {
  if (model_name == "Poisson") {
    fit_poisson(x_boot, xmin, xmax)
  } else if (model_name == "Negative binomial") {
    fit_negative_binomial(x_boot, xmin, xmax)
  } else if (model_name == "Lognormal") {
    fit_lognormal(x_boot, xmin, xmax)
  } else if (model_name == "Weibull") {
    fit_weibull(x_boot, xmin, xmax)
  } else if (model_name == "Gamma") {
    fit_gamma(x_boot, xmin, xmax)
  } else if (model_name == "Pure Pareto") {
    fit_pure_pareto(x_boot, xmin, xmax)
  } else if (model_name == "Truncated Pareto") {
    # If a pure-Pareto bootstrap sample exceeds the observed xmax this function is not used.
    # For finite-support TPL samples, all values are within [xmin, xmax].
    fit_truncated_pareto(x_boot[x_boot <= xmax], xmin, xmax)
  } else {
    stop("Unknown model for refitting: ", model_name)
  }
}

# Parametric-bootstrap KS goodness-of-fit p value.
# p_value = Pr(KS_boot >= KS_observed)
# Because parameters are estimated from the data, this bootstrap p value is preferred over
# a standard one-sample ks.test() p value.
bootstrap_ks_pvalue <- function(x_tail, fit, xmin, xmax, B = 200, seed = 123) {
  if (is.null(B) || is.na(B) || B <= 0) {
    fit$KS_p_value <- NA_real_
    fit$KS_bootstrap_B <- 0
    return(fit)
  }

  set.seed(seed)
  n <- length(x_tail)
  obs_KS <- fit$KS
  boot_KS <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    xb <- tryCatch(
      simulate_from_fit(fit, n = n, xmin = xmin, xmax = xmax),
      error = function(e) NULL
    )

    if (is.null(xb) || length(xb) < 5) next

    fit_b <- tryCatch(
      refit_same_model(fit$model, xb, xmin, xmax),
      error = function(e) NULL
    )

    if (is.null(fit_b)) next

    fit_b <- tryCatch(
      add_model_metrics(fit_b, xb, xmin, xmax),
      error = function(e) NULL
    )

    if (!is.null(fit_b) && is.finite(fit_b$KS)) {
      boot_KS[b] <- fit_b$KS
    }
  }

  valid <- boot_KS[is.finite(boot_KS)]
  if (length(valid) == 0) {
    fit$KS_p_value <- NA_real_
    fit$KS_bootstrap_B <- 0
  } else {
    # Add-one smoothing avoids reporting exactly 0 when B is finite.
    fit$KS_p_value <- (sum(valid >= obs_KS) + 1) / (length(valid) + 1)
    fit$KS_bootstrap_B <- length(valid)
  }

  fit
}


############################
# 3. Model fitting functions
############################

fit_poisson <- function(x, xmin, xmax) {
  support <- xmin:xmax
  n <- length(x)

  nll <- function(par) {
    lambda <- exp(par[1])
    log_obs <- dpois(x, lambda = lambda, log = TRUE)
    log_sup <- dpois(support, lambda = lambda, log = TRUE)
    logZ <- log_sum_exp(log_sup)
    -sum(log_obs - logZ)
  }

  m <- max(mean(x), 1)
  starts <- list(log(m), log(median(x)), log(max(x)))
  opt <- safe_optim_nm(starts, nll)
  lambda <- exp(opt$par[1])
  ll <- -safe_value(nll(log(lambda)))
  list(model = "Poisson", params = list(lambda = lambda), logLik = ll, k = 1, n = n)
}

fit_negative_binomial <- function(x, xmin, xmax) {
  support <- xmin:xmax
  n <- length(x)
  m <- mean(x)
  v <- var(x)
  size0 <- ifelse(is.finite(v) && v > m, m^2 / (v - m), 1000)
  size0 <- max(size0, 1e-3)

  nll <- function(par) {
    size <- exp(par[1])
    mu <- exp(par[2])
    log_obs <- dnbinom(x, size = size, mu = mu, log = TRUE)
    log_sup <- dnbinom(support, size = size, mu = mu, log = TRUE)
    logZ <- log_sum_exp(log_sup)
    -sum(log_obs - logZ)
  }

  starts <- list(
    c(log(size0), log(m)),
    c(log(1), log(m)),
    c(log(10), log(m)),
    c(log(100), log(median(x)))
  )
  opt <- safe_optim_nm(starts, nll)
  size <- exp(opt$par[1])
  mu <- exp(opt$par[2])
  ll <- -safe_value(nll(log(c(size, mu))))
  list(model = "Negative binomial", params = list(size = size, mu = mu), logLik = ll, k = 2, n = n)
}

fit_lognormal <- function(x, xmin, xmax) {
  support <- xmin:xmax
  n <- length(x)
  lx <- log(x)
  meanlog0 <- mean(lx)
  sdlog0 <- max(sd(lx), 0.05)

  nll <- function(par) {
    meanlog <- par[1]
    sdlog <- exp(par[2])
    log_obs <- log(disc_probs(x, "lognormal", c(meanlog, sdlog)))
    log_sup <- log(disc_probs(support, "lognormal", c(meanlog, sdlog)))
    logZ <- log_sum_exp(log_sup)
    -sum(log_obs - logZ)
  }

  starts <- list(
    c(meanlog0, log(sdlog0)),
    c(log(median(x)), log(sdlog0)),
    c(meanlog0, log(0.2)),
    c(meanlog0, log(0.8))
  )
  opt <- safe_optim_nm(starts, nll)
  meanlog <- opt$par[1]
  sdlog <- exp(opt$par[2])
  ll <- -safe_value(nll(c(meanlog, log(sdlog))))
  list(model = "Lognormal", params = list(meanlog = meanlog, sdlog = sdlog), logLik = ll, k = 2, n = n)
}

fit_weibull <- function(x, xmin, xmax) {
  support <- xmin:xmax
  n <- length(x)
  m <- mean(x)

  nll <- function(par) {
    shape <- exp(par[1])
    scale <- exp(par[2])
    log_obs <- log(disc_probs(x, "weibull", c(shape, scale)))
    log_sup <- log(disc_probs(support, "weibull", c(shape, scale)))
    logZ <- log_sum_exp(log_sup)
    -sum(log_obs - logZ)
  }

  starts <- list(
    c(log(1), log(m)),
    c(log(0.7), log(m)),
    c(log(1.5), log(m)),
    c(log(2.5), log(m)),
    c(log(1), log(median(x)))
  )
  opt <- safe_optim_nm(starts, nll)
  shape <- exp(opt$par[1])
  scale <- exp(opt$par[2])
  ll <- -safe_value(nll(log(c(shape, scale))))
  list(model = "Weibull", params = list(shape = shape, scale = scale), logLik = ll, k = 2, n = n)
}

fit_gamma <- function(x, xmin, xmax) {
  support <- xmin:xmax
  n <- length(x)
  m <- mean(x)
  v <- var(x)
  shape0 <- ifelse(is.finite(v) && v > 0, m^2 / v, 1)
  scale0 <- ifelse(is.finite(v) && m > 0, v / m, m)
  shape0 <- max(shape0, 1e-3)
  scale0 <- max(scale0, 1e-3)

  nll <- function(par) {
    shape <- exp(par[1])
    scale <- exp(par[2])
    log_obs <- log(disc_probs(x, "gamma", c(shape, scale)))
    log_sup <- log(disc_probs(support, "gamma", c(shape, scale)))
    logZ <- log_sum_exp(log_sup)
    -sum(log_obs - logZ)
  }

  starts <- list(
    c(log(shape0), log(scale0)),
    c(log(1), log(m)),
    c(log(2), log(m / 2)),
    c(log(5), log(m / 5))
  )
  opt <- safe_optim_nm(starts, nll)
  shape <- exp(opt$par[1])
  scale <- exp(opt$par[2])
  ll <- -safe_value(nll(log(c(shape, scale))))
  list(model = "Gamma", params = list(shape = shape, scale = scale), logLik = ll, k = 2, n = n)
}

fit_truncated_pareto <- function(x, xmin, xmax) {
  support <- xmin:xmax
  n <- length(x)

  nll <- function(mu) {
    if (!is.finite(mu) || mu <= 0) return(PENALTY)
    logZ <- log_sum_exp(-mu * log(support))
    -sum(-mu * log(x) - logZ)
  }

  opt <- safe_optimize_1d(nll, interval = c(0.01, 10))
  mu <- opt$minimum
  ll <- -safe_value(nll(mu))
  list(model = "Truncated Pareto", params = list(mu = mu), logLik = ll, k = 1, n = n)
}

fit_pure_pareto <- function(x, xmin, xmax = NULL) {
  n <- length(x)

  nll <- function(mu) {
    if (!is.finite(mu) || mu <= 1) return(PENALTY)
    z <- power_tail_sum(mu, xmin)
    if (!is.finite(z) || z <= 0) return(PENALTY)
    -sum(-mu * log(x) - log(z))
  }

  opt <- safe_optimize_1d(nll, interval = c(1.0001, 10))
  mu <- opt$minimum
  ll <- -safe_value(nll(mu))
  list(model = "Pure Pareto", params = list(mu = mu), logLik = ll, k = 1, n = n)
}

add_model_metrics <- function(fit, x, xmin, xmax) {
  fit$AIC <- 2 * fit$k - 2 * fit$logLik
  fit$BIC <- log(length(x)) * fit$k - 2 * fit$logLik
  fit$KS <- ks_statistic(x, fit, xmin, xmax)
  fit
}

fit_all_models <- function(x, xmin, xmax) {
  x_tail <- x[x >= xmin & x <= xmax]
  if (length(x_tail) < 5) stop("Too few observations in the selected tail.")

  fits <- list(
    fit_poisson(x_tail, xmin, xmax),
    fit_negative_binomial(x_tail, xmin, xmax),
    fit_lognormal(x_tail, xmin, xmax),
    fit_weibull(x_tail, xmin, xmax),
    fit_gamma(x_tail, xmin, xmax),
    fit_pure_pareto(x_tail, xmin, xmax),
    fit_truncated_pareto(x_tail, xmin, xmax)
  )

  fits <- lapply(fits, add_model_metrics, x = x_tail, xmin = xmin, xmax = xmax)
  fits
}

fit_to_row <- function(fit) {
  par_text <- paste(
    paste(names(fit$params), signif(unlist(fit$params), 6), sep = "="),
    collapse = "; "
  )
  data.frame(
    model = fit$model,
    n_tail = fit$n,
    n_parameters = fit$k,
    logLik = fit$logLik,
    AIC = fit$AIC,
    BIC = fit$BIC,
    KS = fit$KS,
    KS_p_value = if (is.null(fit$KS_p_value)) NA_real_ else fit$KS_p_value,
    KS_bootstrap_B = if (is.null(fit$KS_bootstrap_B)) NA_integer_ else fit$KS_bootstrap_B,
    parameters = par_text,
    stringsAsFactors = FALSE
  )
}

############################
# 4. xmin estimation
############################

estimate_xmin_by_tpl_ks <- function(x, xmax, min_tail_n = 50, candidate_max_quantile = 0.90) {
  x <- x[x > 0 & x <= xmax]
  candidate_upper <- as.numeric(quantile(x, candidate_max_quantile, names = FALSE))
  candidates <- sort(unique(x[x <= candidate_upper]))
  candidates <- candidates[sapply(candidates, function(z) sum(x >= z & x <= xmax) >= min_tail_n)]

  if (length(candidates) == 0) {
    stop("No candidate xmin has at least min_tail_n observations. Reduce min_tail_n.")
  }

  rows <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    xmin <- candidates[i]
    x_tail <- x[x >= xmin & x <= xmax]
    fit <- fit_truncated_pareto(x_tail, xmin, xmax)
    fit <- add_model_metrics(fit, x_tail, xmin, xmax)
    rows[[i]] <- data.frame(
      xmin = xmin,
      xmax = xmax,
      n_tail = length(x_tail),
      mu = fit$params$mu,
      logLik = fit$logLik,
      AIC = fit$AIC,
      KS = fit$KS
    )
  }

  res <- do.call(rbind, rows)
  res <- res[order(res$KS, -res$n_tail), ]
  rownames(res) <- NULL
  res
}

############################
# 5. Main analysis
############################

x <- data_clean

xmax <- ifelse(is.na(xmax_user), max(x), as.integer(xmax_user))
if (xmax < max(x)) {
  warning("xmax_user is smaller than the largest observation. Values above xmax were removed.")
}

if (is.na(xmin_user)) {
  xmin_search <- estimate_xmin_by_tpl_ks(
    x = x,
    xmax = xmax,
    min_tail_n = min_tail_n,
    candidate_max_quantile = candidate_max_quantile
  )
  xmin <- xmin_search$xmin[1]
} else {
  xmin <- as.integer(xmin_user)
  xmin_search <- NULL
}

x_tail <- x[x >= xmin & x <= xmax]

summary_df <- data.frame(
  n_original = length(x),
  min_original = min(x),
  max_original = max(x),
  xmin = xmin,
  xmax = xmax,
  n_tail = length(x_tail),
  stringsAsFactors = FALSE
)

fits <- fit_all_models(x, xmin, xmax)

# Add parametric-bootstrap KS goodness-of-fit p values to each fitted model.
# This step may take some time because each bootstrap sample is refitted.
if (!is.null(ks_bootstrap_B) && ks_bootstrap_B > 0) {
  cat("Calculating parametric-bootstrap KS p values. B =", ks_bootstrap_B, "\n")
  fits <- lapply(seq_along(fits), function(i) {
    bootstrap_ks_pvalue(
      x_tail = x_tail,
      fit = fits[[i]],
      xmin = xmin,
      xmax = xmax,
      B = ks_bootstrap_B,
      seed = ks_bootstrap_seed + i
    )
  })
}

comparison <- do.call(rbind, lapply(fits, fit_to_row))
comparison$delta_AIC <- comparison$AIC - min(comparison$AIC)
comparison$AIC_weight <- exp(-0.5 * comparison$delta_AIC) / sum(exp(-0.5 * comparison$delta_AIC))
comparison <- comparison[order(comparison$AIC), ]
rownames(comparison) <- NULL

write.csv(
  comparison,
  file.path(path_figures, "model_comparison.csv"),
  row.names = FALSE
)

support <- xmin:xmax
emp_x <- sort(unique(x_tail))
emp_ccdf <- sapply(emp_x, function(k) mean(x_tail >= k))

emp_df <- data.frame(
  x = emp_x,
  ccdf = emp_ccdf,
  model = "Simulation"
)

y_min <- max(min(emp_ccdf[emp_ccdf > 0]) / 2, 1e-6)

fit_df <- do.call(rbind, lapply(fits, function(fit) {
  ccdf <- model_ccdf_on_support(fit, support, xmin, xmax)
  ccdf <- pmax(ccdf, y_min / 10)
  
  data.frame(
    x = support,
    ccdf = ccdf,
    model = fit$model
  )
}))

fit_models <- sapply(fits, function(z) z$model)
legend_order <- c("Simulation", fit_models)

plot_cols <- c(
  "Simulation" = "black",
  "Poisson" = "#999999",
  "Negative binomial" = "#7B3294",
  "Lognormal" = "#008837",
  "Weibull" = "#E66101",
  "Gamma" = "#5E3C99",
  "Pure Pareto" = "#1F78B4",
  "Truncated Pareto" = "#D73027"
)

plot_cols <- plot_cols[legend_order]

plot_ltys <- c(
  "Simulation" = "blank",
  "Poisson" = "dotted",
  "Negative binomial" = "dotdash",
  "Lognormal" = "dashed",
  "Weibull" = "longdash",
  "Gamma" = "twodash",
  "Pure Pareto" = "dashed",
  "Truncated Pareto" = "solid"
)

plot_ltys <- plot_ltys[legend_order]

plot_lwds <- c(
  "Simulation" = 0,
  "Poisson" = 0.7,
  "Negative binomial" = 0.7,
  "Lognormal" = 0.7,
  "Weibull" = 0.7,
  "Gamma" = 0.7,
  "Pure Pareto" = 0.8,
  "Truncated Pareto" = 1.3
)

plot_lwds <- plot_lwds[legend_order]



files <- data.frame(
  test = c("test1", "test2", "test3", "test4", "test5"),
  parameter = qqq,
  file = c(
    paste0(path_figures,"/power_law_test1.csv"),
    paste0(path_figures,"/power_law_test2.csv"),
    paste0(path_figures,"/power_law_test3.csv"),
    paste0(path_figures,"/power_law_test4.csv"),
    paste0(path_figures,"/power_law_test5.csv")
  ),
  stringsAsFactors = FALSE
)

# Column name containing integer data. If this column is absent, the first numeric column is used.
value_col <- "Steps"

# xmin selection model:
# "truncated_pareto" = choose xmin by minimizing KS for hard-truncated Pareto
# "pure_pareto"      = choose xmin by minimizing KS for pure discrete Pareto
xmin_selection_model <- "truncated_pareto"

# Minimum tail observations required when searching xmin.
min_tail_n <- 100

# To speed up xmin search, candidate xmin values can be thinned.
# Use Inf to use all unique candidate xmins.
max_xmin_candidates <- 500

# If your experiment/simulation has a known physical upper bound, set it here, e.g. 5000.
# If NULL, the largest observed value across all groups is used.
xmax_common_user <- generation_total

# Bootstrap settings for confidence intervals of truncated Pareto mu.
# Set B_boot <- 0 to skip bootstrap.
B_boot <- 1000
set.seed(123)


############################
# 1. Numerical utility functions
############################

PENALTY <- 1e100

good <- function(z) all(is.finite(z))

clamp01 <- function(z) {
  z[!is.finite(z)] <- NA_real_
  pmin(pmax(z, 0), 1)
}

safe_prob <- function(p) pmax(p, .Machine$double.xmin)

logsumexp <- function(a) {
  m <- max(a)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(a - m)))
}

# log(exp(loga) - exp(logb)), requiring loga >= logb.
logspace_sub <- function(loga, logb) {
  out <- rep(-Inf, length(loga))
  ok <- is.finite(loga) & is.finite(logb) & (loga >= logb)
  out[ok] <- loga[ok] + log1p(-exp(logb[ok] - loga[ok]))
  # If logb is -Inf, subtraction is exp(loga)
  ok2 <- is.finite(loga) & !is.finite(logb)
  out[ok2] <- loga[ok2]
  out
}

thin_candidates <- function(x, max_n = 500) {
  x <- sort(unique(as.integer(x)))
  if (is.infinite(max_n) || length(x) <= max_n) return(x)
  idx <- unique(round(seq(1, length(x), length.out = max_n)))
  x[idx]
}

read_group_data <- function(files, value_col = "Steps") {
  out <- list()
  for (i in seq_len(nrow(files))) {
    f <- files$file[i]
    if (!file.exists(f) && files$test[i] == "test1") {
      ff <- fallback_files[file.exists(fallback_files)]
      if (length(ff) > 0) f <- ff[1]
    }
    if (!file.exists(f)) {
      stop("File not found: ", f,
           "\nPut the CSV files in the working directory or edit the file paths.")
    }

    dat <- read.csv(f)
    if (!is.null(value_col) && value_col %in% names(dat)) {
      x <- dat[[value_col]]
    } else {
      num_cols <- names(dat)[sapply(dat, is.numeric)]
      if (length(num_cols) == 0) stop("No numeric column found in ", f)
      x <- dat[[num_cols[1]]]
      message("Using numeric column '", num_cols[1], "' for file: ", f)
    }

    x <- x[is.finite(x)]
    x <- round(x)
    x <- as.integer(x[x > 0])

    out[[i]] <- data.frame(
      test = files$test[i],
      parameter = files$parameter[i],
      value = x
    )
  }
  do.call(rbind, out)
}

empirical_cdf_at <- function(x, grid) {
  sapply(grid, function(z) mean(x <= z))
}

empirical_ccdf <- function(x) {
  grid <- sort(unique(x))
  data.frame(x = grid, ccdf = sapply(grid, function(z) mean(x >= z)))
}

# Safe objective wrapper. This is the key fix for:
# "L-BFGS-B needs finite values of fn".
# Even if the original likelihood returns Inf, -Inf, NA, NaN, or throws an error,
# the optimizer will only see a large but finite penalty.
make_safe_fn <- function(fn) {
  force(fn)
  function(par) {
    val <- tryCatch(fn(par), error = function(e) PENALTY)
    if (length(val) != 1 || is.na(val) || !is.finite(val)) return(PENALTY)
    as.numeric(min(val, PENALTY))
  }
}

# Multi-start optimizer. All optim() calls receive the safe objective, not the raw nll.
optim_multistart <- function(starts, fn, lower, upper, method = "L-BFGS-B") {
  safe_fn <- make_safe_fn(fn)
  best <- NULL
  best_value <- Inf

  for (i in seq_len(nrow(starts))) {
    st <- as.numeric(starts[i, ])
    st <- pmin(pmax(st, lower), upper)

    fit <- tryCatch(
      optim(st, safe_fn, method = method, lower = lower, upper = upper,
            control = list(maxit = 5000)),
      error = function(e) NULL
    )

    if (!is.null(fit) && is.finite(fit$value) && fit$value < best_value) {
      best <- fit
      best_value <- fit$value
    }
  }

  if (is.null(best)) {
    best <- list(par = as.numeric(pmin(pmax(starts[1, ], lower), upper)),
                 value = PENALTY, convergence = 1)
  }
  best
}

############################
# 2. Discrete Pareto functions
############################

zeta_approx <- function(alpha, xmin) {
  if (!is.finite(alpha) || alpha <= 1) return(Inf)
  K <- max(100000L, as.integer(xmin + 100000L))
  j <- xmin:K
  s <- sum(j^(-alpha))
  tail <- ((K + 0.5)^(1 - alpha)) / (alpha - 1)
  s + tail
}

pure_pareto_ccdf_vec <- function(grid, alpha, xmin) {
  if (alpha <= 1 || !is.finite(alpha)) return(rep(NA_real_, length(grid)))
  grid <- pmax(as.integer(grid), xmin)
  K <- max(max(grid) + 100000L, 100000L)
  support <- xmin:K
  p <- support^(-alpha)
  tail_K <- ((K + 0.5)^(1 - alpha)) / (alpha - 1)
  denom <- sum(p) + tail_K
  rev_cum <- rev(cumsum(rev(p)))
  idx <- grid - xmin + 1L
  clamp01((rev_cum[idx] + tail_K) / denom)
}

pure_pareto_cdf_vec <- function(grid, alpha, xmin) {
  clamp01(1 - pure_pareto_ccdf_vec(grid + 1L, alpha, xmin))
}

truncated_pareto_cdf_vec <- function(grid, mu, xmin, xmax) {
  support <- xmin:xmax
  logp <- -mu * log(support)
  p <- exp(logp - logsumexp(logp))
  cdf <- cumsum(p)
  out <- numeric(length(grid))
  out[grid < xmin] <- 0
  out[grid >= xmax] <- 1
  ok <- grid >= xmin & grid < xmax
  out[ok] <- cdf[grid[ok] - xmin + 1L]
  clamp01(out)
}

truncated_pareto_ccdf_vec <- function(grid, mu, xmin, xmax) {
  support <- xmin:xmax
  logp <- -mu * log(support)
  p <- exp(logp - logsumexp(logp))
  ccdf <- rev(cumsum(rev(p)))
  out <- numeric(length(grid))
  out[grid < xmin] <- 1
  out[grid > xmax] <- 0
  ok <- grid >= xmin & grid <= xmax
  out[ok] <- ccdf[grid[ok] - xmin + 1L]
  clamp01(out)
}

############################
# 3. Generic CCDF/CDF helpers for truncated-at-xmin continuous models
############################

# These are left-truncated at xmin - 0.5 because the data are integers.
# The fitted mass for integer k is F(k+0.5)-F(k-0.5), conditional on X >= xmin-0.5.

continuous_left_trunc_cdf <- function(grid, xmin, log_survival_fun) {
  log_den <- log_survival_fun(pmax(xmin - 0.5, 0))
  log_tail <- log_survival_fun(pmax(grid + 0.5, 0))
  out <- 1 - exp(log_tail - log_den)
  out[grid < xmin] <- 0
  clamp01(out)
}

continuous_left_trunc_ccdf <- function(grid, xmin, log_survival_fun) {
  log_den <- log_survival_fun(pmax(xmin - 0.5, 0))
  log_tail <- log_survival_fun(pmax(grid - 0.5, 0))
  out <- exp(log_tail - log_den)
  out[grid <= xmin] <- 1
  clamp01(out)
}

continuous_interval_logmass <- function(x, log_survival_fun) {
  lo <- pmax(x - 0.5, 0)
  hi <- x + 0.5
  logspace_sub(log_survival_fun(lo), log_survival_fun(hi))
}

############################
# 4. Fitting functions
############################

fit_truncated_pareto <- function(x, xmin, xmax) {
  x <- as.integer(x[x >= xmin & x <= xmax])
  if (length(x) < 2) stop("Too few observations for truncated Pareto.")
  support <- xmin:xmax
  log_support <- log(support)
  nll <- function(mu) {
    if (!is.finite(mu) || mu <= 0) return(PENALTY)
    logZ <- logsumexp(-mu * log_support)
    val <- -sum(-mu * log(x) - logZ)
    ifelse(is.finite(val), val, PENALTY)
  }
  opt <- optimize(nll, interval = c(0.01, 20))
  mu <- opt$minimum
  list(distribution = "Truncated Pareto",
       params = c(mu = mu, xmin = xmin, xmax = xmax),
       k = 1, logLik = -opt$objective, n = length(x),
       cdf = function(grid) truncated_pareto_cdf_vec(grid, mu, xmin, xmax),
       ccdf = function(grid) truncated_pareto_ccdf_vec(grid, mu, xmin, xmax))
}

fit_pure_pareto <- function(x, xmin) {
  x <- as.integer(x[x >= xmin])
  nll <- function(alpha) {
    if (!is.finite(alpha) || alpha <= 1) return(PENALTY)
    logZ <- log(zeta_approx(alpha, xmin))
    val <- -sum(-alpha * log(x) - logZ)
    ifelse(is.finite(val), val, PENALTY)
  }
  opt <- optimize(nll, interval = c(1.0001, 20))
  alpha <- opt$minimum
  list(distribution = "Pure Pareto",
       params = c(alpha = alpha, xmin = xmin),
       k = 1, logLik = -opt$objective, n = length(x),
       cdf = function(grid) pure_pareto_cdf_vec(grid, alpha, xmin),
       ccdf = function(grid) pure_pareto_ccdf_vec(grid, alpha, xmin))
}

fit_poisson <- function(x, xmin) {
  x <- as.integer(x[x >= xmin])
  nll <- function(par) {
    lambda <- exp(par[1])
    log_den <- ppois(xmin - 1, lambda, lower.tail = FALSE, log.p = TRUE)
    if (!is.finite(log_den)) return(PENALTY)
    val <- -sum(dpois(x, lambda, log = TRUE) - log_den)
    ifelse(is.finite(val), val, PENALTY)
  }
  starts <- matrix(log(c(mean(x), median(x), max(mean(x), 1))), ncol = 1)
  opt <- optim_multistart(starts, nll, lower = log(1e-8), upper = log(max(x) * 50))
  lambda <- exp(opt$par[1])
  list(distribution = "Poisson",
       params = c(lambda = lambda, xmin = xmin),
       k = 1, logLik = -opt$value, n = length(x),
       cdf = function(grid) {
         den <- ppois(xmin - 1, lambda, lower.tail = FALSE)
         clamp01((ppois(grid, lambda) - ppois(xmin - 1, lambda)) / den)
       },
       ccdf = function(grid) {
         den <- ppois(xmin - 1, lambda, lower.tail = FALSE)
         clamp01(ppois(grid - 1, lambda, lower.tail = FALSE) / den)
       })
}

fit_negative_binomial <- function(x, xmin) {
  x <- as.integer(x[x >= xmin])
  m <- mean(x); v <- var(x)
  size_start <- ifelse(is.finite(v) && v > m, m^2 / (v - m), 1000)
  size_start <- max(size_start, 1e-3)
  nll <- function(par) {
    size <- exp(par[1]); mu <- exp(par[2])
    log_den <- pnbinom(xmin - 1, size = size, mu = mu, lower.tail = FALSE, log.p = TRUE)
    if (!is.finite(log_den)) return(PENALTY)
    val <- -sum(dnbinom(x, size = size, mu = mu, log = TRUE) - log_den)
    ifelse(is.finite(val), val, PENALTY)
  }
  starts <- rbind(
    c(log(size_start), log(m)),
    c(log(1), log(m)),
    c(log(10), log(m)),
    c(log(1000), log(m))
  )
  opt <- optim_multistart(starts, nll,
                          lower = c(log(1e-8), log(1e-8)),
                          upper = c(log(1e8), log(max(x) * 50)))
  size <- exp(opt$par[1]); mu <- exp(opt$par[2])
  list(distribution = "Negative binomial",
       params = c(size = size, mu = mu, xmin = xmin),
       k = 2, logLik = -opt$value, n = length(x),
       cdf = function(grid) {
         den <- pnbinom(xmin - 1, size = size, mu = mu, lower.tail = FALSE)
         clamp01((pnbinom(grid, size = size, mu = mu) -
                    pnbinom(xmin - 1, size = size, mu = mu)) / den)
       },
       ccdf = function(grid) {
         den <- pnbinom(xmin - 1, size = size, mu = mu, lower.tail = FALSE)
         clamp01(pnbinom(grid - 1, size = size, mu = mu, lower.tail = FALSE) / den)
       })
}

fit_lognormal_discrete <- function(x, xmin) {
  x <- as.integer(x[x >= xmin])
  lx <- log(x)
  sd0 <- max(sd(lx), 1e-3)
  nll <- function(par) {
    meanlog <- par[1]; sdlog <- exp(par[2])
    log_surv <- function(q) plnorm(q, meanlog, sdlog, lower.tail = FALSE, log.p = TRUE)
    log_den <- log_surv(pmax(xmin - 0.5, 0))
    log_mass <- continuous_interval_logmass(x, log_surv)
    if (!is.finite(log_den) || any(!is.finite(log_mass))) return(PENALTY)
    val <- -sum(log_mass - log_den)
    ifelse(is.finite(val), val, PENALTY)
  }
  starts <- rbind(
    c(mean(lx), log(sd0)),
    c(log(mean(x)), log(sd0)),
    c(median(lx), log(sd0 * 1.5)),
    c(mean(lx), log(max(sd0 * 0.5, 1e-3)))
  )
  opt <- optim_multistart(starts, nll,
                          lower = c(log(1e-6), log(1e-5)),
                          upper = c(log(max(x) * 100), log(10)))
  meanlog <- opt$par[1]; sdlog <- exp(opt$par[2])
  log_surv <- function(q) plnorm(q, meanlog, sdlog, lower.tail = FALSE, log.p = TRUE)
  list(distribution = "Lognormal",
       params = c(meanlog = meanlog, sdlog = sdlog, xmin = xmin),
       k = 2, logLik = -opt$value, n = length(x),
       cdf = function(grid) continuous_left_trunc_cdf(grid, xmin, log_surv),
       ccdf = function(grid) continuous_left_trunc_ccdf(grid, xmin, log_surv))
}

fit_gamma_discrete <- function(x, xmin) {
  x <- as.integer(x[x >= xmin])
  m <- mean(x); v <- var(x)
  shape_start <- max(m^2 / v, 1e-3)
  scale_start <- max(v / m, 1e-3)
  nll <- function(par) {
    shape <- exp(par[1]); scale <- exp(par[2])
    log_surv <- function(q) pgamma(q, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
    log_den <- log_surv(pmax(xmin - 0.5, 0))
    log_mass <- continuous_interval_logmass(x, log_surv)
    if (!is.finite(log_den) || any(!is.finite(log_mass))) return(PENALTY)
    val <- -sum(log_mass - log_den)
    ifelse(is.finite(val), val, PENALTY)
  }
  starts <- rbind(
    c(log(shape_start), log(scale_start)),
    c(log(1), log(m)),
    c(log(2), log(m / 2)),
    c(log(10), log(max(m / 10, 1e-6)))
  )
  opt <- optim_multistart(starts, nll,
                          lower = c(log(1e-8), log(1e-8)),
                          upper = c(log(1e8), log(1e8)))
  shape <- exp(opt$par[1]); scale <- exp(opt$par[2])
  log_surv <- function(q) pgamma(q, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
  list(distribution = "Gamma",
       params = c(shape = shape, scale = scale, xmin = xmin),
       k = 2, logLik = -opt$value, n = length(x),
       cdf = function(grid) continuous_left_trunc_cdf(grid, xmin, log_surv),
       ccdf = function(grid) continuous_left_trunc_ccdf(grid, xmin, log_surv))
}

fit_weibull_discrete <- function(x, xmin) {
  x <- as.integer(x[x >= xmin])
  nll <- function(par) {
    shape <- exp(par[1]); scale <- exp(par[2])
    log_surv <- function(q) pweibull(q, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
    log_den <- log_surv(pmax(xmin - 0.5, 0))
    log_mass <- continuous_interval_logmass(x, log_surv)
    if (!is.finite(log_den) || any(!is.finite(log_mass))) return(PENALTY)
    val <- -sum(log_mass - log_den)
    ifelse(is.finite(val), val, PENALTY)
  }
  # Multiple starting points avoid the original error:
  # optim(..., method="L-BFGS-B") cannot start when fn is Inf.
  starts <- rbind(
    c(log(1), log(mean(x))),
    c(log(0.5), log(median(x))),
    c(log(1.5), log(mean(x))),
    c(log(2), log(mean(x))),
    c(log(3), log(max(mean(x), 1)))
  )
  opt <- optim_multistart(starts, nll,
                          lower = c(log(1e-8), log(1e-8)),
                          upper = c(log(1e8), log(max(x) * 100)))
  shape <- exp(opt$par[1]); scale <- exp(opt$par[2])
  log_surv <- function(q) pweibull(q, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
  list(distribution = "Weibull",
       params = c(shape = shape, scale = scale, xmin = xmin),
       k = 2, logLik = -opt$value, n = length(x),
       cdf = function(grid) continuous_left_trunc_cdf(grid, xmin, log_surv),
       ccdf = function(grid) continuous_left_trunc_ccdf(grid, xmin, log_surv))
}

############################
# 5. KS statistic 
############################

ks_statistic <- function(x, fit) {
  grid <- sort(unique(x))
  emp <- empirical_cdf_at(x, grid)
  theo <- fit$cdf(grid)
  max(abs(emp - theo), na.rm = TRUE)
}

# Approximate one-sample KS p-value based on the fitted CDF.
# Note: because parameters are estimated from the same data and data are discrete,
# this p-value is an approximate goodness-of-fit indicator.
ks_pvalue_approx <- function(D, n) {
  if (!is.finite(D) || !is.finite(n) || n <= 0) return(NA_real_)
  en <- sqrt(n)
  z <- (en + 0.12 + 0.11 / en) * D
  terms <- sapply(1:100, function(k) (-1)^(k - 1) * exp(-2 * k^2 * z^2))
  p <- 2 * sum(terms)
  pmin(pmax(p, 0), 1)
}

estimate_xmin <- function(x, xmax, model = "truncated_pareto",
                                 min_tail_n = 100, max_candidates = 500) {
  candidates <- sort(unique(as.integer(x)))
  candidates <- candidates[candidates < xmax]
  candidates <- candidates[sapply(candidates, function(z) sum(x >= z & x <= xmax)) >= min_tail_n]
  candidates <- thin_candidates(candidates, max_candidates)
  if (length(candidates) == 0) stop("No xmin candidates satisfy min_tail_n. Try reducing min_tail_n.")

  res <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    xmin <- candidates[i]
    x_tail <- x[x >= xmin & x <= xmax]
    fit <- if (model == "truncated_pareto") {
      fit_truncated_pareto(x_tail, xmin, xmax)
    } else if (model == "pure_pareto") {
      fit_pure_pareto(x_tail, xmin)
    } else stop("Unknown xmin selection model.")
    D <- ks_statistic(x_tail, fit)
    est_name <- ifelse(model == "truncated_pareto", "mu", "alpha")
    res[[i]] <- data.frame(xmin = xmin, n_tail = length(x_tail),
                           estimate = as.numeric(fit$params[est_name]), KS = D)
  }
  out <- do.call(rbind, res)
  out[order(out$KS), ][1,]
}

############################
# 6. Fit all seven models
############################

fit_all_models <- function(x, xmin, xmax) {
  x_tail <- as.integer(x[x >= xmin & x <= xmax])
  fits <- list(
    fit_poisson(x_tail, xmin),
    fit_negative_binomial(x_tail, xmin),
    fit_lognormal_discrete(x_tail, xmin),
    fit_weibull_discrete(x_tail, xmin),
    fit_gamma_discrete(x_tail, xmin),
    fit_pure_pareto(x_tail, xmin),
    fit_truncated_pareto(x_tail, xmin, xmax)
  )
  tab <- do.call(rbind, lapply(fits, function(fit) {
    D <- ks_statistic(x_tail, fit)
    p_value <- ks_pvalue_approx(D, fit$n)
    aic <- 2 * fit$k - 2 * fit$logLik
    bic <- log(fit$n) * fit$k - 2 * fit$logLik
    data.frame(distribution = fit$distribution, n_tail = fit$n,
               logLik = fit$logLik, k = fit$k, AIC = aic, BIC = bic, KS = D,
               p_value = p_value,
               parameters = paste(names(fit$params), signif(fit$params, 6),
                                  sep = "=", collapse = "; "),
               stringsAsFactors = FALSE)
  }))
  tab$delta_AIC <- tab$AIC - min(tab$AIC)
  tab <- tab[order(tab$AIC), ]
  rownames(tab) <- NULL
  list(fits = fits, table = tab)
}

bootstrap_mu_tpl <- function(x, xmin, xmax, B = 500) {
  x_tail <- as.integer(x[x >= xmin & x <= xmax])
  if (length(x_tail) < 5 || B <= 0) return(c(mu_mean = NA, mu_lower = NA, mu_upper = NA))
  mu_boot <- numeric(B)
  for (b in seq_len(B)) {
    xb <- sample(x_tail, length(x_tail), replace = TRUE)
    mu_boot[b] <- fit_truncated_pareto(xb, xmin, xmax)$params["mu"]
  }
  c(mu_mean = mean(mu_boot),
    mu_lower = as.numeric(quantile(mu_boot, 0.025, na.rm = TRUE)),
    mu_upper = as.numeric(quantile(mu_boot, 0.975, na.rm = TRUE)))
}
############################
# 7. Main analysis
############################

dat <- read_group_data(files, value_col = value_col)
dat$test <- factor(dat$test, levels = files$test)
dat$parameter <- as.numeric(dat$parameter)
all_x <- dat$value

xmax_common <- if (is.null(xmax_common_user)) max(all_x) else xmax_common_user
xmax_common <- as.integer(round(xmax_common))

xmin_search <- NULL
for (q_tem in qqq) {
  xmin_search_tem <- estimate_xmin(
    x = dat$value[dat$parameter == q_tem],
    xmax = xmax_common,
    model = xmin_selection_model,
    min_tail_n = min_tail_n,
    max_candidates = max_xmin_candidates
  ) 
  xmin_search <- rbind(xmin_search, xmin_search_tem)
}



xmin_common <- min(xmin_search$xmin)

common_summary <- data.frame(
  xmin_selection_model = xmin_selection_model,
  xmin_common = xmin_common,
  xmax_common = xmax_common,
  pooled_tail_n = sum(all_x >= xmin_common & all_x <= xmax_common),
  pooled_KS = xmin_search$KS[1],
  stringsAsFactors = FALSE
)

all_results <- list()
mu_results <- list()

for (tt in levels(dat$test)) {
  z <- dat[dat$test == tt, ]
  x <- z$value
  par_value <- unique(z$parameter)

  fit_obj <- fit_all_models(x, xmin_common, xmax_common)
  tab <- fit_obj$table
  tab$test <- tt
  tab$parameter <- par_value
  tab$xmin_common <- xmin_common
  tab$xmax_common <- xmax_common
  tab <- tab[, c("test", "parameter", "distribution", "n_tail", "xmin_common", "xmax_common",
                 "logLik", "k", "AIC", "delta_AIC", "BIC", "KS", "p_value", "parameters")]
  all_results[[tt]] <- tab

  tpl_idx <- which(sapply(fit_obj$fits, function(f) f$distribution) == "Truncated Pareto")
  tpl_fit <- fit_obj$fits[[tpl_idx]]
  mu <- tpl_fit$params["mu"]
  ci <- if (B_boot > 0) bootstrap_mu_tpl(x, xmin_common, xmax_common, B = B_boot) else c(mu_mean = NA, mu_lower = NA, mu_upper = NA)

  mu_results[[tt]] <- data.frame(test = tt, parameter = par_value,
                                 n_tail = sum(x >= xmin_common & x <= xmax_common),
                                 xmin_common = xmin_common, xmax_common = xmax_common,
                                 mu = as.numeric(mu), mu_boot_mean = ci["mu_mean"],
                                 mu_lower = ci["mu_lower"], mu_upper = ci["mu_upper"],
                                 stringsAsFactors = FALSE)

  safe_par <- gsub("\\.", "_", as.character(par_value))
}

comparison_table <- do.call(rbind, all_results)
comparison_table <- comparison_table[order(comparison_table$parameter, comparison_table$AIC), ]
rownames(comparison_table) <- NULL
write.csv(comparison_table, file.path(path_figures, "model_comparison_by_group.csv"), row.names = FALSE)

mu_table <- do.call(rbind, mu_results)
mu_table <- mu_table[order(mu_table$parameter), ]
rownames(mu_table) <- NULL
write.csv(mu_table, file.path(path_figures, "truncated_pareto_mu_by_parameter.csv"), row.names = FALSE)


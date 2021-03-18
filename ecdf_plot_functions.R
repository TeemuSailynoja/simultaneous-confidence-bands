library(bayesplot)
library(ggplot2)
library(parallel)
library(future.apply)
source("alpha_adjustment_functions.R")

get_gamma <- function(n, chains=1, alpha=.95, K=n, nsamples=10000) {
  #' Get Gamma
  #' 
  #' @description This function simulates an adjusted value of gamma to account
  #' for multiplicity when forming an 1-alpha level confidence envelope for
  #' the ECDF of a sample from the uniform distribution on the interval (0,1).
  #' 
  #' @param n integer. Desired sample size.
  #' @param nsamples integer. The number of simualtions to run.
  #' @param chains integer. The number of samples to compare.
  #' @param alpha float. The desired confidence level between 0 and 1.
  #' @param K integer. Default: n. The desired granularity of the uniform partition of the unit interval.
  #' The ecdf will be evaluated at points z_k = k/K, for k= 1,...,K-1.
  #' @return The adjusted gamma value.
  
  if (chains == 1) {
    # Use the optimization-based method
    gamma <- adjust_alpha_optimize(1-alpha, n, K)
  } else if (chains == 2) {
    
    gamma <- adjust_alpha_optimize_chains(1-alpha, n, chains, K = K)
  } else {
    gamma <- adjust_alpha_simulate_chains(1-alpha, n, chains, K, M = nsamples)
  }
  gamma
}

get_lims <- function(n, gamma, chains=1, alpha = 0.95, K = n, nsamples=10000, verbose=TRUE){
  # Simulate limits s.t. 
  # p(lims$lower[i] <= sum(unif(n) < 1/i) <= lims$upper[i] for any i = 1, ..., K-1) = alpha
  #' Get Lims
  #' 
  #' @description This function computes the simultaneous alpha level confidence
  #' bands for 1. uniformity on the interval (0,1) in the case of one chain and 
  #' 2. chains sharing an underlying distribution for multiple chains.
  #' 
  #' @param n integer. Desired sample size.
  #' @param gamma float. Precomputed coverage parameter gamma. Optional.
  #' @param chains integer. The number of samples to compare.
  #' @param alpha float. The desired confidence level between 0 and 1.
  #' @param K integer. Default: n. The desired granularity of the uniform partition of the unit interval.
  #' The ecdf will be evaluated at points z_k = k/K, for k= 1,...,K-1.
  #' @param nsamples integer. The number of simualtions to run.
  #' 
  #' @return List containing the upper and lower simultaneous confidence bands.
  #' Evaluated at z_i = i/K for i = 1, ..., K-1.
  if (missing(gamma)) {
      gamma <- get_gamma(n, chains, alpha, K, nsamples)
    }
  if (chains == 1) {
    lims <- list(
      lower = qbinom(gamma/2, n, (0:K)/K), 
      upper = qbinom(1 - gamma/2, n, (0:K)/K))
  } else {
    m <- n
    n1 <- n * (chains - 1)
    k <- floor( (0:K)/K * n * chains)
    lims <- list(
      lower = qhyper(gamma / 2, m, n1, k),
      upper = qhyper(1 - gamma / 2, m, n1, k)
    )
  }
  lims
}

unif_conf_band <- function(y, yrep, lw, pit, gamma, ..., K=0, chains=0, nsamples=10000, conf=.95, alpha=1, size=1, verbose = FALSE) {
  #' ECDF plot with simultaneous confidence bands for uniformity.
  #' 
  #' @description This function plots the ECDF plot(s) of the given sample(s).
  #' If more than one sample are given, the values are transformed to 
  #' fractional ranks before drawing the ECDF plots.
  #' 
  #'
  #' @param gamma float. Precomputed coverage parameter gamma. Optional.
  #' @param K integer. Default: n. The desired granularity of the uniform partition of the unit interval.
  #' The ecdf will be evaluated at points z_k = k/K, for k= 1,...,K-1.
  #' @param chains integer. The number of samples to compare.
  #' @param nsamples integer. The number of simualtions to run.
  #' @param conf float. The desired confidence level between 0 and 1.
  #' @param alpha float [0,1]. Alpha level of the ecdf plots.
  #' @param size float. Size parameter for the ecdf plots.
  #' 
  #' @return ecdf plot with simultaneous confidence bands.
  bayesplot:::check_ignored_arguments(...)
  if (!missing(pit)) {
    stopifnot(is.numeric(pit))
    message("'pit' specified so ignoring 'y','yrep','lw' if specified.")
  } else {
    # work in progress.
    bayesplot:::suggested_package("rstantools")
    y <- bayesplot:::validate_y(y)
    yrep <- bayesplot:::validate_yrep(yrep, y)
    stopifnot(identical(dim(yrep), dim(lw)))
    pit <- matrix(rstantools::loo_pit(object = yrep, y = y, lw = lw),nrow = 1)
  }
  n <- dim(pit)[2]
  if (chains == 0) {
    chains <- dim(pit)[1]
  }
  if (chains != 1) {
    pit <- u_scale(pit)
  }
  if (K == 0) {
    K <- n
  }
  p <- seq(0,1,length.out = K + 1)
  if (missing(gamma)) {
    gamma <- get_gamma(n, chains, conf, K, nsamples)
  } 
  ecdf_lims <- get_lims(n, gamma = gamma, chains = chains, alpha = conf,
                        K = K, nsamples = nsamples, verbose = verbose
                        )
  ecdfs <- t(apply(pit, 1, function(t) ecdf(t)(p)))
  
  ggplot() +
    geom_step(
      data = bayesplot:::melt_yrep(ecdfs),
      mapping = aes_(
        x = ~ (y_id-1) / K,
        y = ~ value,
        group = ~ rep_id,
        color = ~ factor(rep_id)),
      size = size,
      alpha = alpha
    ) +
    geom_step(
      data = data.frame(lower = ecdf_lims$lower, upper = ecdf_lims$upper),
      mapping = aes_(x = p, y = ~ upper /n)) +
    geom_step(
      data = data.frame(lower = ecdf_lims$lower, upper = ecdf_lims$upper),
      mapping = aes_(x = p, y = ~ lower /n)) +
    labs(
      color = if (chains == 1) element_blank() else"Chain:",
      x = if (chains == 1) element_blank() else element_text("Fractional rank"),
      y = element_blank())
}

unif_conf_band_diff <- function(y, yrep, lw, pit, gamma, ..., K=0, chains=0, nsamples=1000, conf=.95, alpha=1, size=1, verbose = FALSE) {
  #' ECDF plot with simultaneous confidence bands for uniformity.
  #' 
  #' @description This function plots the ECDF difference plot(s) of the given
  #' sample(s). If more than one sample are given, the values are transformed to 
  #' fractional ranks before drawing the ECDF plots.
  #' 
  #'
  #' @param gamma float. Precomputed coverage parameter gamma. Optional.
  #' @param K integer. Default: n. The desired granularity of the uniform partition of the unit interval.
  #' The ecdf will be evaluated at points z_k = k/K, for k= 1,...,K-1.
  #' @param chains integer. The number of samples to compare.
  #' @param nsamples integer. The number of simualtions to run.
  #' @param conf float. The desired confidence level between 0 and 1.
  #' @param alpha float [0,1]. Alpha level of the ecdf plots.
  #' @param size float. Size parameter for the ecdf plots.
  #' 
  #' @return ecdf difference plot with simultaneous confidence bands.
  bayesplot:::check_ignored_arguments(...)
  if (!missing(pit)) {
    stopifnot(is.numeric(pit))
    message("'pit' specified so ignoring 'y','yrep','lw' if specified.")
  } else {
    # work in progress
    bayesplot:::suggested_package("rstantools")
    y <- bayesplot:::validate_y(y)
    yrep <- bayesplot:::validate_yrep(yrep, y)
    stopifnot(identical(dim(yrep), dim(lw)))
    pit <- matrix(rstantools::loo_pit(object = yrep, y = y, lw = lw),nrow = 1)
  }
  n <- dim(pit)[2]
  if (chains == 0) {
    chains <- dim(pit)[1]
  }
  if (chains != 1) {
    pit <- u_scale(pit)
  }
  if (K == 0) {
    K <- n
  }
  
  p <- seq(0,1,length.out = K + 1)
  if (missing(gamma)) {
    gamma <- get_gamma(n, chains, conf, K, nsamples)
  } 
  ecdf_lims <- get_lims(n, gamma = gamma, chains = chains, alpha = conf, K = K, nsamples = nsamples, verbose = verbose)
  ecdfs <- t(apply(pit, 1, function(t) ecdf(t)(p) - p))
  ggplot() +
    geom_step(
      data = bayesplot:::melt_yrep(ecdfs),
      mapping = aes_(x = ~ (y_id-1) / K, y = ~value, group = ~rep_id, color = ~ factor(rep_id)),
      size = size,
      alpha = alpha
    ) +
    geom_step(
      data = data.frame(lower = ecdf_lims$lower, upper = ecdf_lims$upper),
      mapping = aes_(x = p, y = ~ upper / n - p)) +
    geom_step(
      data = data.frame(lower = ecdf_lims$lower, upper = ecdf_lims$upper), 
      mapping = aes_(x = p, y = ~ lower / n - p)) +
    hline_at(
      0,
      size = 0.1,
      linetype = 2
    ) +
    labs(
      color = if (chains == 1) element_blank() else element_text("Chain:"),
      x = if (chains == 1) element_blank() else element_text("Fractional rank"),
      y = element_blank())
}


# Used in plot_notebook-Rmd to streamline drawing histograms with confidence intervals. 
draw_hist <- function(sample, bins=50, alpha=.95, row=1, title_str="", fsize = 21, limit_y = TRUE) {
  y_max <- 0
  if (isTRUE(dim(sample)[1] > 1)) {
    data <- u_scale(sample)
    if (limit_y == TRUE) {
      for (idx in 1:dim(sample)[1]) {
        y_max <- max(
          y_max,
          max(hist(data[idx,], breaks = seq(0,1,length.out = bins + 1), plot = FALSE)$counts)
        )
      }
    }
    data <- data[row, ]
  } else {
    data <- sample
    y_max <- max(
      y_max,
      max(hist(data, breaks = seq(0,1,length.out = bins + 1), plot = FALSE)$counts)
    )
  }
  CI <- qbinom(c((1-alpha)/2, 0.5, 1-(1-alpha)/2),
               size = length(data),
               prob = 1 / bins)
  ggplot(data = data.frame(X1 = data)) + 
    geom_histogram(mapping = aes_(x = ~ X1), breaks = seq(0,1,length.out = bins+1), fill = colors[row], colour="black") +
    geom_polygon(
      data = data.frame(
        x = c(-.05, 0, 1, 0, -.05, 1.05, 1, 1.05,-.05),
        y = c(CI[1], CI[2], CI[2], CI[2],  CI[3], CI[3], CI[2], CI[1], CI[1])
      ),
      aes(x = x, y = y),
      fill = "grey75",
      color = "grey50",
      alpha = 0.6
    ) +
    labs(
      y = element_blank(),
      x = if (isTRUE(dim(sample)[1] > 1)) {"Fractional rank"} else {""}
    ) +
    scale_x_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8, 1.0)) +
    scale_y_continuous(limits = c(NA, max(y_max, CI[3])), n.breaks = 5) +
    theme(
      plot.title = element_text(hjust=.5),
      text = element_text(size=fsize)) +
    ggtitle(title_str)
}
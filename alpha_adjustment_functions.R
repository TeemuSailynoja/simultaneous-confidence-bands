# adjust alpha level via simulations
# @param alpha desired overall type I error rate
# @param N sample size
# @param K granularity size
# @param ncases number of discrete cases. ncases = 0 for continuous values
# @param M number of simulation trials
adjust_alpha_simulate <- function(alpha, N, K = N, ncases=0, M = 10000) {
  z <- get_partition(K)
  gamma <- numeric(M)
  for (i in seq_len(M)) {
    x <- sample_unif(N, ncases)
    ecdfs <- eval_ecdf(x, z)
    # P(X <= x)
    probs1 <- pbinom(ecdfs, N, z)
    # P(X > x - 1) = P(X >= x)
    probs2 <- pbinom(ecdfs - 1, N, z, lower.tail = FALSE)
    # take the minimum across the whole ecdf and consider both sides
    gamma[i] <- 2 * min(probs1, probs2)
  }
  alpha_quantile(gamma, alpha)
}

# adjust alpha level via simulations
# @param alpha desired overall type I error rate
# @param niterations number of iterations
# @param nchains number of chains
# @param K granularity size
# @param ncases number of discrete cases. ncases = 0 for continuous values
# @param M number of simulation trials
# @param ties.method method to rank equal values in discrete data
adjust_alpha_simulate_chains <- function(alpha, niterations, nchains,
                                         K = niterations, ncases=0, M = 10000, ties.method = 'max') {
  z <- get_partition(K)
  ndraws <- niterations * nchains
  # parameters of the hypergeometric distribution
  m <- niterations
  n <- niterations * (nchains - 1)
  k <- floor(z * ndraws)
  gamma <- numeric(M)
  for (i in seq_len(M)) {
    # original distribution of 'x' doesn't matter 
    # as long at it is the samme across all chains
    x <- sample_unif(ndraws, ncases)
    x <- matrix(x, nrow = niterations, ncol = nchains)
    u <- u_scale(x, ties.method)
    ecdfs <- eval_ecdf_matrix(u, z)
    gamma_chains <- numeric(nchains)
    for (j in seq_len(nchains)) {
      # P(X <= x)
      probs1 <- phyper(ecdfs[, j], m, n, k)
      # P(X > x - 1) = P(X >= x)
      probs2 <- phyper(ecdfs[, j] - 1, m, n, k, lower.tail = FALSE)
      # take the minimum across the whole ecdf and consider both sides
      gamma_chains[j] <- 2 * min(probs1, probs2)
    }
    # minimum probability across chains
    gamma[i] <- min(gamma_chains)
  }
  alpha_quantile(gamma, alpha)
}

# compute the ECDF and compare with thresholds
test_uniformity <- function(x, alpha, K = length(x)) {
  N <- length(x)
  z <- get_partition(K)
  thres_lb <- qbinom(alpha / 2, N, z)
  thres_ub <- qbinom(1 - alpha / 2, N, z)
  ecdfs <- eval_ecdf(x, z)
  any(ecdfs < thres_lb | ecdfs > thres_ub)
}

# Test for uniformity 
# @param ties.method method to rank equal values in discrete data
test_uniformity_chains <- function(x, alpha, K = length(x), ties.method = 'max') {
  niterations <- nrow(x)
  nchains <- ncol(x)
  ndraws <- niterations * nchains
  z <- get_partition(K)
  # rank transformations
  u <- u_scale(x, ties.method)
  ecdfs <- eval_ecdf_matrix(u, z)
  # parameters of the hypergeometric distribution
  m <- niterations
  n <- niterations * (nchains - 1)
  k <- floor(z * ndraws)
  # comparison to thresholds
  thres_lb <- qhyper(alpha / 2, m, n, k)
  thres_ub <- qhyper(1 - alpha / 2, m, n, k)
  any(ecdfs < thres_lb | ecdfs > thres_ub)
}

# compute the empirical type I error rate
# @param alpha adjusted type I error rate per test
# @param N sample size
# @param K granularity size
# @param ncases number of discrete cases. ncases = 0 for continuous values
# @param M number of simulation trials
empirical_alpha <- function(alpha, N, K = N, ncases=0, M = 10000) {
  reject <- numeric(M)
  for (m in seq_len(M)) {
    x <- sample_unif(N, ncases)
    reject[m] <- test_uniformity(x, alpha = alpha, K = K)
  }
  mean(reject)
}

# compute the empirical type I error rate
# @param alpha adjusted type I error rate per test
# @param N sample size
# @param K granularity size
# @param ncases number of discrete cases. ncases = 0 for continuous values
# @param M number of simulation trials
# @param ties.method method to rank equal values in discrete data
empirical_alpha_chains <- function(alpha, niterations, nchains,
                                   K = niterations, ncases=0, M = 10000, ties.method = 'max') {
  reject <- numeric(M)
  ndraws <- niterations * nchains
  for (m in seq_len(M)) {
    x <- sample_unif(ndraws, ncases)
    x <- matrix(x, nrow = niterations, ncol = nchains)
    reject[m] <- test_uniformity_chains(x, alpha = alpha, K = K, ties.method)
  }
  mean(reject)
}
  
# quantiles in the ECDF to be investigated
get_partition <- function(K) {
  1:(K - 1) / K 
}

sample_unif <- function(ndraws,ncases) {
  if (ncases == 0) {
    x <- runif(ndraws, 0, 1)
  } else {
    x <- sample(1:ncases / ncases, ndraws, replace = TRUE)
  }
  x
}

# evaluate and scale ecdf at given quantiles 
eval_ecdf <- function(x, z) {
  colSums(outer(x, z, "<="))
}

# evaluate and scale ecdf matrix at given quantiles 
eval_ecdf_matrix <- function(x, z) {
  stopifnot(is.matrix(x))
  out <- matrix(nrow = length(z), ncol = ncol(x))
  for (j in seq_len(ncol(x))) {
    out[, j] <- colSums(outer(x[, j], z, "<=")) 
  }
  out
}

# Transform sample or a matrix of samples into normalized ranks.
u_scale <- function(x, ties.method = 'max') {
  stopifnot(ties.method %in% c('max', 'random'))
  S <- length(x)
  r_max <- rank(x, ties.method = 'max')
  if (ties.method == 'max') {
    r <- r_max
  } else if (ties.method == 'random') {
    r <- numeric(S)
    r_min <- rank(x, ties.method = 'min')
    for (idx in seq_len(S)) {
      r[idx] <- sample(r_min[idx]:r_max[idx], 1)
    }
  }
  # TODO: check which of the two versions of u should be preferred
  # u <- (r - 1 / 2) / S
  u <- r / S
  if (!is.null(dim(x))) {
    # output should have the input dimension
    u <- array(u, dim = dim(x), dimnames = dimnames(x))
  }
  u
}

# alpha percent of the trials are allowed to be rejected 
alpha_quantile <- function(gamma, alpha, tol = 0.001) {
  a <- unname(quantile(gamma, probs = alpha))
  a_tol <- unname(quantile(gamma, probs = alpha + tol))
  if (a == a_tol) {
    # using 'a' will lead to exceedance of alpha by at least tol
    unique_gamma <- sort(unique(gamma))
    if (unique_gamma[1] < a) {
      # choose the largest value smaller than 'a'
      smaller_a <- unique_gamma[unique_gamma < a]
      a <- smaller_a[length(smaller_a)]
    }
  }
  a
}

# adjust the alpha level of the envelope test using optimization
adjust_alpha_optimize <- function(alpha, N, K = N) {
  optimize(target_gamma, c(0, alpha), alpha = alpha, N = N, K = K)$minimum
}

target_gamma <- function(gamma, alpha, N, K) {
  z <- get_partition(K)
  z2 <- c(z, 1)
  z1 <- c(0, z)
  
  # pre-compute quantiles and use symmetry for increased efficiency
  x2_lower <- qbinom(gamma / 2, N, z2)
  x2_upper <- c(N - rev(x2_lower)[seq_len(K)[-1]], N)

  # compute the total probability of trajectories inside the envelope
  # initialize interior set and corresponding probabilities
  # known to be 0 and 1 for the starting value z1 = 0 
  x1 <- 0
  p_int <- 1
  for (i in seq_along(z1)) {
    tmp <- p_interior(
      p_int, x1 = x1, x2 = x2_lower[i]:x2_upper[i], 
      z1 = z1[i], z2 = z2[i], gamma = gamma, N = N
    )
    x1 <- tmp$x1
    p_int <- tmp$p_int
  }
  abs(1 - alpha - sum(p_int))
}

# probability of the ECDF being completely within the envelope as z2
p_interior <- function(p_int, x1, x2, z1, z2, gamma, N) {
  z_tilde <- (z2 - z1) / (1 - z1)
  
  # vectorize quantities for increased efficiency
  N_tilde <- rep(N - x1, each = length(x2))
  p_int <- rep(p_int, each = length(x2))
  x_diff <- outer(x2, x1, "-")
  p_x2_int <- p_int * dbinom(x_diff, N_tilde, z_tilde)
  
  # non-vectorized code below
  # N_tilde <- N - x1
  # p_x2_int <- matrix(NA, length(x2), length(x1))
  # for (j in seq_along(x1)) {
  #   # compute the interior for y2 = x2 - x1 which is binomial given x1
  #   p_x2_given_int <- dbinom(x2 - x1[j], N_tilde[j], z_tilde)
  #   p_x2_int[, j] <- p_int[j] * p_x2_given_int
  # }
  
  list(p_int = rowSums(p_x2_int), x1 = x2)
}

# adjust the alpha level of the envelope test using optimization
adjust_alpha_optimize_chains <- function(alpha, niterations, nchains,
                                         K = niterations) {
  if (nchains <= 1) {
    stop("At least 2 chains are required.")
  } 
  if (nchains == 2) {
    target <- target_gamma_chains_2
  } else {
    target <- target_gamma_chains_multi
  }
  optimize(
    target, c(0, alpha), 
    alpha = alpha, niterations = niterations, 
    nchains = nchains, K = K
  )$minimum
}

target_gamma_chains_2 <- function(gamma, alpha, niterations, nchains, K) {
  stopifnot(nchains == 2)
  
  # pre-compute quantiles
  z <- get_partition(K)
  m <- niterations
  n <- niterations * (nchains - 1)
  s <- floor(z * niterations * nchains)
  s2 <- c(s, niterations * nchains)
  s1 <- c(0, s)
  x2_lower <- qhyper(gamma / 2, m, n, s2)
  x2_upper <- qhyper(1 - gamma / 2, m, n, s2)
  
  # compute the total probability of trajectories inside the envelope
  # initialize interior set and corresponding probabilities
  # known to be 0 and 1 for the starting value z1 = 0 
  x1 <- 0
  p_int <- 1
  for (i in seq_along(s1)) {
    tmp <- p_interior_chains_2(
      p_int, x1 = x1, x2 = x2_lower[i]:x2_upper[i], 
      s1 = s1[i], s2 = s2[i], gamma = gamma, 
      niterations = niterations, nchains = nchains
    )
    x1 <- tmp$x1
    p_int <- tmp$p_int
  }
  abs(1 - alpha - sum(p_int))
}

# probability of the ECDF being completely within the envelope as z2
p_interior_chains_2 <- function(p_int, x1, x2, s1, s2, gamma, 
                              niterations, nchains) {
  stopifnot(nchains == 2)
  
  # vectorize quantities for increased efficiency
  s_tilde <- s2 - s1
  m_tilde <- rep(niterations - x1, each = length(x2))
  n_tilde <- rep(niterations - (s1 - x1), each = length(x2))
  p_int <- rep(p_int, each = length(x2))
  x_diff <- outer(x2, x1, "-")
  p_x2_int <- p_int * dhyper(x_diff, m_tilde, n_tilde, s_tilde)
  list(p_int = rowSums(p_x2_int), x1 = x2)
}

target_gamma_chains_multi <- function(gamma, alpha, niterations, nchains, K) {
  
  # pre-compute quantiles
  z <- get_partition(K)
  m <- niterations
  n <- niterations * (nchains - 1)
  s <- floor(z * niterations * nchains)
  s2 <- c(s, niterations * nchains)
  s1 <- c(0, s)
  x2_lower <- qhyper(gamma / 2, m, n, s2)
  x2_upper <- qhyper(1 - gamma / 2, m, n, s2)
  
  # compute the total probability of trajectories inside the envelope
  # initialize interior set and corresponding probabilities
  # known to be 0 and 1 for the starting value z1 = 0 
  x1 <- matrix(0, ncol = nchains)
  p_int <- 1
  for (i in seq_along(s1)) {
    tmp <- p_interior_chains_multi(
      p_int, x1 = x1, x2 = x2_lower[i]:x2_upper[i], 
      s1 = s1[i], s2 = s2[i], gamma = gamma, 
      niterations = niterations, nchains = nchains
    )
    x1 <- tmp$x1
    p_int <- tmp$p_int
  }
  abs(1 - alpha - sum(p_int))
}

# probability of the ECDF being completely within the envelope as z2
p_interior_chains_multi <- function(p_int, x1, x2, s1, s2, gamma, 
                                    niterations, nchains) {
  
  # vectorize quantities for increased efficiency
  s_tilde <- s2 - s1
  states = get_states(x2, s2, nchains)
  rep_x1 <- rep(seq_len(nrow(x1)), each = nrow(states))
  rep_states <- rep(seq_len(nrow(states)), nrow(x1))
  n_tilde <- (niterations - x1)[rep_x1, , drop = FALSE]
  p_int <- rep(p_int, each = nrow(states))
  x_diff <- states[rep_states, , drop = FALSE] - x1[rep_x1, , drop = FALSE]
  p_x2_int <- matrix(
    p_int * extraDistr::dmvhyper(x_diff, n_tilde, s_tilde), 
    nrow = nrow(states)
  )
  list(p_int = rowSums(p_x2_int), x1 = states)
}

# Get a matrix of transitions from x1 to a new state by adding s_tilde
# ranks and only accepting transitions, where all ranks are between
# x2_lower and x2_upper. 
get_states <- function(x2, s2, nchains) {
  # TODO: improve efficiency
  x2_nchains <- replicate(nchains, x2, simplify = FALSE)
  states <- do.call(expand.grid, args = x2_nchains)
  # Limit to states where ranks add up to s2
  as.matrix(states)[(rowSums(states) == s2), , drop = FALSE]
}

# Probability density of a multivariate hypergeometric distribution
# replaced by extraDistr::dmvhyper
# dmvhypergeom <- function(y, m, fill.na = 0) {
#   # Check, that the number of columns in y and m match.of rows in 
#   stopifnot(isTRUE(ncol(y) == ncol(m)))
#   p <- matrixStats::rowProds(choose(m,y)) / choose(rowSums(m), rowSums(y))
#   p[is.na(p)] = fill.na
#   p
# }

source("alpha_adjustment_functions.R")

# validate alpha adjustment (simple case)
alpha <- 0.05
N <- 1000
K <- 100
M <- 10000

system.time(gamma1 <- adjust_alpha_simulate(alpha, N = N, K = K, M = M))
gamma1
system.time(gamma2 <- adjust_alpha_optimize(alpha, N = N, K = K))
gamma2

# should roughly match original alpha
empirical_alpha(gamma1, N = N, K = K, M = M)
empirical_alpha(gamma2, N = N, K = K, M = M)


# example for power calculation
k <- 1.3
M <- 10000
reject <- numeric(M)
for (m in seq_len(M)) {
  x <- 1 - (1 - runif(N))^k
  reject[m] <- test_uniformity(x, alpha = gamma2, K = K)
}
mean(reject)


# alpha adjustment (discrete case)
alpha <- 0.05
N <- 1000
K <- 5
M <- 10000

system.time(
  gamma1 <- adjust_alpha_simulate(alpha, N = N, K = K, ncases=K, M = M)
)
gamma1
system.time(
  gamma2 <- adjust_alpha_optimize(alpha, N = N, K = K)
)
gamma2

# should roughly match original alpha
empirical_alpha(gamma1, N = N, K = K, ncases = K, M = M)
empirical_alpha(gamma2, N = N, K = K, ncases = K, M = M)


# validate alpha adjustment (multiple chains)
alpha <- 0.05
niterations <- 50
nchains <- 3
K <- 50
M <- 10000

system.time(
  gamma1 <- adjust_alpha_simulate_chains(
    alpha, niterations = niterations, 
    nchains = nchains, K = K, M = M
  )
)
gamma1

system.time(
  gamma2 <- adjust_alpha_optimize_chains(
    alpha, niterations = niterations, 
    nchains = nchains, K = K
  )
)
gamma2

# should roughly match original alpha
empirical_alpha_chains(
  gamma1, niterations = niterations, 
  nchains = nchains, K = K, M = M
)
empirical_alpha_chains(
  gamma2, niterations = niterations, 
  nchains = nchains, K = K, M = M
)

# validate alpha adjustment (multiple chains discrete case)
alpha <- 0.05
niterations <- 100
nchains <- 4
K <- 10
M <- 10000

system.time(
  gamma_max <- adjust_alpha_simulate_chains(
    alpha, niterations = niterations, 
    nchains = nchains, K = K, ncases = K,
    M = M, ties.method = 'max'
  )
)
gamma_max

system.time(
  gamma_rnd <- adjust_alpha_simulate_chains(
    alpha, niterations = niterations, 
    nchains = nchains, ncases = K,
    M = M, ties.method = 'random'
  )
)
gamma_rnd

# should roughly match original alpha
empirical_alpha_chains(
  gamma_max, niterations = niterations, 
  nchains = nchains, K = K, ncases = K,
  M = M, ties.method = 'max'
)

empirical_alpha_chains(
  gamma_rnd, niterations = niterations, 
  nchains = nchains, ncases = K,
  M = M, ties.method = 'random'
)

# example for power calculation
M <- 10000
reject <- numeric(M)
for (m in seq_len(M)) {
  x <- rnorm(niterations * nchains)
  x <- matrix(x, nrow = niterations, ncol = nchains)
  # chain 1 has a slightly different mean
  x[, 1] <- x[, 1] + 0.2
  reject[m] <- test_uniformity_chains(x, alpha = gamma, K = K)
}
mean(reject)

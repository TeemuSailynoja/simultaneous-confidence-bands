source("alpha_adjustment_functions.R")
library(future.apply)
library(future)
plan(multisession)

set.seed(112233)
# Test rejection rate of single sample test.

## Use random assignment to solve ties in ranks.
empirical_pit_discrete <- function(y, fun, ..., m = 999) {
  pit <- numeric(length(y))
  for (k in seq_along(y)) {
    pit[k] <- rank(
      c(y[k], fun(n = m, ...)),
      ties.method = "random"
    )[1] / (m + 1)
  }
  pit
}

n <- 250
lam <- 10
alpha <- adjust_alpha_optimize(.05, n)

mean(
  future_replicate(
    test_uniformity(
      empirical_pit_discrete(
        y = rpois(n = n, lambda = lam),
        fun = rpois,
        lambda = lam
      ),
      alpha = alpha),
    n = 1000
  )
)

# Test discrete multiple sample comparison.
## Use random assignment to solve ties in ranks.

## mask the 'u_scale' function used previously.
u_scale <- function(x, ties.method = 'random') {
  # stopifnot(ties.method %in% c('max', 'random'))
  S <- length(x)
  r <- rank(x, ties.method = ties.method)
  u <- r / S
  if (!is.null(dim(x))) {
    # output should have the input dimension
    u <- array(u, dim = dim(x), dimnames = dimnames(x))
  }
  u
}

l <- 2
n <- 250

alpha_2 <- adjust_alpha_optimize_chains(.05, n, l)

mean(
  future_replicate(
    test_uniformity_chains(
      matrix(rpois(n * l, lambda = lam), ncol = l),
      alpha = alpha_2,
      ties.method = "random"
    ),
    n = 1000
  )
)

l <- 4
alpha_4 <- 0.0005892573 # adjust_alpha_simulate_chains(.05, n, l)

mean(
  future_replicate(
    test_uniformity_chains(
      matrix(rpois(n * l, lambda = lam), ncol = l),
      alpha = alpha_4,
      ties.method = "random"
    ),
    n = 1000
  )
)

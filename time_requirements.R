times = sapply(1:20, function(n) {
  t = mean(replicate(10, {
    ptm <- proc.time()
    adjust_alpha_optimize(.05, n*50)
    t = proc.time() - ptm
    t[["elapsed"]]
  }))
  cat("N:", 50 * n, "t:", t, "\n")
  t
})

times_2chains = sapply(1:4, function(n) {
  t = mean(replicate(10, {
    ptm <- proc.time()
    adjust_alpha_optimize_chains(.05, n*50, nchains = 2)
    t = proc.time() - ptm
    t[["elapsed"]]
  }))
  cat("N:", 50 * n, "t:", t, "\n")
  t
})


times_4chains = sapply(1:4, function(n) {
  t = mean(replicate(10, {
    ptm <- proc.time()
    adjust_alpha_optimize_chains(.05, N=n*10, nchains = 4)
    t = proc.time() - ptm
    t[["elapsed"]]
  }))
  cat("N:", 10 * n, "t:", t, "\n")
  t
})


times_sim = sapply(1:10, function(n) {
  t = mean(replicate(10, {
    ptm <- proc.time()
    adjust_alpha_simulate(.05, n*100)
    t = proc.time() - ptm
    t[["elapsed"]]
  }))
  cat("N:", 50 * n, "t:", t, "\n")
  t
})

times_2chains_sim = sapply(1:10, function(n) {
  t = mean(replicate(10, {
    ptm <- proc.time()
    adjust_alpha_simulate_chains(.05, n*100, nchains = 2)
    t = proc.time() - ptm
    t[["elapsed"]]
  }))
  cat("N:", 50 * n, "t:", t, "\n")
  t
})


times_4chains_sim = sapply(1:10, function(n) {
  t = mean(replicate(10, {
    ptm <- proc.time()
    adjust_alpha_simulate_chains(.05, n*100, nchains = 4)
    t = proc.time() - ptm
    t[["elapsed"]]
  }))
  cat("N:", 50 * n, "t:", t, "\n")
  t
})
Ns = c(10, 20, 50, 70, 100, 200, 500, 700, 1000, 2000)
Ks = Ns
alpha = 0.95

# Test and plot rejection rates of the recommended K = N
gammas = sapply(Ns, function(n) adjust_alpha_optimize(1-alpha, n))

rejs = sapply(seq_along(Ns), function(idx) {
  cat(Ns[idx], ": ", sep="")
  ea = empirical_alpha(gammas[idx], N = Ns[idx])
  cat(ea, "\n", sep="")
  ea
})

# Save the results and adjusted gamma values.
write.csv(gammas, file ="data/gamma_vals.csv", row.names=Ns)
write.csv(rejs, file ="data/empirical_alpha_varying_chain_lengths.csv", row.names=Ns)

# Show results
print(rejs)
plot(Ns, rejs)
abline(h=.05)

# Run the same test for varying K. This takes time.

gammas_K = array(c(rep(0, 100)), dim=c(10,10,1), dimnames = list("N"=Ns, "K"=Ns, "gamma"=1))
# We already ran the diagonal above.
for (idx in seq_along(Ns)){gammas_K[idx,idx,] = gammas[idx]}

rejs_K = array(c(rep(0, 100)), dim=c(10,10,1), dimnames = list("N"=Ns, "K"=Ns, "rejr"=1))
# We already ran the diagonal above.
for (idx in seq_along(Ns)){rejs_K[idx,idx,] = rejs[idx]}

# Estimate the off-diagonal rejection rates.
for (n in Ns) {
  for (k in Ks) {
    if (n != k) {
      cat("N: ", n, " - K: ", k,": ", sep="")
      gammas_K[paste(n),paste(k), paste(1)] = adjust_alpha_optimize(0.05, n, k)
      rejs_K[paste(n),paste(k), paste(1)] = empirical_alpha(gammas_K[paste(n),paste(k), paste(1)], n, k, M=50000)
      print(rejs_K[paste(n),paste(k), paste(1)])
    }
  }
}
# Transform the results to long from
res.df = data.frame(rejs_K)
colnames(res.df) = Ks
res.df["N"] = Ns
res.df = reshape2::melt(res.df, id.vars = "N", variable.name="K", value.name="rejr")

# Plot the results
ggplot(res.df, mapping=aes(x = N, group = K, color=K, y = rejr)) +
  hline_at(0.05, linetype="dashed") +
  geom_point() +
  geom_line(size=.3, alpha=.5) +
  scale_color_discrete()

# Save the results.
write.csv(gammas_K, file ="data/gamma_vals_N_K.csv", row.names=Ns)
write.csv(rejs_K, file ="data/empirical_alpha_varying_N_and_K.csv", row.names=Ns)

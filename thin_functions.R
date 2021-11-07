
thin <- function(x, by) {
  if (!is.integer(by)) {
    by <- max(1, as.integer(ceiling(by)))
  }
  if (any(is.null(dim(x)), length(dim(x)) == 1)) {
    x <- array(x, dim = c(length(x), 1))
  }
  
  S <- nrow(x)
  stopifnot(by < S)
  
  thinned = x[seq(1, S, by), ]
  
  if (length(thinned) < 100) {
    warning(
      paste(
        "Only ",
        length(thinned),
        "samples left after thinning! Chain length = ",
        S,
        ", 'by' = ",
        by,
        sep = ""
      )
    )
  }
  
  thinned
}


by_quantiles <-
  function(x,
           n_quantiles = 20,
           by_fun = geyers_tau,
           ...) {
    by <- 1
    for (x_ in quantile(x, 1:(n_quantiles - 1) / n_quantiles)) {
      by <- max(by, by_fun((x <= x_), ...))
    }
    by
  }

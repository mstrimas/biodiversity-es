calculate_targets <- function(x, minimum_range_size = 1000,
                              proportion_at_minimum_range_size = 1,
                              maximum_range_size = 250000,
                              proportion_at_maximum_range_size = 0.1,
                              cap_range_size = 10000000,
                              cap_target_amount = 1000000) {
  # initialization
  out <- rep(NA_real_, length(x))
  pos <- which(is.finite(x))
  # set targets below threshold
  out[pos] <- stats::approx(
    x = log10(c(minimum_range_size, maximum_range_size)),
    y = c(proportion_at_minimum_range_size,
          proportion_at_maximum_range_size),
    xout = log10(x[pos]), method = "linear", rule = 2)$y
  # set caps
  cap_pos <- which(x > cap_range_size)
  if (length(cap_pos) > 0) {
    out[cap_pos] <- cap_target_amount / x[cap_pos]
  }
  # return values
  return(out)
}
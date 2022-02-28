# source("code.R")

# load packages
library(readr)
library(dplyr)
library(tibble)
library(prioritizr)
library(units)
library(Matrix)

# load data
rij_10km_data <- readRDS("data/rij-matrix_10km.rds")

# define variables
resolution <- 10

# define functions
## https://github.com/mstrimas/biodiversity-es/blob/master/R/calculate-targets.R
calculate_targets <- function(x, minimum_range_size = 1000,
                              proportion_at_minimum_range_size = 1,
                              maximum_range_size = 250000,
                              proportion_at_maximum_range_size = 0.1,
                              cap_range_size = 10000000,
                              cap_target_amount = 1000000) {
  ## initialization
  out <- rep(NA_real_, length(x))
  pos <- which(is.finite(x))
  ## set targets below threshold
  out[pos] <- stats::approx(
    x = log10(c(minimum_range_size, maximum_range_size)),
    y = c(proportion_at_minimum_range_size,
          proportion_at_maximum_range_size),
    xout = log10(x[pos]), method = "linear", rule = 2)$y
  ## set caps
  cap_pos <- which(x > cap_range_size)
  if (length(cap_pos) > 0) {
    out[cap_pos] <- cap_target_amount / x[cap_pos]
  }
  ## return values
  return(out)
}

# re-create minimal version of feature data based only on rij data
feats_data <- tibble(
  ## create id column
  id = seq_len(nrow(rij_10km_data)),
  ## extract feature names
  name = rownames(rij_10km_data),
  ## manually specify which features are ES vs. species
  type = c(rep("es", 12), rep("aoh", nrow(rij_10km_data) - 12))
)

# re-scale ES data
## calculate scaling factors needed to multiply data so max value == 100
es_scaling <- 100 / apply(
  rij_10km_data[which(feats_data$type == "es"), , drop = FALSE],
  MARGIN = 1,
  FUN = max,
  na.rm = TRUE
)

## apply scaling factors
for (i in seq_along(es_scaling)) {
  idx <- which(feats_data$type == "es")[i]
  rij_10km_data[idx, ] <- rij_10km_data[idx, ] * es_scaling[i]
}
rij_10km_data <- drop0(rij_10km_data) # ensure sparsity

## store scaling factors for later
feats_data$es_scaling <- c(es_scaling, rep(1, sum(feats_data$type != "es")))

# calculate total amount of each feature, so we can calculate the targets later
## note that we calculate the row-sums in chunks to avoid memory issues
feats_data$total <- unlist(
  plyr::llply(
    distribute_load(nrow(rij_10km_data), 5),
    function(i) rowSums(rij_10km_data[i, , drop = FALSE], na.rm = TRUE),
    .progress = "text"
  ),
  recursive = TRUE
)

# calculate the targets
## here we will just set the targets for ES as 10% because we are lazy
feats_data <-
  feats_data %>%
  ## calculate targets with Rodrigues approach
  mutate(target_rel = (total / 100) * ((resolution * 1000)^2)) %>%
  mutate(
    target_rel = as.numeric(set_units(set_units(target_rel, m^2), km^2))
  ) %>%
  mutate(target_rel = calculate_targets(target_rel)) %>%
  ## manually over-write targets for ES with 10%
  mutate(target_rel = if_else(type == "es", 0.1, target_rel)) %>%
  ## calculate targets in absolute values
  mutate(target_abs = target_rel * total)

# generate min set solution
## formulate problem
min_set_prob <-
  problem(
    rep(1, ncol(rij_10km_data)),
    feats_data,
    rij_matrix = rij_10km_data,
    run_checks = FALSE
  ) %>%
  add_min_set_objective() %>%
  add_absolute_targets("target_abs") %>%
  add_binary_decisions() %>%
  add_gurobi_solver(
    gap = 0.05, threads = 5, verbose = TRUE
  )

## preview problem
print(min_set_prob)

## solve problem
min_set_sol <- solve(min_set_prob, force = TRUE)

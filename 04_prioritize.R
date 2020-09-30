# run prioritization scenarios

library(Matrix)
library(prioritizr)
library(gurobi)
library(tidyverse)

# define target function
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


# parameters
n_cores <- 10
resolution <- 10
# prioritization scenarios
scenarios <- expand_grid(es = c(0.3, 0.5, 0.9),
                         budget = c(0.3, 0.5, 1)) %>% 
  mutate(scenario = str_glue("es-{100 * es}_budget-{100 * budget}")) %>% 
  select(scenario, es, budget)


# input data ----

# features
features <- str_glue("features_{resolution}km.csv") %>% 
  file.path(DATA_DIR, .) %>% 
  read_csv()
# planning units
pu <- str_glue("pu_{resolution}km.rds") %>% 
  file.path(DATA_DIR, "pu", .) %>% 
  read_rds()
# raster template
r <- str_glue("pu_eck4_{resolution}km.tif") %>% 
  file.path(DATA_DIR, "pu", .) %>% 
  raster() %>% 
  raster()
# rij matrices
rij <- str_glue("rij-matrix_{resolution}km.rds") %>% 
  file.path(DATA_DIR, .) %>% 
  read_rds()
# for testing, remove most of the features
# rij <- rij[1:1000, ]

# total across all planning units
# for biodiversity features, values are % of cell occupied
# aoh = area of planning unit times total representation / 100
features <- rowSums(rij) %>% 
  enframe(value = "total") %>% 
  inner_join(features, ., by = "name")
# set biodiversity targets
max_total <- resolution^2
features <- features %>% 
  mutate(aoh = ifelse(type != "es", total / max_total, NA_real_)) %>% 
  # set target based on aoh
  mutate(prop = calculate_targets(aoh))


# prioritize ---- 

for (i in seq.int(nrow(scenarios))) {
  # ecosystem service targets
  features$prop <- ifelse(features$type == "es", scenarios$es[i], features$prop)
  
  # construct problem, use cost = 1
  p <- problem(rep(1, ncol(rij)),
               features, 
               rij_matrix = rij,
               run_checks = FALSE) %>%
    add_relative_targets("prop") %>%
    add_binary_decisions() %>% 
    add_gurobi_solver(gap = 0.01, threads = n_cores)
  # add objective function based on scenario
  budget <- scenarios$budget[i]
  if (budget == 1) {
    p <- add_min_set_objective(p)
  } else {
    p <- add_min_shortfall_objective(p, budget = budget * ncol(rij))
  }
  
  # solve
  # need to turn off checks due to large range of features
  s <- solve(p, run_checks = FALSE)
  
  # save solution stats
  scenarios[["n_selected"]][i] <- sum(s)
  scenarios[["runtime"]][i] <- attr(s, "runtime")
  scenarios[["solution"]][i] <- str_glue("solution_{scenarios$scenario[i]}",
                                         "_{resolution}km.rds")
  scenarios[["raster"]][i] <- str_glue("solution_{scenarios$scenario[i]}",
                                         "_{resolution}km.tif")
  
  # solution object
  sol <- pu
  sol$selected <- as.logical(s)
  scenarios[["solution"]][i] %>% 
    file.path(OUTPUT_DIR, .) %>% 
    write_rds(sol, ., compress = "gz")
  # raster
  r_sol <- r
  r_sol[sol$cell_id] <- as.integer(sol$selected)
  r_sol <- scenarios[["raster"]][i] %>%
    file.path(OUTPUT_DIR, .) %>% 
    writeRaster(r_sol, ., 
                datatype = "INT1U",
                options = c("COMPRESS=LZW", "TILED=YES"),
                overwrite = TRUE)
  
  # clean up
  rm(s, sol, r_sol)
  co <- capture.output(gc())
  removeTmpFiles(h = 0)
}
str_glue("scenarios_{resolution}km.tif") %>%
  file.path(OUTPUT_DIR, .) %>% 
  write_csv(scenarios, .)
# run prioritization scenarios

library(fs)
library(glue)
library(arrow)
library(Matrix)
library(prioritizr)
library(gurobi)
library(tidyverse)
source("R/calculate-targets.R")
source("R/multi-objective-prioritization.R")

data_dir <- "data/"
output_dir <- "output_4/"

# set all values below this to 0
clamp_value <- 1

# parameters
n_cores <- 30

# 2, 3, 5, 10
resolution <- 10
res_lbl <- paste0(resolution, "km")
res_output_dir <- path(output_dir, res_lbl)
dir_create(res_output_dir)

# prioritization scenarios
# prioritization scenarios
scenarios <- expand_grid(biod = c(0),
                         es = seq(0, 1, 0.05),
                         budget = 1,
                         resolution = resolution,
                         n_selected = NA,
                         runtime = NA,
                         solution = NA,
                         raster = NA) %>% 
  mutate(scenario = str_glue("es-{100 * es}_biod-{biod}_budget-{100 * budget}_resolution-{resolution}km")) %>% 
  select(scenario, es, biod, budget, resolution, n_selected, runtime, solution, raster) %>%
  filter(es > 0 | biod > 0, !(biod == 0 & budget < 1))




# input data ----

# features
features <- path(data_dir, "features.csv") %>% 
  read_csv() %>% 
  filter(res_km == resolution)
# planning units
pu <- glue("pu_{res_lbl}.parquet") %>% 
  path(data_dir, "pu", res_lbl, .) %>% 
  read_parquet()
# raster template
r <- glue("pu_eck4_{res_lbl}.tif") %>% 
  path(data_dir, "pu", res_lbl, .) %>% 
  raster() %>% 
  raster()
# rij matrices
rij <- glue("rij-matrix_{res_lbl}.rds") %>% 
  path(data_dir, .) %>% 
  read_rds()
# protected areas
pa <- glue("protected-areas_{res_lbl}.parquet") %>% 
  path(data_dir, "pu", res_lbl, .) %>% 
  read_parquet() %>%
  mutate(lock = ifelse(coverage_fraction >= 0.5, TRUE, FALSE))

pu <- pu %>% left_join(pa %>% select(pu_id, lock), by = c("id" = "pu_id")) %>%
  mutate(lock = ifelse(is.na(lock), FALSE, lock))

# clamp NCP layers
# only clamp values for layers that have a huge range.
es_features <- features$id[features$type == "es"]
rij_sub <- rij[es_features, ]
for (ii in seq_len(nrow(rij_sub))) {
  if (max(rij_sub[ii, ]) > 1000000) {
    rij_sub[ii, ] <- ifelse(rij_sub[ii, ] < clamp_value, 0, rij_sub[ii, ])
  }
}

rij <- rij_sub
features <- features[features$type == "es",]
# rij[es_features, ] <- rij_sub

# for testing, remove most of the features
# rij <- rij[1:1000, ]

# total across all planning units
# for biodiversity features, values are % of cell occupied
# aoh = area of planning unit times total representation / 100
# features <- rowSums(rij) %>% 
#   enframe(value = "total") %>% 
#   inner_join(features, ., by = "name")
# # set biodiversity targets
# max_total <- resolution^2
# features <- features %>% 
#   mutate(aoh = ifelse(type != "es", total, NA_real_)) %>% 
#   # set target based on aoh
#   mutate(prop0 = calculate_targets(aoh))


# prioritize ---- 

for (i in seq_len(nrow(scenarios))) {
  # ecosystem service targets
  features$prop <- ifelse(features$type == "es", scenarios$es[i],
                          features$prop0 * scenarios$biod[i])
  # 
  # parts <- c(1:10, 11:nrow(features))
  parts <- c(1:nrow(features))
  # construct problem, use cost = 1
  p <- problem(rep(1, ncol(rij)),
               features[parts,], 
               rij_matrix = rij[parts,],
               run_checks = FALSE) %>%
    add_relative_targets("prop") %>%
    # add_binary_decisions() %>% 
    add_proportion_decisions() %>% #FOR TESTING ONLY
    add_locked_in_constraints(pu$lock) %>%
    add_gurobi_solver(gap = 0.05, threads = n_cores)
  # add objective function based on scenario
  
  gc()
  
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
  scenarios[["solution"]][i] <- glue("solution_{scenarios$scenario[i]}",
                                     "_{res_lbl}.rds")
  scenarios[["raster"]][i] <- glue("solution_{scenarios$scenario[i]}",
                                   "_{res_lbl}.tif")
  
  # solution object
  sol <- pu
  sol$selected <- as.numeric(s)
  scenarios[["solution"]][i] %>% 
    path(res_output_dir, .) %>% 
    write_rds(s, ., compress = "gz")
  # raster
  r_sol <- r
  r_sol[sol$cell_id] <- as.numeric(sol$selected)
  r_sol <- scenarios[["raster"]][i] %>%
    file.path(res_output_dir, .) %>% 
    writeRaster(r_sol, ., 
                # datatype = "FLT4S",
                options = c("COMPRESS=LZW", "TILED=YES"),
                overwrite = TRUE)
  
  # clean up
  rm(p, s, sol, r_sol)
  co <- capture.output(gc())
  removeTmpFiles(h = 0)
}
glue("scenarios_{res_lbl}.csv") %>%
  path(res_output_dir, .) %>% 
  write_csv(scenarios, .)

# prepare sparseMatrix format representation matrices for prioritizr

library(Matrix)
library(fs)
library(glue)
library(tidyverse)
library(arrow)
library(tictoc)

data_dir <- "data/"
resolutions <- 10
features <- read_csv(path(data_dir, "features.csv"))

# feature list
for (this_res in resolutions) {
  res_lbl <- paste0(this_res, "km")
  f <- filter(features, res_km == this_res)
  
  # stack layer-specific rij data frames into sparse matrices
  # number of columns, i.e. planning units
  n_features <- nrow(f)
  n_pu <- glue("pu_{res_lbl}.parquet") %>% 
    path(data_dir, "pu", res_lbl,  .) %>% 
    read_parquet() %>% 
    nrow()
  # load representation tables
  stopifnot(all(file.exists(f$rij)))
  tic()
  rij <- NULL
  for (i in seq_len(nrow(f))) {
    rij <- mutate(read_rds(f$rij[i]), species = f$id[i]) %>% 
      bind_rows(rij, .)
  }
  toc()
  
  # convert to sparse matrix
  rij <- sparseMatrix(i = rij$species, j = rij$pu, x = rij$amount,
                      index1 = TRUE, use.last.ij = FALSE,
                      dims = c(n_features, n_pu),
                      dimnames = list(f$layer, NULL))
  glue("rij-matrix_{res_lbl}.rds") %>% 
    path(data_dir, .) %>% 
    write_rds(rij, ., compress = "gz")
}
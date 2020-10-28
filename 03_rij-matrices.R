# prepare sparseMatrix format representation matrices for prioritizr

library(Matrix)
library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel(detectCores() - 2)

# resolution in km
resolution <- 5

# feature list
f <- str_glue("features_{resolution}km.csv") %>% 
  file.path(DATA_DIR, .) %>% 
  read_csv()

# stack layer-specific rij data frames into sparse matrices
# number of columns, i.e. planning units
n_features <- nrow(f)
n_pu <- str_glue("pu_{resolution}km.rds") %>% 
  file.path(DATA_DIR, "pu", .) %>% 
  read_rds() %>% 
  nrow()
# load representation tables
stopifnot(all(file.exists(f$rij)))
rij <- foreach (i = seq.int(nrow(f))) %dopar% {
  mutate(read_rds(f$rij[i]), species = f$id[i])
}
rij <- bind_rows(rij)

# convert to sparse matrix
rij <- sparseMatrix(i = rij$species, j = rij$pu, x = rij$amount,
                    index1 = TRUE, use.last.ij = FALSE,
                    dims = c(n_features, n_pu),
                    dimnames = list(f$name, NULL))
str_glue("rij-matrix_{resolution}km.rds") %>% 
  file.path(DATA_DIR, .) %>% 
  write_rds(rij, ., compress = "gz")
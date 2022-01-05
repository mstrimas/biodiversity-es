library(Matrix)
library(prioritizr)
library(gurobi)
library(tidyverse)
source("R/calculate-targets.R")
source("R/multi-objective-prioritization.R")

resolution <- 10

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

# raster template
pur <- str_glue("pu_eck4_{resolution}km.tif") %>% 
  file.path("data/pu/", .) %>% 
  raster()

LCNA <- raster("data/paper1/A_90_md5_79f5e0d5d5029d90e8f10d5932da93ff.tif") %>% 
  extend(pur) %>% 
  crop(pur)
GCNA <- raster("data/paper1/C_90_md5_bdf604015a7b1c7c78845ad716d568ef.tif") %>% 
  extend(pur) %>% 
  crop(pur)

LCNA_rij <- rij_matrix(pur, LCNA)

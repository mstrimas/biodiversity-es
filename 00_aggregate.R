library(raster)
library(fs)
library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel(cores = 12)
source("R/gdal-aggregate.R")

# set all values below this to 0
clamp_value <- 0.0001

# template at 2 km
r_template <- dir("data/tifs/birds/", "10km.tif$", full.names = TRUE) %>% 
  pluck(1) %>% 
  raster() %>% 
  raster()

# resample at 2 km, then aggregate to 10km, then clamp
f_raw <- list.files("data/raw/", "tif$", full.names = TRUE)
processed_dir <- file.path("data", "tifs", "es")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
assets <- foreach (f = f_raw, .combine = c) %dopar% {
  f_out <- basename(f) %>% 
    str_remove("_(PERAREA|WARP|clamped|md5).*") %>% 
    str_remove("^realized_") %>% 
    str_to_lower() %>% 
    str_replace_all("_", "-") %>% 
    paste0("_eck4_10km.tif") %>% 
    file.path(processed_dir, .)
  r <- gdal_aggregate(f, fact = 5, e = extent(r_template)) %>% 
    reclassify(rcl = c(-Inf, clamp_value, NA_real_),
               filename = f_out, overwrite = TRUE, 
               options = c("COMPRESS=LZW", "TILED=YES"))
}

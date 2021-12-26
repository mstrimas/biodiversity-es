# aggregate the 2-km ncp data to various resolutions for analysis

library(terra)
library(fs)
library(glue)
library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel(cores = 12)

data_dir <- "data"

temp_dir <- tempdir() %>% 
  path("aggregate")
dir_create(temp_dir)

# set all values below this to 0
clamp_value <- 0.0001
# resolutions for analysis
base_res <- 2
resolutions <- c(2, 3, 5, 10)
# eckert iv crs
crs <- path(data_dir, "eck4.wkt") %>% 
  read_lines() %>% 
  trimws() %>% 
  paste(collapse = "")
  


# disaggregate to 1km ----

# some resolutions are not direct aggregations of the raw 2km
# so starting from 1km avoids having to resample
f_raw <- dir_ls(path(data_dir, "raw"), glob = "*.tif")
f_1km <- foreach (f = f_raw, .combine = c) %dopar% {
  f_out <- path(temp_dir, basename(f))
  r <- rast(f)
  crs(r) <- crs
  r <- disagg(r, fact = 2, filename = f_out, overwrite = TRUE,
              datatype = "FLT8S")
  stopifnot(all(res(r) == 1000))
  return(f_out)
}


# aggregate to analysis resolutions ----

for (this_res in resolutions) {
  res_dir <- path(data_dir, "tifs", paste0(this_res, "km"), "es")
  dir_create(res_dir)
  
  layers <- foreach (f = f_1km) %dopar% {
    layer_name <- basename(f) %>% 
      str_remove("_(PERAREA|WARP|clamped|md5).*") %>% 
      str_remove("^realized_") %>% 
      str_to_lower() %>% 
      str_replace_all("_", "-")
    f_out <- glue("{layer_name}_eck4_{this_res}km.tif") %>% 
      path(res_dir, .)
    
    if (this_res == base_res) {
      # use the original rasters where available
      r <- path(data_dir, "raw", basename(f)) %>% 
        rast()
      crs(r) <- crs
    } else {
      # aggregate to coarser resolution
      r <- rast(f) %>% 
        aggregate(fact = this_res, fun = "mean", na.rm = TRUE)
    }
      
    # clamp out small values
    r <- clamp(r, lower = clamp_value, values = FALSE,
               filename = f_out, overwrite = TRUE,
               datatype = "FLT4S",
               gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
    return(f_out)
  }
  # ensure that all layers are aligned
  ncp_stack <- rast(unlist(layers))
}

dir_delete(temp_dir)

# prepare tabular input data for prioritizr

library(terra)
library(sf)
library(fs)
library(glue)
library(tidyverse)
library(arrow)
library(foreach)
library(doParallel)
registerDoParallel(14)
source("R/get-raster-values.R")

data_dir <- "data/"
resolutions <- c(2, 3, 5, 10)

no_rep <- NULL
features <- NULL
for (this_res in resolutions) {
  res_lbl <- paste0(this_res, "km")
  # directories
  tif_dir <- path(data_dir, "tifs", res_lbl)
  pu_dir <- path(data_dir, "pu", res_lbl)
  feature_dir <- path(data_dir, "features", res_lbl)
  dir_create(feature_dir)
  
  # raster template
  pu <- glue("pu_eck4_{this_res}km.tif") %>% 
    path(pu_dir, .) %>% 
    rast()
  # planning unit raster cell numbers
  pu_lookup <- glue("pu_{this_res}km.parquet") %>% 
    path(pu_dir, .) %>% 
    read_parquet()
  
  # define features
  f <- dir_ls(tif_dir, type = "directory") %>% 
    map(dir_ls, regex = "*.tif", recurse = TRUE) %>% 
    unlist() %>% 
    as.character() %>% 
    tibble(tif = .) %>% 
    mutate(res_km = this_res,
           type = basename(dirname(tif)),
           layer = basename(tif) %>% 
             str_remove("_eck4.*") %>% 
             str_remove("_aoh_[0-9]+km.tif") %>% 
             str_to_lower() %>% 
             str_replace_all("(_|\\s)+", "-"),
           rij = glue("rij_{this_res}km_{layer}.rds") %>% 
             file.path(feature_dir, .)) %>% 
    select(type, layer, res_km, tif, rij)
  # check that there are no duplicate rij file names
  stopifnot(anyDuplicated(f$rij) == 0)
  
  # for testing, just run 10 of each type
  # f <- ungroup(slice_sample(group_by(f, type), n = 10))
  
  # build representation matrices for each feature
  f$has_values <- foreach (i = seq.int(nrow(f)), .combine = c) %dopar% {
    x <- f[i, ]
    rkm <- paste0(x$res_km, "km")
    # align raster with planning unit raster
    r <- rast(x$tif) %>% 
      extend(pu) %>% 
      crop(pu)
    if (!compareGeom(r, pu)) {
      return(NA)
    }
    # extract values from valid planning units
    r_cells <- get_raster_values(r, cells = pu_lookup$cell_id)
    if (nrow(r_cells) == 0) {
      return(FALSE)
    }
    # assign a contiguous planning unit id instead of the raster cell id
    r_cells$pu <- pu_lookup$id[match(r_cells$cell_id, pu_lookup$cell_id)]
    r_cells <- r_cells[, c("pu", "value")]
    names(r_cells) <- c("pu", "amount")
    write_rds(r_cells, x$rij, compress = "gz")
    rm(r_cells)
    co <- capture.output(gc())
    return(TRUE)
  }
  
  # errors resulting from mis-aligned rasters
  if (any(is.na(f$has_values))) {
    stop("These features have tifs not aligning with planning units:\n\t",
         paste(f$tif[is.na(f$has_values)], collapse = "\n"))
  }
  # features with zero representation in planning units
  if (any(is.na(f$has_values))) {
    message("These features are not represented in any planning units:\n\t",
            paste(f$tif[!f$has_values], collapse = "\n\t"))
  }
  no_rep <- f %>% 
    filter(!has_values) %>% 
    bind_rows(no_rep, .)
  
  # assign features ids accounting for layers that have no values associated
  features <- f %>% 
    filter(has_values) %>% 
    # put ecosystem services first
    arrange(res_km, fct_relevel(type, "es"), layer) %>% 
    mutate(id = row_number(),
           name = paste(type, layer, sep = "_")) %>% 
    select(id, name, type, layer, res_km, tif, rij) %>% 
    bind_rows(features, .)
}
write_csv(no_rep, path(data_dir, "no-representation.csv"), na = "")
write_csv(features, path(data_dir, "features.csv"), na = "")
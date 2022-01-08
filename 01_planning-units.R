# identify the set of planning units to use in the prioritization
# only include rasters cells with at least one non-zero feature

library(terra)
library(fasterize)
library(sf)
library(fs)
library(glue)
library(tidyverse)
library(arrow)
library(foreach)
library(doParallel)
registerDoParallel(cores = 12)
source("R/get-raster-values.R")

data_dir <- "data"
resolutions <- c(2, 3, 5, 10)

# directories
tif_dir <- path(data_dir, "tifs")
pu_dir <- path(data_dir, "pu")
dir_create(pu_dir)

# land
crs <- path(data_dir, "eck4.wkt") %>% 
  read_lines() %>% 
  trimws() %>% 
  paste(collapse = "") %>% 
  st_crs()
land <- path(data_dir, "esri-countries.gpkg") %>% 
  read_sf()  %>% 
  filter(NAME != "Antarctica") %>% 
  st_transform(crs = crs) %>% 
  st_buffer(dist = 10000) %>% 
  st_cast("MULTIPOLYGON") %>% 
  transmute(land = 1) %>% 
  vect()

# identify planning units with non-NA values for each feature type
for (this_res in resolutions) {
  res_lbl <- paste0(this_res, "km")
  tif_res_dir <- path(tif_dir, res_lbl)
  pu_res_dir <- path(pu_dir, res_lbl)
  dir_create(pu_res_dir)
  
  feature_types <- dir_ls(tif_res_dir, type = "directory") %>% 
    basename()
  
  # for each feature type make rasters with number of non-zero features
  for (ft in feature_types) {
    tifs <- path(tif_res_dir, ft) %>% 
      dir_ls(glob = "*.tif")
    f_counts <- glue("n-features_{ft}_eck4_{this_res}km.tif") %>% 
      path(tif_res_dir, .)
    temp_dir <- tempdir()
    if (length(tifs) >= 1000) {
      n_groups <- ceiling(length(tifs) / 1000)
      tif_groups <- split(tifs, cut(seq_along(tifs), n_groups, labels = FALSE))
      tifs <- glue("{ft}_group-{seq_len(n_groups)}.tif") %>% 
        path(temp_dir, .)
      s <- foreach (g = seq_len(n_groups)) %dopar% {
        cmd <- glue("gdal-summarize.py -q -f count -o {tifs[g]} ",
                    "{paste(tif_groups[[g]], collapse = ' ')} ",
                    "--overwrite")
        system(cmd)
      }
      delete_tifs <- TRUE
      summary_stat <- "sum"
    } else {
      delete_tifs <- FALSE
      summary_stat <- "count"
    }
    cmd <- glue("gdal-summarize.py -q -f {summary_stat} -o {f_counts} ",
                "--co 'COMPRESS=DEFLATE' {paste(tifs, collapse = ' ')} ",
                "--overwrite")
    system(cmd)
    
    if (delete_tifs) {
      file_delete(tifs)
    }
  }
  
  # total count across all feature types
  f_counts <- dir_ls(path(tif_dir, res_lbl), 
                     regexp = "n-features_[a-z]+_eck4", 
                     recurse = FALSE)
  f_total <- glue("n-features_eck4_{this_res}km.tif") %>% 
    path(tif_res_dir, .)
  cmd <- glue("gdal-summarize.py -q -f sum -o {f_total} --overwrite ",
              "--co 'COMPRESS=DEFLATE' {paste(f_counts, collapse = ' ')}")
  system(cmd)
  r_counts <- rast(f_total)
  global(r_counts, range, na.rm = TRUE)
  
  # create mask/template based on cells that have data
  # only consider cells within 1 grid cell of of land, exclude antarctica
  # land mask
  f_lm <- glue("land-mask_eck4_{this_res}km.tif") %>% 
    path(pu_res_dir, .)
  if (!file_exists(f_lm)) {
    land_mask <- rasterize(land, r_counts, field = 1) %>% 
      writeRaster(filename = f_lm, 
                  datatype = "INT1U", overwrite = T,
                  gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  }
  land_mask <- rast(f_lm)
  # planning unit mask
  f_mask <- glue("pu-mask_eck4_{this_res}km.tif") %>% 
    path(pu_res_dir, .)
  r_mask <- r_counts %>% 
    # remove zeros
    subst(0, NA) %>% 
    mask(land_mask, 
         filename = f_mask, 
         datatype = "INT4U", 
         gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  global(r_mask, range, na.rm = TRUE)
  
  # create raster will cells equal to planning unit id
  r <- rast(r_mask)
  # assign cell id
  stopifnot(ncell(r) < 2^31)
  values(r) <- seq_len(ncell(r))
  # mask to non-zero asset cells
  f <- glue("pu_eck4_{this_res}km.tif") %>% 
    path(pu_res_dir, .)
  r <- mask(r, r_mask,
            filename = f, overwrite = TRUE, 
            datatype = "INT4U",
            gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  
  # save the non-masked cell ids for each template
  v <- get_raster_values(r) %>% 
    sort()
  stopifnot(anyDuplicated(v) == 0)
  v <- tibble(id = as.integer(seq_along(v)),
              cell_id = as.integer(v))
  glue("pu_{this_res}km.parquet") %>% 
    path(pu_res_dir, .) %>% 
    write_parquet(v, ., compression = "gzip")
  rm(v)
  co <- capture.output(gc())
}


# protected areas
pa <- read_sf("data/protected-areas/global-2021-10-30.shp") %>% 
  st_transform(crs = crs) %>% 
  st_geometry() %>% 
  st_combine()
for (this_res in resolutions) {
  res_lbl <- paste0(this_res, "km")
  pu_res_dir <- path(pu_dir, res_lbl)
  
  r <- glue("pu-mask_eck4_{this_res}km.tif") %>% 
    path(pu_res_dir, .) %>% 
    rast()
  
  pu <- glue("pu_{this_res}km.parquet") %>% 
    path(pu_res_dir, .) %>% 
    read_parquet()
  r_pa <- exact_extract(r, pa, include_cell = TRUE)[[1]] %>% 
    filter(!is.na(value)) %>% 
    select(cell_id = cell, coverage_fraction) %>% 
    inner_join(pu, by = "cell_id") %>% 
    select(pu_id = id, coverage_fraction) %>% 
    arrange(pu_id)
  
  glue("protected-areas_{this_res}km.parquet") %>% 
    path(pu_res_dir, .) %>% 
    write_parquet(r_pa, ., compression = "gzip")
  rm(r_pa)
}

library(raster)
library(fasterize)
library(sf)
library(fs)
library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel(cores = 12)
source("R/get-raster-values.R")

# directories
pu_dir <- file.path(DATA_DIR, "pu")
dir.create(pu_dir, showWarnings = FALSE)

# target resolutions in km
res <- 10

# full set of planning units
tif_dir <- "data/tifs"
feature_types <- dir_ls(tif_dir, type = "directory") %>% 
  basename()
for (ft in feature_types) {
  tifs <- dir(path(tif_dir, ft), "10km.tif$", full.names = TRUE)
  f_counts <- str_glue("n-features_{ft}_eck4_{res}km.tif") %>% 
    file.path(tif_dir, .)
  if (file.exists(f_counts)) {
    next()
  }
  if (length(tifs) >= 1000) {
    n_groups <- ceiling(length(tifs) / 1000)
    tif_groups <- split(tifs, cut(seq_along(tifs), n_groups, labels = FALSE))
    td <- tempdir()
    tifs <- str_glue("{ft}_group-{seq_len(n_groups)}.tif") %>% 
      path(td, .)
    s <- foreach (g = seq_len(n_groups)) %dopar% {
      cmd <- str_glue("gdal-summarize.py -q -f count -o {tifs[g]} ",
                      "{paste(tif_groups[[g]], collapse = ' ')}")
      system(cmd)
    }
    summary_stat <- "sum"
  } else {
    summary_stat <- "count"
  }
  cmd <- str_glue("gdal-summarize.py -q -f {summary_stat} -o {f_counts} ",
                  "--co 'COMPRESS=DEFLATE' {paste(tifs, collapse = ' ')}")
  system(cmd)
  unlink(td)
}
f_counts <- dir(tif_dir, "n-features_[a-z]+_eck4", full.names = TRUE)
f_total <- str_glue("n-features_eck4_{res}km.tif") %>% 
  file.path(tif_dir, .)
cmd <- str_glue("gdal-summarize.py -q -f sum -o {f_total} --overwrite ",
                "--co 'COMPRESS=DEFLATE' {paste(f_counts, collapse = ' ')}")
system(cmd)
r_counts <- raster(f_total)
cellStats(r_counts, range)

# create mask/template based on cells that have data
# only consider cells within 10 km of land, exclude antarctica
# land mask
f_lm <- str_glue("land-mask_eck4_{res}km.tif") %>% 
  file.path(pu_dir, .)
if (!file.exists(f_lm)) {
  land_mask <- file.path(DATA_DIR, "esri-countries.gpkg") %>% 
    read_sf()  %>% 
    st_transform(crs = projection(r_counts)) %>% 
    filter(NAME != "Antarctica") %>% 
    st_buffer(dist = 10000) %>% 
    st_cast("MULTIPOLYGON") %>% 
    fasterize(r_counts) %>% 
    writeRaster(f_lm, overwrite = TRUE,
                options = c("COMPRESS=LZW", "TILED=YES"))
}
land_mask <- raster(f_lm)
# planning unit mask
f_mask <- str_glue("pu-mask_eck4_{res}km.tif") %>% 
  file.path(pu_dir, .)
r_mask <- r_counts %>% 
  # remove zeros
  subs(data.frame(id = 0, v = NA), subsWithNA = FALSE) %>% 
  mask(land_mask, filename = f_mask, options = c("COMPRESS=LZW", "TILED=YES"))
cellStats(r_mask, range)

# create raster will cells equal to planning unit id
r <- raster(r_mask)
dataType(r) <- "INT4S"
# assign cell id
stopifnot(ncell(r) < 2^31)
values(r) <- seq.int(ncell(r))
# mask to non-zero asset cells
f <- str_glue("pu_eck4_{res}km.tif") %>% 
  file.path(pu_dir, .)
r <- mask(r, r_mask,
          filename = f, overwrite = TRUE, datatype = "INT4S",
          options = c("COMPRESS=LZW", "TILED=YES"))

# save the non-masked cell ids for each template
f_cells <- str_glue("pu_{res}km.rds") %>% 
  file.path(pu_dir, .)
v <- get_raster_values(r) %>% 
  sort()
stopifnot(anyDuplicated(v) == 0)
v <- tibble(id = seq_along(v), cell_id = v)
saveRDS(v, f_cells)
rm(v)
co <- capture.output(gc())

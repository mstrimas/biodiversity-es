# create template raster for planning units are different scales
# cells within land area are assigned cell ids, all other cells get na

library(raster)
library(fasterize)
library(sf)
library(tidyverse)
walk(list.files("R", full.names = TRUE), source)

# directories
pu_dir <- file.path(DATA_DIR, "pu")
dir.create(pu_dir, showWarnings = FALSE)

# target resolutions in km
resolutions <- c(3, 5, 10)

# define raster properties
proj <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
r_template <- extent(c(xmin = -180, xmax = 180, ymin = -90, ymax = 90)) %>% 
  raster(crs = st_crs(4326)$proj4string) %>% 
  projectExtent(crs = proj)
res(r_template) <- resolutions[1] * 1000
dataType(r_template) <- "INT4S"
# assign cell id
stopifnot(ncell(r_template) < 2^31)
values(r_template) <- seq.int(ncell(r_template))

# land
land <- file.path(DATA_DIR, "esri-countries.gpkg") %>% 
  read_sf()  %>% 
  st_transform(crs = proj)

# rasterize
f <- str_glue("pu_eck4_{resolutions[1]}km.tif") %>% 
  file.path(pu_dir, .)
r <- fasterize(land, r_template, field = NULL) %>% 
  mask(r_template, .,
       filename = f, overwrite = TRUE, datatype = "INT4S",
       options = c("COMPRESS=LZW", "TILED=YES"))
rm(r_template)

# rasterize at lower resolutions
for (res in resolutions[-1]) {
  f_agg <- str_glue("pu_eck4_{res}km.tif") %>% 
    file.path(pu_dir, .)
  fact <- res / resolutions[1]
  r_agg <- gdal_aggregate(f, fact = fact)
  # assign cell id
  cell_id <- seq.int(ncell(r_agg))
  cell_id[is.na(r_agg[])] <- NA_integer_
  values(r_agg) <- cell_id
  r_agg <- writeRaster(r_agg, f_agg, overwrite = TRUE, datatype = "INT4S",
                       options = c("COMPRESS=LZW", "TILED=YES"))
}

# save the non-masked cell ids for each template
r_pu <- str_glue("pu_eck4_{resolutions}km.tif") %>% 
  file.path(pu_dir, .) %>% 
  map(raster)
f_cells <- str_glue("pu_{resolutions}km.rds") %>% 
  file.path(pu_dir, .)
for (i in seq_along(r_pu)) {
  v <- get_raster_values(r_pu[[i]]) %>% 
    sort()
  stopifnot(anyDuplicated(v) == 0)
  v <- tibble(id = seq_along(v), cell_id = v)
  saveRDS(v, f_cells[i])
  rm(v)
  co <- capture.output(gc())
}

# prepare tabular input data for prioritizr

library(raster)
library(sf)
library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel(detectCores() - 2)
walk(list.files("R", full.names = TRUE), source)

# resolution in km
resolution <- 5

# directories
tif_dir <- file.path(DATA_DIR, "tifs", paste0(resolution, "km"))
feature_dir <- file.path(DATA_DIR, "features")
pu_dir <- file.path(DATA_DIR, "pu")
dir.create(feature_dir, showWarnings = FALSE)

# raster template
pu <- str_glue("pu_eck4_{resolution}km.tif") %>% 
  file.path(pu_dir, .) %>% 
  raster()
# planning unit raster cell numbers
pu_lookup <- str_glue("pu_{resolution}km.rds") %>% 
  file.path(pu_dir, .) %>% 
  read_rds()

# define features
f <- list.files(tif_dir, 
                str_glue("_{resolution}km.*tif$"), 
                recursive = TRUE, 
                full.names = TRUE) %>% 
  tibble(tif = .) %>% 
  mutate(res_km = resolution,
         type = basename(dirname(tif)),
         layer = basename(tif) %>% 
           str_remove("_[0-9]+km_.*tif$") %>% 
           str_remove("_aoh_[0-9]+km.tif") %>% 
           str_to_lower() %>% 
           str_replace_all("(_|\\s)+", "-"),
         rij = str_glue("rij_{res_km}km_{layer}.rds") %>% 
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
  r <- raster(x$tif) %>% 
    extend(pu) %>% 
    crop(pu)
  if (!compareRaster(r, pu, stopiffalse = FALSE)) {
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
f %>% 
  filter(!has_values) %>% 
  write_csv(file.path(DATA_DIR, 
                      str_glue("no-representation_{resolution}km.csv")))

# assign features ids accounting for layers that have no values associated
f <- f %>% 
  filter(has_values) %>% 
  # put ecosystem services first
  arrange(res_km, fct_relevel(type, "es"), layer) %>% 
  mutate(id = row_number(),
         name = paste(type, layer, sep = "_")) %>% 
  select(id, name, type, layer, res_km, tif, rij)
str_glue("features_{resolution}km.csv") %>% 
  file.path(DATA_DIR, .) %>% 
  write_csv(f, .)

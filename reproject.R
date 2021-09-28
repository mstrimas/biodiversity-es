library(fs)
library(terra)
library(foreach)
library(doParallel)
library(magrittr)
library(glue)
registerDoParallel(cores = 14)

crs <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# generate a raster template
template <- rast("data/tifs_grs80/amphibians/Abavorana-luctuosa_aoh_10km.tif")
names(template) <- "cell_id"
values(template) <- seq_len(ncell(template))
tf <- tempfile(fileext = ".tif")
writeRaster(template, filename = tf)
f_template <- "data/template_eck4_10km.tif"
glue("gdalwarp -t_srs '{crs}' -tr 10000 10000 -r near -co 'COMPRESS=DEFLATE' ",
     "{tf} {f_template}") %>% 
  system()
file_delete(tf)
template <- rast(f_template)

files <- dir("data/tifs_grs80", recursive = TRUE)
projected <- foreach (f = files) %dopar% {
  f_out <- path("data/tifs/", f)
  f %>% 
    path("data/tifs_grs80", .) %>% 
    rast() %>% 
    project(y = template, method = "ngb", filename = f_out)
}
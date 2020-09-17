# aggregate and crop a raster layer using gdal
#
# x: character or Raster object; input raster
# f_out: character; output raster
# fact: integer; aggregation factor
# e: Extent object; crop extent in units of output crs
# overwrite: logical; overwrite output file
#
# return:
#   RasterLayer or RasterStack object
gdal_aggregate <- function(x, f_out, fact, e, overwrite = FALSE) {
  stopifnot(is.numeric(fact), length(fact) == 1)
  
  # input
  if (is.character(x)) {
    stopifnot(length(x) == 1, file.exists(x))
  } else if (inherits(x, "Raster")) {
    x <- raster::filename(x)
    if (!is.character(x) || !file.exists(x)) {
      stop("Can't find input raster on disk.")
    }
  } else {
    stop("Invalid input raster.")
  }
  
  # default to a temp file
  if (missing(f_out)) {
    f_out <- paste0(raster::rasterTmpFile(), ".tif")
  }
  stopifnot(is.character(f_out), length(f_out) == 1)
  if (isTRUE(overwrite)) {
    unlink(f_out)
  }
  
  # output resolution
  tr <- paste(fact * raster::res(raster::raster(x)), collapse = " ")
  tr <- paste("-tr", tr)
  
  # extent to crop to
  if (missing(e)) {
    e <- raster::extent(raster::raster(x))
  }
  stopifnot(inherits(e, c("Extent", "Raster")))
  e <- raster::extent(e)
  ext_arg <- stringr::str_glue("-te {e@xmin} {e@ymin} {e@xmax} {e@ymax}")
  
  # build command
  cmd <- stringr::str_glue("gdalwarp {tr} -r average {ext_arg} ",
                           "-co 'COMPRESS=LZW' {x} {f_out}")
  system(cmd, ignore.stdout = TRUE)
  r <- raster::brick(f_out)
  if (raster::nlayers(r) == 1) {
    r <- raster::raster(f_out)
  }
  return(r)
}
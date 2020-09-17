# extract non-missing, positive values from a raster
#
# x: Raster object to extract from
# cells: cell IDs of values to extract, provided as an integer vector or 
#   filename linking to an rds file containing the integers
# drop_na: should missing values be dropped
# drop_zero: should zero and negative values be dropped
#
# return:
#   numeric vector if cells is missing, data frame if cells if provided
get_raster_values <- function(x, cells = NULL, drop_na = TRUE, 
                              drop_zero = TRUE) {
  stopifnot(inherits(x, "Raster"))
  stopifnot(is.logical(drop_na), length(drop_na) == 1)
  stopifnot(is.logical(drop_zero), length(drop_zero) == 1)
  
  v <- raster::values(x)
  co <- capture.output(gc())
  if (is.null(cells)) {
    if (isTRUE(drop_na)) {
      v <- as.vector(na.omit(v))
    }
    if (isTRUE(drop_zero)) {
      v <- v[v > 0]
    }
  } else {
    if (is.character(cells)) {
      stopifnot(length(cells) == 1, file.exists(cells))
      cells <- readRDS(cells)
    }
    v <- data.frame(cell_id = cells, value = v[cells])
    rm(cells)
    co <- capture.output(gc())
    if (isTRUE(drop_na)) {
      v <- v[!is.na(v$value), ]
      co <- capture.output(gc())
    }
    if (isTRUE(drop_zero)) {
      v <- v[v$value > 0, ]
      co <- capture.output(gc())
    }
  }
  return(v)
}

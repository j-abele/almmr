# Shared helpers for the test suite -------------------------------------------

# Skip a test if the bundled example DEM is not available (e.g. odd install)
skip_without_dem <- function() {
  testthat::skip_if_not(
    nzchar(system.file("extdata", "dem.tif", package = "almmr")),
    "example DEM (inst/extdata/dem.tif) not available"
  )
}

# The example DEM as a SpatRaster
get_dem <- function() almmr::load_dem()

# n random points located on non-NA cells of r (guaranteed inside the surface)
valid_points <- function(r, n, seed = 1) {
  set.seed(seed)
  terra::spatSample(
    r, size = n, method = "random",
    na.rm = TRUE, as.points = TRUE
  )
}

# A coordinate at fractional position (fx, fy) of the DEM extent
xy_at <- function(r, fx, fy) {
  e <- as.vector(terra::ext(r))
  cbind(
    e[["xmin"]] + fx * (e[["xmax"]] - e[["xmin"]]),
    e[["ymin"]] + fy * (e[["ymax"]] - e[["ymin"]])
  )
}

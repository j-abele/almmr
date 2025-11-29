#' Create Hillshdae from SpatRaster-DEM
#'
#' @param dem Digital elevation modell of class SpatRaster
#' @param angle Angle of shading
#' @param direction Light direction
#' @export
#'
create_hillshade <- function(dem,  angle = 45, direction = 315 ) {
  # check if SpatRast
  if (class(dem) != "SpatRaster") {
    dem <- terra::rast(dem)
  }
  slope <- terra::terrain(dem, "slope", unit="radians")
  aspect <- terra::terrain(dem, "aspect", unit="radians")
  hill <- terra::shade(slope, aspect, angle=angle, direction)
}


#' Load the example DEM
#'
#' Loads the example digital elevation model (DEM) that comes with the package.
#' The data file itself is stored as a GeoTIFF in `inst/extdata/` but is
#' accessed through this helper so that users can simply call `data(dem)`.
#'
#' @return A [`SpatRaster`][terra::SpatRaster] object.
#' @examples
#' \dontrun{
#' data(dem)
#' plot(dem)
#' }
#' @export
load_dem <- function() {
  terra::rast(system.file("extdata", "dem.tif", package = "almmr", mustWork = TRUE))
}

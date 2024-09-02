#' Koordinatenbezugssystem vereinheitlichen und Überführung von geo-Objekten in terra-clasen
#'
#' @param raster_object Geo-Raster, e.g. of raster/sp/sf-class
#' @param vector_object Geo-Vector, e.g. of raster/sp/sf-class
#' #export


#' Create Hillshdae from SpatRaster-DEM
#'
#' @param dem Digital elevation modell of class SpatRaster
#' @param angle Angle of shading
#' @param direction Light direction
#' @export
#'
create_hillshade <- function(dem,  angle = 45, direction = 315 ) {
  # check if SpatRast
  slope <- terra::terrain(dem, "slope", unit="radians")
  aspect <- terra::terrain(dem, "aspect", unit="radians")
  hill <- terra::shade(slope, aspect, angle=angle, direction)
}



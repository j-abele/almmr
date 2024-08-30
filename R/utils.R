#' Koordinatenbezugssystem vereinheitlichen und Überführung von geo-Objekten in terra-clasen
#'
#' @param raster_object Geo-Raster, e.g. of raster/sp/sf-class
#' @param vector_object Geo-Vector, e.g. of raster/sp/sf-class
#' #export


#' Create Hillshdae from SpatRaster-DEM
#'
#' @param dem Digital elevation modell of class SpatRaster
#' @param vector_object Geo-Vector, e.g. of raster/sp/sf-class
#' @export
#'
create_hillshade <- function(dem,  angle = 45, driection = 315 ) {
  # check if SpatRast
  slope <- terrain(dem, "slope", unit="radians")
  aspect <- terrain(dem, "aspect", unit="radians")
  hill <- shade(slope, aspect,angle, driection)
}




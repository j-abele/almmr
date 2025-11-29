#' Example Digital Elevation Model
#'
#' Loads the example DEM (Heuneburg region, 25 m resolution).
#' Use `data(dem)` to access the raster via a lazy-loading wrapper.
#'
#' @format A function that returns a [`SpatRaster`][terra::SpatRaster].
#' @source Landesamt für Geoinformation und Landentwicklung Baden-Württemberg
#' @examples
#' \dontrun{
#' data(dem)
#' r <- dem()
#' plot(r)
#' }
"dem"


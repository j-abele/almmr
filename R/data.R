#' Digital elevation model
#'
#' Digital elevation model of the Heuneburg region (5 km rectangle from Heuneburg)
#' with a resolution of 25 m, aggregated from high resolution LiDAR files from the
#' Landesamt für Geoinformation und Landentwicklung Baden-Württemberg.
#'
#' @format A [`SpatRaster`][terra::SpatRaster] object named `dem`.
#' The DEM is projected in the ETRS89 / UTM Zone 32N coordinate system
#' (`EPSG:25832`).

#' @source Landesamt für Geoinformation und Landentwicklung Baden-Württemberg (LGL),
#'  \url{https://www.lgl-bw.de}
#'  #' @examples
#' data(dem)
#' plot(dem)
"dem"

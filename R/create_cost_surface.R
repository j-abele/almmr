#' Create Cost Surface
#'
#' Creates a cost surface based on a digital elevation model (DEM) using
#' different cost functions and optional barrier or wetland factors.
#'
#' @param dem Raster. Digital elevation data. Best as class SpatRaster; sp, sf and raster objects may also work.
#' @param epsg EPSG-Code as string (e.g. "EPSG:25832"). Currently not used internally.
#' @param numberOfNeighbors Whole number. Number of neighbours to be considered when calculating slope (16 is default).
#' @param numberOfDirections Whole number. Number of neighbours for cost surface calculation (4, 8 or 16).
#' @param slopeGainFactor Logical. Default FALSE. If TRUE, a quadratic function starting at \code{slopeGainStart} will be used.
#' @param slopeGainStart Numeric. Start slope value of the quadratic function in slope m (e.g. 0.1 for 10% slope).
#' @param slopeBarrier Logical. Default FALSE. If TRUE, slopes higher than \code{slopeBarrierValue} get NA values (not passable).
#' @param slopeBarrierValue Numeric. Start slope value for barriers (e.g. 0.07 for 7% slope).
#' @param barriers Polygon of class SpatVector. Area of barriers (e.g. rivers, lakes, political boundaries).
#' @param wetlands Polygon of class SpatVector. Area of wetlands/fens.
#' @param wetlandsFactor Numeric. Factor of wetlands friction (e.g. 1.78); applied directly on speed (1 / wetlandsFactor).
#' @param costFunction Character. Cost function. Implemented: "ToblersHikingFunction", "Irmischer-Clarke's", "Wheeled-Vehicels".
#' @param critical_slope Numeric. Critical slope in percent (only relevant for "Wheeled-Vehicels"). Default is 10.
#'
#' @return A \code{TransitionLayer} object representing the cost surface.
#' @export
create_cost_surface <- function(
    dem,
    epsg = NULL, # e.g. "EPSG:25832"
    numberOfNeighbors = 16,
    numberOfDirections = 16,
    slopeGainFactor = FALSE,
    slopeGainStart = 0.10,
    slopeBarrier = FALSE,
    slopeBarrierValue = 0.07,
    barriers = NULL,
    wetlands = NULL,
    wetlandsFactor = 1.78,
    costFunction = "ToblersHikingFunction",
    critical_slope = 10
) {

  # --- Input validation -------------------------------------------------------
  stopifnot(inherits(dem, "SpatRaster"))
  stopifnot(numberOfDirections %in% c(4, 8, 16))
  stopifnot(costFunction %in% c("ToblersHikingFunction", "Irmischer-Clarke's", "Wheeled-Vehicels"))

  # --- Barrier and wetland transitions ---------------------------------------
  if (!is.null(barriers)) {
    barriers_r <- terra::rasterize(barriers, dem, field = 0, background = 1)
    barriers_TR <- gdistance::transition(
      raster::raster(barriers_r),
      transitionFunction = min,
      directions = numberOfDirections
    )
  }

  if (!is.null(wetlands)) {
    wetlands_r <- terra::rasterize(wetlands, dem, field = 1 / wetlandsFactor, background = 1)
    wetlands_TR <- gdistance::transition(
      raster::raster(wetlands_r),
      transitionFunction = min,
      directions = numberOfDirections
    )
  }

  # --- Slope transition object ------------------------------------------------
  # Convert SpatRaster (terra) to RasterLayer (raster) for gdistance
  r_dem <- raster::raster(dem)

  # Altitude difference function for transition
  alt_diff <- function(x) x[2] - x[1]

  # Base slope transition and geo-correction
  slope_trans <- gdistance::transition(r_dem, alt_diff, numberOfNeighbors, symm = FALSE, na.rm = FALSE)
  slope <- gdistance::geoCorrection(slope_trans)

  # --- Apply slope gain / barrier modifications ------------------------------
  m_slope <- gdistance::transitionMatrix(slope)

  # Apply slope gain factor (quadratic increase beyond threshold)
  if (slopeGainFactor) {
    m_slope@x <- ifelse(
      abs(m_slope@x) > slopeGainStart,
      sign(m_slope@x) * ((abs(m_slope@x) / slopeGainStart)^2) * slopeGainStart,
      m_slope@x
    )
  }

  # Apply slope barrier (cells above threshold become impassable)
  if (slopeBarrier) {
    m_slope@x <- ifelse(abs(m_slope@x) > slopeBarrierValue, NA, m_slope@x)
  }

  # Put the modified matrix back into the existing TransitionLayer
  slope@transitionMatrix <- m_slope

  # --- Cost surface calculation ----------------------------------------------
  # Start from the modified slope transition
  speed <- slope

  if (costFunction == "ToblersHikingFunction") {
    # Tobler's Hiking Function: speed in km/h
    m_speed <- gdistance::transitionMatrix(slope)
    m_speed@x <- 6 * exp(-3.5 * abs(m_speed@x + 0.05))

    # Convert to m/s
    m_speed_ms <- m_speed
    m_speed_ms@x <- m_speed@x * 0.278

    speed@transitionMatrix <- m_speed_ms
    cost_surface <- gdistance::geoCorrection(speed)
  }

  if (costFunction == "Irmischer-Clarke's") {
    # Irmischer/Clarke on-path male in km/h
    m_speed <- gdistance::transitionMatrix(slope)
    m_speed@x <- (0.11 + exp(-(abs(m_speed@x) * 100 + 5)^2 / (2 * 30^2))) * 3.6

    # Convert to m/s
    m_speed_ms <- m_speed
    m_speed_ms@x <- m_speed@x * 0.278

    speed@transitionMatrix <- m_speed_ms
    cost_surface <- gdistance::geoCorrection(speed)
  }

  if (costFunction == "Wheeled-Vehicels") {
    # Wheeled vehicles: cost factor based on critical slope
    m_speed <- gdistance::transitionMatrix(slope)
    m_speed@x <- 1 / (1 + ((abs(m_slope@x) * 100) / critical_slope)^2)

    speed@transitionMatrix <- m_speed
    cost_surface <- gdistance::geoCorrection(speed)
  }

  # --- Apply wetlands and barriers -------------------------------------------
  if (!is.null(wetlands)) {
    wetlands_TR <- raster::resample(wetlands_TR, cost_surface)
    cost_surface <- cost_surface * wetlands_TR
  }

  if (!is.null(barriers)) {
    barriers_TR <- raster::resample(barriers_TR, cost_surface)
    cost_surface <- barriers_TR * cost_surface
  }

  # --- Return final cost surface ---------------------------------------------
  return(cost_surface)
}

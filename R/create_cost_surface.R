#' Create Cost Surface
#'
#' Creates a cost surface based on a digital elevation model (DEM) using
#' different cost functions and optional barrier or wetland factors.
#'
#' @param dem Raster. Digital elevation data. Best as class SpatRaster; sp, sf and raster objects may also work.
#' @param epsg EPSG-Code as string (e. g. "EPSG:25832").
#' @param numberOfNeighbors Whole number. Number of neighbours to be considered when calculating slope (16 is default).
#' @param numberOfDirections Whole number. Number of neighbours for cost surface calculation (4, 8 or 16).
#' @param slopeGainFactor TRUE/FALSE. Default FALSE. If TRUE, a quadratic function starting slopeGainStart value will be used.
#' @param slopeGainStart A number.  Start slope value of the quadratic function in slope m (e.g. 0.1 for 10% slope).
#' @param slopeBarrier TRUE/FALSE. Default FALSE. If TRUE, slope higher than given slopeBarrierValue get NA values (not passable).
#' @param slopeBarrierValue Start slope value for barriers (e.g. 0.07 for 7% slope).
#' @param barriers Polygon of class SpatVector. Area of barriers (e.g. a river, lake, political boundaries).
#' @param wetlands Polygon of class SpatVector. Area of wetlands/fens.
#' @param wetlandsFactor Factor of wetlands-friction (e.g. 0.9 for factor 10); calculated directly on speed.
#' @param costFunction Cost function. Implemented: "ToblersHikingFunction",  "Irmischer-Clarke's", "Wheeled-Vehicels".
#' @param critical_slope Critical Slope, only relevant for Wheeled Vehicles function. Default is 10%.
#' @return  A \code{TransitionLayer} object representing the cost surface.
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
    critical_slope = NULL
) {

  # --- Input validation -------------------------------------------------------
  stopifnot(inherits(dem, "SpatRaster"))
  stopifnot(numberOfDirections %in% c(4, 8, 16))
  stopifnot(costFunction %in% c("ToblersHikingFunction", "Irmischer-Clarke's", "Wheeled-Vehicels"))

  # --- Barrier and wetland transitions ---------------------------------------
  if (!is.null(barriers)) {
    barriers_r <- terra::rasterize(barriers, dem, field = 0, background = 1)
    barriers_TR <- gdistance::transition(raster::raster(barriers_r), transitionFunction = min, directions = numberOfDirections)
  }

  if (!is.null(wetlands)) {
    wetlands_r <- terra::rasterize(wetlands, dem, field = 1 / wetlandsFactor, background = 1)
    wetlands_TR <- gdistance::transition(raster::raster(wetlands_r), transitionFunction = min, directions = numberOfDirections)
  }

  # --- Slope transition object ------------------------------------------------
  alt_diff <- function(x) x[2] - x[1]
  slope_trans <- gdistance::transition(dem_r, alt_diff, numberOfNeighbors, symm = FALSE)
  slope <- gdistance::geoCorrection(slope_trans)

  adj <- terra::adjacent(dem, cells = 1:terra::ncell(dem), pairs = TRUE, directions = numberOfDirections)

  # --- Apply slope gain / barrier modifications ------------------------------
  if (slopeGainFactor) {
    slope[adj] <- ifelse(
      abs(slope[adj]) > slopeGainStart,
      sign(slope[adj]) * ((abs(slope[adj]) / slopeGainStart)^2) * slopeGainStart,
      slope[adj]
    )
  }

  if (slopeBarrier) {
    slope[adj] <- ifelse(abs(slope[adj]) > slopeBarrierValue, NA, slope[adj])
  }

  ### Cost Surface calculation ----
  speed <- slope
  # Toblers
  if (costFunction == "ToblersHikingFunction") {
    speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
    speed.ms <- speed
    speed.ms[adj] <- speed[adj] * 0.278
    cost_surface <- gdistance::geoCorrection(speed.ms)
    }

  # Irmischer/Clarke clark on-path male in kmh
  if (costFunction == "Irmischer-Clarke's") {
    speed[adj] <- ((0.11 + exp(-(abs(slope[adj])*100 + 5)^2 / (2 * 30^2))) * 3.6)
    # irmsch/clark on-path male in kmh
    speed.ms <- speed
    speed.ms[adj] <- speed[adj] * 0.278
    cost_surface <- gdistance::geoCorrection(speed.ms)
  }
  ## wheeled-vehicels result is time! 1 + (ŝ / š)²
  if (costFunction == "Wheeled-Vehicels") {
    # Critical slope default 10; +1 could be used for off-terrain factors
    speed[adj] <- (1 / ((1 + ((abs(slope[adj])*100) / critical_slope)^2) * 1))
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

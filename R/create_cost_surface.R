#' Create Cost Surface
#' @param dem Raster. Digital elevation data. Best as class SpatRaster; sp, sf and raster objects may also work. 
#' @param epsg EPSG-Code as string (e. g. "EPSG:25832").
#' @param numberOfNeighbors Whole number. Number of neighbours to be considered when calculating slope (16 is default). 
#' @param numberOfDirections Whole number. Number of neighbours for cost surface calculation (4, 8 or 16).
#' @param slopeGainFactor TRUE/FALSE. Default FALSE. If TRUE, a quadratic function starting slopeGainStart value will be used
#' @param slopeGainStart A number.  Start slope value of the quadratic function in slope m (e.g. 0.1 for 10% slope)
#' @param slopeBarrier TRUE/FALSE. Default FALSE. If TRUE, slope higher than given slopeBarrierValue get NA values (not passable)
#' @param slopeBarrierValue Start slope value for barriers (e.g. 0.07 for 7% slope). 
#' @param barriers Polygon of class SpatVector. Area of barriers (e.g. a river, lake, political boundaries)
#' @param wetlands Polygon of class SpatVector. Area of wetlands/fens. 
#' @param wetlandsFactor Factor of wetlands-friction (e.g. 0.9 for factor 10); calculated directly on speed
#' @param costFunction Cost function. Implemented: "ToblersHikingFunction",  "Irmischer-Clarke's", "Wheeled-Vehicels" 
#' @param critical_slope Critical Slope, only relevant for Wheeled Vehicles function. Default is 10%.
#' @return Cost surface

create_cost_surface <- function (dem, 
                                     epsg = NULL, # urm32n: "EPSG:25832"; GK3: "EPSG:31467" 
                                     numberOfNeighbors = 16,
                                     numberOfDirections = 16,
                                     slopeGainFactor = FALSE,
                                     slopeGainStart = .10, 
                                     slopeBarrier = FALSE,
                                     slopeBarrierValue = .07,
                                     barriers = NULL,
                                     wetlands = NULL,
                                     wetlandsFactor = 1.78,
                                     costFunction = "ToblersHikingFunction",
                                     critical_slope = 10 
) {


  ### Barrier and wetland transition object, if present ----
  if (!is.null(barriers)) {
    barriers.r <- terra::rasterize(barriers, dem, field = 0, background = 1)
    barriers.TR <- gdistance::transition(raster(barriers.r), transitionFunction = min, directions = numberOfDirections)
  }
  
  if (!is.null(wetlands)) {
    wetlands.r <- terra::rasterize(wetlands, dem, field = 1 / wetlandsFactor, background = 1)
    wetlands.TR <- gdistance::transition(raster(wetlands.r), transitionFunction = min, directions = numberOfDirections)
  }
  
  ### Calculation of the slope and transition object of the slope ----
  # Slope for friction surface creation:
  alt.diff <- function(x) {x[2] - x[1]} 
  slope.trans <- gdistance::transition(raster::raster(dem), alt.diff, numberOfNeighbors, symm = F)
  #  Divide by distance between cells:: 
  slope <- gdistance::geoCorrection(slope.trans) 
  adj <- adjacent(dem, cells = 1:ncell(dem), pairs = TRUE, directions = numberOfDirections)
  
  # include slope-barrier-factor and slope-gain factor (calculation is done on transitionMatrix, 
  ## every cell-connection/ edge gets its own value)
  if (slopeGainFactor) {
    slope[adj] <- ifelse(
      abs(slope[adj]) > slopeGainStart,
      sign(slope[adj]) * ((abs(slope[adj]) / slopeGainStart) ^ 2) * slopeGainStart,
      slope[adj]
    )
  }
  # Apply slope barrier factor
  if (slopeBarrier) {
    slope[adj] <- ifelse(
      abs(slope[adj]) > slopeBarrierValue,
      NA,
      slope[adj]
    )
  }
  
  ### Cost Surface calculation ----
  speed <- slope
  # Toblers
  if (costFunction == "ToblersHikingFunction") {
    speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05)) # alternativ 08.01.2022
    speed.ms <- speed
    speed.ms[adj] <- speed[adj] * 0.278 
    cost_surface <- gdistance::geoCorrection(speed.ms)
  }
  # Irmischer/Clarke clark on-path male in kmh
  if (costFunction == "Irmischer-Clarke's") {
    speed[adj] <- ((0.11 + exp(-(abs(slope[adj])*100 + 5)^2 / (2 * 30^2))) * 3.6) # irmsch/clark on-path male in kmh
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
  
  # Apply barriers and wetlands (by multiplying) ----
  if (!is.null(wetlands)) {
    res(wetlands.TR) <- res(cost_surface) 
    cost_surface <- cost_surface * wetlands.TR
  }
  
  if (!is.null(barriers)) {
    res(barriers.TR) <- res(cost_surface) 
    cost_surface <- barriers.TR * cost_surface
  }
  # Return final cost-surface 
  return(cost_surface)
}

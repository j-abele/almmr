#' Create Cost Surface
#'
#' Creates a cost surface based on a digital elevation model (DEM) using
#' different cost functions and optional barrier or wetland factors.
#'
#' @param dem SpatRaster. Digital elevation data.
#' @param numberOfNeighbors Integer. Number of neighbours for slope calculation (4, 8, or 16).
#' @param numberOfDirections Integer. Number of neighbours for graph connectivity (4, 8, 16, 32, or 48).
#' @param slopeGainFactor Logical. Default FALSE. If TRUE, a quadratic penalty is
#'   applied to slopes steeper than \code{slopeGainStart}.
#' @param slopeGainStart Numeric. Slope threshold (rise/run) for quadratic penalty (e.g. 0.1 = 10 percent).
#' @param slopeBarrier Logical. Default FALSE. If TRUE, slopes steeper than
#'   \code{slopeBarrierValue} become impassable.
#' @param slopeBarrierValue Numeric. Slope threshold as ratio (e.g. 0.07 = 7 percent).
#' @param barriers SpatVector (Polygon). Impassable areas (e.g. rivers, lakes).
#' @param wetlands SpatVector (Polygon). Areas with reduced movement speed.
#' @param wetlandsFactor Numeric. Speed reduction factor for wetlands (e.g. 1.78).
#' @param costFunction Character. One of \code{"ToblersHikingFunction"},
#'   \code{"Irmischer-Clarke's"}, \code{"Wheeled-Vehicels"}.
#' @param critical_slope Numeric. Critical slope in percent for \code{"Wheeled-Vehicels"}. Default 10.
#' @param lazy Logical. Default FALSE. If TRUE, edge weights are not pre-computed.
#'   Instead, \code{.build_graph()} computes them on demand for the relevant extent only.
#'   Recommended for large DEMs and functions like \code{perform_tpla()} that operate
#'   on a spatial subset.
#'
#' @return An \code{almmr_cs} object containing the DEM and all parameters.
#'   In eager mode (\code{lazy = FALSE}), also contains pre-computed edge pairs
#'   and weights. Pass to \code{perform_tpla()}, \code{lcsc_territory()},
#'   and \code{sbr_network()}.
#' @export
create_cost_surface <- function(
    dem,
    numberOfNeighbors   = 16,
    numberOfDirections  = 16,
    slopeGainFactor     = FALSE,
    slopeGainStart      = 0.10,
    slopeBarrier        = FALSE,
    slopeBarrierValue   = 0.07,
    barriers            = NULL,
    wetlands            = NULL,
    wetlandsFactor      = 1.78,
    costFunction        = "ToblersHikingFunction",
    critical_slope      = 10,
    lazy                = FALSE
) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!inherits(dem, "SpatRaster"))
    stop("'dem' must be a SpatRaster. Convert with terra::rast() first.")

  if (terra::is.lonlat(dem))
    stop(
      "Geographic (lon/lat) CRS detected. ",
      "Please provide a projected DEM (e.g. UTM). ",
      "Use terra::project() to reproject first."
    )

  if (!numberOfDirections %in% c(4, 8, 16, 32, 48))
    stop("'numberOfDirections' must be one of: 4, 8, 16, 32, 48.")

  if (!numberOfNeighbors %in% c(4, 8, 16))
    stop("'numberOfNeighbors' must be one of: 4, 8, 16.")

  valid_funs <- c("ToblersHikingFunction", "Irmischer-Clarke's", "Wheeled-Vehicels")
  if (!costFunction %in% valid_funs)
    stop(paste("'costFunction' must be one of:", paste(valid_funs, collapse = ", ")))

  if (!is.null(barriers) && !inherits(barriers, "SpatVector"))
    stop("'barriers' must be a SpatVector.")

  if (!is.null(wetlands) && !inherits(wetlands, "SpatVector"))
    stop("'wetlands' must be a SpatVector.")

  if (slopeBarrier && is.null(slopeBarrierValue))
    stop("'slopeBarrierValue' must be set when 'slopeBarrier = TRUE'.")

  if (slopeGainFactor && is.null(slopeGainStart))
    stop("'slopeGainStart' must be set when 'slopeGainFactor = TRUE'.")

  # ---------------------------------------------------------------------------
  # Rasterize barriers and wetlands
  # In lazy mode these are stored as SpatVectors in params and rasterized
  # later in .build_graph() on the clipped DEM extent.
  # In eager mode they are rasterized now on the full DEM.
  # ---------------------------------------------------------------------------
  if (!lazy) {
    if (!is.null(barriers)) {
      if (!terra::same.crs(barriers, dem)) {
        message("'barriers' CRS differs from DEM - reprojecting automatically.")
        barriers <- terra::project(barriers, dem)
      }
      barriers_r <- terra::rasterize(barriers, dem, field = 1L, background = 0L)
    }

    if (!is.null(wetlands)) {
      if (!terra::same.crs(wetlands, dem)) {
        message("'wetlands' CRS differs from DEM - reprojecting automatically.")
        wetlands <- terra::project(wetlands, dem)
      }
      wetlands_r <- terra::rasterize(wetlands, dem, field = 1L, background = 0L)
    }
  }

  # ---------------------------------------------------------------------------
  # Eager mode: compute adj and weights on the full DEM now
  # ---------------------------------------------------------------------------
  if (!lazy) {

    # Valid cells (excluding NA)
    vals  <- terra::values(dem, mat = FALSE)
    valid <- which(!is.na(vals))

    # Neighbour pairs for slope calculation
    adj_slope <- terra::adjacent(dem, cells = valid,
                                 directions = numberOfNeighbors,
                                 pairs = TRUE)

    # Keep only pairs where both cells have values
    adj_slope <- adj_slope[!is.na(vals[adj_slope[, 2]]), ]

    # Euclidean distance between cell centroids (geo-correction)
    xy   <- terra::xyFromCell(dem, seq_len(terra::ncell(dem)))
    dx   <- xy[adj_slope[, 2], 1] - xy[adj_slope[, 1], 1]
    dy   <- xy[adj_slope[, 2], 2] - xy[adj_slope[, 1], 2]
    dist <- sqrt(dx^2 + dy^2)

    # Slope as rise/run (directed: uphill and downhill differ)
    dz    <- vals[adj_slope[, 2]] - vals[adj_slope[, 1]]
    slope <- dz / dist

    # Apply quadratic penalty to slopes steeper than threshold
    if (slopeGainFactor) {
      steep        <- abs(slope) > slopeGainStart
      slope[steep] <- sign(slope[steep]) *
        ((abs(slope[steep]) / slopeGainStart)^2) * slopeGainStart
    }

    # Set slopes steeper than threshold to NA (impassable)
    if (slopeBarrier) {
      slope[abs(slope) > slopeBarrierValue] <- NA
    }

    # Convert slope to travel speed (m/s) using the chosen cost function
    speed_ms <- switch(costFunction,
                       "ToblersHikingFunction" = {
                         # Tobler's Hiking Function: speed in km/h, converted to m/s
                         6 * exp(-3.5 * abs(slope + 0.05)) * (1000 / 3600)
                       },
                       "Irmischer-Clarke's" = {
                         # Irmischer-Clarke on-path male: speed in km/h, converted to m/s
                         (0.11 + exp(-(abs(slope) * 100 + 5)^2 / (2 * 30^2))) * 3.6 * (1000 / 3600)
                       },
                       "Wheeled-Vehicels" = {
                         # Wheeled vehicles: dimensionless cost factor based on critical slope
                         1 / (1 + ((abs(slope) * 100) / critical_slope)^2)
                       }
    )

    # Edge weight = travel time in seconds (distance / speed)
    # NA slopes (from slopeBarrier) propagate naturally to NA weights
    weight <- dist / speed_ms

    # Increase travel time in wetland cells (reduce speed by wetlandsFactor)
    if (!is.null(wetlands)) {
      wetland_cells    <- terra::values(wetlands_r, mat = FALSE) == 1L
      from_wet         <- wetland_cells[adj_slope[, 1]]
      weight[from_wet] <- weight[from_wet] * wetlandsFactor
    }

    # Set travel time to NA in barrier cells (impassable)
    if (!is.null(barriers)) {
      barrier_cells      <- terra::values(barriers_r, mat = FALSE) == 1L
      on_barrier         <- barrier_cells[adj_slope[, 1]] | barrier_cells[adj_slope[, 2]]
      weight[on_barrier] <- NA
    }

  } else {
    # Lazy mode: no computation yet
    adj_slope <- NULL
    weight    <- NULL
  }

  # ---------------------------------------------------------------------------
  # Return almmr_cs object
  # In eager mode: contains adj and weights ready for .build_graph()
  # In lazy mode: contains only DEM and params - .build_graph() computes
  # adj and weights on demand for the clipped extent
  # ---------------------------------------------------------------------------
  structure(
    list(
      dem     = dem,
      adj     = adj_slope,
      weights = weight,
      params  = list(
        cost_function       = costFunction,
        directions          = numberOfDirections,
        slope_neighbors     = numberOfNeighbors,
        slope_gain          = slopeGainFactor,
        slope_gain_start    = slopeGainStart,
        slope_barrier       = slopeBarrier,
        slope_barrier_value = slopeBarrierValue,
        wetlands            = wetlands,
        wetlands_factor     = wetlandsFactor,
        barriers            = barriers,
        critical_slope      = critical_slope,
        lazy                = lazy
      )
    ),
    class = "almmr_cs"
  )
}

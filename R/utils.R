
#' Create Hillshade from SpatRaster-DEM
#'
#' @param dem Digital elevation modell of class SpatRaster
#' @param angle Angle of shading
#' @param direction Light direction
#' @export
#'
  create_hillshade <- function(dem, angle = 45, direction = 315) {
    if (!inherits(dem, "SpatRaster")) {
      dem <- terra::rast(dem)
    }
    slope  <- terra::terrain(dem, "slope",  unit = "radians")
    aspect <- terra::terrain(dem, "aspect", unit = "radians")
    hill   <- terra::shade(slope, aspect, angle = angle, direction = direction)
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

# Internal function - not exported
# Builds a directed igraph object from an almmr_cs object,
# optionally clipped to a given extent.
# Called internally e.g by perform_tpla(), lcsc_territory(), and sbr_network().
# The graph is built on demand and discarded after use.
#
# Supports two modes depending on how create_cost_surface() was called:
#   Eager (lazy = FALSE): adj and weights are pre-computed, only clipped here.
#   Lazy  (lazy = TRUE):  adj and weights are computed here on the clipped DEM.
#
# @param cs   An almmr_cs object from create_cost_surface().
# @param ext  Optional SpatExtent from terra::ext(). If NULL the full DEM is used.
# @return Named list: graph (igraph), cell_to_node (integer), node_to_cell (integer).
.build_graph <- function(cs, ext = NULL) {

  # ---------------------------------------------------------------------------
  # Lazy mode: compute adj and weights on the clipped DEM now
  # ---------------------------------------------------------------------------
  if (isTRUE(cs$params$lazy)) {

    # Clip DEM to extent if provided
    dem <- if (!is.null(ext)) terra::crop(cs$dem, ext) else cs$dem
    p   <- cs$params

    # Rasterize barriers and wetlands on the clipped DEM
    if (!is.null(p$barriers)) {
      if (!terra::same.crs(p$barriers, dem)) {
        p$barriers <- terra::project(p$barriers, dem)
      }
      barriers_r <- terra::rasterize(p$barriers, dem, field = 1L, background = 0L)
    }
    if (!is.null(p$wetlands)) {
      if (!terra::same.crs(p$wetlands, dem)) {
        p$wetlands <- terra::project(p$wetlands, dem)
      }
      wetlands_r <- terra::rasterize(p$wetlands, dem, field = 1L, background = 0L)
    }

    # Valid cells (excluding NA)
    vals  <- terra::values(dem, mat = FALSE)
    valid <- which(!is.na(vals))

    # Neighbour pairs for slope calculation
    adj_slope <- terra::adjacent(dem, cells = valid,
                                 directions = p$slope_neighbors,
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
    if (p$slope_gain) {
      steep        <- abs(slope) > p$slope_gain_start
      slope[steep] <- sign(slope[steep]) *
        ((abs(slope[steep]) / p$slope_gain_start)^2) * p$slope_gain_start
    }

    # Set slopes steeper than threshold to NA (impassable)
    if (p$slope_barrier) {
      slope[abs(slope) > p$slope_barrier_value] <- NA
    }

    # Convert slope to travel speed (m/s) using the chosen cost function
    speed_ms <- switch(p$cost_function,
                       "ToblersHikingFunction" = {
                         6 * exp(-3.5 * abs(slope + 0.05)) * (1000 / 3600)
                       },
                       "Irmischer-Clarke's" = {
                         (0.11 + exp(-(abs(slope) * 100 + 5)^2 / (2 * 30^2))) * 3.6 * (1000 / 3600)
                       },
                       "Wheeled-Vehicels" = {
                         1 / (1 + ((abs(slope) * 100) / p$critical_slope)^2)
                       }
    )

    # Edge weight = travel time in seconds
    weights <- dist / speed_ms

    # Increase travel time in wetland cells
    if (!is.null(p$wetlands)) {
      wetland_cells     <- terra::values(wetlands_r, mat = FALSE) == 1L
      from_wet          <- wetland_cells[adj_slope[, 1]]
      weights[from_wet] <- weights[from_wet] * p$wetlands_factor
    }

    # Set travel time to NA in barrier cells
    if (!is.null(p$barriers)) {
      barrier_cells       <- terra::values(barriers_r, mat = FALSE) == 1L
      on_barrier          <- barrier_cells[adj_slope[, 1]] | barrier_cells[adj_slope[, 2]]
      weights[on_barrier] <- NA
    }

    adj       <- adj_slope
    ncell_ref <- terra::ncell(dem)

  } else {

    # -------------------------------------------------------------------------
    # Eager mode: adj and weights already computed, clip if needed
    # -------------------------------------------------------------------------
    if (!is.null(ext)) {

      # Identify which cells fall within the extent
      cells_in_ext <- terra::cells(cs$dem, terra::vect(ext))

      # Keep only edges where both endpoints are within the extent
      in_ext  <- cs$adj[, 1] %in% cells_in_ext &
        cs$adj[, 2] %in% cells_in_ext

      adj     <- cs$adj[in_ext, ]
      weights <- cs$weights[in_ext]

    } else {

      adj     <- cs$adj
      weights <- cs$weights

    }

    ncell_ref <- terra::ncell(cs$dem)

  }

  # ---------------------------------------------------------------------------
  # Remove edges with NA weights (barriers or impassable slopes)
  # ---------------------------------------------------------------------------
  valid_edges <- !is.na(weights)
  adj         <- adj[valid_edges, ]
  weights     <- weights[valid_edges]

  if (nrow(adj) == 0)
    stop("No valid edges remaining after clip and NA removal. ",
         "Check barriers, slopeBarrier, or reduce the extent.")

  # ---------------------------------------------------------------------------
  # Remap cell indices to contiguous node IDs (1, 2, 3, ...)
  # Required because igraph expects consecutive IDs starting at 1,
  # but cell indices may have gaps after extent clipping.
  # cell_to_node spans the full cell range so all valid indices can be looked up.
  # ---------------------------------------------------------------------------
  unique_cells <- unique(as.integer(adj))
  cell_to_node <- integer(ncell_ref)
  cell_to_node[unique_cells] <- seq_along(unique_cells)

  # Translate edge list from cell indices to node IDs
  adj_remapped <- matrix(
    cell_to_node[as.integer(adj)],
    ncol = 2
  )

  # Build directed weighted graph
  g <- igraph::graph_from_edgelist(adj_remapped, directed = TRUE)
  igraph::E(g)$weight <- weights

  # Return graph + lookup tables so callers can translate
  # between node IDs, cell indices, and coordinates
  list(
    graph        = g,
    cell_to_node = cell_to_node,  # cell index -> node ID
    node_to_cell = unique_cells   # node ID    -> cell index
  )
}

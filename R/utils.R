
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
# Called internally by perform_tpla(), lcsc_territory(), and sbr_network().
# The graph is built on demand and discarded after use (lazy approach).
#
# @param cs   An almmr_cs object from create_cost_surface().
# @param ext  Optional SpatExtent from terra::ext(). If NULL the full DEM is used.
# @return A directed igraph object with edge attribute 'weight' (travel time in seconds).
.build_graph <- function(cs, ext = NULL) {

  # ---------------------------------------------------------------------------
  # Clip edges to extent if provided
  # ---------------------------------------------------------------------------
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

  # ---------------------------------------------------------------------------
  # Remove edges with NA weights (barriers or impassable slopes)
  # ---------------------------------------------------------------------------
  valid   <- !is.na(weights)
  adj     <- adj[valid, ]
  weights <- weights[valid]

  if (nrow(adj) == 0)
    stop("No valid edges remaining after clip and NA removal. ",
         "Check barriers, slopeBarrier, or reduce the extent.")

  # ---------------------------------------------------------------------------
  # Build directed weighted graph from edge list
  # ---------------------------------------------------------------------------

  # Remap cell indices to contiguous node IDs (1, 2, 3, ...)
  # Required because igraph expects consecutive IDs starting at 1,
  # but cell indices may have gaps after extent clipping
  unique_cells <- unique(as.integer(adj))
  # Ensure cell_to_node covers the full DEM cell range
  cell_to_node <- integer(terra::ncell(cs$dem))  # volle DEM-Größe!
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

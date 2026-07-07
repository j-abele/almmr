#' Compute Least-Cost Path
#'
#' Computes a least-cost path between two points using a cost surface
#' configuration. Optionally uses hierarchical resolution refinement
#' for improved performance on large DEMs.
#'
#' Two modes are available - provide exactly one of \code{cs_params}
#' or \code{cs}:
#'
#' \itemize{
#'   \item \strong{Lazy mode} (\code{cs_params}): DEM clipping and graph
#'     construction are handled internally for each path. No pre-computed
#'     cost surface is needed. Recommended for single paths or large DEMs.
#'   \item \strong{Eager mode} (\code{cs}): Uses a pre-computed
#'     \code{almmr_cs} object. Recommended when the same cost surface
#'     is also used by other functions in the same workflow.
#' }
#'
#' @param dem SpatRaster. Digital elevation model.
#' @param origin SpatVector. Start point (single point).
#' @param destination SpatVector. End point (single point).
#' @param cs_params List. Cost surface parameters that define how the
#'   cost surface is computed. Obtain via \code{create_cost_surface()$params}.
#'   In this mode \code{compute_lcp()} handles all DEM clipping and
#'   graph construction internally for each path - no pre-computed
#'   cost surface is needed. This is the recommended approach for
#'   single paths or large DEMs where computing a full cost surface
#'   would be wasteful.
#'   Mutually exclusive with \code{cs}: provide either \code{cs_params}
#'   or \code{cs}, not both.
#' @param cs Optional. An \code{almmr_cs} object from
#'   \code{create_cost_surface(lazy = FALSE)}. If provided, the
#'   pre-computed edge weights are used directly and clipped to the
#'   relevant extent. Recommended when the same cost surface is also
#'   passed to \code{perform_tpla()}, \code{lcsc_territory()}, or
#'   \code{sbr_network()} in the same workflow, avoiding redundant
#'   computation. Mutually exclusive with \code{cs_params}.
#' @param initial_buffer Numeric. Optional absolute buffer in meters
#'   around origin and destination for the initial DEM clip. If set,
#'   takes precedence over \code{initial_buffer_factor}. Recommended
#'   for direct calls where the expected detour distance is known.
#' @param initial_buffer_factor Numeric. Buffer as a fraction of the
#'   straight-line distance between origin and destination. Default 0.3
#'   (30 percent). Used when \code{initial_buffer} is NULL. Recommended for
#'   \code{perform_tpla()} and other functions where detour distance
#'   scales with path length.
#' @param resolutions Numeric vector. DEM resolutions in meters for
#'   hierarchical refinement, ordered coarse to fine
#'   (e.g. \code{c(500, 250)}). The original DEM resolution is always
#'   used as the final step. \code{NULL} skips hierarchical refinement.
#' @param corridor_factor Numeric. Buffer multiplier for corridor
#'   between iterations. Corridor width = \code{corridor_factor *
#'   resolution}. Default 5.
#' @param bidirectional Logical. Default FALSE. If TRUE, both the path
#'   from origin to destination and the return path from destination to
#'   origin are computed. Relevant for anisotropic cost functions such
#'   as \code{"ToblersHikingFunction"} and \code{"Irmischer-Clarke's"}
#'   where uphill and downhill travel speeds differ, meaning the forward
#'   and return paths may diverge. If TRUE, the result contains two
#'   features with a \code{direction} column (\code{"there"} and
#'   \code{"back"}).
#' @param output Character vector. Controls which outputs are returned.
#'   \code{"path"} is always included regardless of this argument.
#'   Optional values:
#'   \describe{
#'     \item{\code{"travel_time_s"}}{Total travel time along the path
#'       in seconds, derived from the sum of edge weights (distance /
#'       speed).}
#'     \item{\code{"distance_m"}}{Total path length in meters, measured
#'       along the actual least-cost path geometry.}
#'     \item{\code{"straight_m"}}{Straight-line (Euclidean) distance
#'       between origin and destination in meters.}
#'     \item{\code{"detour_index"}}{Ratio of path length to straight-line
#'       distance (\code{distance_m / straight_m}). Values above 1
#'       indicate detours caused by terrain. A value of 1 means the
#'       least-cost path follows the direct line.}
#'     \item{\code{"resolution"}}{The DEM resolution in meters at which
#'       the final path was computed. Useful when hierarchical
#'       refinement is used (\code{resolutions} argument) to verify
#'       which resolution level produced the result.}
#'   }
#'
#' @return An \code{sf} object with one feature (or two if
#'   \code{bidirectional = TRUE}). Requested outputs are stored as
#'   attribute columns for direct GIS export via
#'   \code{sf::st_write()}.
#'
#' @examples
#' \dontrun{
#' r      <- load_dem()
#' origin <- terra::vect(cbind(530657, 5326988), crs = terra::crs(r))
#' dest   <- terra::vect(cbind(534500, 5323000), crs = terra::crs(r))
#'
#' # Lazy mode: create cost surface first to obtain params, then compute LCP
#' # No full cost surface is pre-computed - graph built on demand per path
#' cs_lazy <- create_cost_surface(r, lazy = TRUE)
#' path    <- compute_lcp(r, origin, dest, cs_params = cs_lazy$params)
#'
#' # Eager mode: pre-computed cost surface is reused directly
#' # Recommended when cs is also used by perform_tpla() or lcsc_territory()
#' cs_eager <- create_cost_surface(r, lazy = FALSE)
#' path     <- compute_lcp(r, origin, dest, cs = cs_eager)
#'
#' # Hierarchical refinement: coarse to fine for large DEMs
#' # Graph is built at 200m and 100m before final pass at native resolution
#' path_hier <- compute_lcp(
#'   r, origin, dest,
#'   cs_params       = cs_lazy$params,
#'   resolutions     = c(200, 100),
#'   corridor_factor = 5,
#'   initial_buffer  = 5000
#' )
#'
#' # With all outputs and bidirectional
#' path_bi <- compute_lcp(
#'   r, origin, dest,
#'   cs_params     = cs_lazy$params,
#'   bidirectional = TRUE,
#'   output        = c("path", "travel_time_s", "distance_m",
#'                     "straight_m", "detour_index", "resolution")
#' )
#'
#' # Export to GeoPackage
#' sf::st_write(path_bi, "lcp.gpkg")
#'
#' # Visualisation with hillshade
#' hill <- create_hillshade(r)
#' par(mar = c(5, 2, 2, 2))  # extra bottom margin for legend
#' plot(hill, col = grey(0:100 / 100), legend = FALSE, axes = FALSE)
#' plot(r, col = terrain.colors(50, alpha = 0.4),
#'      legend = FALSE, axes = FALSE, add = TRUE)
#' plot(sf::st_geometry(path_bi[path_bi$direction == "there", ]),
#'      add = TRUE, col = "#8b0000", lwd = 2)
#' plot(sf::st_geometry(path_bi[path_bi$direction == "back", ]),
#'      add = TRUE, col = "#00008b", lwd = 2, lty = 2)
#' plot(origin, add = TRUE, col = "black", pch = 21, bg = "#666666", cex = 1.2)
#' plot(dest,   add = TRUE, col = "black", pch = 21, bg = "#666666", cex = 1.2)
#' }
#' @export
compute_lcp <- function(
    dem,
    origin,
    destination,
    cs_params             = NULL,
    cs                    = NULL,
    initial_buffer        = NULL,
    initial_buffer_factor = 0.3,
    resolutions           = NULL,
    corridor_factor       = 5,
    bidirectional         = FALSE,
    output                = "path"
) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  # Exactly one of cs_params or cs must be provided
  if (is.null(cs_params) && is.null(cs))
    stop("Either 'cs_params' or 'cs' must be provided.")
  if (!is.null(cs_params) && !is.null(cs))
    stop("'cs_params' and 'cs' are mutually exclusive. Provide one or the other.")
  if (!is.null(cs) && isTRUE(cs$params$lazy))
    stop(
      "'cs' must be an eager cost surface (lazy = FALSE). ",
      "Use create_cost_surface(lazy = FALSE) or provide 'cs_params' instead."
    )
  if (!is.null(cs) && !inherits(cs, "almmr_cs"))
    stop("'cs' must be an almmr_cs object from create_cost_surface().")

  # Derive cs_params from cs if in eager mode
  if (!is.null(cs))
    cs_params <- cs$params

  # Convert origin and destination to SpatVector if needed
  if (!inherits(origin, "SpatVector"))
    origin <- terra::vect(origin)
  if (!inherits(destination, "SpatVector"))
    destination <- terra::vect(destination)

  # Both must be single points
  if (nrow(origin) != 1)
    stop("'origin' must be a single point.")
  if (nrow(destination) != 1)
    stop("'destination' must be a single point.")

  # DEM must be a projected SpatRaster
  if (!inherits(dem, "SpatRaster"))
    stop("'dem' must be a SpatRaster. Convert with terra::rast() first.")
  if (terra::is.lonlat(dem))
    stop(
      "Geographic (lon/lat) CRS detected. ",
      "Please provide a projected DEM (e.g. UTM). ",
      "Use terra::project() to reproject first."
    )

  # cs_params must contain at least cost_function
  if (!is.list(cs_params) || is.null(cs_params$cost_function))
    stop(
      "'cs_params' must be a named list with at least 'cost_function'. ",
      "Use create_cost_surface()$params or build the list manually."
    )

  # resolutions must be in descending order (coarse to fine)
  if (!is.null(resolutions)) {
    if (!is.numeric(resolutions) || any(resolutions <= 0))
      stop("'resolutions' must be a numeric vector of positive values.")
    if (is.unsorted(rev(resolutions)))
      stop("'resolutions' must be ordered coarse to fine (e.g. c(500, 250)).")
    native_res <- terra::res(dem)[1]
    if (any(resolutions <= native_res))
      stop(sprintf(
        "'resolutions' must be coarser than the native DEM resolution (%.1fm).",
        native_res
      ))
  }

  # output must contain only valid values
  valid_outputs <- c("path", "travel_time_s", "distance_m",
                     "straight_m", "detour_index", "resolution")
  unknown <- setdiff(output, valid_outputs)
  if (length(unknown) > 0)
    stop(paste("Unknown output value(s):", paste(unknown, collapse = ", ")))

  # Reproject origin and destination to DEM CRS if needed
  if (!terra::same.crs(origin, dem)) {
    if (!identical(terra::crs(origin, describe = TRUE)$code,
                   terra::crs(dem,    describe = TRUE)$code)) {
      message("'origin' CRS differs from DEM - reprojecting automatically.")
      origin <- terra::project(origin, dem)
    }
  }
  if (!terra::same.crs(destination, dem)) {
    if (!identical(terra::crs(destination, describe = TRUE)$code,
                   terra::crs(dem,         describe = TRUE)$code)) {
      message("'destination' CRS differs from DEM - reprojecting automatically.")
      destination <- terra::project(destination, dem)
    }
  }

  # ---------------------------------------------------------------------------
  # Clip DEM to study area
  # Buffer is applied around both points individually, then merged via
  # convex hull. This ensures that potential detours behind either point
  # are captured - the optimal path may initially move away from the
  # destination before reaching it.
  # ---------------------------------------------------------------------------
  # Compute initial buffer - absolute value takes precedence over factor
  straight_dist  <- as.numeric(sf::st_distance(
    sf::st_as_sf(origin),
    sf::st_as_sf(destination)
  ))
  buf <- if (!is.null(initial_buffer)) initial_buffer else straight_dist * initial_buffer_factor

  clip_extent <- terra::convHull(
    terra::buffer(
      terra::vect(
        rbind(terra::crds(origin), terra::crds(destination)),
        type = "points",
        crs  = terra::crs(dem)
      ),
      buf
    )
  )
  dem_clip <- terra::crop(dem, clip_extent)

  # Check that both points fall within the clipped DEM
  if (any(is.na(terra::extract(dem_clip, origin)[, 2])) ||
      any(is.na(terra::extract(dem_clip, destination)[, 2])))
    stop(
      "Origin or destination falls outside the DEM extent. ",
      "Check coordinates or increase 'initial_buffer'."
    )

  # ---------------------------------------------------------------------------
  # Internal helper: compute single LCP on a given DEM
  # In lazy mode: builds a temporary almmr_cs and calls .build_graph()
  # In eager mode: clips pre-computed cs to the DEM extent
  # ---------------------------------------------------------------------------
  .compute_single_lcp <- function(dem_input, from, to) {

    if (!is.null(cs) && !isTRUE(cs$params$lazy)) {
      # Eager mode: clip pre-computed adj/weights to dem_input extent
      graph_obj <- .build_graph(cs, ext = terra::ext(dem_input))
    } else {
      # Lazy mode: build graph from scratch on clipped DEM
      params_lazy        <- cs_params
      params_lazy$lazy   <- TRUE   # make sure lazy
      cs_clip <- structure(
        list(
          dem     = dem_input,
          adj     = NULL,
          weights = NULL,
          params  =  params_lazy
          ),
        class = "almmr_cs"
      )
      graph_obj <- .build_graph(cs_clip)
    }

    # Translate points to node IDs
    # In eager mode cell indices must reference the full DEM since
    # cell_to_node is built from full DEM indices in .build_graph()
    # In lazy mode cell indices reference the clipped DEM
    dem_ref   <- if (!is.null(cs)) cs$dem else dem_input
    from_cell <- terra::cellFromXY(dem_ref, terra::crds(from))
    to_cell   <- terra::cellFromXY(dem_ref, terra::crds(to))
    from_node <- graph_obj$cell_to_node[from_cell]
    to_node   <- graph_obj$cell_to_node[to_cell]

    if (is.na(from_node) || from_node == 0)
      stop("Origin does not fall on a valid graph node. ",
           "Increase 'initial_buffer' or check barriers.")
    if (is.na(to_node) || to_node == 0)
      stop("Destination does not fall on a valid graph node. ",
           "Increase 'initial_buffer' or check barriers.")

    # Compute shortest path
    path <- igraph::shortest_paths(
      graph_obj$graph,
      from    = from_node,
      to      = to_node,
      weights = igraph::E(graph_obj$graph)$weight,
      output  = "vpath"
    )$vpath[[1]]

    if (length(path) < 2)
      stop("No path found between origin and destination. ",
           "Check barriers or increase 'initial_buffer'.")

    # Sum edge weights along path for travel time
    path_int   <- as.integer(path)
    edge_times <- igraph::E(graph_obj$graph, path = path_int)$weight

    # Convert node path to coordinates
    # In eager mode node_to_cell references full DEM cell indices
    # In lazy mode it references clipped DEM cell indices
    dem_ref_coords <- if (!is.null(cs)) cs$dem else dem_input
    cells  <- graph_obj$node_to_cell[path_int]
    coords <- terra::xyFromCell(dem_ref_coords, cells)

    list(
      coords      = coords,
      travel_time = sum(edge_times)
    )
  }

  # ---------------------------------------------------------------------------
  # Hierarchical refinement
  # If resolutions is NULL: single pass on clipped DEM at native resolution
  # If resolutions is provided: coarse to fine, each iteration clips the DEM
  # to a corridor around the previous path
  # ---------------------------------------------------------------------------
  if (is.null(resolutions)) {

    # Single pass - no hierarchical refinement
    result_there <- .compute_single_lcp(dem_clip, origin, destination)
    final_res    <- terra::res(dem_clip)[1]

    if (isTRUE(bidirectional))
      result_back <- .compute_single_lcp(dem_clip, destination, origin)

  } else {

    # Build all aggregated DEMs upfront and store for reuse across iterations
    native_res <- terra::res(dem)[1]
    dem_levels <- lapply(resolutions, function(res) {
      fact <- round(res / native_res)
      terra::aggregate(dem_clip, fact = fact, fun = "mean")
    })
    # Always add native resolution as final level
    dem_levels <- c(dem_levels, list(dem_clip))
    res_levels <- c(resolutions, native_res)

    # First iteration: full clipped extent at coarsest resolution
    current_result <- .compute_single_lcp(dem_levels[[1]], origin, destination)
    if (isTRUE(bidirectional))
      current_back <- .compute_single_lcp(dem_levels[[1]], destination, origin)

    # Subsequent iterations: refine within corridor around previous path
    for (i in seq(2, length(dem_levels))) {

      # Build corridor around previous path (there)
      prev_path_sf   <- sf::st_sfc(
        sf::st_linestring(current_result$coords),
        crs = terra::crs(dem)
      )
      # Corridor based on next iteration resolution, not previous
      # e.g. 250m path buffered by 5 * 50m = 250m for the 50m iteration
      corridor_width <- corridor_factor * res_levels[i]
      corridor       <- sf::st_buffer(prev_path_sf, corridor_width)

      # Clip DEM to corridor and refine path
      dem_iter       <- terra::crop(dem_levels[[i]],
                                    terra::vect(corridor), mask = TRUE)
      current_result <- .compute_single_lcp(dem_iter, origin, destination)

      if (isTRUE(bidirectional)) {
        prev_back_sf  <- sf::st_sfc(
          sf::st_linestring(current_back$coords),
          crs = terra::crs(dem)
        )
        corridor_back <- sf::st_buffer(prev_back_sf, corridor_width)
        dem_iter_back <- terra::crop(dem_levels[[i]],
                                     terra::vect(corridor_back), mask = TRUE)
        current_back  <- .compute_single_lcp(dem_iter_back, destination, origin)
      }
    }

    result_there <- current_result
    final_res    <- native_res
    if (isTRUE(bidirectional))
      result_back <- current_back
  }

  # ---------------------------------------------------------------------------
  # Internal helper: build output sf object with requested attribute columns
  # ---------------------------------------------------------------------------
  .build_output_sf <- function(result, from, to, res, direction = NULL) {

    # Path geometry is always included
    geom <- sf::st_sfc(
      sf::st_linestring(result$coords),
      crs = terra::crs(dem)
    )

    # Start with empty data frame
    attrs <- data.frame(row.names = 1)

    # Add direction column if bidirectional
    if (!is.null(direction))
      attrs$direction <- direction

    if ("travel_time_s" %in% output)
      attrs$travel_time_s <- result$travel_time

    if ("distance_m" %in% output)
      attrs$distance_m <- as.numeric(sf::st_length(geom))

    if ("straight_m" %in% output)
      attrs$straight_m <- as.numeric(
        sf::st_distance(sf::st_as_sf(from), sf::st_as_sf(to))
      )

    if ("detour_index" %in% output) {
      dist_m   <- as.numeric(sf::st_length(geom))
      straight <- as.numeric(
        sf::st_distance(sf::st_as_sf(from), sf::st_as_sf(to))
      )
      attrs$detour_index <- dist_m / straight
    }

    if ("resolution" %in% output)
      attrs$resolution <- res

    sf::st_sf(attrs, geometry = geom)
  }

  # ---------------------------------------------------------------------------
  # Assemble and return final result
  # ---------------------------------------------------------------------------
  if (isTRUE(bidirectional)) {
    rbind(
      .build_output_sf(result_there, origin,      destination, final_res, "there"),
      .build_output_sf(result_back,  destination, origin,      final_res, "back")
    )
  } else {
    .build_output_sf(result_there, origin, destination, final_res)
  }
}

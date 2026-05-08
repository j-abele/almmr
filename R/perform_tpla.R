#' Total Passability Landscape Analysis (TPLA)
#'
#' TPLA stands for Total Passability Landscape Analysis. With the help of least-cost paths and a density calculation,
#' regions with high route or movement potential are calculated within a circular analysis area.
#'
#' @param cost_surface almmr_cs object from \code{create_cost_surface()}.
#'   Supports both eager (\code{lazy = FALSE}) and lazy (\code{lazy = TRUE})
#'   modes. In lazy mode the graph is built on demand for the buffer extent of the LCP
#'   only, recommended for large DEMs.
#' @param center_point Center point of the study area. SpatVector (terra).
#'   Will be reprojected automatically if CRS differs from cost surface.
#' @param radius_tpla Numeric. Radius in meters for the circle of starting
#'   points around the center point.
#' @param number_of_points Integer. Number of starting points sampled
#'   regularly along the circle. At least 25 recommended.
#' @param sigma_density_calc Numeric. Standard deviation in meters for the
#'   kernel density estimation applied to the least-cost paths.
#' @param keep_lines Logical. Default FALSE. If TRUE, the least-cost paths
#'   are included in the result as an sf object.
#' @return SpatRaster with kernel density values, or a named list with
#'   elements \code{density} (SpatRaster) and \code{lines} (sf) if
#'   \code{keep_lines = TRUE}.
#' @param resolutions Numeric vector. DEM resolutions in meters for
#'   hierarchical refinement, ordered coarse to fine
#'   (e.g. \code{c(500, 250)}). The original DEM resolution is always
#'   used as the final step. \code{NULL} skips hierarchical refinement.
#'   Passed directly to \code{compute_lcp()}.
#' @param corridor_factor Numeric. Buffer multiplier for corridor between
#'   hierarchical iterations. Corridor width = \code{corridor_factor *
#'   resolution}. Default 5. Passed directly to \code{compute_lcp()}.
#' @export

perform_tpla <- function(cost_surface,
                         center_point,
                         radius_tpla,
                         number_of_points,
                         sigma_density_calc,
                         keep_lines      = FALSE,
                         resolutions     = NULL,
                         corridor_factor = 5) {

  # Convert center_point to SpatVector if needed
  if (!inherits(center_point, "SpatVector")) {
    center_point <- terra::vect(center_point)
  }

  # Check if cost_surface is an almmr_cs object
  if (!inherits(cost_surface, "almmr_cs")) {
    stop("'cost_surface' must be an almmr_cs object from create_cost_surface().")
  }

  # Check if center_point has a CRS
  if (is.na(terra::crs(center_point)) || terra::crs(center_point) == "") {
    stop("'center_point' does not have a valid CRS.")
  }

  # Reproject center_point to DEM CRS if needed
  if (!terra::same.crs(center_point, cost_surface$dem)) {
    epsg_point <- terra::crs(center_point,    describe = TRUE)$code
    epsg_dem   <- terra::crs(cost_surface$dem, describe = TRUE)$code
    if (!identical(epsg_point, epsg_dem)) {
      message("'center_point' CRS differs from cost surface - reprojecting automatically.")
      center_point <- terra::project(center_point, cost_surface$dem)
    }
  }

  # Generate circle for starting points
  # orig buffer_center <- as(buffer(center, radius), "SpatialLines")
  buffer_center <- terra::as.lines(terra::buffer(center_point, radius_tpla))

  # Check if buffer circle is within the cost surface extent
  extent_cost_surface <- terra::vect(terra::ext(cost_surface$dem))
  terra::crs(extent_cost_surface) <- terra::crs(cost_surface$dem)
  if (!terra::relate(extent_cost_surface, buffer_center, "contains")) {
    stop("The buffer circle extends outside the cost surface extent. ",
         "Choose a smaller radius or use a larger DEM.")
  }

  # Warning if less than 25 points (25 because 2 are removed for each set during calculation)
  if (length(1:number_of_points) < 25) {
    warning('Few starting points selected. At least 25 are recommended.')
  }

  # Generate starting points on the circle
  start_points <- terra::vect(sf::st_cast(
    sf::st_line_sample(sf::st_as_sf(buffer_center), number_of_points, type = "regular"), "POINT"
    ))

  # Check for starting points in NA regions (barriers or edge of DEM)
  if (isTRUE(cost_surface$params$lazy)) {
    # Lazy mode: use DEM elevation as proxy - NA = outside DEM or no-data
    passable <- cost_surface$dem
  } else {
    # Eager mode: use mean outgoing weight per cell as passability proxy
    avg_weights <- tapply(cost_surface$weights, cost_surface$adj[, 1], mean, na.rm = FALSE)
    passable    <- terra::rast(cost_surface$dem)
    terra::values(passable) <- NA
    terra::values(passable)[as.integer(names(avg_weights))] <- avg_weights
  }

  count_start_points       <- nrow(start_points)
  start_points$cost_value  <- terra::extract(passable, start_points)[, 2]
  start_points             <- start_points[!is.na(start_points$cost_value), ]

  if (nrow(start_points) < count_start_points) {
    warning(sprintf(
      "%d starting point(s) removed: located in NA region (barrier or DEM edge).",
      count_start_points - nrow(start_points)
    ))
  }

  # Compute least-cost paths between all starting points via compute_lcp()
  # compute_lcp() handles lazy/eager mode, hierarchical refinement and
  # DEM clipping internally - no manual graph management needed here
  sp_paths <- lapply(seq_len(nrow(start_points)), function(i) {
    start         <- start_points[i]
    target_points <- start_points[-i]

    # Remove immediate neighbours to reduce edge effect
    min_dist      <- min(terra::distance(start, target_points))
    target_points <- terra::erase(target_points,
                                  terra::buffer(start, min_dist + min_dist / 2))

    if (nrow(target_points) == 0) return(NULL)

    # Compute LCP from start to each target point
    lapply(seq_len(nrow(target_points)), function(j) {
      tryCatch(
        compute_lcp(
          dem             = cost_surface$dem,
          origin          = start,
          destination     = target_points[j],
          cs_params       = cost_surface$params,
          initial_buffer  = radius_tpla,
          resolutions     = resolutions,
          corridor_factor = corridor_factor
        ),
        error = function(e) NULL  # skip unreachable targets silently
      )
    })
  })

  # Flatten, remove NULLs and combine into single sf object
  all_paths    <- Filter(Negate(is.null), unlist(sp_paths, recursive = FALSE))
  merged_lines <- do.call(rbind, all_paths)

  # Extract segment endpoints from sf lines
  # merged_lines is now a multi-row sf object - one feature per path
  coords_list <- lapply(sf::st_geometry(merged_lines), function(line) {
    coords <- sf::st_coordinates(line)[, 1:2]
    list(
      x0 = coords[-nrow(coords), 1],
      y0 = coords[-nrow(coords), 2],
      x1 = coords[-1,            1],
      y1 = coords[-1,            2]
    )
  })

  # Combine all segments into vectors
  x0 <- unlist(lapply(coords_list, `[[`, "x0"))
  y0 <- unlist(lapply(coords_list, `[[`, "y0"))
  x1 <- unlist(lapply(coords_list, `[[`, "x1"))
  y1 <- unlist(lapply(coords_list, `[[`, "y1"))

  # Build observation window and PSP object for kernel density
  xrange    <- range(c(x0, x1))
  yrange    <- range(c(y0, y1))
  window    <- spatstat.geom::owin(xrange = xrange, yrange = yrange)
  paths_psp <- spatstat.geom::psp(x0 = x0, y0 = y0,
                                  x1 = x1, y1 = y1,
                                  window = window)

  # Kernel density estimation on path segments
  # Resolution matches the original DEM cell size
  res_value    <- terra::res(cost_surface$dem)[1]
  paths_kernel <- spatstat.explore::density.psp(
    paths_psp,
    sigma = sigma_density_calc,
    eps   = res_value
  )
  paths_kernel <- terra::rast(paths_kernel)

  # Return density raster, optionally with sf lines
  if (keep_lines) {
    return(list(density = paths_kernel, lines = merged_lines))
  } else {
    return(paths_kernel)
  }
}


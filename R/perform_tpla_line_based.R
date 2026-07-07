#' Total Passability Landscape Analysis (TPLA) line based
#'
#' TPLA stands for Total Passability Landscape Analysis. Using least-cost paths
#' and a density calculation, it calculates regions with high route or movement
#' potential between two spatial lines.
#'
#' @param cost_surface almmr_cs object from \code{create_cost_surface()}.
#'   Supports both eager (\code{lazy = FALSE}) and lazy (\code{lazy = TRUE})
#'   modes. In lazy mode the graph is built on demand per path, recommended
#'   for large DEMs.
#' @param first_line LINESTRING (sf/sfc). Source line for the start points.
#'   Reprojected to the cost surface CRS automatically if needed.
#' @param second_line LINESTRING (sf/sfc). Target line for the end points.
#' @param number_of_points Integer. Number of points sampled regularly along
#'   each line. At least 20 recommended.
#' @param sigma_density_calc Numeric. Standard deviation in meters for the
#'   kernel density estimation applied to the least-cost paths.
#' @param project_crs Optional target CRS (any input accepted by
#'   \code{sf::st_transform}) applied to both lines before computation.
#' @param keep_lines Logical. Default FALSE. If TRUE, the least-cost paths are
#'   included in the result as an sf object.
#' @param resolutions Numeric vector. DEM resolutions in meters for
#'   hierarchical refinement, ordered coarse to fine (e.g. \code{c(500, 250)}).
#'   \code{NULL} skips hierarchical refinement. Passed to \code{compute_lcp()}.
#' @param corridor_factor Numeric. Buffer multiplier for the corridor between
#'   hierarchical iterations. Default 5. Passed to \code{compute_lcp()}.
#' @return SpatRaster with kernel density values, or a named list with elements
#'   \code{density} (SpatRaster) and \code{lines} (sf) if \code{keep_lines = TRUE}.
#' @export
perform_tpla_line_based <- function(cost_surface,
                                    first_line,
                                    second_line,
                                    number_of_points,
                                    sigma_density_calc,
                                    project_crs     = NULL,
                                    keep_lines      = FALSE,
                                    resolutions     = NULL,
                                    corridor_factor = 5) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!inherits(cost_surface, "almmr_cs"))
    stop("'cost_surface' must be an almmr_cs object from create_cost_surface().")

  # Accept sf or sfc LINESTRING input
  first_line  <- sf::st_geometry(sf::st_as_sf(first_line))
  second_line <- sf::st_geometry(sf::st_as_sf(second_line))

  if (!all(sf::st_geometry_type(first_line)  == "LINESTRING") ||
      !all(sf::st_geometry_type(second_line) == "LINESTRING"))
    stop("Both 'first_line' and 'second_line' must be LINESTRING geometries.")

  if (!is.null(project_crs)) {
    first_line  <- sf::st_transform(first_line,  project_crs)
    second_line <- sf::st_transform(second_line, project_crs)
  }

  if (number_of_points < 20)
    warning("Few starting points selected. At least 20 are recommended.")

  # ---------------------------------------------------------------------------
  # Sample regularly spaced start and end points along the two lines
  # ---------------------------------------------------------------------------
  set.seed(21)
  start_points <- terra::vect(sf::st_cast(
    sf::st_line_sample(first_line,  n = number_of_points, type = "regular"), "POINT"))
  end_points   <- terra::vect(sf::st_cast(
    sf::st_line_sample(second_line, n = number_of_points, type = "regular"), "POINT"))

  # ---------------------------------------------------------------------------
  # Drop points that fall in NA regions (barriers or DEM edge)
  # Passability proxy: DEM in lazy mode, mean outgoing edge weight in eager mode
  # ---------------------------------------------------------------------------
  if (isTRUE(cost_surface$params$lazy)) {
    passable <- cost_surface$dem
  } else {
    avg_weights <- tapply(cost_surface$weights, cost_surface$adj[, 1], mean, na.rm = FALSE)
    passable    <- terra::rast(cost_surface$dem)
    terra::values(passable) <- NA
    terra::values(passable)[as.integer(names(avg_weights))] <- avg_weights
  }

  n_start_in <- nrow(start_points)
  n_end_in   <- nrow(end_points)
  start_points$cost_value <- terra::extract(passable, start_points)[, 2]
  end_points$cost_value   <- terra::extract(passable, end_points)[, 2]
  start_points <- start_points[!is.na(start_points$cost_value), ]
  end_points   <- end_points[!is.na(end_points$cost_value), ]

  if (nrow(start_points) < n_start_in)
    warning(sprintf("%d start point(s) removed: located in NA region (barrier or DEM edge).",
                    n_start_in - nrow(start_points)))
  if (nrow(end_points) < n_end_in)
    warning(sprintf("%d end point(s) removed: located in NA region (barrier or DEM edge).",
                    n_end_in - nrow(end_points)))

  if (nrow(start_points) == 0 || nrow(end_points) == 0)
    stop("No valid start or end points remain after NA removal.")

  # ---------------------------------------------------------------------------
  # Least-cost path from every start point to every end point via compute_lcp()
  # compute_lcp() handles lazy/eager mode, hierarchical refinement and clipping
  # ---------------------------------------------------------------------------
  sp_paths <- lapply(seq_len(nrow(start_points)), function(i) {
    start <- start_points[i]
    lapply(seq_len(nrow(end_points)), function(j) {
      tryCatch(
        compute_lcp(
          dem             = cost_surface$dem,
          origin          = start,
          destination     = end_points[j],
          cs_params       = cost_surface$params,
          resolutions     = resolutions,
          corridor_factor = corridor_factor
        ),
        error = function(e) NULL  # skip unreachable targets silently
      )
    })
  })

  all_paths    <- Filter(Negate(is.null), unlist(sp_paths, recursive = FALSE))
  if (length(all_paths) == 0)
    stop("No least-cost paths could be computed between the two lines.")
  merged_lines <- do.call(rbind, all_paths)

  # ---------------------------------------------------------------------------
  # Convert paths to line segments for the PSP kernel density estimation
  # ---------------------------------------------------------------------------
  coords_list <- lapply(sf::st_geometry(merged_lines), function(line) {
    coords <- sf::st_coordinates(line)[, 1:2]
    list(
      x0 = coords[-nrow(coords), 1],
      y0 = coords[-nrow(coords), 2],
      x1 = coords[-1,            1],
      y1 = coords[-1,            2]
    )
  })

  x0 <- unlist(lapply(coords_list, `[[`, "x0"))
  y0 <- unlist(lapply(coords_list, `[[`, "y0"))
  x1 <- unlist(lapply(coords_list, `[[`, "x1"))
  y1 <- unlist(lapply(coords_list, `[[`, "y1"))

  xrange    <- range(c(x0, x1))
  yrange    <- range(c(y0, y1))
  window    <- spatstat.geom::owin(xrange = xrange, yrange = yrange)
  paths_psp <- spatstat.geom::psp(x0 = x0, y0 = y0,
                                  x1 = x1, y1 = y1,
                                  window = window)

  # Kernel density estimation on path segments (resolution = DEM cell size)
  res_value    <- terra::res(cost_surface$dem)[1]
  paths_kernel <- spatstat.explore::density.psp(
    paths_psp,
    sigma = sigma_density_calc,
    eps   = res_value
  )
  paths_kernel <- terra::rast(paths_kernel)

  if (keep_lines) {
    return(list(density = paths_kernel, lines = merged_lines))
  } else {
    return(paths_kernel)
  }
}

#' Total Passability Landscape Analysis (TPLA)
#'
#' TPLA stands for Total Passability Landscape Analysis. With the help of least-cost paths and a density calculation,
#' regions with high route or movement potential are calculated within a circular analysis area.
#'
#' @param cost_surface Transition Object. Cost surface (Class: Transition, calculated with the gdistance package)
#' @param center_point Center point of the study area. Class of SpatVector (terra).
#' @param radius_tpla Integer. Radius for the circle of starting points (in meters), in wich tpla will take place.
#' @param number_of_points Integer. Number of starting points on the circle.
#' @param sigma_density_calc Number. Standard deviation for the kernel density estimation.
#' @param keep_lines TRUE or FALSE. Default is FALSE. If TRUE, the cost-optimal paths will be included in the result object.
#' @return List or raster. If keep_lines = TRUE, a list object containing the result raster of the kernel density estimation and the cost-optimal paths will be returned.
#' @export

perform_tpla <- function(cost_surface,
                         center_point,
                         radius_tpla,
                         number_of_points,
                         sigma_density_calc,
                         keep_lines=FALSE) {

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
  # Use mean outgoing weight per cell as proxy for passability
  avg_weights <- tapply(cost_surface$weights, cost_surface$adj[, 1], mean, na.rm = FALSE)
  passable    <- terra::rast(cost_surface$dem)
  terra::values(passable) <- NA
  terra::values(passable)[as.integer(names(avg_weights))] <- avg_weights

  count_start_points       <- nrow(start_points)
  start_points$cost_value  <- terra::extract(passable, start_points)[, 2]
  start_points             <- start_points[!is.na(start_points$cost_value), ]

  if (nrow(start_points) < count_start_points) {
    warning(sprintf(
      "%d starting point(s) removed: located in NA region (barrier or DEM edge).",
      count_start_points - nrow(start_points)
    ))
  }

  # Build graph once for the full buffer extent (lazy)
  # Use center point buffer as extent - guarantees all circle points are included
  graph_obj <- .build_graph(
    cost_surface,
    ext = terra::ext(terra::buffer(center_point, radius_tpla + terra::res(cost_surface$dem)[1] * 3))
  )

  # build spiders-web :)
  sp_paths <- lapply(seq_len(nrow(start_points)), function(i) {
    start         <- start_points[i]
    target_points <- start_points[-i]

    # Remove immediate neighbours to reduce edge effect
    min_dist      <- min(terra::distance(start, target_points))
    target_points <- terra::erase(target_points,
                                  terra::buffer(start, min_dist + min_dist / 2))

    # Translate start cell to node ID
    start_cell <- terra::cellFromXY(cost_surface$dem, terra::crds(start))
    start_node <- graph_obj$cell_to_node[start_cell]

    # Translate target cells to node IDs
    target_cells <- terra::cellFromXY(cost_surface$dem, terra::crds(target_points))
    target_nodes <- graph_obj$cell_to_node[target_cells]

    # Remove targets that fall outside the clipped graph (node ID = 0)
    target_nodes <- target_nodes[target_nodes > 0]
    if (length(target_nodes) == 0) return(NULL)

    # Compute shortest paths from start to all targets
    paths <- igraph::shortest_paths(
      graph_obj$graph,
      from    = start_node,
      to      = target_nodes,
      weights = igraph::E(graph_obj$graph)$weight,
      output  = "vpath"
    )

    # returns coordinate matrix, no sf object
    lapply(paths$vpath, function(path) {
      if (length(path) < 2) return(NULL)
      cells  <- graph_obj$node_to_cell[as.integer(path)]
      coords <- terra::xyFromCell(cost_surface$dem, cells)
      # Return plain matrix - sf objects built later to avoid nested list issues
      coords
    })
  })

  # Collect all coordinate matrices and remove NULLs
  all_coords <- Filter(Negate(is.null),
                       do.call(c, sp_paths))

  # Build sf object from coordinate matrices
  merged_lines <- sf::st_sf(
    geometry = sf::st_sfc(
      lapply(all_coords, sf::st_linestring),
      crs = terra::crs(cost_surface$dem)
    )
  )

  # Extract segment endpoints from sf lines directly
  # sf geometries are already coordinate matrices - no S4 slot access needed
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


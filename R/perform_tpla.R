#' Total Passability Landscape Analysis (TPLA)
#'
#' TPLA stands for Total Passability Landscape Analysis. With the help of least-cost paths and a density calculation,
#' regions with high route or movement potential are calculated within a circular analysis area.
#' Placeholder: Description to follow (suggestions for research-area, number of points, Spatial-classes etc.)
#'
#' @param cost_surface Transition Object. Cost surface (Class: Transition, calculated with the gdistance package)
#' @param center_point Center point of the study area. Class of SpatVector (terra).
#' @param radius_tpla Integer. Radius for the circle of starting points (in meters), in wich tpla will take place.
#' @param number_of_points Integer. Number of starting points on the circle.
#' @param sigma_density_calc Number. Standard deviation for the kernel density estimation.
#' @param keep_lines TRUE or FALSE. Default is FALSE. If TRUE, the cost-optimal paths will be included in the result object.
#' @return S4-Class or Raster. If keep_lines = TRUE, an S4 object containing the result raster of the kernel density estimation and the cost-optimal paths is returned.
#' @export
#'
#' @examples
#' # fun_TPLA(cost_surface, center_point, radius_tpla, number_of_points, sigma_density_calc, keep_lines=FALSE)

perform_tpla <- function(cost_surface,
                         center_point,
                         radius_tpla,
                         number_of_points,
                         sigma_density_calc,
                         project_crs=NULL,
                         keep_lines=FALSE) {


  # center_point als terra::vect
  if (class(center_point) != "SpatVector") {
    center_point <- terra::vect(center_point)
  }

  # Check if project CRS is given, otherwise use center_point point CRS. Warning and stop if no CRS is given.
  if (!is.null(project_crs)) {
    crs(center_point) <- project_crs
  } else if (is.null(project_crs) & !is.null(crs(center_point))) {
    project_crs <- crs(center_point)
  } else {
    warning('No CRS provided as an argument; the focal point has no CRS. Please specify a coordinate reference system (project_crs = "")')
    stop('No valid CRS provided.')
  }

  # Generate circle for starting points
  # orig buffer_center <- as(buffer(center, radius), "SpatialLines")
  buffer_center <- terra::as.lines(terra::buffer(center_point, radius_tpla))

  # Check if buffer is outside of cost_surface extent, stop process if yes
  # original: extent_cost_surface <- as(extent(raster(cost_surface)), "SpatialPolygons")
  extent_cost_surface <- terra::vect(terra::ext(cost_surface@extent))
  crs(extent_cost_surface) <- crs(center_point)
  crs(buffer_center) <- crs(center_point)

  # Return warning if not:
  if (!relate(extent_cost_surface, buffer_center, relation="contains")) {
    warning('The radius might be outside the cost surface. Choose a smaller radius or use a larger cost surface as input.')
    stop()
  }

  # Warning if less than 25 points (25 because 2 are removed for each set during calculation)
  if (length(number_of_points) < 25) {
    warning('Few starting points selected. At least 25 are recommended.')
  }

  # Generate starting points on the circle
  start_points <- terra::vect(sf::st_cast(sf::st_line_sample(sf::st_as_sf(buffer_center), number_of_points, type = "regular"), "POINT"))

  # exclude start_points in NA-regions (barriers in cost surface) and give warning
  start_points$cost_value <- terra::extract(terra::rast(raster::raster(cost_surface)), start_points)[,2]
  check_vect <- is.na(start_points$cost_value)
  count_start_points <- length(start_points)
  start_points <- na.omit(start_points, "cost_value")
  if (!count_start_points == length(start_points)) {
    warning(paste("Some starting points (total: ", count_start_points - length(start_points), ") are located within a region of the cost surface that contains no values (NA). This is likely due to a barrier or an error in the data. This will be excluded from the calculation.", sep=""))
  }

  # Create shortest paths
  sp_paths <-  lapply(seq_along(start_points), function(i) {
    start <- start_points[i]
    target_points <- start_points[-i]
    crs(target_points) <- project_crs
    crs(start) <- project_crs
    # Remove neighborhood points to reduce edge-effect
    # select neigbouring point by Buffer of closest point (may be better with distance-function of terra)
    min_dist <- min(terra::distance(start, target_points))
    target_points <- terra::erase(target_points, terra::buffer(start, min_dist+min_dist/2)) # + half of mininmum distance to make sure, closest points are eraesed
    # start and target to sp for gdistance-package
    target_points <- as(target_points , "Spatial")
    start <- as(start, "Spatial")
    gdistance::shortestPath(cost_surface, start, target_points, output = "SpatialLines")
  })
  ### merge Lines
  merged_lines <- do.call(rbind, sp_paths)

  ###########
  ### SpatialLines into PSP-Object (since maptools is autdated there is no working solution e.g. with sf)

  # Extrahieren aller Linienkoordinaten in einer Liste
  all_coords <- lapply(merged_lines@lines, function(line) {
    do.call(rbind, lapply(line@Lines, function(segment) {
      segment@coords
    }))
  })

  # Anzahl der Linien
  num_lines <- length(all_coords)

  # Initialisierung der Vektoren für die Start- und Endpunkte
  x0 <- y0 <- x1 <- y1 <- numeric(0)

  # Vektorisierte Verarbeitung aller Linien
  for (i in seq_len(num_lines)) {
    coords <- all_coords[[i]]
    n <- nrow(coords)
    x0 <- c(x0, coords[-n, 1])
    y0 <- c(y0, coords[-n, 2])
    x1 <- c(x1, coords[-1, 1])
    y1 <- c(y1, coords[-1, 2])
  }

  # Bestimmen der Bounding Box für das Beobachtungsfenster
  xrange <- range(c(x0, x1))
  yrange <- range(c(y0, y1))
  window <- spatstat.geom::owin(xrange = xrange, yrange = yrange)

  # Erstellen des psp-Objekts
  paths_psp <- spatstat.geom::psp(x0 = x0, y0 = y0, x1 = x1, y1 = y1, window = window)


  # Kerndichteberchnung
  paths_kernel <- terra::rast(density(paths_psp, sigma=sigma_density_calc))

  ########
  ## create list if keep_lines = TRUE (better S4-Object?)
  if(keep_lines) {
    return(list(merged_lines, paths_kernel))
  } else {
    return(paths_kernel)
  }
}


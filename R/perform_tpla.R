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
#' @return List or raster. If keep_lines = TRUE, a list object containing the result raster of the kernel density estimation and the cost-optimal paths will be returned.
#' @examples
#'
#' # Create cost surface
#' cost_surface <- almmr::create_cost_surface(dem, slopeGainFactor=TRUE)
#' site <- terra::vect(cbind(530681, 5326949))
#' terra::crs(site) <- cost_surface@srs
#'
#' # Perform TPLA
#' tpla_result <- almmr::perform_tpla(cost_surface, site, 4500, 50, 200)
#' plot(tpla_result)
#' plot(site, add=T)
#'
#'
#' Optional:
#' dem_hill <- almmr::create_hillshade(dem)
#' plot(dem_hill, col = gray.colors(100), main = "TPLA Results with Hillshade", legend = FALSE)
#' plot(tpla_result, add = TRUE, col = hcl.colors(100, "YlOrRd"))
#' plot(site, add=T)

perform_tpla <- function(cost_surface,
                         center_point,
                         radius_tpla,
                         number_of_points,
                         sigma_density_calc,
                         keep_lines=FALSE) {

  # center_point as terra::vect
  if (class(center_point) != "SpatVector") {
    center_point <- terra::vect(center_point)
  }

  # Check if cost_surface is a Transition-Objekt
  if (!inherits(cost_surface, "Transition")) {
    stop("The cost_surface must be a Transition object.")
  }

  # EPSG-Code from cost_surface
  cost_surface_epsg <- cost_surface@srs
  if (is.null(cost_surface_epsg) || cost_surface_epsg == "") {
    stop("The Transition object does not have a valid EPSG code.")
  }

  # Check if center_point has CRS
  center_point_epsg <- terra::crs(center_point, proj = TRUE)
  if (is.null(center_point_epsg) || center_point_epsg == "") {
    stop("The center_point does not have a valid CRS or EPSG code.")
  }


  # check for same CRS
  if (cost_surface_epsg != center_point_epsg) {
    stop(sprintf("CRS mismatch: cost_surface EPSG (%s) and center_point EPSG (%s) do not match.",
                 cost_surface_epsg, center_point_epsg))
  }

  # Generate circle for starting points
  # orig buffer_center <- as(buffer(center, radius), "SpatialLines")
  buffer_center <- terra::as.lines(terra::buffer(center_point, radius_tpla))
  terra::crs(buffer_center) <- terra::crs(center_point)

  # Check if buffer is outside of cost_surface extent, stop process if yes
  extent_cost_surface <- terra::vect(terra::ext(cost_surface@extent))
  terra::crs(extent_cost_surface) <- terra::crs(center_point)
  if (!terra::relate(extent_cost_surface, buffer_center, "contains")) {
    stop("The radius might be outside the cost surface. Choose a smaller radius or use a larger cost surface.")
  }

  # Warning if less than 25 points (25 because 2 are removed for each set during calculation)
  if (length(1:number_of_points) < 25) {
    warning('Few starting points selected. At least 25 are recommended.')
  }

  # Generate starting points on the circle
  start_points <- terra::vect(sf::st_cast(
    sf::st_line_sample(sf::st_as_sf(buffer_center), number_of_points, type = "regular"), "POINT"
    ))

  # exclude start_points in NA-regions (barriers in cost surface) and give warning
  start_points$cost_value <- terra::extract(terra::rast(raster::raster(cost_surface)), start_points)[,2]
  count_start_points <- length(start_points)
  start_points <- na.omit(start_points, "cost_value")
  if (!count_start_points == length(start_points)) {
    warning(paste("Some starting points (total: ", count_start_points - length(start_points), ") are located within a region of the cost surface that contains no values (NA). This is likely due to a barrier or an error in the data. This will be excluded from the calculation.", sep=""))
  }

  # Create shortest paths
  sp_paths <-  lapply(1:nrow(start_points), function(i) {
    start <- start_points[i]
    target_points <- start_points[-i]
    crs(target_points) <- terra::crs(center_point)
    crs(start) <- terra::crs(center_point)

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
  paths_kernel <- spatstat.explore::density.psp(paths_psp, sigma=sigma_density_calc, eps=(terra::res(terra::rast(dem))))
  paths_kernel <- terra::rast(paths_kernel)

  ########
  ## create list if keep_lines = TRUE (better S4-Object?)
  if(keep_lines) {
    return(list(merged_lines, paths_kernel))
  } else {
    return(paths_kernel)
  }
}


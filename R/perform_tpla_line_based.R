#' Total Passability Landscape Analysis (TPLA) line based
#'
#' TPLA stands for Total Passability Landscape Analysis. Using least-cost paths and a density calculation,
#' it calculates regions with high routes or movement potential between two spatial lines.
#'
#' @param cost_surface Transition Object. Cost surface (Class: Transition, calculated with the gdistance package)
#' @param first_line Spatial Line. Class of LINESTRING (sf).
#' @param second_line  Spatial Line. Class of LINESTRING (sf).
#' @param number_of_points Integer. Number of starting points on the circle.
#' @param sigma_density_calc Number. Standard deviation for the kernel density estimation.
#' @param keep_lines TRUE or FALSE. Default is FALSE. If TRUE, the cost-optimal paths will be included in the result object.
#' @return List or Raster. If keep_lines = TRUE, an S4 object containing the result raster of the kernel density estimation and the cost-optimal paths is returned.

perform_tpla_line_based <- function(cost_surface,
                                    first_line,
                                    second_line,
                                    number_of_points,
                                    sigma_density_calc,
                                    project_crs=NULL,
                                    keep_lines=FALSE) {


  # Überprüfen und Transformieren der CRS für Linien; prüfen ob linien außerhalb von cost-surface; warning wenn relativ na an grenze?
  if (!inherits(first_line, "sfc_LINESTRING") || !inherits(second_line, "sfc_LINESTRING")) {
    stop("Both first_line and second_line must be LINESTRING objects of class sf.")
  }

  if (!is.null(project_crs)) {
    first_line <- sf::st_transform(first_line, project_crs)
    second_line <- sf::st_transform(second_line, project_crs)
  }

  # Warning if less than 20 points
  if (length(number_of_points) < 20) {
    warning('Few starting points selected. At least 20 are recommended.')
  }

  # to-do: Prüfroutine für CRS und Object-Classe implementieren; ggfls. Transformaiton
  set.seed(21)
  start_points <- terra::vect(sf::st_line_sample(first_line, n = number_of_points, type = "regular"))
  end_points <-   terra::vect(sf::st_line_sample(second_line, n = number_of_points, type = "regular"))


  # exclude start_points and end_points in NA-regions (barriers in cost surface) and give warning
  start_points$cost_value <- terra::extract(terra::rast(raster::raster(cost_surface)), start_points)[,2]
  end_points$cost_value <- terra::extract(terra::rast(raster::raster(cost_surface)), end_points)[,2]



  # Count original Points (for checking) and remove points with NA-values
  count_start_points <- length(sp::SpatialPoints(terra::crds(start_points)))
  count_end_points <- length(sp::SpatialPoints(terra::crds(end_points)))
  start_points <- na.omit(start_points, "cost_value")
  end_points <- na.omit(end_points, "cost_value")
  # Convert Points to sp-SpatialClass for gdistance::shortestPath
  start_points <- sp::SpatialPoints(terra::crds(start_points), proj4string = sp::CRS(terra::crs(start_points)))
  end_points <- sp::SpatialPoints(terra::crds(end_points), proj4string = sp::CRS(terra::crs(end_points)))

  if (!count_start_points == length(start_points)) {
    warning(paste("Some starting points (total: ", length(start_points) - count_start_points, ") are located within a region of the cost surface that contains no values (NA). This is likely due to a barrier or an error in the data. This will be excluded from the calculation.", sep=""))
  }

  if (!count_end_points == length(start_points)) {
    warning(paste("Some target points (total: ", length(end_points) - count_end_points, ") are located within a region of the cost surface that contains no values (NA). This is likely due to a barrier or an error in the data. This will be excluded from the calculation.", sep=""))
  }


  # Create shortest paths
  sp_paths <-  lapply(seq_along(start_points), function(i) {
    start <- start_points[i]
    # start and target to sp for gdistance-package
    gdistance::shortestPath(cost_surface, start,  end_points, output = "SpatialLines")
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
  paths_kernel <- density(paths_psp, sigma=sigma_density_calc)
  paths_kernel <- terra::rast(paths_kernel)

  ########
  ## create list if keep_lines = TRUE (better S4-Object?)
  if(keep_lines) {
    return(list(merged_lines, paths_kernel))
  } else {
    return(paths_kernel)
  }
}

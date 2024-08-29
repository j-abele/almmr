#' Site-based (theoretical) route network (sbr_network)
#'
#' Generates a theoretical route network based on site locations, a given route, based on conductance values. 
#' Strongly inspired by the least-cost network-to-builder process described by David Waugh in 2000 (Waugh 2000, 615, see also Herzog 2013, 239).
#' 
#' Literature:
#' D. Waugh, Geography. An Integrated Approach. Third edition (Cheltenham 2000).
#' 
#' I. Herzog, Least-cost Networks. In: G. Earl/T. Sly/A. Chrysanthi/P. Murrieta-Flores/C.
#' Papadopoulos/I. Romanowska/D. Wheatley (Hrsg.), Archaeology in the Digital Era II. e-Papers
#' from the 40th Conference on Computer Applications and Quantitative Methods in Archaeology
#' (Amsterdam 2013) 237-248.
#'
#' @param sites SpatialPointsDataFrame or compatible object. The site locations.
#' @param lines SpatialLinesDataFrame. The base lines to start from (known route).
#' @param conductance Transition Layer. Cost surface (gidstance).
#' @param steps_points Integer. The number of steps (normaly in Meter) for sampling points along the lines.
#' @param max_speed Numeric. The maximum speed (optional, can be refined based on use case).
#' @return A SpatialLinesDataFrame object with the calculated routes.
#' @export
#'
#' @examples
#' 

sbr_network <- function(sites, lines, conductance, steps_points, max_speed) { 
  # Überprüfen und ggf. Konvertieren der Eingabedaten
  if (!inherits(sites, "SpatialPointsDataFrame")) {
    stop("The 'sites' object must be a SpatialPointsDataFrame.")
  }
  
  if (!inherits(lines, "SpatialLinesDataFrame")) {
    lines <- as(lines, "SpatialLinesDataFrame")
  }
  
  # Erstellen der Ausdehnung für die Kostenoberfläche
  extent.area <- as(extent(raster::raster(conductance)), "SpatialPolygons")
  
  # Zuschneiden der Punkte und Linien auf die Kostenoberfläche
  sites <- raster::crop(sites, extent.area)
  lines <- raster::crop(lines, extent.area)
  
  # Initialisieren der zusätzlichen Attribute für die Linien
  lines@data <- data.frame(ID = 1:length(lines), dist_m = 0, time_minutes = 0, hub_dist_next_site = 0, time_to_next_site = 0)
  
  # Ursprüngliche Linien und Standorte sichern
  line_start <- lines
  sites_original <- sites
  
  stop_it <- FALSE
  while (!stop_it) {
    # Bestimmen der Anzahl der Punkte entlang der Linie
    numOfPoints <- round(rgeos::gLength(lines) / steps_points)
    
    # Erstellen von Wegpunkten in regelmäßigen Abständen entlang der Linie
    way.points <- sp::spsample(lines, n = numOfPoints, type = "regular")
    way.points <- as(way.points, "SpatialPointsDataFrame")
    way.points@data <- data.frame(dist = rep(1, length(way.points)))
    
    # Berechnung der minimalen Distanz zu den Standorten
    acc_way_points <- gdistance::accCost(conductance, way.points) / 60
    sites$min_dist <- raster::extract(acc_way_points, sites)
    
    # NA-Werte entfernen und eine Warnung ausgeben, falls vorhanden
    initial_start_count <- nrow(sites)
    sites <- sites[!is.na(sites$min_dist), ]
    removed_start_count <- initial_start_count - nrow(sites)
    if (removed_start_count > 0) {
      warning(paste("Removed", removed_start_count, "sites due to NA values in cost surface extraction."))
    }
    
    # Standort für den nächsten Durchlauf auswählen
    origin_new <- sites[sites$used == "no" & sites$min_dist == min(sites$min_dist),] 
    
    # Standort für die zukünftige Berechnung ausschließen
    sites$used[sites$min_dist == origin_new$min_dist & sites$used == "no"] <- "yes"
    
    # Ausschneiden und Berechnen der Distanzen
    if (terra::nrow(way.points) > 1) {
      way.points$dist <- gdistance::costDistance(conductance, way.points, origin_new)
    }
    
    # Neuen Zielpunkt bestimmen
    min_dist <- min(way.points$dist, na.rm = TRUE)
    target_new <- way.points[way.points$dist == min_dist, ][1, ]
    
    # Kürzesten Pfad berechnen und Attribute zuweisen
    route_new <- gdistance::shortestPath(conductance, origin_new, target_new, "SpatialLines")
    route_new$dist_m <- min_dist
    route_new$time_minutes <- min_dist / 60
    route_new$hub_dist_next_site <- 0
    
    # Minimale Entfernung zum nächsten Standort berechnen
    if (!any(terra::relate(terra::vect(sites_original), terra::vect(target_new), "touches"))) {
      sites_without_origin <- sites_original[!sites_original$id %in% origin_new$id, ]
      route_new$hub_dist_next_site <- min(gdistance::costDistance(conductance, sites_without_origin, target_new)) / 60
    }
    
    route_new$time_to_next_site <- route_new$time_minutes + route_new$hub_dist_next_site
    lines <- rbind(lines, route_new)
    
    # Entfernen der 'min_dist'-Spalte und Filterung
    sites$min_dist <- NULL
    sites <- sites[sites$used == "no", ]
    if (terra::nrow(sites) == 0) {
      stop_it <- TRUE
    }
  }
  
  return(lines)
}
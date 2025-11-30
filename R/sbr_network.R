#' Site-based (theoretical) route network (sbr_network)
#'
#' Generates a theoretical route network based on site locations and a given route,
#' based on conductance values. Each iteration connects the next unlinked site
#' with the existing network through the least-cost path.
#'
#' Strongly inspired by the least-cost network-to-builder process described by
#' David Waugh (2000, 615), see also Herzog (2013, 239).
#'
#' Literature:
#' D. Waugh, Geography. An Integrated Approach. Third edition (Cheltenham 2000).
#' I. Herzog, Least-cost Networks. In: G. Earl/T. Sly/A. Chrysanthi/P. Murrieta-Flores/C. Papadopoulos/I. Romanowska/D. Wheatley (Hrsg.), Archaeology in the Digital Era II. e-Papers from the 40th Conference on Computer Applications and Quantitative Methods in Archaeology (Amsterdam 2013) 237-248.
#'
#' @param sites SpatialPoints, SpatialPointsDataFrame, sf, or SpatVector. The site locations.
#' @param lines SpatialLines, SpatialLinesDataFrame, sf, or SpatVector. The base lines to start from (known route).
#' @param conductance Transition Layer. Cost surface (gdistance).
#' @param steps_points Integer. Number of steps (normally in meters) for sampling points along the lines.
#' @param max_speed Numeric. Maximum speed (optional, for travel time calculation).
#' @return A SpatialLinesDataFrame object with the calculated routes.
#' @export
sbr_network <- function(sites, lines, conductance, steps_points = 500, max_speed = 6) {

  # --- Ensure 'sites' is SpatialPointsDataFrame ------------------------------
  if (inherits(sites, "sf")) sites <- as(sites, "Spatial")
  if (inherits(sites, "SpatVector")) sites <- as(sites, "Spatial")

  if (inherits(sites, "SpatialPoints") && !inherits(sites, "SpatialPointsDataFrame")) {
    df <- data.frame(id = seq_along(sites))
    rownames(df) <- df$id
    sites <- sp::SpatialPointsDataFrame(sites, data = df)
  }

  if (!inherits(sites, "SpatialPointsDataFrame")) {
    stop("The 'sites' object must be or convertible to a SpatialPointsDataFrame.")
  }

  if (!"id" %in% names(sites@data)) sites@data$id <- seq_len(nrow(sites))
  if (!"used" %in% names(sites@data)) sites@data$used <- "no"

  # --- Ensure 'lines' is SpatialLinesDataFrame -------------------------------
  if (inherits(lines, "sf")) lines <- as(lines, "Spatial")
  if (inherits(lines, "SpatVector")) lines <- as(lines, "Spatial")

  if (inherits(lines, "SpatialLines") && !inherits(lines, "SpatialLinesDataFrame")) {
    df <- data.frame(ID = seq_along(lines))
    rownames(df) <- df$ID
    lines <- sp::SpatialLinesDataFrame(lines, data = df)
  }

  if (!inherits(lines, "SpatialLinesDataFrame")) {
    stop("The 'lines' object must be or convertible to a SpatialLinesDataFrame.")
  }

  # --- Crop inputs to conductance extent -------------------------------------
  extent.area <- as(raster::extent(raster::raster(conductance)), "SpatialPolygons")
  raster::crs(extent.area) <- raster::crs(raster::raster(conductance))
  sites <- raster::crop(sites, extent.area)
  lines <- raster::crop(lines, extent.area)

  if (nrow(sites) == 0) stop("No valid sites within the conductance extent.")
  if (length(lines) == 0) stop("No valid lines within the conductance extent.")

  # --- Initialize line attributes --------------------------------------------
  required_cols <- c("ID", "dist_m", "time_minutes", "hub_dist_next_site", "time_to_next_site")
  for (col in required_cols) {
    if (!col %in% names(lines@data)) lines@data[[col]] <- NA_real_
  }

  # Convert ID to numeric if needed
  suppressWarnings({
    lines@data$ID <- as.numeric(as.character(lines@data$ID))
  })
  lines@data$ID[is.na(lines@data$ID)] <- seq_len(nrow(lines))

  # Keep originals for hub-distance calculation
  sites_original <- sites
  stop_it <- FALSE

  # --- Iterative linking process --------------------------------------------
  while (!stop_it) {

    # Determine number of sampling points
    total_len_m <- as.numeric(sum(sf::st_length(sf::st_as_sf(lines))))
    numOfPoints <- max(2, round(total_len_m / steps_points))

    way.points <- sp::spsample(lines, n = numOfPoints, type = "regular")
    if (is.null(way.points) || length(way.points) == 0) {
      warning("No waypoints generated. Adjust 'steps_points' or check line geometry.")
      break
    }

    way.points <- sp::SpatialPointsDataFrame(
      way.points,
      data = data.frame(dist = rep(1, length(way.points)))
    )

    # Calculate accumulated cost surface
    acc_way_points <- gdistance::accCost(conductance, way.points) / 60
    sites$min_dist <- raster::extract(acc_way_points, sites)

    # NA-Werte entfernen und eine Warnung ausgeben, falls vorhanden
    initial_start_count <- nrow(sites)
    sites <- sites[!is.na(sites$min_dist), ]
    removed_start_count <- initial_start_count - nrow(sites)
    if (removed_start_count > 0) {
      warning(paste("Removed", removed_start_count, "sites due to NA values in cost surface extraction."))
    }

    # Select next site (unused, minimal cost)
    origin_new <- sites[sites$used == "no" & sites$min_dist == min(sites$min_dist[sites$used == "no"]), , drop = FALSE]

    # Standort für die zukünftige Berechnung ausschließen
    sites$used[sites$min_dist == origin_new$min_dist & sites$used == "no"] <- "yes"

    # Ausschneiden und Berechnen der Distanzen
    if (terra::nrow(way.points) > 1) {
      way.points$dist <- gdistance::costDistance(conductance, way.points, origin_new)
    }

    # Neuen Zielpunkt bestimmen
    idx_target <- which.min(way.points$dist)
    target_new <- way.points[idx_target, , drop = FALSE]
    min_dist <- way.points$dist[idx_target]


    # Compute least-cost path
    route_sl <- gdistance::shortestPath(conductance, origin_new, target_new, output = "SpatialLines")

    # Create new route line with attributes
    id_vals <- suppressWarnings(as.numeric(as.character(lines@data$ID)))
    new_id <- if (any(is.finite(id_vals))) max(id_vals, na.rm = TRUE) + 1 else nrow(lines) + 1

    route_new <- sp::SpatialLinesDataFrame(
      route_sl,
      data = data.frame(
        ID = new_id,
        dist_m = min_dist,
        time_minutes = min_dist / (max_speed * (1000 / 60 / 60)),  # m / (m/min)
        hub_dist_next_site = 0,
        time_to_next_site = min_dist / (max_speed * (1000 / 60 / 60))
      )
    )

    # Append route to existing network
    lines <- raster::bind(lines, route_new)

    # Prepare next iteration
    sites <- sites[sites$used == "no", , drop = FALSE]
    if (nrow(sites) == 0) stop_it <- TRUE
  }

  message("SBR network construction completed successfully.")
  return(lines)
}

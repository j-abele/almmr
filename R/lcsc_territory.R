#' LCSC Territory - Potential settlement area
#'
#' Generates theoretical territories or catchment areas based on a cost surface
#' for one or more locations. For each site the accumulated travel time from the
#' site outward is computed, and all cells reachable within \code{movement_time}
#' form the territory polygon.
#'
#' The accumulated cost is computed as a one-to-all shortest-path problem on a
#' directed, weighted graph (\code{igraph}), replacing the former
#' \code{gdistance::accCost()}. Travel time is measured \emph{outward} from the
#' site (\code{mode = "out"}), matching the semantics of \code{accCost()}. For
#' anisotropic cost functions (e.g. Tobler) the outward time differs from the
#' time required to reach the site.
#'
#' @param dem SpatRaster. Digital elevation model with a projected CRS.
#'   RasterLayer input is converted automatically.
#' @param sites SpatVector or sf. Site locations as points. Reprojected to the
#'   DEM CRS automatically if needed.
#' @param movement_time Numeric. Time frame in minutes.
#' @param max_speed Numeric. Maximum speed of the cost function in km/h
#'   (e.g. 6 for \code{"ToblersHikingFunction"}). Used only to bound the DEM
#'   clip per site, which speeds up the calculation.
#' @param slopeBarrier Logical. Default FALSE. Passed to
#'   \code{create_cost_surface()}. If TRUE, slopes steeper than
#'   \code{slopeBarrierValue} become impassable.
#' @param slopeBarrierValue Numeric. Slope threshold as ratio (e.g. 0.07 = 7\%).
#' @param write_polygon Logical. Default FALSE. If TRUE and \code{wd} is set,
#'   the combined territory polygons are written as a GeoPackage.
#' @param write_raster Logical. Default FALSE. If TRUE and \code{wd} is set, the
#'   accumulated-cost raster (in minutes) is written per site as a GeoTIFF.
#' @param name_index Column name or index holding the site names. If NULL, names
#'   are generated automatically as \code{"LCSC_territory_site_id<i>"}.
#' @param wd Character. Output directory for optional file writing.
#' @param costFunction Character. Cost function passed to
#'   \code{create_cost_surface()}. Designed for \code{"ToblersHikingFunction"}.
#' @param numberOfNeighbors Integer. Neighbours for slope calculation (4, 8, 16).
#'   Passed to \code{create_cost_surface()}.
#' @param numberOfDirections Integer. Neighbours for graph connectivity
#'   (4, 8, 16, 32, 48). Passed to \code{create_cost_surface()}.
#' @return A SpatVector of polygons, one feature per site, with attribute
#'   columns \code{name} and \code{area} (territory size in hectares).
#' @export

lcsc_territory <- function(dem,
                           sites,
                           movement_time,
                           max_speed,
                           slopeBarrier       = FALSE,
                           slopeBarrierValue  = NULL,
                           write_polygon      = FALSE,
                           write_raster       = FALSE,
                           name_index         = NULL,
                           wd                 = NULL,
                           costFunction       = "ToblersHikingFunction",
                           numberOfNeighbors  = 16,
                           numberOfDirections = 16) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  stopifnot(is.numeric(movement_time), is.numeric(max_speed))

  if (!inherits(dem, "SpatRaster")) {
    if (inherits(dem, "RasterLayer")) {
      dem <- terra::rast(dem)
    } else {
      stop("'dem' must be a SpatRaster (or RasterLayer).")
    }
  }
  if (terra::is.lonlat(dem))
    stop("Geographic (lon/lat) CRS detected. ",
         "Please provide a projected DEM (e.g. UTM).")

  if (inherits(sites, "sf")) sites <- terra::vect(sites)
  if (!inherits(sites, "SpatVector"))
    stop("'sites' must be of class SpatVector or sf.")

  # Reproject sites to DEM CRS if needed
  if (!terra::same.crs(sites, dem)) {
    message("'sites' CRS differs from DEM - reprojecting automatically.")
    sites <- terra::project(sites, dem)
  }

  if (is.null(wd) && (write_raster || write_polygon))
    warning("No working directory (wd = ) specified. ",
            "The process continues without saving files.")

  # ---------------------------------------------------------------------------
  # Setup
  # Maximum reachable distance = v_max * t, plus a 3-cell margin for the clip.
  # (s = v * t; used only to bound the DEM clip per site.)
  # ---------------------------------------------------------------------------
  res_dem               <- terra::res(dem)[1]
  max_possible_distance <- (max_speed * 1000 / 3600) * (movement_time * 60) +
    res_dem * 3

  n_sites        <- nrow(sites)
  territory_list <- vector("list", n_sites)

  # ---------------------------------------------------------------------------
  # Per-site loop
  # ---------------------------------------------------------------------------
  for (i in seq_len(n_sites)) {

    site <- sites[i, ]

    # Site name
    name <- if (!is.null(name_index)) {
      paste0("LCSC_territory_", as.data.frame(site)[, name_index], "_", i)
    } else {
      paste0("LCSC_territory_site_id", i)
    }

    # Clip DEM to the site's reachable buffer (performance)
    site_buffer <- terra::buffer(site, max_possible_distance)
    dem_clip    <- terra::crop(dem, site_buffer)

    # Cost surface on the clipped DEM.
    # Eager: accCost is one-to-all, so the full graph on the clip is needed
    # anyway - building it once eagerly is fastest for many small areas.
    cs <- create_cost_surface(
      dem                = dem_clip,
      slopeBarrier       = slopeBarrier,
      slopeBarrierValue  = slopeBarrierValue,
      costFunction       = costFunction,
      numberOfNeighbors  = numberOfNeighbors,
      numberOfDirections = numberOfDirections,
      lazy               = FALSE
    )

    # Build directed weighted graph
    graph_obj <- .build_graph(cs)

    # Locate the site on the graph
    site_cell <- terra::cellFromXY(dem_clip, terra::crds(site))
    site_node <- graph_obj$cell_to_node[site_cell]

    if (is.na(site_node) || site_node == 0) {
      warning(sprintf(
        "Site %d falls on an NA or barrier cell - skipped.", i))
      next
    }

    # accCost equivalent: one-to-all shortest path.
    # Edge weights are travel time in seconds; mode = "out" gives travel time
    # FROM the site outward (matching gdistance::accCost).
    d_sec <- igraph::distances(
      graph_obj$graph,
      v       = site_node,
      mode    = "out",
      weights = igraph::E(graph_obj$graph)$weight
    )[1, ]
    d_min <- d_sec / 60

    # Accumulated-cost raster (minutes)
    v                          <- rep(NA_real_, terra::ncell(dem_clip))
    v[graph_obj$node_to_cell]  <- d_min
    acc_rast                   <- terra::rast(dem_clip)
    terra::values(acc_rast)    <- v
    names(acc_rast)            <- "cost_min"

    # Reachable cells within movement_time (Inf/unreachable -> NA)
    reachable <- terra::ifel(acc_rast < movement_time, 1L, NA)

    if (all(is.na(terra::values(reachable)))) {
      warning(sprintf(
        "Site %d: no cells reachable within %g minutes - skipped.",
        i, movement_time))
      next
    }

    # Polygonise and dissolve into a single territory feature
    poly      <- terra::as.polygons(reachable, dissolve = TRUE)
    poly$name <- name
    poly$area <- as.numeric(terra::expanse(poly, unit = "ha"))
    poly      <- poly[, c("name", "area")]

    territory_list[[i]] <- poly

    # Optional per-site raster output
    if (write_raster && !is.null(wd)) {
      terra::writeRaster(
        acc_rast,
        file.path(wd, paste0(name, "_", movement_time,
                             "_min_LCSC_potential.tif")),
        overwrite = TRUE
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Combine all site territories
  # ---------------------------------------------------------------------------
  territory_list <- Filter(Negate(is.null), territory_list)
  if (length(territory_list) == 0)
    stop("No territories could be computed for any site. ",
         "Check barriers, movement_time, or site locations.")

  LCSC_potential <- do.call(rbind, territory_list)

  # Optional combined polygon output (written after the loop -
  # fixes the previous reference to an undefined object inside the loop)
  if (write_polygon && !is.null(wd)) {
    terra::writeVector(
      LCSC_potential,
      file.path(wd, paste0("sites_LCSC_potential_", movement_time, "min_",
                           slopeBarrierValue, "_perc.gpkg")),
      overwrite = TRUE
    )
  }

  LCSC_potential
}

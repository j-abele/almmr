#' LCSC Territory - Potential settlement area
#'
#' The function generates theoretical territories or catchment areas based on a cost surface for one or more locations.
#' The catchment areas are calculated on the basis of a certain period of time, which, starting from a central location, the movement can take place in all directions depending on the costs.
#'
#' @param dem Raster. Digital Elevation Data.
#' @param sites Sites to analyse as points.
#' @param movement_time Time-Frame in minutes
#' @param max_speed Maximum possible speed of the cost function (e.g. 6 km/h for ToblersHickingFunction. Optional, but good for speeding up the cost_surface calculation)
#' @param slopeBarrier TRUE/FALSE. Default FALSE. If TRUE, slope higher than given slopeBarrierValue get a slope value of NA
#' @param slopeBarrierValue Start slope value for barriers.
#' @param write_polygon Write polygon of LCSC territory per site to local directory (specified with wd = )
#' @param write_raster Write accumulate cost of LCSC territory per site to local directory (specified with wd = )
#' @param name_index Index of a column with the names of the sites. If not specified, a name is automatically generated according to the pattern site_ [consecutive number].
#' @param wd Directory in which the polygons and any raster data can be saved. Optional.
#' @param costFunction Cost-Function. Possible: "ToblersHikingFunction"
#' @return Polygon-Shapefile. Optional: Raster TIFF

lcsc_territory <- function(dem,
                           sites,
                           movement_time,
                           max_speed,
                           slopeBarrier = FALSE,
                           slopeBarrierValue = NULL,
                           write_polygon = FALSE,
                           write_raster=FALSE,
                           name_index=NULL,
                           wd=NULL,
                           epsg = "EPSG:25832",
                           costFunction = "ToblersHikingFunction") {

  # Überprüfung der Eingaben
  stopifnot(is.numeric(movement_time), is.numeric(max_speed))
  if (!inherits(dem, c("RasterLayer", "SpatRaster"))) {
    stop("The 'dem' must be a RasterLayer or SpatRaster.")
  }
  if (is.null(wd) && (write_raster || write_polygon)) {
    warning("No working directory specified (wd = ). The process will continue without saving files locally.")
  }

  # Initialisierung
  sites$area <- 0
  site_acc_poly_list <- list()
  ## Calculat maximum possible distance for given movement-time and add tribble of raster-cell resolution. (s = v x t)
  # max_possible_distance is used for dem clipping before cost-surface calculaiton to speed up process
  max_possible_distance <- ((max_speed * 1000) / 3600) * (movement_time * 60) + (terra::res(rast(dem))[1] * 3)

  for (i in 1:nrow(sites)) {
    #buffer
    site <- sites[i,]
    # Generierung eines Namens: Entweder basierend auf dem name_index oder fortlaufend nummeriert
    if (!is.null(name_index )) {
      name <- paste0("LCSC_territory_",data.frame(site)[, name_index], "_", i)
    } else {
      name <- paste0("LCSC_territory_site_id", i)
    }

    ##datierung <- site$Datierung braucht es das??
    site_buffer <- terra::buffer(site, max_possible_distance)
    #clip
    rast <- dem %>% # in raster ueberfuehren falls es nicht schon ist
      #raster::mask(site_buffer) %>%
      raster::crop(site_buffer) #%>%
      raster::raster()
      #conductance
    conductance <- create_cost_surface(dem=rast,
                                       epsg = epsg,
                                       slopeBarrier=slopeBarrier,
                                       slopeBarrierValue=slopeBarrierValue,
                                       costFunction = costFunction)




    # accCost ACHTUNG: Noch Warnung einbauen, falls slopeBarrier zu NA WErten an site führt
    site_acc <- gdistance::accCost(conductance, as(site, "Spatial")) / 60

    #
    # Raster zu Polygon konvertieren und aggregieren
    site_acc_poly <- raster::rasterToPolygons(site_acc,
                                              fun = function(x) {x < movement_time},
                                              n = 4, na.rm = TRUE,
                                              digits = 2, dissolve = FALSE)
    site_acc_poly <- raster::aggregate(site_acc_poly, fun = mean)
    site_acc_poly$name <- name
    site_acc_poly$area <- raster::area(site_acc_poly) / 10000  # in Hektar
    site_acc_poly_list[[i]] <- site_acc_poly

    # Optionales Schreiben des Rasters
    if (write_raster && !is.null(wd)) {
      raster::writeRaster(site_acc,
                          file.path(wd, paste0(name, "_", movement_time, "_min_LCSC_potential.tiff")),
                          overwrite = TRUE)
    }
    # Optionales Schreiben der Polygon-Datei
    if (write_polygon && !is.null(wd)) {
      sf::st_write(LCSC_potential,
                   file.path(wd, paste0("sites_LCSC_potential_", movement_time, "min_", slopeBarrierValue, "_perc_.shp")),
                   delete_dsn = TRUE)
    }

  }
  #write contour polygons
  LCSC_potential <- do.call(bind, site_acc_poly_list)
  return(LCSC_potential)
}

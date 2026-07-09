# Internal helpers for river / waterbody travel weights.
# Not exported. Called from create_cost_surface() (eager) and .build_graph()
# (lazy) so the water logic lives in exactly one place.

# Oriented along-river position for every river cell.
#
# Internal sub-helper of .apply_water_weights(), which passes in the rivers
# already cropped and projected to `dem`. Not called from anywhere else.
#
# Returns a named numeric vector: names are river cell indices (as character),
# values are the position along the river, oriented so that a HIGHER value means
# further downstream. Direction is derived robustly: local order comes from the
# river geometry (position along the line), the overall downstream orientation
# from the elevation trend over the whole river - so noise in the DEM between
# adjacent cells cannot flip it.
#
# Note: in very flat regions or for very short rivers the elevation trend can be
# weak or ambiguous, so the downstream orientation may be determined
# incorrectly. In flat terrain this matters little in practice, since the up-
# and downstream speeds would typically be chosen closer together there anyway,
# making the direction largely irrelevant.
#
# @param rivers SpatVector (lines), already cropped/projected to `dem`.
# @param dem    SpatRaster.
.river_flow_position <- function(rivers, dem) {

  res_half <- terra::res(dem)[1] / 2
  vals     <- terra::values(dem, mat = FALSE)

  # Group river features: by 'name' if present, otherwise each feature on its own
  groups <- if ("name" %in% names(rivers)) {
    split(seq_len(nrow(rivers)), rivers$name)
  } else {
    as.list(seq_len(nrow(rivers)))
  }

  cell_s <- numeric(0)

  for (idx in groups) {
    grp <- rivers[idx, ]

    # Cells covered by this river group
    grp_r <- terra::rasterize(grp, dem, field = 1L, touches = TRUE)
    cells <- which(!is.na(terra::values(grp_r, mat = FALSE)))
    if (length(cells) < 2L) next

    # Merge the group's lines into a single geometry and densify into ordered
    # points with cumulative along-line distance. st_line_merge() only accepts
    # MULTILINESTRING; a single segment is already a LINESTRING and used as is.
    geom  <- sf::st_union(sf::st_geometry(sf::st_as_sf(grp)))
    line  <- if (inherits(geom, "sfc_MULTILINESTRING")) sf::st_line_merge(geom) else geom
    dense <- sf::st_cast(sf::st_segmentize(line, dfMaxLength = res_half), "POINT")
    dco   <- sf::st_coordinates(dense)[, 1:2, drop = FALSE]
    if (nrow(dco) < 2L) next
    cumd  <- c(0, cumsum(sqrt(rowSums(diff(dco)^2))))

    # Along-line position of each river cell = position of its nearest dense point
    cxy      <- terra::xyFromCell(dem, cells)
    cell_pts <- sf::st_as_sf(data.frame(x = cxy[, 1], y = cxy[, 2]),
                             coords = c("x", "y"), crs = sf::st_crs(line))
    nn <- sf::st_nearest_feature(cell_pts, dense)
    s  <- cumd[nn]

    # Robust downstream orientation: if elevation rises with s, flip s so that
    # larger s is always downstream.
    elev <- vals[cells]
    if (sum(is.finite(elev)) >= 2L && stats::sd(s) > 0) {
      b <- tryCatch(stats::coef(stats::lm(elev ~ s))[[2]], error = function(e) NA)
      if (is.finite(b) && b > 0) s <- max(s) - s
    }

    names(s) <- as.character(cells)
    cell_s   <- c(cell_s, s)
  }

  # On confluences a cell may appear twice; keep the last assignment
  cell_s[!duplicated(names(cell_s), fromLast = TRUE)]
}


# Apply river / waterbody travel weights to an edge list.
#
# @param adj     2-column matrix [from, to] of cell indices in `dem`.
# @param weights Numeric vector, parallel to `adj` rows (travel time in seconds).
# @param dem     SpatRaster (full DEM in eager mode, clipped DEM in lazy mode).
# @param rivers,waterbodies  SpatVector or NULL.
# @param *_speed_kmh  Travel speeds in km/h.
# @return The modified `weights` vector.
.apply_water_weights <- function(adj, weights, dem,
                                 rivers = NULL, waterbodies = NULL,
                                 downstream_speed_kmh  = 12,
                                 upstream_speed_kmh    = 3,
                                 waterbodies_speed_kmh = 5) {

  ms <- function(kmh) kmh * 1000 / 3600  # km/h -> m/s

  # Euclidean distance between the two cells of the given edge rows
  edge_dist <- function(rows) {
    a <- terra::xyFromCell(dem, adj[rows, 1])
    b <- terra::xyFromCell(dem, adj[rows, 2])
    sqrt((b[, 1] - a[, 1])^2 + (b[, 2] - a[, 2])^2)
  }

  # --- Waterbodies first (isotropic); rivers are applied afterwards and win on
  #     any overlapping edge.
  if (!is.null(waterbodies)) {
    wb <- waterbodies
    if (!terra::same.crs(wb, dem)) wb <- terra::project(wb, terra::crs(dem))
    wb <- suppressWarnings(terra::crop(wb, dem))
    if (!is.null(wb) && terra::nrow(wb) > 0) {
      wb_r     <- terra::rasterize(wb, dem, field = 1L, touches = TRUE)
      wb_cells <- which(!is.na(terra::values(wb_r, mat = FALSE)))
      if (length(wb_cells) >= 2L) {
        both <- which(adj[, 1] %in% wb_cells & adj[, 2] %in% wb_cells)
        if (length(both)) weights[both] <- edge_dist(both) / ms(waterbodies_speed_kmh)
      } else {
        warning("Too few waterbody cells in the DEM extent - waterbodies skipped.")
      }
    }
  }

  # --- Rivers (directed by along-river flow position)
  if (!is.null(rivers)) {
    rv <- rivers
    if (!terra::same.crs(rv, dem)) rv <- terra::project(rv, terra::crs(dem))
    rv <- suppressWarnings(terra::crop(rv, dem))
    if (!is.null(rv) && terra::nrow(rv) > 0) {
      cell_s <- .river_flow_position(rv, dem)
      if (length(cell_s) >= 2L) {
        river_cells <- as.integer(names(cell_s))
        both <- which(adj[, 1] %in% river_cells & adj[, 2] %in% river_cells)
        if (length(both)) {
          s_from     <- cell_s[as.character(adj[both, 1])]
          s_to       <- cell_s[as.character(adj[both, 2])]
          downstream <- s_to >= s_from  # ties treated as downstream
          d          <- edge_dist(both)
          weights[both] <- ifelse(downstream,
                                  d / ms(downstream_speed_kmh),
                                  d / ms(upstream_speed_kmh))
        }
      } else {
        warning("Too few river cells in the DEM extent - rivers skipped.")
      }
    }
  }

  weights
}

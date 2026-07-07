#' Site-based (theoretical) route network (sbr_network)
#'
#' Generates a theoretical route network based on site locations and a given
#' route, using a cost surface. Each iteration connects the next unlinked site
#' to the existing network via the least-cost path.
#'
#' Strongly inspired by the least-cost network-to-builder process described by
#' David Waugh (2000, 615), see also Herzog (2013, 239).
#'
#' The former gdistance operators are replaced by graph operations on a single
#' directed, weighted \code{igraph} built once from the cost surface:
#' the accumulated cost from the network to each site and the cost to the
#' nearest waypoint both use \code{igraph::distances()}, and each connecting
#' route uses \code{igraph::shortest_paths()}.
#'
#' Literature:
#' D. Waugh, Geography. An Integrated Approach. Third edition (Cheltenham 2000).
#' I. Herzog, Least-cost Networks. In: G. Earl/T. Sly/A. Chrysanthi/P.
#' Murrieta-Flores/C. Papadopoulos/I. Romanowska/D. Wheatley (Hrsg.),
#' Archaeology in the Digital Era II (Amsterdam 2013), 237-248.
#'
#' @param sites sf or SpatVector. Point locations of the sites. Reprojected to
#'   the cost surface CRS automatically if needed.
#' @param lines sf or SpatVector. Base line(s) of the known route to start from.
#' @param cost_surface almmr_cs object from \code{create_cost_surface(lazy = FALSE)}.
#'   An eager cost surface is required, as network-wide distances are computed
#'   repeatedly on a single graph.
#' @param steps_points Numeric. Spacing in meters for sampling waypoints along
#'   the current network in each iteration.
#' @param max_speed Numeric. Retained for backward compatibility. Travel time is
#'   now derived directly from the cost surface, so this argument is currently
#'   informational only.
#' @return An sf object of LINESTRING geometries: the base line(s) plus one
#'   feature per connecting route, with columns \code{ID}, \code{dist_m},
#'   \code{time_minutes} and \code{connected_site}.
#' @export
sbr_network <- function(sites, lines, cost_surface, steps_points = 500, max_speed = 6) {

  # ---------------------------------------------------------------------------
  # Cost surface
  # ---------------------------------------------------------------------------
  if (!inherits(cost_surface, "almmr_cs"))
    stop("'cost_surface' must be an almmr_cs object from create_cost_surface().")
  if (isTRUE(cost_surface$params$lazy))
    stop("sbr_network() requires an eager cost surface (lazy = FALSE), because ",
         "network-wide distances are computed repeatedly on one graph.")

  dem <- cost_surface$dem

  # ---------------------------------------------------------------------------
  # Normalise sites to SpatVector points in the DEM CRS
  # ---------------------------------------------------------------------------
  if (inherits(sites, "sf")) sites <- terra::vect(sites)
  if (!inherits(sites, "SpatVector"))
    stop("'sites' must be an sf or SpatVector point layer.")
  if (terra::geomtype(sites) != "points")
    stop("'sites' must have point geometries.")
  if (!terra::same.crs(sites, dem)) {
    message("'sites' CRS differs from cost surface - reprojecting automatically.")
    sites <- terra::project(sites, dem)
  }

  # ---------------------------------------------------------------------------
  # Normalise lines to sf LINESTRING in the DEM CRS
  # ---------------------------------------------------------------------------
  if (inherits(lines, "SpatVector")) lines <- sf::st_as_sf(lines)
  lines <- sf::st_as_sf(lines)
  if (is.na(sf::st_crs(lines))) sf::st_crs(lines) <- sf::st_crs(dem)
  if (sf::st_crs(lines) != sf::st_crs(dem))
    lines <- sf::st_transform(lines, sf::st_crs(dem))
  lines <- sf::st_cast(lines, "LINESTRING", warn = FALSE)

  base_n <- nrow(lines)
  if (base_n == 0) stop("'lines' contains no line geometries to start from.")

  # ---------------------------------------------------------------------------
  # Build the directed weighted graph once (cost surface is static)
  # ---------------------------------------------------------------------------
  graph_obj    <- .build_graph(cost_surface)
  g            <- graph_obj$graph
  cell_to_node <- graph_obj$cell_to_node
  node_to_cell <- graph_obj$node_to_cell
  w            <- igraph::E(g)$weight

  # Map point geometries to graph node IDs (NA if off-graph)
  node_of <- function(v) {
    nd <- cell_to_node[terra::cellFromXY(dem, terra::crds(v))]
    nd[is.na(nd) | nd == 0] <- NA_integer_
    nd
  }

  # Shortest-path geometry + travel time (seconds) between two nodes
  path_between <- function(from_node, to_node) {
    res <- igraph::shortest_paths(g, from = from_node, to = to_node,
                                  mode = "out", weights = w, output = "both")
    vp <- as.integer(res$vpath[[1]])
    if (length(vp) < 2) return(NULL)
    ttime  <- sum(w[as.integer(res$epath[[1]])])
    coords <- terra::xyFromCell(dem, node_to_cell[vp])
    list(geometry = sf::st_linestring(coords), travel_time_s = ttime)
  }

  # ---------------------------------------------------------------------------
  # Drop sites that fall outside the graph (NA cells, barriers, out of extent)
  # ---------------------------------------------------------------------------
  site_nodes_all <- node_of(sites)
  outside        <- is.na(site_nodes_all)
  if (any(outside)) {
    warning(sprintf("%d site(s) outside the cost surface or on NA cells were dropped.",
                    sum(outside)))
    sites          <- sites[!outside]
    site_nodes_all <- site_nodes_all[!outside]
  }
  if (nrow(sites) == 0) stop("No sites fall on the cost surface graph.")

  # ---------------------------------------------------------------------------
  # Initialise the network as sf with route attributes
  # ---------------------------------------------------------------------------
  network <- sf::st_sf(
    ID             = seq_len(base_n),
    dist_m         = NA_real_,
    time_minutes   = NA_real_,
    connected_site = NA_integer_,
    geometry       = sf::st_geometry(lines)
  )

  used <- rep(FALSE, nrow(sites))

  # ---------------------------------------------------------------------------
  # Iterative linking
  # ---------------------------------------------------------------------------
  repeat {
    unused <- which(!used)
    if (length(unused) == 0) break

    # Sample waypoints along the current network.
    # Each network feature is a LINESTRING; sample per feature (avoid st_union,
    # which merges them into a MULTILINESTRING that st_line_sample rejects).
    net_lines <- sf::st_geometry(network)
    if (any(sf::st_geometry_type(net_lines) == "MULTILINESTRING"))
      net_lines <- sf::st_cast(net_lines, "LINESTRING")

    wp_sfc <- sf::st_line_sample(net_lines, density = 1 / steps_points,
                                 type = "regular")
    wp_sfc <- sf::st_cast(wp_sfc, "POINT")

    # Fall back to line vertices if the network is too short for the spacing
    if (length(wp_sfc) == 0)
      wp_sfc <- sf::st_cast(net_lines, "POINT")

    if (length(wp_sfc) == 0) {
      warning("No waypoints generated; stopping. Adjust 'steps_points'.")
      break
    }

    wp_nodes <- node_of(terra::vect(wp_sfc))
    wp_nodes <- unique(wp_nodes[!is.na(wp_nodes)])
    if (length(wp_nodes) == 0) {
      warning("No network waypoint falls on the graph; stopping.")
      break
    }

    # Accumulated cost from the network to each unused site (one-to-all)
    site_nodes <- site_nodes_all[unused]
    D          <- igraph::distances(g, v = wp_nodes, to = site_nodes,
                                    mode = "out", weights = w)
    site_min   <- apply(D, 2, min)          # minimum over waypoints per site

    if (all(!is.finite(site_min))) {
      warning("No remaining site is reachable from the network; stopping.")
      break
    }

    sel         <- which.min(site_min)      # position within 'unused'
    origin_glob <- unused[sel]
    origin_node <- site_nodes[sel]

    # Nearest network waypoint to the selected site
    d_to_site   <- igraph::distances(g, v = wp_nodes, to = origin_node,
                                     mode = "out", weights = w)[, 1]
    target_node <- wp_nodes[which.min(d_to_site)]

    # Least-cost route from the site to the nearest waypoint
    route <- path_between(origin_node, target_node)
    used[origin_glob] <- TRUE
    if (is.null(route)) {
      warning(sprintf("Site %d could not be connected to the network; skipped.",
                      origin_glob))
      next
    }

    geom      <- sf::st_sfc(route$geometry, crs = sf::st_crs(dem))
    route_row <- sf::st_sf(
      ID             = max(network$ID, na.rm = TRUE) + 1,
      dist_m         = as.numeric(sf::st_length(geom)),
      time_minutes   = route$travel_time_s / 60,
      connected_site = origin_glob,
      geometry       = geom
    )

    network <- rbind(network, route_row)
  }

  message("SBR network construction completed.")
  network
}

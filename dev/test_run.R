# =============================================================================
# dev/test_run.R  -  interactive sanity check for the migrated almmr functions
#
# Not part of the package build (dev/ is .Rbuildignore'd). Run it line by line
# in RStudio after devtools::load_all() to eyeball every function's result.
# =============================================================================

# devtools::load_all(".")
library(almmr)
library(terra)
library(sf)

r  <- load_dem()
hs <- create_hillshade(r)

# Reproducible non-NA sample points inside the DEM
set.seed(1)
pts <- terra::spatSample(r, 6, method = "random", na.rm = TRUE, as.points = TRUE)

op <- par(no.readonly = TRUE)

# -----------------------------------------------------------------------------
# 1. Cost surface (eager + lazy) -----------------------------------------------
# -----------------------------------------------------------------------------
cs      <- create_cost_surface(r, costFunction = "ToblersHikingFunction")
cs_lazy <- create_cost_surface(r, lazy = TRUE)

stopifnot(inherits(cs, "almmr_cs"), inherits(cs_lazy, "almmr_cs"))
cat("cost surface OK  | eager lazy =", isTRUE(cs$params$lazy),
    "| lazy lazy =", isTRUE(cs_lazy$params$lazy), "\n")

plot(hs, col = grey(0:100 / 100), legend = FALSE, main = "DEM over hillshade")
plot(cs$dem, add = TRUE, alpha = 0.5)

# -----------------------------------------------------------------------------
# 2. compute_lcp: single path, metrics, hierarchical, bidirectional ------------
# -----------------------------------------------------------------------------
path <- compute_lcp(
  r, pts[1], pts[2], cs = cs,
  output = c("path", "travel_time_s", "distance_m", "straight_m", "detour_index")
)
print(sf::st_drop_geometry(path))
stopifnot(nrow(path) == 1, path$detour_index > 0.95)

path_hier <- compute_lcp(r, pts[1], pts[2], cs_params = cs$params,
                         resolutions = c(200, 100))
path_bi   <- compute_lcp(r, pts[1], pts[2], cs = cs, bidirectional = TRUE)
stopifnot(nrow(path_bi) == 2)

plot(hs, col = grey(0:100 / 100), legend = FALSE, main = "compute_lcp")
plot(sf::st_geometry(path),      add = TRUE, col = "#8b0000", lwd = 3)
plot(sf::st_geometry(path_hier), add = TRUE, col = "#1f78b4", lwd = 2, lty = 2)
plot(pts[1:2], add = TRUE, pch = 21, bg = "white")
cat("compute_lcp OK\n")

# -----------------------------------------------------------------------------
# 3. perform_tpla --------------------------------------------------------------
# -----------------------------------------------------------------------------
ctr  <- terra::vect(cbind(mean(as.vector(ext(r))[1:2]),
                          mean(as.vector(ext(r))[3:4])), crs = crs(r))
tpla <- perform_tpla(cs, center_point = ctr, radius_tpla = 4500,
                     number_of_points = 25, sigma_density_calc = 90)
stopifnot(inherits(tpla, "SpatRaster"))

plot(hs, col = grey(0:100 / 100), legend = FALSE, main = "perform_tpla")
plot(terra::mask(tpla, tpla > 0.04, maskvalues = FALSE), add = TRUE, alpha = 0.7)
plot(ctr, add = TRUE, col = "red", pch = 16)
cat("perform_tpla OK\n")

# -----------------------------------------------------------------------------
# 4. lcsc_territory ------------------------------------------------------------
# -----------------------------------------------------------------------------
terr <- lcsc_territory(r, sites = pts[1:4], movement_time = 5, max_speed = 6)
print(as.data.frame(terr))
stopifnot(inherits(terr, "SpatVector"), all(terr$area > 0))

plot(hs, col = grey(0:100 / 100), legend = FALSE, main = "lcsc_territory")
plot(terr, col = hcl.colors(nrow(terr), "Spectral", alpha = 0.5), add = TRUE)
plot(pts[1:4], add = TRUE, pch = 16)
cat("lcsc_territory OK\n")

# -----------------------------------------------------------------------------
# 5. perform_tpla_line_based ---------------------------------------------------
# -----------------------------------------------------------------------------
frac <- function(fx, fy) cbind(as.vector(ext(r))[1] + fx * diff(as.vector(ext(r))[1:2]),
                               as.vector(ext(r))[3] + fy * diff(as.vector(ext(r))[3:4]))
l1 <- sf::st_sfc(sf::st_linestring(rbind(frac(.3, .35), frac(.3, .55))), crs = sf::st_crs(crs(r)))
l2 <- sf::st_sfc(sf::st_linestring(rbind(frac(.65, .35), frac(.65, .55))), crs = sf::st_crs(crs(r)))
tpla_lb <- perform_tpla_line_based(cs, l1, l2, number_of_points = 12,
                                   sigma_density_calc = 90)
stopifnot(inherits(tpla_lb, "SpatRaster"))

plot(hs, col = grey(0:100 / 100), legend = FALSE, main = "perform_tpla_line_based")
plot(terra::mask(tpla_lb, tpla_lb > 0.04, maskvalues = FALSE), add = TRUE, alpha = 0.7)
plot(l1, add = TRUE, col = "blue", lwd = 2); plot(l2, add = TRUE, col = "blue", lwd = 2)
cat("perform_tpla_line_based OK\n")

# -----------------------------------------------------------------------------
# 6. sbr_network (slower) ------------------------------------------------------
# -----------------------------------------------------------------------------
ref <- compute_lcp(r, pts[1], pts[2], cs = cs)
sbr <- sbr_network(sites = pts[3:6], lines = ref, cost_surface = cs,
                   steps_points = 150)
print(sf::st_drop_geometry(sbr))
stopifnot(inherits(sbr, "sf"), nrow(sbr) >= 2)

plot(hs, col = grey(0:100 / 100), legend = FALSE, main = "sbr_network")
plot(sf::st_geometry(sbr), add = TRUE, col = "darkgreen", lwd = 2)
plot(pts[3:6], add = TRUE, pch = 16)
cat("sbr_network OK\n")

par(op)
cat("\nAll functions ran without error. Inspect the plots for plausibility.\n")

# =============================================================================
# OPTIONAL: compare new (branch) vs old (main) results
# -----------------------------------------------------------------------------
# You cannot load two versions of 'almmr' at once, so compare across a git
# checkout. The two branches also use different APIs, so save comparable,
# API-neutral outputs (geometries / rasters) on each side, then diff.
#
# On main (old gdistance code): run the old equivalent of a single LCP, then
#   saveRDS(sf::st_coordinates(<old_path_as_sf>), "dev/lcp_main.rds")
# On this branch:
#   saveRDS(sf::st_coordinates(path), "dev/lcp_branch.rds")
# Then, in either session:
#   a <- readRDS("dev/lcp_main.rds"); b <- readRDS("dev/lcp_branch.rds")
#   # visual overlay or summary of coordinate differences
#   plot(a[, 1:2], type = "l"); lines(b[, 1:2], col = "red")
#
# Expect the paths to be very close (identical cost model, new solver). Small
# differences are normal where several cells share an equal least cost.
# =============================================================================


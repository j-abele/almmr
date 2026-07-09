# Helpers assume the shared helper-almmr.R (get_dem, valid_points, skip_without_dem)

# A straight river line down the vertical centre of the DEM
make_river <- function(r) {
  e  <- as.vector(terra::ext(r))
  cx <- mean(e[1:2])
  terra::vect(sf::st_sf(
    name     = "r1",
    geometry = sf::st_sfc(
      sf::st_linestring(rbind(
        c(cx, e[["ymin"]] + 0.2 * (e[["ymax"]] - e[["ymin"]])),
        c(cx, e[["ymin"]] + 0.8 * (e[["ymax"]] - e[["ymin"]]))
      )),
      crs = sf::st_crs(terra::crs(r))
    )
  ))
}

# A rectangular waterbody near the DEM centre
make_waterbody <- function(r) {
  e  <- as.vector(terra::ext(r))
  x0 <- e[["xmin"]] + 0.45 * (e[["xmax"]] - e[["xmin"]])
  x1 <- e[["xmin"]] + 0.55 * (e[["xmax"]] - e[["xmin"]])
  y0 <- e[["ymin"]] + 0.45 * (e[["ymax"]] - e[["ymin"]])
  y1 <- e[["ymin"]] + 0.55 * (e[["ymax"]] - e[["ymin"]])
  terra::vect(sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0)))),
      crs = sf::st_crs(terra::crs(r))
    )
  ))
}

test_that("create_cost_surface accepts rivers and changes weights (eager)", {
  skip_without_dem()
  r     <- get_dem()
  base  <- almmr::create_cost_surface(r)
  water <- almmr::create_cost_surface(r, rivers = make_river(r))

  expect_s3_class(water, "almmr_cs")
  expect_equal(nrow(water$adj), nrow(base$adj))                 # same edge set
  expect_false(isTRUE(all.equal(water$weights, base$weights)))  # some weights changed
})

test_that("create_cost_surface accepts waterbodies and changes weights (eager)", {
  skip_without_dem()
  r     <- get_dem()
  base  <- almmr::create_cost_surface(r)
  water <- almmr::create_cost_surface(r, waterbodies = make_waterbody(r))

  expect_s3_class(water, "almmr_cs")
  expect_false(isTRUE(all.equal(water$weights, base$weights)))
})

test_that("river cost surface also works in lazy mode", {
  skip_without_dem()
  r  <- get_dem()
  cs <- almmr::create_cost_surface(r, rivers = make_river(r), lazy = TRUE)

  expect_true(isTRUE(cs$params$lazy))
  expect_false(is.null(cs$params$rivers))

  p    <- valid_points(r, 2, seed = 41)
  path <- almmr::compute_lcp(r, p[1], p[2], cs_params = cs$params)
  expect_s3_class(path, "sf")
})

test_that("river flow direction is cheaper downstream than upstream", {
  # Synthetic DEM tilted so elevation decreases with x (west -> east downhill)
  r <- terra::rast(nrows = 5, ncols = 40, xmin = 0, xmax = 4000,
                   ymin = 0, ymax = 500, crs = "EPSG:25832")
  xcell <- terra::xyFromCell(r, seq_len(terra::ncell(r)))[, 1]
  terra::values(r) <- 200 - (xcell / 4000) * 180

  ymid  <- 250
  river <- terra::vect(sf::st_sf(
    name = "r1",
    geometry = sf::st_sfc(
      sf::st_linestring(rbind(c(100, ymid), c(3900, ymid))),
      crs = sf::st_crs(terra::crs(r))
    )
  ))

  cs <- almmr::create_cost_surface(r, rivers = river,
                                   numberOfNeighbors = 4, numberOfDirections = 4)

  # Identify horizontal river edges along the middle row and compare directions
  rr    <- terra::rasterize(river, r, field = 1L, touches = TRUE)
  rcell <- which(!is.na(terra::values(rr, mat = FALSE)))
  is_r  <- cs$adj[, 1] %in% rcell & cs$adj[, 2] %in% rcell
  xf    <- terra::xyFromCell(r, cs$adj[is_r, 1])[, 1]
  xt    <- terra::xyFromCell(r, cs$adj[is_r, 2])[, 1]
  w     <- cs$weights[is_r]

  down <- w[xt > xf]   # moving east = downhill = downstream (fast)
  up   <- w[xt < xf]   # moving west = upstream (slow)

  expect_gt(length(down), 0)
  expect_gt(length(up), 0)
  expect_lt(mean(down), mean(up))
})

test_that("invalid water inputs are rejected", {
  skip_without_dem()
  r <- get_dem()
  expect_error(almmr::create_cost_surface(r, rivers = "not a spatvector"))
  expect_error(almmr::create_cost_surface(r, waterbodies = 42))
})

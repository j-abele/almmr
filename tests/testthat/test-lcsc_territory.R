test_that("lcsc_territory returns a SpatVector of territories", {
  skip_without_dem()
  r     <- get_dem()
  sites <- valid_points(r, 3, seed = 21)

  terr <- almmr::lcsc_territory(
    r,
    sites         = sites,
    movement_time = 5,
    max_speed     = 6
  )

  expect_s4_class(terr, "SpatVector")
  expect_identical(terra::geomtype(terr), "polygons")
  expect_true(all(c("name", "area") %in% names(terr)))
  expect_true(all(terr$area > 0))
})

test_that("lcsc_territory accepts sf input for sites", {
  skip_without_dem()
  r        <- get_dem()
  sites_sf <- sf::st_as_sf(valid_points(r, 2, seed = 22))

  terr <- almmr::lcsc_territory(
    r,
    sites         = sites_sf,
    movement_time = 5,
    max_speed     = 6
  )

  expect_s4_class(terr, "SpatVector")
  expect_gte(nrow(terr), 1)
})

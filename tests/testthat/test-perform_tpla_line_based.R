test_that("perform_tpla_line_based returns a density SpatRaster", {
  skip_without_dem()
  r  <- get_dem()
  cs <- almmr::create_cost_surface(r)

  crs_r <- sf::st_crs(terra::crs(r))
  l1 <- sf::st_sfc(
    sf::st_linestring(rbind(xy_at(r, 0.30, 0.35), xy_at(r, 0.30, 0.55))),
    crs = crs_r
  )
  l2 <- sf::st_sfc(
    sf::st_linestring(rbind(xy_at(r, 0.65, 0.35), xy_at(r, 0.65, 0.55))),
    crs = crs_r
  )

  dens <- suppressWarnings(
    almmr::perform_tpla_line_based(
      cs,
      first_line         = l1,
      second_line        = l2,
      number_of_points   = 6,
      sigma_density_calc = 90
    )
  )

  expect_s4_class(dens, "SpatRaster")
})

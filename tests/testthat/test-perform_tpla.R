test_that("perform_tpla returns a density SpatRaster", {
  skip_without_dem()
  r   <- get_dem()
  cs  <- almmr::create_cost_surface(r)
  ctr <- terra::vect(xy_at(r, 0.5, 0.5), crs = terra::crs(r))

  dens <- suppressWarnings(
    almmr::perform_tpla(
      cs,
      center_point       = ctr,
      radius_tpla        = 1500,
      number_of_points   = 10,
      sigma_density_calc = 90
    )
  )

  expect_s4_class(dens, "SpatRaster")
})

test_that("perform_tpla can additionally return the least-cost lines", {
  skip_without_dem()
  r   <- get_dem()
  cs  <- almmr::create_cost_surface(r)
  ctr <- terra::vect(xy_at(r, 0.5, 0.5), crs = terra::crs(r))

  res <- suppressWarnings(
    almmr::perform_tpla(
      cs,
      center_point       = ctr,
      radius_tpla        = 1500,
      number_of_points   = 10,
      sigma_density_calc = 90,
      keep_lines         = TRUE
    )
  )

  expect_type(res, "list")
  expect_true(all(c("density", "lines") %in% names(res)))
  expect_s4_class(res$density, "SpatRaster")
  expect_s3_class(res$lines, "sf")
})

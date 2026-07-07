test_that("compute_lcp returns a single LINESTRING path (eager cost surface)", {
  skip_without_dem()
  r  <- get_dem()
  cs <- almmr::create_cost_surface(r)
  p  <- valid_points(r, 2, seed = 11)

  path <- almmr::compute_lcp(r, p[1], p[2], cs = cs)

  expect_s3_class(path, "sf")
  expect_equal(nrow(path), 1)
  expect_identical(as.character(sf::st_geometry_type(path)), "LINESTRING")
})

test_that("compute_lcp works in lazy mode via cs_params", {
  skip_without_dem()
  r  <- get_dem()
  cs <- almmr::create_cost_surface(r, lazy = TRUE)
  p  <- valid_points(r, 2, seed = 11)

  path <- almmr::compute_lcp(r, p[1], p[2], cs_params = cs$params)

  expect_s3_class(path, "sf")
  expect_equal(nrow(path), 1)
})

test_that("compute_lcp returns requested metrics with sensible values", {
  skip_without_dem()
  r  <- get_dem()
  cs <- almmr::create_cost_surface(r)
  p  <- valid_points(r, 2, seed = 11)

  path <- almmr::compute_lcp(
    r, p[1], p[2], cs = cs,
    output = c("path", "travel_time_s", "distance_m", "straight_m", "detour_index")
  )

  expect_true(all(c("travel_time_s", "distance_m", "straight_m", "detour_index")
                  %in% names(path)))
  expect_gt(path$travel_time_s, 0)
  expect_gt(path$distance_m, 0)
  # The least-cost path is never meaningfully shorter than the straight line
  expect_gt(path$detour_index, 0.95)
})

test_that("compute_lcp bidirectional returns two directed features", {
  skip_without_dem()
  r  <- get_dem()
  cs <- almmr::create_cost_surface(r)
  p  <- valid_points(r, 2, seed = 11)

  path <- almmr::compute_lcp(r, p[1], p[2], cs = cs, bidirectional = TRUE)

  expect_equal(nrow(path), 2)
  expect_true("direction" %in% names(path))
  expect_setequal(path$direction, c("there", "back"))
})

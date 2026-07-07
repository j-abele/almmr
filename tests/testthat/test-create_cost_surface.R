test_that("create_cost_surface (eager) returns an almmr_cs object", {
  skip_without_dem()
  r  <- get_dem()
  cs <- almmr::create_cost_surface(r, costFunction = "ToblersHikingFunction")

  expect_s3_class(cs, "almmr_cs")
  expect_s4_class(cs$dem, "SpatRaster")
  expect_false(isTRUE(cs$params$lazy))
})

test_that("create_cost_surface supports lazy mode", {
  skip_without_dem()
  r  <- get_dem()
  cs <- almmr::create_cost_surface(r, lazy = TRUE)

  expect_s3_class(cs, "almmr_cs")
  expect_true(isTRUE(cs$params$lazy))
})

test_that("create_cost_surface validates neighbour/direction arguments", {
  skip_without_dem()
  r <- get_dem()

  expect_error(almmr::create_cost_surface(r, numberOfDirections = 7))
  expect_error(almmr::create_cost_surface(r, numberOfNeighbors = 5))
})

test_that("create_cost_surface accepts alternative cost functions", {
  skip_without_dem()
  r <- get_dem()

  expect_s3_class(
    almmr::create_cost_surface(r, costFunction = "Irmischer-Clarke's"),
    "almmr_cs"
  )
})

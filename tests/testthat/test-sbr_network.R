test_that("sbr_network builds an sf route network (reduced run)", {
  skip_without_dem()
  skip_on_cran()  # iterative graph search - keep off CRAN for runtime

  r   <- get_dem()
  cs  <- almmr::create_cost_surface(r)          # eager cost surface required
  pts <- valid_points(r, 4, seed = 31)

  # Reference line between two of the points (replaces gdistance::shortestPath)
  ref <- almmr::compute_lcp(r, pts[1], pts[2], cs = cs)

  net <- almmr::sbr_network(
    sites        = pts[3:4],
    lines        = ref,
    cost_surface = cs,
    steps_points = 200
  )

  expect_s3_class(net, "sf")
  expect_true(all(c("ID", "dist_m", "time_minutes", "connected_site")
                  %in% names(net)))
  # base line + at least one connecting route
  expect_gte(nrow(net), 2)
})

test_that("sbr_network rejects a lazy cost surface", {
  skip_without_dem()
  r   <- get_dem()
  cs  <- almmr::create_cost_surface(r, lazy = TRUE)
  pts <- valid_points(r, 3, seed = 32)
  ref <- almmr::compute_lcp(r, pts[1], pts[2], cs_params = cs$params)

  expect_error(
    almmr::sbr_network(sites = pts[3], lines = ref, cost_surface = cs),
    regexp = "eager"
  )
})

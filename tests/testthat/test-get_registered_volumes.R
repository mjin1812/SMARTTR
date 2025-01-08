test_that("calculating registered volumes for a slice object works", {

  skip_if_not_installed("wholebrain")
  path <- file.path(test_object_path("slice"), "mapped_slice.RDS")
  s <- readRDS(file = path)
  s$volumes <- NULL
  expect_error(get_registered_volumes(s),
               regexp = NA)
  s$volumes <- NULL
  expect_error(get_registered_volumes(s, simplify_regions = FALSE),
               regexp = NA)
})

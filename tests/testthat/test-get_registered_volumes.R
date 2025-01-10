test_that("calculating registered volumes for a mouse object works", {

  skip_if_not_installed("wholebrain")
  skip_on_cran()
  skip_on_ci()

  m <- make_test_mapped_mouse_light()
  m$slices$`1_4`$volumes <- NULL

  expect_error(get_registered_volumes(m,
                       slice_ID = "1_4"),
             regexp = NA)

  expect_error(get_registered_volumes(m, slice_ID = "1_4", simplify_regions = FALSE),
  regexp = NA)

})

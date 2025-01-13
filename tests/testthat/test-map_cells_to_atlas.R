test_that("mapping a cells to atlas works on mouse and slice objects", {

  skip_if_not_installed("wholebrain")
  skip_on_cran()
  skip_on_ci()
  skip("This check takes too long. Manually comment to test forward warping")
  m <- make_test_mapped_mouse_light()
  m$slices$`1_4`$forward_warped_data <- NULL

   expect_error(m <- map_cells_to_atlas(m,
                      slice_ID = "1_4",
                      clean =  TRUE,
                      display = FALSE,
                      channels = NULL),
                regexp = NA)
})






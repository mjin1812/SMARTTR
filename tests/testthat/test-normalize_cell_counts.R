test_that("normalizing cell counts works on a standard mouse slice", {

  skip("Skipping")  # uncomment when testing this function
  skip_on_ci()
  skip_on_cran()

  path <- file.path(test_object_path("mouse"), "mapped_mouse.RDS")
  m <- readRDS(file = path)

  m <- get_cell_table(m)
  expect_error(normalize_cell_counts(m,
                                    combine_hemispheres = TRUE,
                                    simplify_regions = TRUE,
                                    split_hipp_DV = TRUE),
               regexp = NA)

  m <- get_cell_table(m)
  expect_error(normalize_cell_counts(m,
                                     combine_hemispheres = TRUE,
                                     simplify_regions = FALSE,
                                     split_hipp_DV = FALSE),
               regexp = NA)

  m <- get_cell_table(m)
  expect_error(normalize_cell_counts(m,
                                     combine_hemispheres = FALSE,
                                     simplify_regions = FALSE,
                                     split_hipp_DV = FALSE),
               regexp = NA)

})

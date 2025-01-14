test_that("Check enough_mice_per_group for testing", {

  skip_on_cran()
  skip_on_ci()

  e <- make_test_experiment_light()
  expect_message(e <- enough_mice_per_group(e, by = c("group"), min_n = 5, remove = TRUE, log = TRUE),
                 regexp = "There were regions below the minimum n:")

  file_path <- file.path(attr(e, "info")$output_path, "regions_below_N_thresh.csv")
  expect_true(file.exists(file_path))
})

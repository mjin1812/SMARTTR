test_that("Flexible combining of cell counts", {
  skip_on_cran()
  skip_on_ci()

  path <- file.path(test_object_path("experiment"), "experiment.RDS")
  e <- readRDS(file = path)

  attr(e, "info")$channels  <- c("cfos", "eyfp", "colabel")
  expect_error(e <- combine_cell_counts(e, by = c("group", "sex", "age")),
               regexp = NA)
  e <- combine_cell_counts(e, by = c("group", "sex", "age"))

  expect_snapshot(e$combined_counts_per_slice_list$cfos %>% head())
  expect_snapshot(e$combined_normalized_counts$cfos %>% head())
})




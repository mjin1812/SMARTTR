test_that("finding outliers works", {
  path <- file.path(test_object_path("experiment"), "experiment.RDS")
  e <- readRDS(file = path)
  expect_error(normalize_colabel_counts(e, denominator_channel = "eyfp"),
               regexp = NA)
  expect_error(normalize_colabel_counts(e, denominator_channel = "cfos"),
               regexp = NA)
})

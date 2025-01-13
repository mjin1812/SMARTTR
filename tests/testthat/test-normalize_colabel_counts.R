test_that("normalizing colabelled counts works", {
  skip_on_cran()
  skip_on_ci()

  e <- make_test_experiment_light()

  e$combined_normalized_counts$colabel_over_cfos <- NULL
  e$combined_normalized_counts$colabel_over_eyfp <- NULL

  e <- normalize_colabel_counts(e, denominator_channel = "eyfp")
  e <- normalize_colabel_counts(e, denominator_channel = "cfos")
  expect_setequal(names(e$combined_normalized_counts), c("eyfp", "cfos", "colabel", "colabel_over_eyfp", "colabel_over_cfos"))
  expect_snapshot(e$combined_normalized_counts$colabel_over_eyfp %>% head())
  expect_snapshot(e$combined_normalized_counts$colabel_over_cfos %>% head())
})

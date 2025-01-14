
test_that("check that the values vector should match the by vector in get_correlations", {
  skip_on_cran()
  skip_on_ci()

  e <- make_test_experiment_light()
  e$correlation_list <- NULL

  expect_warning(get_correlations(e, by = "group", values = c("Context", "Shock"),
                                  channels = "cfos", , method = "pearson"),
                 regexp = "Your values vector")
})


test_that("spearman and pearson options work for get_correlation",{
  skip_on_cran()
  skip_on_ci()
  e <- make_test_experiment_light()
  e$correlation_list <- NULL

  suppressWarnings(e <- get_correlations(e, by = "group", values = c("Context"),
                        channels = "cfos", method = "spearman"))

  expect_snapshot(e$correlation_list$Context$cfos$r)

  e$correlation_list <- NULL
  suppressWarnings(e <- get_correlations(e, by = "group", values = c("Shock"),
                        channels = "cfos", , method = "pearson"))
  expect_snapshot(e$correlation_list$Context$cfos$r)

})

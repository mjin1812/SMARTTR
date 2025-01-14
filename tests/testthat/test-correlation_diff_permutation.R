test_that("permutation and saving works as expected", {
  skip_on_cran()
  skip_on_ci()
  e <- make_test_experiment_light()

  expect_warning(e <- correlation_diff_permutation(e,
                               correlation_list_name_1 = "Context",
                               correlation_list_name_2 = "Shock", channels = "cfos", n_shuffle = 100,
                               p_adjust_method = "none", alpha = 0.1), regexp = "NaNs produced")
  expect_snapshot(e$permutation_p_matrix$Context_vs_Shock$cfos$test_statistic)
  expect_silent(export_permutation_results(e, channels = "cfos", filter_significant = TRUE))

})

test_that("saving to output directory works", {

  e <- make_test_experiment_light()
  out_path <- attr(e, "info")$output_path
  expect_message(get_percent_colabel(e, "group", colabel_channel = "colabel",
                      channel = "eyfp", save_table = TRUE), regexp = "Saved colabel count percentages at location:")

  file_path <- file.path(out_path, "tables/colabel_percentage_eyfp_individual.csv")
  expect_snapshot_file(file_path)
})


test_that("entering non-existent grouping variable throws error", {
  e <- make_test_experiment_light()
  out_path <- attr(e, "info")$output_path

  expect_error(get_percent_colabel(e, "sex", colabel_channel = "colabel",
                                  channel = "eyfp"), regexp = "doesn't exist")

})

test_that("You can analyze subsets of rois", {

  e <- make_test_experiment_light()
  out_path <- attr(e, "info")$output_path

  expect_message(e <- get_percent_colabel(e, "group", colabel_channel = "colabel",
                                     channel = "eyfp", save_table = FALSE, rois = "HPF"),
                 regexp = "Checking also for child regions")
  expect_snapshot( e$colabel_percent$eyfp$individual %>% head())
})

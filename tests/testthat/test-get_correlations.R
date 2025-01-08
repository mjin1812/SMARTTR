test_that("getting correlations works", {

  path <- file.path(test_object_path("experiment"), "experiment.RDS")
  e <- readRDS(file = path)

  expect_error(e <- get_correlations(e,
                       by = c("group"),
                       values = c("Context"),
                       channels = c("eyfp"),
                       p_adjust_method = "none",
                       alpha = 0.01),
               regexp = NA)

  expect_warning(e <- get_correlations(e,
                           by = c("group"),
                           values = c("Shock"),
                           channels = c("eyfp"),
                           p_adjust_method = "none",
                           alpha = 0.01),
               regexp = "NaNs produced")
})

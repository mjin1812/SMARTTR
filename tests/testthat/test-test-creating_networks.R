test_that("Test different parameters to create single networks", {
  skip_on_cran()
  skip_on_ci()
  e <- make_test_experiment_light()

  expect_silent(e <- create_networks(e,
                       correlation_list_name = "Shock",
                       channels = "eyfp",
                       alpha = 0.05))

  expect_snapshot(e$networks$Shock$eyfp %>% tidygraph::activate(nodes) %>% tibble::as_tibble())
})

test_that("Test creating joined networks works", {
  skip_on_cran()
  skip_on_ci()
  e <- make_test_experiment_light()
  expect_error(e <- create_joined_networks(e,
                                             correlation_list_names = c("Context", "Shock"),
                                             channels = "eyfp",
                                             ontology = "allen",
                                             alpha =  0.01,
                                             pearson_thresh = 0.9,
                                             proportional_thresh = NULL,
                                             alpha2 = NULL,
                                             pearson_thresh2 = NULL,
                                             proportional_thresh2 = NULL,
                                             export_overlapping_edges = FALSE), regexp = NA)
})


test_that("proportional thresholding works", {
  skip_on_cran()
  skip_on_ci()
  e <- make_test_experiment_light()
  expect_error(e <- create_joined_networks(e,
                              correlation_list_names = c("Context", "Shock"),
                              channels = "eyfp",
                              ontology = "allen",
                              alpha =  1,
                              pearson_thresh = 0,
                              proportional_thresh = 0.001,
                              alpha2 = NULL,
                              pearson_thresh2 = NULL,
                              proportional_thresh2 = NULL,
                              export_overlapping_edges = FALSE,
                              filter_isolates =  TRUE), regexp = NA)
})



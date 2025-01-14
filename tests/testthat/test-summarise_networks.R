test_that("summarizing networks works", {
  skip_on_cran()
  skip_on_ci()

  e <- make_analyzed_experiment()
  e$networks$Context_Shock <- NULL
  e$networks_summaries <- NULL

  expect_error(e <-  summarise_networks(e,
                           network_names = c("Context", "Shock"),
                           channels = "eyfp",
                           save_stats = FALSE,
                           save_degree_distribution = FALSE,
                           save_betweenness_distribution = FALSE,
                           save_efficiency_distribution = FALSE),
               regexp = NA)

  expect_snapshot(e$networks_summaries$eyfp$networks_stats)
})

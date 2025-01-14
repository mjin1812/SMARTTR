test_that("test that rewire_network() and summarizing networks works as expected", {
  skip_on_cran()
  skip_on_ci()
  e <- make_analyzed_experiment()


  expect_error(
    summarytable <- rewire_network(e,
                                   network_name = "Context",
                                   channels = "eyfp",
                                   n_rewires = 1000,
                                   n_networks = 10,
                                   return_graphs = FALSE,
                                   seed = 5),
    regexp = NA
    )

  # expect_snapshot(summarytable$eyfp)

  expect_error(summary_list <- summarize_null_networks(list("Context"=summarytable),
                                         network_names = "Context",
                                         channel = "eyfp"),
                regexp = NA)

  # expect_snapshot(summary_list$global_summary)
})

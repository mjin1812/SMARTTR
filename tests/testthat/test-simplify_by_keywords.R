test_that("Test that simplifying keywords works for all ontologys", {
  # Test Allen Atlas
  df <- dplyr::tibble(acronym = c("LGd", "GU4", "dCA1so" ),
                      name = c("Dorsal part of the lateral geniculate complex",
                               "Gustatory areas, layer 4", "Field CA1, stratum oriens (dorsal)"))
  df2 <- dplyr::tibble(acronym = c("LGd", "GU", "dCA1" ),
                      name = c("Dorsal part of the lateral geniculate complex",
                               "Gustatory areas", "Field CA1 (dorsal)"))
  expect_equal(simplify_by_keywords(df), df2)

  # Test Kim Unified Atlas
  df <- dplyr::tibble(acronym = c("LPBD", "MPBE", "CPre"),
                      name = c("Lateral parabrachial nucleus, dorsal part",
                                "Medial parabrachial nucleus, external part",
                                "Caudoputamen- rostral extreme"))
  df2 <- dplyr::tibble(acronym = c("LPB", "MPB", "CP"),
                      name = c("Lateral parabrachial nucleus",
                                "Medial parabrachial nucleus",
                                "Caudate Putamen"))
  expect_equal(simplify_by_keywords(df,
                                    keywords = c("dorsal part", "external part", "Caudoputamen-"),
                                    ontology = "unified"),
               df2)
})




test_that("Test that simplifying vector of acronyms by keywords works for all ontologys", {
  # Test Allen Atlas
  df <- dplyr::tibble(name = c("Dorsal part of the lateral geniculate complex",
                               "Gustatory areas", "Field CA1 (dorsal)"), acronym = c("LGd", "GU", "dCA1" ))
  acronyms <- c("LGd", "GU4", "dCA1so" )
  expect_equal(simplify_vec_by_keywords(acronyms), df)

  # Test Kim Unified Atlas
  df <- dplyr::tibble(name = c("Lateral parabrachial nucleus",
                               "Medial parabrachial nucleus",
                               "Caudate Putamen"), acronym = c("LPB", "MPB", "CP"))
  acronyms <- c("LPBD", "MPBE", "CPre")
  expect_equal(simplify_vec_by_keywords(acronyms,
                                        keywords = c("dorsal part", "external part", "Caudoputamen-"),
                                        ontology = "unified"), df)
})





test_that("test simplify by keywords for experiment object", {
    skip_on_ci()
    skip_on_cran()

    simplify_keywords <- readRDS(file.path(test_object_path("external"), "simplify_keywords.RDS"))
    e <- experiment(experiment_name = "external unified dataset",
                    experimenters = c(""),
                    channels = "cfos",
                    experiment_groups = c("Control", "Isoflurane"),
                    sex_groups = "male",
                    output_path = tempdir())
    path <- file.path(test_object_path("external"), "reformatted_ontologychecked.csv")
    e <- import_mapped_datasets(e, normalized_count_paths = path, show_col_types = FALSE)
    expect_error(e <- simplify_cell_count(e, ontology = "unified", simplify_keywords = simplify_keywords),
                 regexp = NA)

})



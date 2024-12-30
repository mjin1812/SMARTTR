test_that("Test that simplifying vector of acronyms by keywords works for all ontologys", {
  # Test Allen Atlas
  df <- dplyr::tibble(name = c("Dorsal part of the lateral geniculate complex",
                                "Gustatory areas", "Field CA1 (dorsal)"), acronym = c("LGd", "GU", "dCA1" ))
  acronyms <- c("LGd", "GU4", "dCA1so" )
  expect_equal(simplify_vec_by_keywords(acronyms), df)

  # Test Kim Unified Atlas
  df <- dplyr::tibble(name = c("Lateral parabrachial nucleus",
                                "Medial parabrachial nucleus",
                                "Caudat/e Putamen"), acronym = c("LPB", "MPB", "CP"))
  acronyms <- c("LPBD", "MPBE", "CPre")
  expect_equal(simplify_vec_by_keywords(acronyms,
                                        keywords = c("dorsal part", "external part", "Caudoputamen-"),
                                        ontology = "unified"), df)
})

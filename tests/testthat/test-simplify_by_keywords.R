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

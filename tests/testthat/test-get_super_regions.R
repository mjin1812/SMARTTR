test_that("Finding parent super regions works for all ontologyies", {
  # Testing Allen Atlas
  out <- get.super.regions(acronym = c("DG-sg", "ACA5", "CLA", "ACB", "BLA", "PERI" , "PB"))
  expect_equal(out, c("HPF", "Isocortex", "CTXsp", "CNU", "CTXsp", "Isocortex", "HB"))

  # Testing Unified Atlas
  out <- c("CNU", "TH", "HY", "CNU", "CTXsp")
  expect_equal(get.super.regions(acronym = c("CPi", "MHb", "ZI", "Acb", "BL"), ontology = "unified"),
               out)
})

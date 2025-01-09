test_that("normalizing colabelled counts works", {

  skip_on_cran()
  skip_on_ci()

  path <- file.path(test_object_path("experiment"), "experiment.RDS")
  e <- readRDS(file = path)

})

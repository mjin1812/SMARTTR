test_that("import_mapped_datasets works", {
  skip_on_ci()
  skip_on_cran()

  e <- experiment(experiment_name = "external unified dataset",
                  experimenters = c(""),
                  channels = "cfos",
                  experiment_groups = c("Control", "Isoflurane"),
                  sex_groups = "male",
                  output_path = tempdir())
  path <- file.path(test_object_path("external"), "reformatted_ontologyunchecked.csv")

  expect_silent(import_mapped_datasets(e, normalized_count_paths = c("cfos" = path), show_col_types = FALSE))
  expect_silent(import_mapped_datasets(e, normalized_count_paths = path, show_col_types = FALSE))
  expect_error(import_mapped_datasets(e, normalized_count_paths = c("eyfp" = path), show_col_types = FALSE),
               regexp = "does not match")
})





test_that("Checking ontology coding works",{
  skip_on_cran()
  skip_on_ci()

  e <- experiment(experiment_name = "external unified dataset",
                  experimenters = c(""),
                  channels = "cfos",
                  experiment_groups = c("Control", "Isoflurane"),
                  sex_groups = "male",
                  output_path = tempdir())
  path <- file.path(test_object_path("external"), "reformatted_ontologyunchecked.csv")
  e <- import_mapped_datasets(e, normalized_count_paths = path,  show_col_types = FALSE)

  expect_message(e <- check_ontology_coding(e, ontology = "unified"),
                 regexp = "Finished")
})







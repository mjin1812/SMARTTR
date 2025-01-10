test_that("filtering the regions of a dataset works", {

  skip_on_ci()
  skip_on_cran()

  path <- file.path(test_object_path("experiment"), "experiment_light.RDS")
  e <- readRDS(path)

  e <- filter_regions(e,
                      base_regions = c("HY"),
                       ontology = "allen")
  expect_equal(e$combined_normalized_counts$cfos$acronym[1], "LHA")
  expect_snapshot(e$combined_normalized_counts$cfos %>% head())
})





test_that("redundant parent regions are appropriately excluded",{

  skip_on_ci()
  skip_on_cran()

  e <- experiment(experiment_name = "external unified dataset",
                  experimenters = c(""),
                  channels = "cfos",
                  experiment_groups = c("Control", "Isoflurane"),
                  sex_groups = "male",
                  output_path = tempdir())
  path <- file.path(test_object_path("external"), "reformatted_ontologyunchecked.csv")
  e <- import_mapped_datasets(e, normalized_count_paths = path, show_col_types = FALSE)
  prev_acro <- e$combined_normalized_counts$cfos$acronym %>% unique()
  e <- exclude_redundant_regions(e,
                                 ontology = "unified")
  new_acro <- e$combined_normalized_counts$cfos$acronym %>% unique()
  expect_equal(setdiff(prev_acro, new_acro), c("GENd", "MED", "PVA", "MM", "SCm", "APT", "Gi", "PGRNl"))
})

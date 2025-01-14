## usage:
## test_object_path("slice")

test_object_path <- function(object_name) testthat::test_path("testdata", object_name)


make_test_mouse <- function(){
  m <- mouse(mouse_ID = "733",
             sex = "female",
             strain ="129s",
             experiment = "learned helplessness",
             group = "context",
             cohort = "group_A",
             output_path = test_object_path("mouse"))
  s <- slice(slice_ID = "1_4",
             coordinate = -2.14,
             conversion_factor = 1.0833,
             bin = 1,
             z_width = 9,
             hemisphere = NULL,
             registration_path = file.path(test_object_path("slice"), "MAX_733_1_4.tif"),
             slice_directory = test_object_path("slice"))
  m <- add_slice(m, s)
  return(m)
}


make_test_mapped_mouse_light <- function(){
  path <- file.path(test_object_path("slice"), "mapped_slice.RDS")
  s <- readRDS(file = path)
  s$registration_obj$outputfile  <- file.path(test_object_path("slice"), "registrations",
                                              "Registration_MAX_733_1_4")

  m <- mouse(mouse_ID = "733",
             sex = "female",
             strain ="129s",
             experiment = "learned helplessness",
             group = "context",
             cohort = "group_A",
             output_path = test_object_path("mouse"))
  m <- add_slice(m, s)
  return(m)
}




make_test_experiment_light <- function(){
  path <- file.path(test_object_path("experiment"), "experiment_light.RDS")
  e <- readRDS(file = path)
  attr(e, "info")$output_path <- tempdir()
  return(e)

}


make_analyzed_experiment <- function(){
  path <- file.path(test_object_path("experiment"), "analyzed_experiment.RDS")
  e <- readRDS(file = path)
  attr(e, "info")$output_path <- tempdir()
  return(e)

}


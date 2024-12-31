## Testing slice objects
test_that("reading segmentation info from a slice object works", {
  s <- slice(slice_ID = "1_4",
             slice_directory = test_object_path("slice"))
  expect_message(
    s <- import_segmentation_ij(s,
                                mouse_ID = "733",
                                channels = c("eyfp", "cfos", "colabel")),
    regexp = NA,
    class = "simpleMessage/message/condition"
  )
})

test_that("Reading a non-existent channel from a slice object results in error", {
  s <- slice(slice_ID = "1_4",
             slice_directory = test_object_path("slice"))

  expect_error(
    suppressWarnings(s <- import_segmentation_ij(s,
                                                  mouse_ID = "733",
                                                  channels = "fakechannel")),
    "cannot open the connection",
    # NULL,
    inherit = TRUE,
    fixed = TRUE
  )
})

# Testing mouse objects
test_that("reading segmentation data from a mouse object works", {

 m <- make_test_mouse()
  expect_message(
  m <- import_segmentation_ij(m, slice_ID = "1_4", channels = c("eyfp", "cfos", "colabel")),
  regexp = NA,
  class = "simpleMessage/message/condition",
  inherit = TRUE
  )
})


test_that("Reading a non-existent channel from a mouse object results in error", {
  m <- make_test_mouse()
  expect_error(
    suppressWarnings(m <- import_segmentation_ij(m, slice_ID = "1_4",
                                                 channels = "fakechannel")),
    "cannot open the connection",
    # NULL,
    inherit = TRUE,
    fixed = TRUE
  )
})

test_that("Test for existing segmentation data in an object", {
  m <-  make_test_mouse()
  m <- import_segmentation_ij(m, slice_ID = "1_4", channels = c("eyfp", "cfos", "colabel"))

  expect_error(import_segmentation_ij(m, slice_ID = "1_4", channels = c("eyfp", "cfos", "colabel")),
              "There is existing segmentation data",
              inherit = TRUE,
              fixed = TRUE)

})


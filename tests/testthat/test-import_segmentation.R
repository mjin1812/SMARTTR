## Testing ___________________ import_segmentation_ij __________________

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


## _____________________ testing custom import functions ___________________________



## Testing slice objects
test_that("reading custom segmentation info from a slice object works", {


  s <- slice(slice_ID = "1_4",
             slice_directory = test_object_path("slice"))

  # Testing that auto search function works if x_col and y_col are NULL
  expect_message(
    import_segmentation_custom(s,
                               channel = "eYFP",
                               x_col = 5,
                               y_col = 6,
                               meas_path = NULL,
                               quant_path = NULL),
    regexp = NA,
    class = "simpleMessage/message/condition"
  )

  # Testing that import works with supplied import paths
  m_path <-  file.path( test_object_path("slice"), "M_C1_eYFP_733_1_4.txt")
  q_path <-  file.path( test_object_path("slice"), "Q_C1_eYFP_733_1_4_eYFP.txt")

  expect_message(
    import_segmentation_custom(s,
                               channel = "eYFP",
                               x_col = 5,
                               y_col = 6,
                               meas_path = m_path,
                               quant_path = q_path),
    regexp = NA,
    class = "simpleMessage/message/condition"
  )

  # Test mismatch between channel and path parameters
  expect_warning(
    s <- import_segmentation_custom(s,
                                    channel = "nonsense",
                                    x_col = 5,
                                    y_col = 6,
                                    meas_path = m_path,
                                    quant_path = q_path),
    regexp = "does not match",
  )
})


test_that("reading custom segmentation info from a mouse object works", {
  m <- make_test_mouse()

  expect_message(
    import_segmentation_custom(m,
                               channel = "eYFP",
                               slice_ID = "1_4",
                               x_col = 5,
                               y_col = 6,
                               meas_path = NULL,
                               quant_path = NULL),
    regexp = NA,
    class = "simpleMessage/message/condition"
  )

  m_path <-  file.path( test_object_path("slice"), "M_C1_eYFP_733_1_4.txt")
  q_path <-  file.path( test_object_path("slice"), "Q_C1_eYFP_733_1_4_eYFP.txt")

  expect_message(
    import_segmentation_custom(m,
                               channel = "eYFP",
                               slice_ID = "1_4",
                               x_col = 5,
                               y_col = 6,
                               meas_path = m_path,
                               quant_path = q_path),
    regexp = NA,
    class = "simpleMessage/message/condition"
  )

  # Test can't find slice in mouse object
  expect_error(import_segmentation_custom(m,
                                          channel = "eYFP",
                                          slice_ID = "1_1",
                                          x_col = 5,
                                          y_col = 6,
                                          meas_path = NULL,
                                          quant_path = NULL),
               regexp = " no slices matching")

})


# _____________________ testing make segmentation object ________________________

test_that("reading segmentation info from a mouse and slice object works", {
  s <- slice(slice_ID = "1_4",
             slice_directory = test_object_path("slice"))
  channels <- c("eyfp", "cfos", "colabel")

  s <- import_segmentation_ij(s,
                              mouse_ID = "733",
                              channels = channels)
  m <- mouse(mouse_ID = "733")
  m <- add_slice(m, s)

  expect_silent(
    s <- make_segmentation_object(s,
                                  mouse_ID = "733",
                                  channels = channels,
                                  use_filter = FALSE))

  expect_silent(
  m <- make_segmentation_object(m,
                                slice_ID = "1_4",
                                hemisphere = NULL,
                                channels = channels))

  expect_error(
    m <- make_segmentation_object(m,
                                  slice_ID = "1_4",
                                  hemisphere = NULL,
                                  channels = channels), regexp = "existing segmentation")


})



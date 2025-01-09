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

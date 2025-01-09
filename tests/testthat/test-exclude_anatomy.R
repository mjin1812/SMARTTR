test_that("Excluding anatomy works for mouse and slice objects", {

  m <- make_test_mapped_mouse_light()

  expect_silent(m <- exclude_anatomy(m,
                  slice_ID = "1_4",
                  hemisphere = NULL,
                  channels = NULL,
                  clean = TRUE,
                  exclude_right_regions = NULL,
                  exclude_left_regions = NULL,
                  exclude_hemisphere = FALSE,
                  exclude_layer_1 = TRUE,
                  include_right_regions = NULL,
                  include_left_regions = NULL,
                  simplify_regions = TRUE,
                  plot_filtered = FALSE))


  expect_silent(m <- exclude_anatomy(m,
                                     slice_ID = "1_4",
                                     channels = NULL,
                                     clean = TRUE,
                                     exclude_right_regions = c("PTLp", "AUD", "ENT"),
                                     exclude_left_regions = c("PERI"),
                                     exclude_hemisphere = FALSE,
                                     exclude_layer_1 = TRUE,
                                     include_right_regions = NULL,
                                     include_left_regions = NULL,
                                     simplify_regions = TRUE,
                                     plot_filtered = FALSE))


  m <- exclude_anatomy(m,
                       slice_ID = "1_4",
                       channels = NULL,
                       clean = TRUE,
                       exclude_right_regions = c("PTLp", "AUD", "ENT"),
                       exclude_left_regions = c("PERI"),
                       exclude_hemisphere = FALSE,
                       exclude_layer_1 = TRUE,
                       include_right_regions = c("BLA"),
                       include_left_regions = c("BLA"),
                       simplify_regions = TRUE,
                       plot_filtered = FALSE)

  out <- m$slices$`1_4`$forward_warped_data$cfos$acronym %>% unique()
  expect_setequal(out, c("BLAa", "BLAp", "BLAv"))


})

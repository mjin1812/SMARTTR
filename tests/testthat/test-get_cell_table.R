test_that("Getting cell table works", {

  skip_on_ci()
  skip_on_cran()

  m <- make_test_mapped_mouse_light()

  m <- get_cell_table(m)

  expect_s3_class(m$cell_table$cfos, "tbl_df")
  expect_snapshot(m$cell_table$cfos %>% head())
})

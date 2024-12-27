test_that("add_slice returns a mouse object", {
  expect_equal(class(add_slice(mouse(), slice(slice_ID = "1_1"))), "mouse")
})

test_that("Check for an existing slice and user-optional replacement", {
  expect_error(add_slice(add_slice(mouse(), slice(slice_ID = "1_1")), slice(slice_ID = "1_1"), replace = FALSE),
                "If you want to replace a previous slice object")
  expect_equal(add_slice(mouse(), slice(slice_ID = "1_1")),
               add_slice(add_slice(mouse(), slice(slice_ID = "1_1")), slice(slice_ID = "1_1"), replace = TRUE))
})



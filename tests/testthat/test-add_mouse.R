test_that("add_mouse returns an experiment object", {
  expect_equal(class(add_mouse(experiment(), mouse(mouse_ID = "0"))), "experiment")
})

test_that("Check for an existing mouse and user-optional replacement", {
  expect_error(add_mouse(add_mouse(experiment(), mouse(mouse_ID = "0")), mouse(mouse_ID = "0"), replace = FALSE),
               "If you want to replace a previous mouse object")
  expect_equal(add_mouse(experiment(), mouse(mouse_ID = "0")),
               add_mouse(add_mouse(experiment(), mouse(mouse_ID = "0")), mouse(mouse_ID = "0"), replace = TRUE))
})



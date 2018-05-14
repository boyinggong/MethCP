
context("Constructor")

params1 <- list(test = "myTest",
                group1 = paste0("Treatment", 1:3),
                group2 = paste0("Control", 1:3),
                chr = rep("Chr01", 5),
                pos = c(2, 5, 9, 10, 18),
                effect.size = c(1, -1, NA, 9, Inf))

test_that("contructor works", {
  expect_message(obj <- do.call("MethCP", c(params1, list(pvals = c(0, 0.1, 0.9, NA, 0.02)))))
  expect_equal(length(obj@stat[[1]]), 3)
  expect_equal(length(obj@segmentation), 0)
})

test_that("contructor throw errors for invalid inputs", {
  expect_error(do.call("MethCP", c(params1, list(pvals = c(0, 0.1, 0.9, NA, Inf)))))
  expect_error(do.call("MethCP", c(params1, list(pvals = c(0, 0.1, 0.9, 0.7)))))
})


context("methcpFromStat")

params1 <- list(
    testName = "myTest",
    chr = rep("Chr01", 5),
    pos = c(2, 5, 9, 10, 18),
    effect.size = c(1, -1, NA, 9, Inf))

test_that("function works", {
    expect_message(obj <- do.call(
        "methcpFromStat", c(params1, list(pvals = c(0, 0.1, 0.9, NA, 0.02)))))
    expect_equal(length(obj@stat[[1]]), 3)
    expect_equal(length(obj@segmentation), 0)
})

test_that("function throw errors for invalid inputs", {
    expect_error(
        do.call(
            "methcpFromStat",
            c(params1, list(pvals = c(0, 0.1, 0.9, NA, Inf)))))
    expect_error(
        do.call(
            "methcpFromStat",
            c(params1, list(pvals = c(0, 0.1, 0.9, 0.7)))))
})



context("MethCPFromStat")

library(GenomicRanges)

data <- data.frame(
    chr = rep("Chr01", 5),
    pos = c(2, 5, 9, 10, 18),
    effect.size = c(1,-1, NA, 9, Inf),
    pvals = c(0, 0.1, 0.9, NA, 0.02))
test_that("function works for data frame", {
    expect_message(obj <- MethCPFromStat(
        data, test.name="myTest",
        pvals.field = "pvals",
        effect.size.field="effect.size",
        seqnames.field="chr",
        pos.field="pos"))
    expect_equal(length(obj@stat[[1]]), 3)
    expect_equal(length(obj@segmentation), 0)
})

data <- GRanges(
    "Chr01", IRanges(c(2, 5, 9, 10, 18), c(2, 5, 9, 10, 18)),
    pvals=c(0, 0.1, 0.9, NA, 0.02), effect.size = c(1,-1, NA, 9, Inf))
test_that("function works for GRanges object", {
    expect_message(obj <- MethCPFromStat(
        data, test.name="myTest",
        pvals.field = "pvals",
        effect.size.field="effect.size"
    ))
})

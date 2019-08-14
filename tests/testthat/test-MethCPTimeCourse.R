
context("MethCP time-course data")

test_that("methylKit", {
    expect_silent(
        obj_ts <- calcLociStatTimeCourse(bs_object_ts, meta)
    )
    expect_silent(
        obj_ts <- segmentMethCP(
            obj_ts, bs_object_ts,
            region.test = "stouffer", mc.cores = 1)
    )
    expect_silent(res_ts <- getSigRegion(obj_ts))
    pos <- unlist(lapply(1:nrow(res_ts), function(x)
        ((res_ts[x, 2]/10):(res_ts[x, 3]/10))*10))
    jaccard_ts <- .calc_jaccard(pos, DMRs_pos)
    expect_true(jaccard_ts >= 0.8)
})

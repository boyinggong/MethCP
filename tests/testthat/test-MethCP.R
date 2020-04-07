
context("MethCP two population groups")

test_that("DSS", {
    expect_silent(methcp_obj1 <- calcLociStat(
        bs_object,
        group1 = paste0("treatment", 1:3),
        group2 = paste0("control", 1:3),
        test = "DSS"))
    expect_silent(
        methcp_obj1 <- segmentMethCP(methcp_obj1, bs_object, region.test = "weighted-coverage"))
    expect_silent(methcp_res1 <- getSigRegion(methcp_obj1))
    pos <- unlist(lapply(1:nrow(methcp_res1), function(x)
        ((methcp_res1[x, 2]/10):(methcp_res1[x, 3]/10))*10))
    jaccard1 <- .calc_jaccard(pos, DMRs_pos)
    expect_true(jaccard1 >= 0.8)
})

test_that("methylKit", {
    expect_message(methcp_obj2 <- calcLociStat(
        bs_object,
        group1 = paste0("treatment", 1:3),
        group2 = paste0("control", 1:3),
        test = "methylKit"))
    expect_silent(
        methcp_obj2 <- segmentMethCP(methcp_obj2, bs_object, region.test = "fisher"))
    expect_silent(methcp_res2 <- getSigRegion(methcp_obj2))
    pos <- unlist(lapply(1:nrow(methcp_res2), function(x)
        ((methcp_res2[x, 2]/10):(methcp_res2[x, 3]/10))*10))
    jaccard2 <- .calc_jaccard(pos, DMRs_pos)
    expect_true(jaccard2 >= 0.8)
})


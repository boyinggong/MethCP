
context("MethCP two population groups")

test_that("DSS", {
  expect_silent(
    methcp_obj1 <- calcLociStat(bs_object,
                                group1 = paste0("treatment", 1:3),
                                group2 = paste0("control", 1:3),
                                test = "DSS"))
  expect_silent(
    methcp_obj1 <- segmentMethCP(methcp_obj1, bs_object, region.test = "weighted-coverage", mc.cores = 1))
  expect_silent(methcp_res1 <- getSigRegion(methcp_obj1))
})

test_that("methylKit", {
  expect_message(
    methcp_obj2 <- calcLociStat(bs_object,
                                group1 = paste0("treatment", 1:3),
                                group2 = paste0("control", 1:3),
                                test = "methylKit"))
  expect_silent(
    methcp_obj2 <- segmentMethCP(methcp_obj2, bs_object, region.test = "fisher", mc.cores = 1))
  expect_silent(methcp_res2 <- getSigRegion(methcp_obj2))
})

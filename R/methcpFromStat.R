
#' @title Create a \code{MethCP} object given the per-cytosine based tests.
#'
#' @description
#' Given the per-cytosine based \code{p}-values and effect sizes,
#' \code{methcpFromStat} create a \code{MethCP} object for segmentation.
#'
#' @details
#' For each cytosine, \code{calcLociStatTimeCourse} calculates a statistics
#' using either package \code{DSS} or \code{methylKit} to test the
#' differences between two groups, and returns a \code{MethCP} object.
#' For customized per-cytosine statistics, please use the function
#' \code{methcpFromStat}.
#' The input \code{bs.object} is a \code{BSseq} object from the \code{bsseq}
#' package which contains the raw data including coverges, methylated
#' counts and position infomation for every cytosine in the dataset.
#'
#' @param chr a character vector containing the cytosine chromosome infomation.
#' @param pos a numeric vector containing the cytosine positions.
#' @param pvals a numeric vector containing the \code{p}-values for each cytosine.
#' @param effect.size a numeric vector containing the effect sizes for each cytosine.
#' @param testName a character string containing the names of the test to be
#' performed per cytosine.
#'
#' @return a \code{MethCP} object that is not segmented.
#'
#' @examples
#' obj <- methcpFromStat(testName = "myTest",
#'                       chr = rep("Chr01", 5),
#'                       pos = c(2, 5, 9, 10, 18),
#'                       effect.size = c(1, -1, NA, 9, Inf),
#'                       pvals = c(0, 0.1, 0.9, NA, 0.02))
#' obj
#'
#' @import parallel
#'
#' @export
methcpFromStat <- function(chr, pos, pvals, effect.size, testName){
  return(MethCP(test = testName,
                group1 = "notApplicable",
                group2 = "notApplicable",
                chr = chr, pos = pos, pvals = pvals, effect.size = effect.size))
}

#' @title Obtain the significant DMRs.
#'
#' @usage
#' getSigRegion(
#'     object, sig.level = 0.01, mean.coverage = 1,
#'     mean.diff = 0.1, nC.valid = 10)
#'
#' @description
#' \code{getSigRegion} returns the significant DMRs giving the segmented
#' \code{MethCP} object.
#'
#' @param object a \code{MethCP} object that is segmented using function
#' \code{segmentMethCP}.
#' @param sig.level significance level to call a region DMR.
#' @param mean.coverage The minimum average coverage required for the
#' reported DMRs.
#' @param mean.diff The minimum differences between groups required for
#' the reported DMRs.
#' @param nC.valid number of valid cytosines required for the reported DMRs.
#'
#' @return a \code{data.frame} containing the DMRs.
#'
#' @examples
#' library(bsseq)
#' # Simulate a small dataset with 2000 cyotsine and 6 samples,
#' # 3 in the treatment group and 3 in the control group. The
#' # methylation ratio are generated using Binomial distribution
#' # with probability 0.3.
#' nC <- 2000
#' sim_cov <- rnbinom(6*nC, 5, 0.5) + 5
#' sim_M <- vapply(
#'     sim_cov, function(x) rbinom(1, x, 0.3), FUN.VALUE = numeric(1))
#' sim_cov <- matrix(sim_cov, ncol = 6)
#' sim_M <- matrix(sim_M, ncol = 6)
#' # methylation ratios in the DMRs in the treatment group are
#' # generated using Binomial(0.7)
#' DMRs <- c(600:622, 1089:1103, 1698:1750)
#' sim_M[DMRs, 1:3] <- vapply(
#'     sim_cov[DMRs, 1:3], function(x) rbinom(1, x, 0.7),
#'     FUN.VALUE = numeric(1))
#' # sample names
#' sample_names <- c(paste0("treatment", 1:3), paste0("control", 1:3))
#' colnames(sim_cov) <- sample_names
#' colnames(sim_M) <- sample_names
#'
#' # create a bs.object
#' bs_object <- BSseq(gr = GRanges(
#'     seqnames = "Chr01",
#'     IRanges(start = (1:nC)*10, width = 1)),
#'     Cov = sim_cov, M = sim_M, sampleNames = sample_names)
#' DMRs_pos <- DMRs*10
#' methcp_obj1 <- calcLociStat(
#'     bs_object,
#'     group1 = paste0("treatment", 1:3),
#'     group2 = paste0("control", 1:3),
#'     test = "DSS")
#' methcp_obj1 <- segmentMethCP(
#'     methcp_obj1, bs_object,
#'     region.test = "weighted-coverage",
#'     mc.cores = 1)
#' methcp_res1 <- getSigRegion(methcp_obj1)
#'
#' @export
#' @aliases getSigRegion
getSigRegion <- function(
    object, sig.level = 0.01, mean.coverage = 1,
    mean.diff = 0.1, nC.valid = 10)
{
    if (!is(object, "MethCP")){
        stop("Input must be an object of class \"MethCP\".")
    }
    res <- as.data.frame(segmentation(object))
    res$strand <- NULL
    res$width <- NULL
    res <- res[
        res$region.pval <= sig.level & res$mean.cov >= mean.coverage &
            abs(res$mean.diff) >= mean.diff & res$nC.valid >= nC.valid, ]
    res$mean.diff <- round(res$mean.diff, 4)
    res$mean.cov <- round(res$mean.cov, 4)
    return(res)
}

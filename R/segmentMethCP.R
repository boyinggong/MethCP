#' @title Perform segmentation on a \code{MethCP} object.
#'
#' @usage
#' segmentMethCP(
#'     methcp.object, bs.object,
#'     region.test = c(
#'         "fisher", "stouffer", "weighted-variance", "weighted-coverage"),
#'     min.width = 2, sig.level = 0.01,
#'     presegment_dist = 600, ...)
#'
#' @description Perform CBS algorithm that segments the genome into
#' similar levels of sigficance.
#'
#' @details
#' The \code{MethCP} object \code{methcp.object} can be generated from
#' functions \code{calcLociStat}, \code{calcLociStatTimeCourse}, or
#' \code{methcpFromStat}.
#'
#' If \code{region.test = "fisher"}, Fisher's combined probability test is used.
#'
#' If \code{region.test = stouffer} Stouffer's test is applied.
#'
#' If \code{region.test = "weighted-variance"} we use the variance of the
#' test to combine per-cytosine based statistcis into a region-based statistic.
#'
#' If \code{region.test = "weighted-coverage"} we use the coverage of the
#' test to combine per-cytosine based statistcis into a region-based statistic.
#'
#' @param methcp.object a \code{MethCP} object.
#' @param bs.object a \code{BSseq} object from the \code{bsseq} package.
#' @param region.test The meta-analysis method used to create region-based
#' test statistics.
#' @param min.width the minimum width for the segments, which is used as
#' termination rule for the segmentation algorithm.
#' @param sig.level the significance level of the segments, which is used as
#' termination rule for the segmentation algorithm.
#' @param presegment_dist the maximum distance between cytosines for the
#' presegmentation.
#' @param ... argument to be passed to segment function in DNAcopy package
#'
#' @return a \code{MethCP} object that is not segmented.
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
#'     seqnames = "Chr01", IRanges(start = (1:nC)*10, width = 1)),
#'     Cov = sim_cov, M = sim_M,
#'     sampleNames = sample_names)
#' DMRs_pos <- DMRs*10
#' methcp_obj1 <- calcLociStat(
#'     bs_object,
#'     group1 = paste0("treatment", 1:3),
#'     group2 = paste0("control", 1:3),
#'     test = "DSS")
#' methcp_obj1 <- segmentMethCP(
#'     methcp_obj1, bs_object,
#'     region.test = "weighted-coverage")
#'
#' @import parallel
#' @importFrom DNAcopy CNA segment
#' @importFrom bsseq getCoverage
#' @importFrom BiocParallel bplapply MulticoreParam
#'
#' @export
#' @aliases segmentMethCP
segmentMethCP <- function(
    methcp.object, bs.object,
    region.test = c(
        "fisher", "stouffer", "weighted-variance", "weighted-coverage"),
    min.width = 2, sig.level = 0.01,
    presegment_dist = 600, ...)
{
    object <- methcp.object
    if (!is(object, "MethCP")){
        stop("Input must be an object of class \"MethCP\".")
    }
    if (test(object) == "methylKit" &
        region.test %in% c("weighted-variance", "weigxhted-coverage")){
        stop(paste(
            "can not apply weighted effect size method",
            "with methylKit."))
    }
    if (length(segmentation(object)) != 0){
        a <- readline(paste(
            "Object has been segmented.",
            "Remove previous segmentation? (y/n) > "))
        if (a == "y") {
            message("Removed. Start running new segmentation ...")
            segmentation(object) <- GRanges()
        } else {
            stop("Stopped.")
        }
    }
    if (test(object) == "methylKit") {
        stat(object) <- GRangesList(lapply(names(stat(object)), function(o){
            tmp <- stat(object)[[o]]
            tmp$stat <- .pvalToStat(tmp$pval, tmp$methDiff)
            tmp
        }))
    }
    # calculate total coverage and methylated counts for each loci
    stat(object) <- GRangesList(
        lapply(seq_len(length(stat(object))), function(o){
            tmp <- stat(object)[[o]]
            ovrlp <- findOverlaps(granges(bs.object), tmp)
            tmp$CovGroup1 <- rowSums(
                as.data.frame(
                    getCoverage(bs.object)[from(ovrlp), group1(object)]))
            tmp$CovGroup2 <- rowSums(
                as.data.frame(
                    getCoverage(bs.object)[from(ovrlp), group2(object)]))
            tmp
        }))
    segments <- list()
    for (chr in seq_len(length(stat(object)))) {
        o <- stat(object)[[chr]]
        pos <- start(o)
        presegments <- split(
            seq_len(length(pos)),
            cumsum(c(TRUE, abs(diff(pos)) >= presegment_dist)))
        res <- BiocParallel::bplapply(
            presegments,
            function(idx){
                cp.object <- CNA(
                    o[idx]$stat,
                    chrom=as.vector(seqnames(o[idx])),
                    maploc=start(o[idx]),
                    data.type="logratio")
                invisible(capture.output(segment.cp.object <- segment(
                    cp.object, verbose = 1, min.width = min.width,
                    alpha = sig.level, ...)))
                return(segment.cp.object$output)
            }, BPPARAM=BiocParallel::MulticoreParam())
        res <- do.call("rbind", res)
        res$ID <- NULL
        res$seg.mean <- NULL
        colnames(res)[4] <- "nC.valid"
        segments[[chr]] <- res
    }
    segments <- as.data.frame(do.call("rbind", segments))
    segments <- GRanges(segments)
    
    # calculate region summary
    ovrlp <- findOverlaps(granges(bs.object), segments)
    segments$nC <- as.numeric(table(to(ovrlp)))
    M1 <- as.numeric(by(getCoverage(
        bs.object, type = "M")[from(ovrlp), group1(object)], to(ovrlp), sum))
    M2 <- as.numeric(by(getCoverage(
        bs.object, type = "M")[from(ovrlp), group2(object)], to(ovrlp), sum))
    Cov1 <- as.numeric(by(getCoverage(bs.object)[
        from(ovrlp), group1(object)], to(ovrlp), sum))
    Cov2 <- as.numeric(by(getCoverage(bs.object)[
        from(ovrlp), group2(object)], to(ovrlp), sum))
    segments$mean.diff <- M1/Cov1 - M2/Cov2
    segments$mean.cov <- (Cov1 + Cov2)/length(
        c(group1(object), group2(object)))/segments$nC.valid
    
    # calculate region statistics
    ovrlp <- findOverlaps(unlist(stat(object)), segments)
    if (region.test == "fisher"){
        segments$region.pval <- as.numeric(tapply(
            unlist(stat(object))$pval[from(ovrlp)],
            to(ovrlp), .calcFisherPval))
    } else if (region.test == "stouffer"){
        segments$region.pval <- as.numeric(tapply(
            unlist(stat(object))$pval[from(ovrlp)],
            to(ovrlp), .calcStoufferPvalOneSided))
    } else if (region.test == "weighted-variance"){
        if (is.null(stat(object)[[1]]$mu)){
            stop(paste(
                "weighted test not applicable,",
                "consider Stouffer's test or Fisher's test."))
        }
        segments$region.pval <- as.numeric(by(
            elementMetadata(unlist(stat(object)))[from(ovrlp), ], to(ovrlp),
            function(x) .calcWeightedPval(x$mu, x$se, 1/x$se)))
    } else if (region.test == "weighted-coverage"){
        if (is.null(stat(object)[[1]]$mu)){
            stop(paste(
                "weighted test not applicable,",
                "consider Stouffer's test or Fisher's test."))
        }
        segments$region.pval <- as.numeric(by(
            elementMetadata(unlist(stat(object)))[from(ovrlp), ], to(ovrlp),
            function(x)
                .calcWeightedPval(x$mu, x$se, x$CovGroup1 + x$CovGroup2)))
    }
    segmentation(methcp.object) <- segments
    return(methcp.object)
}


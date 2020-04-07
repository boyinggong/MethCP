
#' @title Calculate the per-cytosine statistics for simple two-group
#' comparisons.
#'
#' @usage
#' calcLociStat(
#'     bs.object, group1, group2, test = c("DSS", "methylKit"))
#'
#' @description
#' \code{calcLociStat} calculates per-cytosine based statistics
#' between two population groups.
#'
#' @details
#' For each cytosine, \code{calcLociStat} calculates a statistics
#' using either package \code{DSS} or \code{methylKit} to test the
#' differences between two groups, and returns a \code{MethCP} object.
#' For customized per-cytosine statistics, please use the function
#' \code{methcpFromStat}.
#' The input \code{bs.object} is a \code{BSseq} object from the \code{bsseq}
#' package which contains the raw data including coverges, methylated
#' counts and position infomation for every cytosine in the dataset.
#'
#' @param bs.object a \code{BSseq} object from the \code{bsseq} package.
#' @param group1 a character vector containing the sample names of the
#' treatment group.
#' @param group2 a character vector containing the sample names of the
#' control group.
#' @param test a character string containing the names of the test to be
#' performed per cytosine.
#'
#' @return a \code{MethCP} object that is not segmented.
#'
#' @examples
#' library(bsseq)
#' library(GenomicRanges)
#' library(IRanges)
#'
#' set.seed(0286374)
#'
#' # Similate a small dataset with 11 cyotsine and 6 samples,
#' # 3 in the treatment group and 3 in the control group. The
#' # methylation ratio are generated using Binomial distribution
#' # with probability 0.3.
#' nC <- 2000
#' sim_cov <- rnbinom(6*nC, 5, 0.5) + 5
#' sim_M <- vapply(
#'     sim_cov, function(x) rbinom(1, x, 0.3),
#'     FUN.VALUE = numeric(1))
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
#'     seqnames = "Chr01", IRanges(
#'         start = (1:nC)*10, width = 1)),
#'     Cov = sim_cov, M = sim_M, sampleNames = sample_names)
#' methcp_obj1 <- calcLociStat(
#'     bs_object,
#'     group1 = paste0("treatment", 1:3),
#'     group2 = paste0("control", 1:3),
#'     test = "DSS")
#' methcp_obj2 <- calcLociStat(
#'     bs_object,
#'     group1 = paste0("treatment", 1:3),
#'     group2 = paste0("control", 1:3),
#'     test = "methylKit")
#'
#' @importFrom methylKit calculateDiffMeth
#' @importClassesFrom methylKit methylBase
#' @importFrom DSS DMLtest
#' @importFrom bsseq sampleNames
#' @importFrom GenomeInfoDb sortSeqlevels
#' @import GenomicRanges
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom utils capture.output
#'
#' @export
calcLociStat <- function(
    bs.object, group1, group2,
    test = c("DSS", "methylKit")){

    if (!is(bs.object, "BSseq")){
        stop(paste0(
            "Input must be an object of class",
            "\"BSseq\" from bsseq package."))
    }
    if (!all(c(group1, group2) %in% sampleNames(bs.object))){
        stop(paste0(
            "Group sample names must be a subset of",
            "`sampleNames(bs.object)`."))
    }
    if(!(test %in% c("DSS", "methylKit"))){
        stop("Test type is not supported.")
    }
    object_list <- list()
    for (chr in unique(seqnames(bs.object))){
        object_list[[chr]] <- bs.object[seqnames(bs.object) == chr]
    }
    if(test == "DSS"){
        stat <- GenomicRanges::GRangesList(
            BiocParallel::bplapply(object_list, function(o){
            invisible(capture.output(tmp <- DMLtest(
                o, group1, group2, equal.disp = FALSE, smoothing = FALSE)))
            GRanges(
                seqnames = tmp$chr,
                IRanges(start = tmp$pos, end = tmp$pos),
                mu = tmp$diff,
                se = tmp$diff.se,
                stat = tmp$stat,
                pval = tmp$pval)
        }, BPPARAM=BiocParallel::MulticoreParam()))
    } else if(test == "methylKit"){
        nsample <- length(group1) + length(group2)
        stat <- GenomicRanges::GRangesList(
            BiocParallel::bplapply(object_list, function(o){
            df <- cbind(
                as.data.frame(granges(o)),
                getCoverage(o, type = "Cov")[, c(group1, group2)],
                getCoverage(o, type = "M")[, c(group1, group2)],
                (getCoverage(o, type = "Cov") - getCoverage(
                    o, type = "M"))[, c(group1, group2)])
            df$width <- NULL
            colnames(df) <- c(
                "chr", "start", "end", "strand",
                paste0("coverage", seq_len(nsample)),
                paste0("numCs", seq_len(nsample)),
                paste0("numTs", seq_len(nsample)))
            df <- df[
                , c("chr", "start", "end", "strand",
                    vapply(seq_len(nsample), function(x)
                        paste0(c("coverage", "numCs", "numTs"), x),
                        FUN.VALUE = character(3)))]
            coverage.ind <- which(vapply(colnames(df), function(x)
                grepl("coverage", x), FUN.VALUE = logical(1)))
            names(coverage.ind) <- NULL
            obj <- new(
                "methylBase", (df),
                sample.ids = c(group1, group2),
                assembly = "-",
                context = "-",
                treatment = c(
                    rep(1, length(group1)),
                    rep(0, length(group2))),
                coverage.index = coverage.ind,
                numCs.index = coverage.ind+1,
                numTs.index = coverage.ind+2,
                destranded = TRUE,
                resolution = "base" )
            filter <- (rowSums(as.data.frame(getCoverage(o)[, group1])) >= 1) &
                (rowSums(as.data.frame(getCoverage(o)[, group2])) >= 1)
            errorInd <- rowSums(as.data.frame(getCoverage(o, type = "M"))) /
                rowSums(as.data.frame(getCoverage(o)))
            errorInd <- (errorInd == 0) | (errorInd == 1)
            errorInd[is.na(errorInd)] = FALSE
            obj_filtered <- obj[filter & !errorInd, ]
            invisible(capture.output(tmp <- calculateDiffMeth(obj_filtered)))
            gr <- GRanges(
                seqnames = tmp$chr,
                IRanges(start = tmp$start, end = tmp$start),
                pval = tmp$pval, methDiff = tmp$meth.diff/100)
            if(sum(errorInd) != 0){
                gr0 <- granges(o[which(errorInd)])
                gr0$pval = 1
                gr0$methDiff = 0
                gr <- unlist(GRangesList(gr, gr0))
                gr <- sortSeqlevels(gr)
                gr <- sort(gr)
            }
            return(gr)
        }, BPPARAM=BiocParallel::MulticoreParam()))
    }
    methcpObj <- new(
        "MethCP",
        group1 = group1,
        group2 = group2,
        test = test,
        stat = stat)
    return(methcpObj)
}

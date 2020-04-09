
#' @title Calculate the per-cytosine statistics for time-course data.
#'
#' @usage
#' calcLociStatTimeCourse(
#'     bs.object, meta, force.slope = FALSE,
#'     BPPARAM = bpparam())
#'
#' @description
#' For each cytosine, \code{calcLociStatTimeCourse} fits a linear model
#' on the arcsin-tranformed methylation ratios, and test the differences
#' of the slope between the treatment and the control group.
#'
#' @details
#'
#' \code{bs.object} is a \code{BSseq} object from the \code{bsseq} package,
#' which contains the raw data including coverges, methylated counts and
#' position infomation for every cytosine in the dataset.
#' \code{meta} must contain columns named \code{Condition}, \code{Time}
#' and \code{SampleName} in the dataframe. They are used to fit the linear
#' model.
#'
#' @param bs.object a \code{BSseq} object from the \code{bsseq} package.
#' @param meta a \code{data.frame}. See details.
#' @param force.slope if \code{TRUE}, we force the slope in the linear
#' model to be the same between two conditions. Otherwise, the slopes are
#' fitted separately but not tested.
#' @param BPPARAM An optional BiocParallelParam instance determining the 
#' parallel back-end to be used during evaluation, or a list of 
#' BiocParallelParam instances, to be applied in sequence for nested calls 
#' to BiocParallel functions. Default bpparam().
#'
#' @return A \code{MethCP} object that is not segmented.
#'
#' @examples
#' library(bsseq)
#' # Simulate a small dataset with 2000 cyotsine and 10 samples,
#' # 5 in the treatment group and 5 in the control group. The
#' # methylation ratio are generated using Binomial distribution
#' # with probability 0.3, 0.4, 0.5, 0.6 and 0.7 for 5 time points.
#' nC <- 2000
#' nsamples <- 5
#' sim_cov <- rnbinom(10*nC, 5, 0.5) + 5
#' sim_cov <- matrix(sim_cov, ncol = 10)
#' time_point <- rep(1:nsamples, 2)
#' ratios <- time_point/10 + 0.2
#' sim_M <- sapply(1:(2*nsamples), function(i){
#'     sapply(sim_cov[, i], function(j) rbinom(1, j, ratios[i]))
#' })
#' sim_M <- matrix(sim_M, ncol = 2*nsamples)
#' # methylation ratios in the DMRs in the treatment group are
#' # generated using Binomial(0.3)
#' DMRs <- c(600:622, 1089:1103, 1698:1750)
#' sim_M[DMRs, 1:5] <- sapply(
#'     sim_cov[DMRs, 1:5], function(x) rbinom(1, x, 0.3))
#' # sample names
#' sample_names <- c(paste0("treatment", 1:nsamples),
#' paste0("control", 1:nsamples))
#' colnames(sim_cov) <- sample_names
#' colnames(sim_M) <- sample_names
#'
#' # create a bs.object
#' bs_object_ts <- BSseq(gr = GRanges(
#'     seqnames = "Chr01", IRanges(
#'         start = (1:nC)*10, width = 1)),
#'     Cov = sim_cov, M = sim_M, sampleNames = sample_names)
#' DMRs_pos_ts <- DMRs*10
#' meta <- data.frame(
#'     Condition = rep(
#'         c("treatment", "control"),
#'         each = nsamples),
#'     SampleName = sample_names,
#'     Time = time_point)
#' obj_ts <- calcLociStatTimeCourse(bs_object_ts, meta)
#' obj_ts
#'
#' @import BiocParallel
#'
#' @export
calcLociStatTimeCourse <- function(
    bs.object, meta, force.slope = FALSE,
    BPPARAM = bpparam()){

    if (!is(bs.object, "BSseq")){
        stop(paste0(
            "ERROR: Input must be an object of class",
            "\"BSseq\" from bsseq package."))
    }
    if (!is(meta, "data.frame")){
        stop("ERROR: Meta data must be a data frame.")
    }
    if (!all(c("Condition", "Time", "SampleName") %in% colnames(meta))){
        stop("ERROR: meta data must contain Condition, Time and SampleName")
    }
    coverage <- as.data.frame(getCoverage(
        bs.object, type = "Cov"))[, meta$SampleName]
    M <- as.data.frame(getCoverage(bs.object, type = "M"))[, meta$SampleName]
    ratios <- .asinTransform(M/coverage)
    res <- BiocParallel::bplapply(seq_len(nrow(coverage)), function(i){
        suppressWarnings({
            data <- meta
            data$r <- t(ratios[i, ])[, 1]
            data$cov <- t(coverage[i, ])[, 1]
            data <- na.omit(data)
            if (force.slope){
                summary(lm(
                    r~Time+Condition:Time, data=data,
                    weights=cov))$coefficients[3, 3:4]
            } else {
                summary(lm(
                    r~Condition+Time+Condition:Time, data=data,
                    weights=cov))$coefficients[4, 3:4]
            }
        })
    }, BPPARAM=BPPARAM)
    res <- do.call("rbind", res)
    stat <- granges(bs.object)
    stat$stat <- qnorm(1-res[, 2]/2)*sign(res[, 1])
    stat$pval <-res[, 2]
    stat <- GenomicRanges::GRangesList(list(stat))
    methcpObj <- new(
        "MethCP",
        group1 = meta$SampleName[meta$Condition == unique(meta$Condition)[1]],
        group2 = meta$SampleName[meta$Condition == unique(meta$Condition)[2]],
        test = "TimeCourse", stat = stat)
    return(methcpObj)
}

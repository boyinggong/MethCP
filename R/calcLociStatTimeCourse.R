
#' @title Calculate the per-cytosine statistics for time-course data.
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
#' @param meta a \code{data.frame}.
#' @param force.slope if \code{TRUE}, we force the slope in the linear
#' model to be the same between two conditions. Otherwise, the slopes are
#' fitted separately but not tested.
#' @param mc.cores number of cores used for the parallelization.
#'
#' @return A \code{MethCP} object that is not segmented.
#'
#' @examples
#'
#' @import parallel
#'
#' @export
calcLociStatTimeCourse <- function(bs.object, meta, force.slope = FALSE, mc.cores = 1){

  if (!is(bs.object, "BSseq")){
    stop("ERROR: Input must be an object of class \"BSseq\" from bsseq package.")
  }
  if (class(meta) != "data.frame"){
    stop("ERROR: Meta data must be a data frame.")
  }
  if (!all(c("Condition", "Time", "SampleName") %in% colnames(meta))){
    stop("ERROR: meta data must contain Condition, Time and SampleName")
  }
  coverage <- as.data.frame(getCoverage(bs.object, type = "Cov"))[, meta$SampleName]
  M <- as.data.frame(getCoverage(bs.object, type = "M"))[, meta$SampleName]
  ratios <- .asinTransform(M/coverage)
  res <- mclapply(1:nrow(coverage), function(i){
    suppressWarnings({
      data <- meta
      data$r <- t(ratios[i, ])[, 1]
      data$cov <- t(coverage[i, ])[, 1]
      data <- na.omit(data)
      if (force.slope){
        summary(lm(r~Time+Condition:Time, data=data, weights=cov))$coefficients[3, 3:4]
      } else {
        summary(lm(r~Condition+Time+Condition:Time, data=data, weights=cov))$coefficients[4, 3:4]
      }
    })
  }, mc.cores = mc.cores)
  res <- do.call("rbind", res)
  stat <- granges(bs.object)
  stat$stat <- qnorm(1-res[, 2]/2)*sign(res[, 1])
  stat$pval <-res[, 2]
  stat <- GenomicRanges::GRangesList(list(stat))
  methcpObj <- new("MethCP",
                   group1 = meta$SampleName[meta$Condition == unique(meta$Condition)[1]],
                   group2 = meta$SampleName[meta$Condition == unique(meta$Condition)[2]],
                   test = "TimeCourse",
                   stat = stat)
  return(methcpObj)
}

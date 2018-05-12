
#' @title Sum of vector elements.
#'
#' @description
#' \code{sum} returns the sum of all the values present in its arguments.
#'
#' @details
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#'
#' @param ... Numeric, complex, or logical vectors.
#' @param na.rm A logical scalar. Should missing values (including NaN)
#'   be removed?
#' @return If all inputs are integer and logical, then the output
#'   will be an integer. If integer overflow
#'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#'   will be NA with a warning. Otherwise it will be a length-one numeric or
#'   complex vector.
#'
#'   Zero-length vectors have sum 0 by definition. See
#'   \url{http://en.wikipedia.org/wiki/Empty_sum} for more details.
#' @examples
#' sum(1:10)
#' sum(1:5, 6:10)
#' sum(F, F, F, T, T)
#'
#' sum(.Machine$integer.max, 1L)
#' sum(.Machine$integer.max, 1)
#'
#' \dontrun{
#' sum("a")
#' }
#'
#'
#' @export
calcLociStat <- function(bs.object, group1, group2,
                         test = c("DSS", "methylKit"), mc.cores = 1){

  if (!is(bs.object, "BSseq")){
    stop("ERROR: Input must be an object of class \"BSseq\" from bsseq package.")
  }
  if (!all(c(group1, group2) %in% sampleNames(bs.object))){
    stop("ERROR: Group sample names must be a subset of `sampleNames(bs.object)`.")
  }
  if(!(test %in% c("DSS", "methylKit"))){
    stop("ERROR: Test type is not supported.")
  }
  object_list <- list()
  for (chr in unique(seqnames(bs.object))){
    object_list[[chr]] <- bs.object[seqnames(bs.object) == chr]
  }
  if(test == "DSS"){
    stat <- GenomicRanges::GRangesList(mclapply(object_list, function(o){
      invisible(capture.output(tmp <- DMLtest(o,
                     group1,
                     group2,
                     equal.disp = FALSE,
                     smoothing = FALSE)))
      GRanges(seqnames = tmp$chr,
              IRanges(start = tmp$pos, end = tmp$pos),
              mu = tmp$diff,
              se = tmp$diff.se,
              stat = tmp$stat,
              pval = tmp$pval)
    }, mc.cores = mc.cores))
  } else if(test == "methylKit"){
    nsample <- length(group1) + length(group2)
    stat <- GenomicRanges::GRangesList(mclapply(object_list, function(o){
      df <- cbind(
        as.data.frame(granges(o)),
        getCoverage(o, type = "Cov")[, c(group1, group2)],
        getCoverage(o, type = "M")[, c(group1, group2)],
        (getCoverage(o, type = "Cov") - getCoverage(o, type = "M"))[, c(group1, group2)])
      df$width <- NULL
      colnames(df) <- c("chr", "start", "end", "strand",
                        paste0("coverage", 1:nsample),
                        paste0("numCs", 1:nsample),
                        paste0("numTs", 1:nsample))
      df <- df[, c("chr", "start", "end", "strand",
                   sapply(1:nsample, function(x) paste0(c("coverage", "numCs", "numTs"), x)))]
      coverage.ind <- which(sapply(colnames(df), function(x) grepl("coverage", x)))
      names(coverage.ind) <- NULL
      obj <- new("methylBase", (df),
                 sample.ids = c(group1, group2),
                 assembly = "-",
                 context = "-",
                 treatment = c(rep(1, length(group1)), rep(0, length(group2))),
                 coverage.index = coverage.ind,
                 numCs.index = coverage.ind+1,
                 numTs.index = coverage.ind+2,
                 destranded = TRUE,
                 resolution = "base" )
      filter <- (rowSums(getCoverage(o)[, group1]) >= 1) & (rowSums(getCoverage(o)[, group2]) >= 1)
      errorInd <- rowSums(getCoverage(o, type = "M")) / rowSums(getCoverage(o))
      errorInd <- (errorInd == 0) | (errorInd == 1)
      obj_filtered <- obj[filter & !errorInd, ]
      invisible(capture.output(tmp <- calculateDiffMeth(obj_filtered)))
      gr <- GRanges(seqnames = tmp$chr,
                    IRanges(start = tmp$start, end = tmp$start),
                    pval = tmp$pval, methDiff = tmp$meth.diff/100)
      gr0 <- granges(o[which(errorInd & filter)])
      gr0$pval = 1
      gr0$methDiff = 0
      gr <- unlist(GRangesList(gr, gr0))
      gr <- sortSeqlevels(gr)
      gr <- sort(gr)
      gr
    }, mc.cores = mc.cores))
  }
  methcpObj <- new("MethCP",
                    group1 = group1,
                    group2 = group2,
                    test = test,
                    stat = stat)
  return(methcpObj)
}

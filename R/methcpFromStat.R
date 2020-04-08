
#' @title Create a \code{MethCP} object given the per-cytosine based tests.
#'
#' @description
#' Given the per-cytosine based \code{p}-values and effect sizes,
#' \code{MethCPFromStat} create a \code{MethCP} object for segmentation.
#'
#' @param data a data.frame or GPos or GRanges object.
#' @param test.name a character string containing the name of the test to be
#' performed per cytosine.
#' @param pvals.field A character vector of recognized names for the column 
#' (if `data` is a data.frame) or meta data column (is `data` is GPos or 
#' GRanges object) in `data` that contains the p-value.
#' @param effect.size.field A character vector of recognized names for the 
#' column (if `data` is a data.frame) or meta data column (is `data` is GPos or 
#' GRanges object) in `data` that contains the effect size.
#' @param seqnames.field A character vector of recognized names for the column
#'  in `data` that contains the chromosome name (sequence name) associated 
#'  with each position. Only the first name in seqnames.field that is 
#'  found in colnames(data) is used. If no one is found, then an error is 
#'  raised. This column is only used when `data` is a data.frame. 
#'  Otherwise, chromosome name is obtained from GPos or GRanges object.
#' @param pos.field A character vector of recognized names for the column
#'  in df that contains the position integer associated 
#'  with each position. Only the first name in pos.field that is 
#'  found in colnames(data) is used. If no one is found, then an error is 
#'  raised. This column is only used when `data` is a data.frame. 
#'  Otherwise, position is obtained from GPos or GRanges object.
#'
#' @return a \code{MethCP} object that is not segmented.
#'
#' @examples 
#' # ====== construct using data frame
#' data <- data.frame(
#'     chr = rep("Chr01", 5),
#'     pos = c(2, 5, 9, 10, 18),
#'     effect.size = c(1,-1, NA, 9, Inf),
#'     pvals = c(0, 0.1, 0.9, NA, 0.02))
#' obj <- MethCPFromStat(
#'     data, test.name="myTest", 
#'     pvals.field = "pvals",
#'     effect.size.field="effect.size",
#'     seqnames.field="chr",
#'     pos.field="pos"
#' )
#' # ====== construct using GRanges
#' library(GenomicRanges)
#' data <- GRanges(
#'     "Chr01", IRanges(c(2, 5, 9, 10, 18), c(2, 5, 9, 10, 18)),
#'     pvals=c(0, 0.1, 0.9, NA, 0.02), effect.size = c(1,-1, NA, 9, Inf))
#' obj <- MethCPFromStat(
#'     data, test.name="myTest", 
#'     pvals.field = "pvals",
#'     effect.size.field="effect.size"
#' )
#'
#' @export
MethCPFromStat <- function(
    data, test.name, pvals.field="pvals", effect.size.field="effect.size",
    seqnames.field=c(
        "seqnames", "seqname", "chromosome", "chrom",
        "chr", "chromosome_name", "seqid"), 
    pos.field="pos"){
    if (length(pvals.field) >= 1) {
        pvals.field = pvals.field[1]
    }
    if (length(effect.size.field) >= 1) {
        pvals.field = pvals.field[1]
    }
    if (length(seqnames.field) >= 1) {
        pvals.field = pvals.field[1]
    }
    if (length(pos.field) >= 1) {
        pvals.field = pvals.field[1]
    }
    if (is(data, 'GPos')){
        return(MethCP(
            test = test.name,
            group1 = "notApplicable",
            group2 = "notApplicable",
            chr = as.character(seqnames(data)), pos = as.integer(pos(data)), 
            pvals = data$pvals, effect.size = data$effect.size))
    }
    if (is(data, 'GRanges')){
        return(MethCP(
            test = test.name,
            group1 = "notApplicable",
            group2 = "notApplicable",
            chr = as.character(seqnames(data)), pos = as.integer(start(data)), 
            pvals = data$pvals, effect.size = data$effect.size))
    }
    if (is(data, 'data.frame')){
        if (seqnames.field == "seqnames"){
            seqname_idx = which(c(
                "seqnames", "seqname", "chromosome", "chrom",
                "chr", "chromosome_name", "seqid") %in% colnames(data))
            if (length(seqname_idx) == 0){
                stop("Chromosome name column not found.")
            } else {
                seqname_idx = seqname_idx[1]
            }
            seqnames.field = c(
                "seqnames", "seqname", "chromosome", "chrom",
                "chr", "chromosome_name", "seqid")[seqname_idx]
        }
        return(MethCP(
            test = test.name,
            group1 = "notApplicable",
            group2 = "notApplicable",
            chr = as.character(data[, seqnames.field]), pos = data$pos, 
            pvals = data$pvals, effect.size = data$effect.size))
    }
}

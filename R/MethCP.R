
setClassUnion("characterORnumeric", c("character", "numeric"))

#' @export
setClass("MethCP", representation(test = "character",
                                  stat = "GRangesList",
                                  group1 = "characterORnumeric",
                                  group2 = "characterORnumeric",
                                  segmentation = "GRanges"),
         prototype(test = NA_character_,
                   group1 = NA,
                   group2 = NA,
                   stat = GenomicRanges::GRangesList(list()),
                   segmentation = GenomicRanges::GRanges()))

#' @export
setMethod("show", signature("MethCP"), function(object){
  cat(paste0("MethCP object with ", length(object@stat), " chromosomes, ",
             sum(sapply(object@stat, length)),
             " methylation loci\n"))
  cat(paste0("test: ", object@test, "\n"))
  cat(paste0("group1: ", paste(object@group1, collapse = " "),
             "\ngroup2: ", paste(object@group2, collapse = " ")))
  if (length(object@segmentation) == 0){
    cat("\nhas not been segmented")
  } else {
    cat("\nhas been segmented")
  }
})

#' @export
setGeneric("segmentMethCP",
           function(methcp.object, bs.object,
                    region.test = c("fisher", "stouffer",
                                    "weighted-variance", "weighted-coverage"),
                    mc.cores = 1, min.width = 2, sig.level = 0.01, ...)
             standardGeneric("segmentMethCP"))
setMethod(
  "segmentMethCP", signature(methcp.object = "MethCP"),
  function(methcp.object, bs.object,
           region.test = c("fisher", "stouffer",
                           "weighted-variance", "weighted-coverage"),
           mc.cores = 1, min.width = 2, sig.level = 0.01, ...){

    object <- methcp.object
    if (!is(object, "MethCP")){
      stop("ERROR: Input must be an object of class \"MethCP\".")
    }
    if (object@test == "methylKit" &
        region.test %in% c("weighted-variance", "weighted-coverage")){
      stop("ERROR: can not apply weighted effect size method with methylKit.")
    }
    if (length(object@segmentation) != 0){
      a <- readline("Object has been segmented. Remove previous segmentation? (y/n) > ")
      if (a == "y") {
        message("Removed. Start running new segmentation ...")
        object@segmentation <- GRanges()
      } else {
        message("Stopped.")
        return(methcp.object)
      }
    }
    if (object@test == "methylKit") {
      object@stat <- GRangesList(lapply(names(object@stat), function(o){
        tmp <- object@stat[[o]]
        tmp$stat <- .pvalToStat(tmp$pval, tmp$methDiff)
        tmp
      }))
    }
    # calculate total coverage and methylated counts for each loci
    object@stat <- GRangesList(lapply(1:length(object@stat), function(o){
      tmp <- object@stat[[o]]
      ovrlp <- findOverlaps(granges(bs.object), tmp)
      tmp$CovGroup1 <- rowSums(
        getCoverage(bs.object)[ovrlp@from, object@group1])
      tmp$CovGroup2 <- rowSums(
        getCoverage(bs.object)[ovrlp@from, object@group2])
      tmp
    }))
    # segmentation
    segments <- mclapply(object@stat, function(o){
      cp.object <- CNA(o$stat,
                       chrom=as.vector(o@seqnames),
                       maploc=start(o),
                       data.type="logratio")
      invisible(capture.output(segment.cp.object <- segment(
        cp.object, verbose = 1, min.width = min.width, alpha = sig.level, ...)))
      return(segment.cp.object$output)
    }, mc.cores = mc.cores)
    segments <- do.call("rbind", segments)
    segments$ID <- NULL
    segments$seg.mean <- NULL
    colnames(segments)[4] <- "nC.valid"
    segments <- GRanges(segments)

    # calculate region summary
    ovrlp <- findOverlaps(granges(bs.object), segments)
    segments$nC <- as.numeric(table(ovrlp@to))
    M1 <- as.numeric(by(getCoverage(bs.object, type = "M")[
      ovrlp@from, object@group1], ovrlp@to, sum))
    M2 <- as.numeric(by(getCoverage(bs.object, type = "M")[
      ovrlp@from, object@group2], ovrlp@to, sum))
    Cov1 <- as.numeric(by(getCoverage(bs.object)[
      ovrlp@from, object@group1], ovrlp@to, sum))
    Cov2 <- as.numeric(by(getCoverage(bs.object)[
      ovrlp@from, object@group2], ovrlp@to, sum))
    segments$mean.diff <- M1/Cov1 - M2/Cov2
    segments$mean.cov <- (Cov1 + Cov2)/length(c(object@group1, object@group2))/segments$nC.valid

    # calculate region statistics
    ovrlp <- findOverlaps(unlist(object@stat), segments)
    if (region.test == "fisher"){
      segments$region.pval <- tapply(
        unlist(object@stat)$pval, ovrlp@to, .calcFisherPval)
    } else if (region.test == "stouffer"){
      # segments$region.pval <- as.numeric(by(
      #   unlist(object@stat)[, c("mu", "pval")], ovrlp@to,
      #   function(x) .calcStoufferPval(x$pval, x$mu)))
      segments$region.pval <- tapply(
        unlist(object@stat)$pval, ovrlp@to, .calcStoufferPvalOneSided)
    } else if (region.test == "weighted-variance"){
      segments$region.pval <- as.numeric(by(
        unlist(object@stat)@elementMetadata, ovrlp@to,
        function(x) .calcWeightedPval(x$mu, x$se, 1/x$se)))
    } else if (region.test == "weighted-coverage"){
      segments$region.pval <- as.numeric(by(
        unlist(object@stat)@elementMetadata, ovrlp@to,
        function(x) .calcWeightedPval(x$mu, x$se, x$CovGroup1 + x$CovGroup2)))
    }
    methcp.object@segmentation <- segments
    return(methcp.object)
  }
)

#' @export
setGeneric("getSigRegion",
           function(object, sig.level = 0.01, mean.coverage = 1,
                    mean.diff = 0.1, nC.valid = 10)
             standardGeneric("getSigRegion"))
setMethod(
  "getSigRegion", "MethCP",
  function(object, sig.level = 0.01, mean.coverage = 1,
           mean.diff = 0.1, nC.valid = 10){

    if (!is(object, "MethCP")){
      stop("ERROR: Input must be an object of class \"MethCP\".")
    }
    res <- as.data.frame(object@segmentation)
    res$strand <- NULL
    res$width <- NULL
    res <- res[res$region.pval <= sig.level & res$mean.cov >= mean.coverage &
                 abs(res$mean.diff) >= mean.diff & res$nC.valid >= nC.valid, ]
    res$mean.diff <- round(res$mean.diff, 4)
    res$mean.cov <- round(res$mean.cov, 4)
    return(res)
  })

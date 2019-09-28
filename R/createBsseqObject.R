#' @title Helper function to read text files and create a bsseq object.
#'
#' @usage
#' createBsseqObject(
#'     files, sample_names, 
#'     chr_col, pos_col, m_col, cov_col, header = TRUE)
#'
#' @description Create a bsseq object when the data for each sample is 
#' stored in a separate text file. 
#'
#' @param files a charactor vector of file names with full path to the file.
#' @param sample_names a charactor vector of sample names. It should have the
#' same length as files vector.
#' @param chr_col name or index of the chromosome column in data files.
#' @param pos_col name or index of the position column in data files.
#' @param m_col name or index of the methylated counts column.
#' @param cov_col name or index of the coverage counts column.
#' @param header a logical value indicating whether the file contains the 
#' names of the variables as its first line. dedault TRUE.
#'
#' @return a \code{MethCP} object that is not segmented.
#'
#' @examples
#' library(bsseq)
#' # The dataset is consist of 6 samples. 3 samples are H2A.Z mutant 
#' # plants, and 3 samples are controls.
#' sample_names <- c(
#'     paste0("control", seq_len(3)), 
#'     paste0("treatment", seq_len(3))
#' )
#' 
#' # Get the vector of file path and names 
#' raw_files <- system.file(
#'     "extdata", paste0(sample_names, ".txt"), package = "methcp")
#' 
#' # load the data
#' bs_object <- createBsseqObject(
#'     files = raw_files, sample_names = sample_names, 
#'     chr_col = 'Chr', pos_col = 'Pos', m_col = "M", cov_col = 'Cov')
#' 
#' @importFrom bsseq BSseq combine
#' @importFrom utils read.table
#' @importFrom S4Vectors to from
#'
#' @export
createBsseqObject <- function(
    files, sample_names, chr_col, pos_col, m_col, cov_col, header = TRUE)
{
    n_sample <- length(sample_names)
    bs_object_list <- list()
    for (i in seq_len(n_sample)) {
        dt <- read.table(
            files[i], stringsAsFactors = FALSE, header = header)
        bs_object_list[[i]] <- BSseq(
            chr = dt[, chr_col], 
            pos = as.numeric(dt[, pos_col]), 
            M = as.matrix(dt[, m_col]), 
            Cov = as.matrix(dt[, cov_col]), 
            sampleNames = sample_names[i])
    }
    
    bs_object <- do.call("combine", c(bs_object_list))
    return(bs_object)
}

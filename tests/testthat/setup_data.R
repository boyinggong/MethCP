
library(bsseq)
library(GenomicRanges)
library(IRanges)

.calc_jaccard <- function(set1, set2){
    return(length(intersect(set1, set2))/length(union(set1, set2)))
}

set.seed(0286374)

############################
# simulate two-group comparison data
############################

# Simulate a small dataset with 2000 cyotsine and 6 samples,
# 3 in the treatment group and 3 in the control group. The
# methylation ratio are generated using Binomial distribution
# with probability 0.3.
nC <- 2000
sim_cov <- rnbinom(6*nC, 5, 0.5) + 5
sim_M <- sapply(sim_cov, function(x) rbinom(1, x, 0.3))
sim_cov <- matrix(sim_cov, ncol = 6)
sim_M <- matrix(sim_M, ncol = 6)
# methylation ratios in the DMRs in the treatment group are
# generated using Binomial(0.7)
DMRs <- c(600:622, 1089:1103, 1698:1750)
sim_M[DMRs, 1:3] <- sapply(sim_cov[DMRs, 1:3], function(x) rbinom(1, x, 0.7))
# sample names
sample_names <- c(paste0("treatment", 1:3), paste0("control", 1:3))
colnames(sim_cov) <- sample_names
colnames(sim_M) <- sample_names

# create a bs.object
bs_object <- BSseq(gr = GRanges(
    seqnames = "Chr01",
    IRanges(start = (1:nC)*10, width = 1)),
    Cov = sim_cov, M = sim_M, sampleNames = sample_names)
DMRs_pos <- DMRs*10

############################
# simulate time-course data
############################


# Simulate a small dataset with 2000 cyotsine and 10 samples,
# 5 in the treatment group and 5 in the control group. The
# methylation ratio are generated using Binomial distribution
# with probability 0.3, 0.4, 0.5, 0.6 and 0.7 for 5 time points.
nC <- 2000
nsamples <- 5
sim_cov <- rnbinom(10*nC, 5, 0.5) + 5
sim_cov <- matrix(sim_cov, ncol = 10)
time_point <- rep(1:nsamples, 2)
ratios <- time_point/10 + 0.2
sim_M <- sapply(1:(2*nsamples), function(i){
    sapply(sim_cov[, i], function(j) rbinom(1, j, ratios[i]))
})
sim_M <- matrix(sim_M, ncol = 2*nsamples)
# methylation ratios in the DMRs in the treatment group are
# generated using Binomial(0.3)
DMRs <- c(600:622, 1089:1103, 1698:1750)
sim_M[DMRs, 1:5] <- sapply(sim_cov[DMRs, 1:5], function(x) rbinom(1, x, 0.3))
# sample names
sample_names <- c(
    paste0("treatment", 1:nsamples), paste0("control", 1:nsamples))
colnames(sim_cov) <- sample_names
colnames(sim_M) <- sample_names

# create a bs.object
bs_object_ts <- BSseq(
    gr = GRanges(seqnames = "Chr01", IRanges(
        start = (1:nC)*10, width = 1)),
    Cov = sim_cov, M = sim_M, sampleNames = sample_names)
DMRs_pos_ts <- DMRs*10
meta <- data.frame(
    Condition = rep(c("treatment", "control"), each = nsamples),
    SampleName = sample_names,
    Time = time_point)

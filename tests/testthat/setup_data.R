
library(bsseq)
library(GenomicRanges)
library(IRanges)

set.seed(0286374)

############################
# simulate two-group comparison data
############################

# Similate a small dataset with 11 cyotsine and 6 samples,
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
bs_object <- BSseq(gr = GRanges(seqnames = "Chr01",
                                IRanges(start = (1:nC)*10,
                                        width = 1)),
                   Cov = sim_cov, M = sim_M, sampleNames = sample_names)

############################
# simulate time-course data
############################

#
# # Similate a small dataset with 11 cyotsine and 6 samples,
# # 3 in the treatment group and 3 in the control group. The
# # methylation ratio are generated using Binomial distribution
# # with probability 0.3.
# nC <- 2000
# sim_cov <- rnbinom(10*nC, 5, 0.5) + 5
# sim_M <- sapply(sim_cov, function(x) rbinom(1, x, 0.3))
# sim_cov <- matrix(sim_cov, ncol = 10)
# sim_M <- matrix(sim_M, ncol = 6)
# # methylation ratios in the DMRs in the treatment group are
# # generated using Binomial(0.7)
# DMRs <- c(600:622, 1089:1103, 1698:1750)
# sim_M[DMRs, 1:3] <- sapply(sim_cov[DMRs, 1:3], function(x) rbinom(1, x, 0.7))
# # sample names
# sample_names <- c(paste0("treatment", 1:3), paste0("control", 1:3))
# colnames(sim_cov) <- sample_names
# colnames(sim_M) <- sample_names
#
# # create a bs.object
# bs_object <- BSseq(gr = GRanges(seqnames = "Chr01",
#                                 IRanges(start = (1:nC)*10,
#                                         width = 1)),
#                    Cov = sim_cov, M = sim_M, sampleNames = sample_names)


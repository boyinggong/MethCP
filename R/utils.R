
.pvalToStat <- function(pval, effect.size, jitter = 1e-16){
  pval[pval < jitter] <- 2*jitter
  stats <- qnorm(1-pval/2) * sign(effect.size)
  stats[stats == Inf] = max(stats)
  stats[stats == -Inf] = min(stats)
  return(stats)
}

.calcFisherStat <- function(pval, jitter = 1e-16){
  if (any(pval==0)){
    pval[pval==0] = jitter
  }
  return(sum(-2*log(pval)))
}

.calcFisherPval <- function(pval, ...){
  stat <- .calcFisherStat(pval, ...)
  return(1-pchisq(stat, 2*length(pval)))
}

.calcStoufferPval <- function(pval, effect.size, jitter = 1e-16){
  pval[pval==1] <- 1-jitter
  pval[pval==0] <- jitter
  tmp <- qnorm(1-pval/2)*sign(effect.size)
  stat <- sum(tmp)/sqrt(length(tmp))
  return((1-pnorm(abs(stat)))*2)
}

.calcStoufferPvalOneSided <- function(pval, jitter = 1e-16){
  pval[pval==1] <- 1-jitter
  tmp <- qnorm(1-pval)
  stat <- sum(tmp)/sqrt(length(tmp))
  return(1-pnorm(stat))
}

.calcWeightedStat <- function(effect.size, se, weights){
  stat.mean <- sum(effect.size * weights)
  stat.var <- sum(se^2 * weights^2)
  return(stat.mean / sqrt(stat.var))
}

.calcWeightedPval <- function(effect.size, se, weights){
  stat <- .calcWeightedStat(effect.size, se, weights)
  return((1-pnorm(abs(stat)))*2)
}

.asinTransform <- function(p) { asin(sqrt(p)) }


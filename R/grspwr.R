## Genetic Risk Score power simulator with support for imputed genotypes.
# The motivating idea was to check whether inclusion of (badly) imputed SNPs still improves power.
# Christian Theil Have, 2018.

# Generate permutation of genotypes that minimizes fn
# TODO: This is a candidate for Rcpp or similar because of R's copy-on-write semantics this is heavy
genetic.optim <- function(genotypes,fn) {
  randomswap <- function(genotypes) {
    i <- sample(length(genotypes),1)
    j <- sample(which(genotypes!=genotypes[i]),1)
    tmp <- genotypes[i]
    genotypes[i] <- genotypes[j]
    genotypes[j] <- tmp
    genotypes
  }
  current.score <- fn(genotypes)
  max.iter=length(genotypes)*2
  iter <- 0
  repeat {
    permuted <- randomswap(genotypes)
    permuted.score <- fn(permuted)
    if(permuted.score < current.score) {
      current.score <- permuted.score
      genotypes <- permuted
    }
    iter <- iter + 1
    if (iter >= max.iter)
      break
  }
  genotypes
}

#' pairs diploid-alleles assorted according hardy-weinberg equibrilium principles iff alleles are in random order
alleles2genotypes <- function(alleles)
  alleles[seq(1,length(alleles)-1,2)] + alleles[seq(2,length(alleles),2)]

sample.genotypes.normal <- function(phenotypes, beta, n, maf) {
  num.variant.alleles <- round(maf*2*n)
  allele.pool <- sample(c(rep(T,num.variant.alleles), rep(F, (2*n)-num.variant.alleles)), 2*n)
  genotypes <- alleles2genotypes(allele.pool)
  optim.func <- function(x) abs(beta-lm(phenotypes~x)$coefficients[2])
  genetic.optim(genotypes, optim.func)
}

#' Scale a numeric vector to values between 0 and 2.
#' @param unscaled The unscaled input
scale.to.dosages <- function(unscaled) {
  minima <- min(unscaled)
  maxima <- max(unscaled)
  2 * (unscaled + abs(minima)) / (maxima+abs(minima) - (minima+abs(minima)))
}

#' Sample dosages from genotypes using fast binary search for r-squared
#' It will result in  approximately correct r-squared, but have a tendency to "look" to nice..
sample.dosages <- function(genotypes,info,epsilon=0.01) {
  jitter.amount = 0
  dosages <- genotypes
  repeat {
    rsq <- summary(lm(genotypes ~ dosages))$r.squared
    diff <- abs(rsq - info)
    if(diff < epsilon)
      break
    else if (rsq > info)
      jitter.amount <- jitter.amount + (diff / 2)
    else
      jitter.amount <- jitter.amount - (diff / 2)
    dosages <- scale.to.dosages(jitter(genotypes,amount=jitter.amount))
  }
  list(genotypes = genotypes, dosages = dosages,r.squared = rsq)
}

#' pwr.grs calculates the power to detect of a GRS to detect a signal at given significance level
#' using a GRS  weights in a certain sample size and with given imputatation quality.
#' Note: It assumes that the weights used are accurate estimates of SNP effects.
#'
#' @param snps: A dataframe with the following columns:
#' \itemize{
#' \item SNP: Unique name of SNP
#' \item Weight: a weight associated with the _minor_ allele of SNP. Weights should be normalized   relative to a standard normal distribution (mean 0, variance 1).
#' \item EAF: Effect allele frequency
#' \item Estimated R-squared correlation of imputation with actual genotype.For imputed genotypes this will be [0-1] and for directly genotyped SNPs it should be 1.
#' }
#' @param n The number of individuals to include in the GRS
#' @param max.iter  The number of iterations to use for the power calculation. This affects the precision of the calculated power. Should be at least 100.
#' @param popsize The size of the population sample from which is sampled. Should be larger than N. Affects accuracy of power estimate.
################################################################################################
pwr.grs <- function(snps, n, alpha = 0.05, max.iter=1000, popsize=n*2) {
  # Sample a population of individuals with normally distributed phenotypes
  phenotypes <- rnorm(popsize)
  genotypes <- list()
  dosages <- list()

  print("Sampling population data")
  for(snp_idx in 1:nrow(snps)) {
    genotypes[[snp_idx]] <- sample.genotypes.normal(phenotypes,snps[snp_idx,]$Weight,popsize, snps[snp_idx,]$EAF)
    dosages[[snp_idx]] <- sample.dosages(genotypes[[snp_idx]], snps[snp_idx,]$INFO)$dosages
  }

  print("Sampling test pvalues")
  pvalues <- sapply(1:max.iter, function(i) {
    popsample <- sample(1:popsize,n)
    sample_phenotypes <- phenotypes[popsample]
    grs <- rep(0,n)
    for(snp_idx in 1:nrow(snps)) {
      w <- snps[snp_idx,]$Weight
      grs <- grs + dosages[[snp_idx]][popsample] * w
    }
    summary(lm(sample_phenotypes~grs))$coefficients[2,4]
  })

  # Power is the fraction of succesfull tests
  sum(pvalues<alpha) / max.iter
}

example.grs.data<- data.frame(
  SNP = c("rs1", "rs2", "rs3"),
  Weight = c(0.01, 0.01, 0.2),
  EAF = c(0.42, 0.42, 0.42),
  INFO= c(0.9, 0.9, 0.2))

#pwr.grs(example.grs.data,n=400, alpha = 0.05)

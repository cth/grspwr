## Genetic Risk Score power simulator with support for imputed genotypes.
# The motivating idea was to check whether inclusion of (badly) imputed SNPs still improves power.
# Christian Theil Have, 2018.

#' Swap two genotypes
#'
#' Swaps a random pair of indices in a vector of genotypes
#' @param genotypes A vector of genotypes
randomSwap <- function(genotypes) {
  # TODO: This is a candidate for Rcpp or similar because of R's copy-on-write semantics this is heavy
  i <- sample(length(genotypes),1)
  j <- sample(which(genotypes!=genotypes[i]),1)
  tmp <- genotypes[i]
  genotypes[i] <- genotypes[j]
  genotypes[j] <- tmp
  genotypes
}

#' Permute genotypes subject to opmization
#'
#' Generates a permutation of genotypes that minimizes fn.
#'
#' @param genotypes A vector of genotypes
#' @param fn A function which takes a vector of genotypes as first argument and return a score.
optimizedPermutation <- function(genotypes,fn) {
  current.score <- fn(genotypes)
  max.iter=length(genotypes)*2
  iter <- 0
  repeat {
    permuted <- randomSwap(genotypes)
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

#' Combine single alleles into genotypes
#'
#' pairs diploid-alleles assorted according hardy-weinberg equibrilium principles iff alleles are in random order
allelesToGenotypes <- function(alleles)
  alleles[seq(1,length(alleles)-1,2)] + alleles[seq(2,length(alleles),2)]


#' Generate genotypes with a given effect on phetypes
#'
#' Make a vector of genotypes with a certain frequency and such that the
#' per-allele effect one the given phenotypes is as specified.
#'
#' @param phenotypes a vector of phenotypes (assumed to be normally distributed)
#' @param beta desired effect size
#' @param n sample size
#' @param maf minor (effect) allele frequency
sampleGenotypes <- function(phenotypes, beta, n, maf) {
  numberOfVariantAlleles <- round(maf*2*n)
  allelePool <- sample(c(rep(T,numberOfVariantAlleles), rep(F, (2*n)-numberOfVariantAlleles)), 2*n)
  genotypes <- allelesToGenotypes(allelePool)
  optimizationFunction <- function(x) abs(beta-lm(phenotypes~x)$coefficients[2])
  optimizedPermutation(genotypes, optimizationFunction)
}

#' Scale a numeric vector to values between 0 and 2.
#' @param unscaled The unscaled input
scaleDosages <- function(unscaled) {
  minima <- min(unscaled)
  maxima <- max(unscaled)
  2 * (unscaled + abs(minima)) / (maxima+abs(minima) - (minima+abs(minima)))
}

#' Sample dosages from genotypes
#'
#' Creates random dosages from genotypes using fast binary search for r-squared.
#' It will result in  approximately correct r-squared.
#'
#' @param genotypes a vector of numeric genotypes (0,1,2)
#' @param info Desired \eqn{R^2} correlation between dosages and genotypes
#' @param epsilon A halting parameter: When the total absolute difference
#' between desired and obtained \eqn{R^2} reaches epsilon the dosages are returned.
sampleDosages <- function(genotypes,info,epsilon=0.01) {
  jitterAmount = 0
  dosages <- genotypes
  if (info < 1) {
    repeat {
      rsq <- summary(lm(genotypes ~ dosages))$r.squared
      diff <- abs(rsq - info)
      if(diff < epsilon)
        break
      else if (rsq > info)
        jitterAmount <- jitterAmount + (diff / 2)
      else
        jitterAmount <- jitterAmount - (diff / 2)
      dosages <- scaleDosages(jitter(genotypes,amount=jitterAmount))
    }
  }
  list(genotypes = genotypes, dosages = dosages,r.squared = rsq)
}

#' Sample a population with dosages for each SNP related to a phenotype with effects as indicated by SNP weights.
#' This resulting object can be used as input to the grspwr calculation.
#'
#' @param snps: A dataframe with the following columns:
#' \itemize{
#' \item SNP: Unique name of SNP
#' \item Weight: a weight associated with the _minor_ allele of SNP. Weights should be normalized   relative to a standard normal distribution (mean 0, variance 1).
#' \item EAF: Effect allele frequency
#' \item Estimated R-squared correlation of imputation with actual genotype.For imputed genotypes this will be [0-1] and for directly genotyped SNPs it should be 1.
#' }
#' @param n The number of individuals to include in the GRS
samplePopulationDosages <- function(snps, n) {
  phenotypes <- rnorm(n)
  genotypes <- list()
  dosages <- list()

  # We should do this with mclapply instead
  for(snp_idx in 1:nrow(snps)) {
    genotypes[[snp_idx]] <- sampleGenotypes(phenotypes,snps[snp_idx,]$Weight,n, snps[snp_idx,]$EAF)
    dosages[[snp_idx]] <- sampleDosages(genotypes[[snp_idx]], snps[snp_idx,]$INFO)$dosages
  }
  dat <- list(phenotypes=phenotypes, snps=snps$SNP, genotypes=genotypes, dosages=dosages)
  class(dat) <- "population"
  dat
}


#' GRS power calculation
#'
#' Calculates the power to detect of a GRS to detect a signal at given significance level
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
#'
#' @examples
#' exampleGRS <- data.frame(SNP = c("rs1", "rs2", "rs3"),
#'                               Weight = c(0.01, 0.01, 0.2),
#'                               EAF = c(0.42, 0.42, 0.42),
#'                               INFO= c(0.9, 0.9, 0.2))
#' grspwr(exampleGRS,n=400, alpha = 0.05)
################################################################################################
grspwr <- function(snps, n, alpha = 0.05, max.iter=1000, popsize=n*2) {
  attach(samplePopulationDosages(snps,popsize))

  pvalues <- rep(NA,max.iter)
  betas <- rep(NA,max.iter)

  for(i in 1:max.iter) {
    popsample <- sample(1:popsize,n)
    sample_phenotypes <- phenotypes[popsample]
    grs <- rep(0,n)
    for(snp_idx in 1:nrow(snps)) {
      w <- snps[snp_idx,]$Weight
      grs <- grs + dosages[[snp_idx]][popsample] * w
    }
    model<-lm(sample_phenotypes~grs)
    betas[i]<- unname(model$coef[2])
    pvalues[i] <- summary(model)$coefficients[2,4]
  }

  # Power is the fraction of succesfull tests
  list(power = sum(pvalues<alpha) / max.iter,
       pvalues = pvalues,
       betas = betas)
}

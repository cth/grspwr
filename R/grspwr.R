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

#' Generate a genotype vector of length n with a given minor allele frequency (maf)
sampleGenotypes <- function(n,maf) {
  numberOfVariantAlleles <- round(maf*2*n)
  allelePool <- sample(c(rep(T,numberOfVariantAlleles), rep(F, (2*n)-numberOfVariantAlleles)), 2*n)
  allelesToGenotypes(allelePool)
}

#' Generate genotypes with a given effect on phetypes
#'
#' Make a vector of genotypes with a certain frequency and such that the
#' per-allele effect one the given phenotypes is as specified.
#'
#' @param phenotypes a vector of phenotypes (assumed to be normally distributed)
#' @param beta desired effect size
#' @param n sample size
#' @param maf minor (effect) allele frequency
sampleCorrelatedGenotypes <- function(phenotypes, beta, n, maf) {
  optimizationFunction <- function(x) abs(beta-lm(phenotypes~x)$coefficients[2])
  optimizedPermutation(sampleGenotypes(n,maf), optimizationFunction)
}

#' Generate genotypes with a given effect on phetypes (C version)
#'
#' Make a vector of genotypes with a certain frequency and such that the
#' per-allele effect one the given phenotypes is as specified.
#'
#' @param phenotypes a vector of phenotypes (assumed to be normally distributed)
#' @param beta desired effect size
#' @param n sample size
#' @param maf minor (effect) allele frequency
sampleCorrelatedGenotypesC <- function(phenotypes,beta,n,maf) {
  permute_genotypes_C(beta,sampleGenotypes(n,maf),phenotypes)
}

#' Scale a numeric vector to values between 0 and 2.
#' @param unscaled The unscaled input
#' @param unscaled The unscaled input
#' @param unscaled The unscaled input
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
      # This could be done inline:
      dosages <- scaleDosages(jitter(genotypes,amount=jitterAmount))
    }
    dosages
  } else {
    dosages
  }
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
samplePopulationDosages <- function(snps, n, phenotypes = rnorm(n)) {
  genotypes <- list()
  dosages <- list()

  # We should do this with mclapply instead
  for(snp_idx in 1:nrow(snps)) {
    genotypes[[snp_idx]] <- sampleCorrelatedGenotypesC(phenotypes,snps[snp_idx,]$Weight,n, snps[snp_idx,]$EAF)
    dosages[[snp_idx]] <- sampleDosages(genotypes[[snp_idx]], snps[snp_idx,]$INFO)
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
#' \item WEIGHT: a weight associated with the _minor_ allele of SNP. Weights should be normalized   relative to a standard normal distribution (mean 0, variance 1).
#' \item EAF: Effect allele frequency
#' \item Estimated R-squared correlation of imputation with actual genotype.For imputed genotypes this will be [0-1] and for directly genotyped SNPs it should be 1.
#' }
#' @param n The number of individuals to include in the GRS
#' @param max.iter  The number of iterations to use for the power calculation. This affects the precision of the calculated power. Should be at least 100.
#' @param popsize The size of the population sample from which is sampled. Should be larger than N. Affects accuracy of power estimate.
#'
#' @examples
#' exampleGRS <- data.frame(SNP = c("rs1", "rs2", "rs3"),
#'                               WEIGHT = c(0.01, 0.01, 0.2),
#'                               EAF = c(0.42, 0.42, 0.42),
#'                               INFO= c(0.9, 0.9, 0.2))
#' grspwr(exampleGRS,n=400, alpha = 0.05)
################################################################################################
grspwr <- function(snps, n, alpha = 0.05, max.iter=1000, popsize=max(n)^2) {
  print(n)
  cat("Generating population phenotype and genotypes..\n")
  pop <- samplePopulationDosages(snps,popsize)

  cat("Subsampling population calculate power..\n")

  powers <- lapply(n,function(n_i) {
    pvalues <- rep(NA,max.iter)
    betas <- rep(NA,max.iter)
    for(i in 1:max.iter) {
      popsample <- sample(1:popsize,n_i)
      sample_phenotypes <- pop$phenotypes[popsample]
      grs <- rep(0,n_i)
      for(snp_idx in 1:nrow(snps)) {
        w <- snps[snp_idx,]$Weight
        grs <- grs + pop$dosages[[snp_idx]][popsample] * w
      }
      model<-lm(sample_phenotypes~grs)
      betas[i]<- unname(model$coef[2])
      pvalues[i] <- summary(model)$coefficients[2,4]
    }

    # Power is the fraction of succesfull tests
    pwr <- list(
         n = n_i,
         power = sum(pvalues<alpha) / max.iter,
         pvalues = pvalues,
         betas = betas)
    class(pwr) <- "grspwr"
    pwr
  })

  if (length(powers)==1) {
      powers[[1]]
  } else {
    powers
  }
}

#' Generates simulated pvalues for ranges of the exponent e in the equation:
#'   $w_i = w * (R^2)^e$ for all i
#' where $w_i$ is used as weights in such that the average p-value is minimized when testing the GRS
#' on simulated data.
#'
#' @param snps: A dataframe with the following columns:
#' \itemize{
#' \item SNP: Unique name of SNP
#' \item WEIGHT: a weight associated with the _minor_ allele of SNP. Weights should be normalized   relative to a standard normal distribution (mean 0, variance 1).
#' \item EAF: Effect allele frequency
#' \item Estimated R-squared correlation of imputation with actual genotype.For imputed genotypes this will be [0-1] and for directly genotyped SNPs it should be 1.
#' }
#' @param n The number of individuals to include in the GRS
#' @param popsize The size of the population sample from which is sampled. Should be larger than N. Affects accuracy of power estimate.
#' @param exponents A vector of exponent values to be tried. By default 20 different exponents in the range 0.01-2.0 are tried.
#' If the range does not include the optimal exponent, then it can be extended. Similarly, it can be narrowed
#' to achieve a finer resolution.
testAdjustmentExponents <- function(snps, n, popsize=n^2, exponents=10*1:20/100) {
  cat("Generating population phenotype and genotypes..\n")
  pop <- samplePopulationDosages(snps,popsize)

  betas <- list()
  pvalues <- list()

  cat("Generating powers..\n")
  for(exp_idx in 1:length(exponents)) {
    exponent <- exponents[exp_idx]
    cat(paste("exponent: ", exponent, "\n"))
    betas[[exp_idx]] <- c()
    pvalues[[exp_idx]] <- c()
    for(i in 1:1000) {
      popsample <- sample(1:popsize,n)
      sample_phenotypes <- pop$phenotypes[popsample]
      grs <- rep(0,n)
      for(snp_idx in 1:nrow(snps)) {
        w <- snps[snp_idx,]$Weight * snps[snp_idx,]$INFO^exponent
        grs <- grs + pop$dosages[[snp_idx]][popsample] * w
      }
      model<-lm(sample_phenotypes~grs)
      if(i==1) {
        betas[[exp_idx]] <- unname(model$coef[2])
        pvalues[[exp_idx]] <- summary(model)$coefficients[2,4]
      } else {
        betas[[exp_idx]] <- c(betas[[exp_idx]], unname(model$coef[2]))
        pvalues[[exp_idx]] <- c(pvalues[[exp_idx]], summary(model)$coefficients[2,4])
      }
   }
  }

  pvalues_all_points <- unlist(pvalues)
  exponents_all_points <- as.vector(matrix(rep(exponents,1000),nrow=1000,byrow=T))

  data.frame(exponents = exponents_all_points, pvalues = pvalues_all_points)
}

adjustWeights <- function(snps, adj.exp=0) {
  snps$Weight <- snps$Weight * snps$INFO^adj.exp
  snps
}

plot.AdjustWeights <- function(adjw) {
  spline=smooth.spline(x = adjw$exponents, y = unlist(adjw$pvalues))
  plot(spline, xlab="exponent c", ylab="mean p-value", main="Imputation adjustment", type="l")
}


#' autoAdjustWeights finds an optimal adjustment the weights to
#' compensate for markers with suboptimal imputation quality.
#' The adjustment is the exponent `e` in the model:
#'   $w_i = w * (R^2)^e$ for all i
#' The model is used to adjust each of the i weights. Note that
#' weight adjustment occurs when the $R^2$ is less than 1.
#' @param snps: A dataframe with the following columns:
#' \itemize{
#' \item SNP: Unique name of SNP
#' \item WEIGHT: a weight associated with the _minor_ allele of SNP. Weights should be normalized   relative to a standard normal distribution (mean 0, variance 1).
#' \item EAF: Effect allele frequency
#' \item Estimated R-squared correlation of imputation with actual genotype.For imputed genotypes this will be [0-1] and for directly genotyped SNPs it should be 1.
#' }'
#' @param n The number of individuals to include in the GRS
#' @param popsize The size of the population sample from which is sampled. Should be larger than N. Affects accuracy of power estimate.
#' The function returns a new data frame identicals to snps except that the weights $w_i$ will be
#' replaced with the updated weights
autoAdjustWeights <- function(snps, n, popsize=n^2) {
  max_exponent <- 1.0
  # TODO: We could test whether p-vlaue means are significantly non-random
  e_vs_p <- data.frame(exponents=c(),pvalues=c())
  stepsize <- 0.05
  minima <- NA
  repeat {
    exponents=seq(stepsize,max_exponent,stepsize)
    e_vs_p <- rbind(e_vs_p, testAdjustmentExponents(snps, n, popsize, exponents))
    # Try fit a spline and find minima
    spline <- smooth.spline(x = e_vs_p$exponents, y = e_vs_p$pvalues)
    minima <- spline$fit$knot[which(spline$fit$coef == min(spline$fit$coef))]
    # Did we get a non-random fit? Mann-Kendall Rank Test of randomness..
    p_nonrandomfit=cor.test(spline$fit$coef,order(spline$fit$coef), method="kendall")$pvalue
    # Minima should not be in either end of the spline. That would likely mean that
    # we have not sampled the range in which the true minima lies. Hence, expand range and re-iterate..
    valid_minima <- (split$fit$knot[1] != minima) & (split$fit$knot[length(split$fit$knot)] != minima)

    if (p_nonrandomfit < 0.01 && valid_minima)
      break

    # We need to check that
    stepsize <- stepsize * 2
    max_exponent <- max_exponent * 2
  }

  # Adjust the weights
  snps$WEIGHT <- snps$WEIGHT * snps$INFO^minima

  snps
}

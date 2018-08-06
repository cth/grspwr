# grspwr

grspwr is a power calculator for Genetic Risk Scores (GRS) constructed from imputed genotypes.
It calculates power given a set of SNPs, their effect allele frequency, their presumed effect and the  $R^2$ measure of imputation quality.

To install
```{r}
devtools::install_github("cth/grspwr")
```


## Short tutorial introduction

In this tutorial introduction we showcase the two major applications of this package. Power calculation and Weight adjustment based on imputation quality.

### Power Calculation 

The power calculation is based on a simulation that can be broken into several steps. First, we sample a population of individuals with a normally distribution random phenotype, genotypes for each SNP in the wGRS such that the per-allele effect on the phenotype mimics the weight specified in  the wGRS. In addition, for each SNP we sample dosages with a specific R<sup>2</sup> correlation to the generated genotypes, mimicking the uncertainty of the genotype imputation. 

Given this population we can subsample to a specific sample size and test whether an alpha-level significant correlation exists between GRS constructed from the dosages and the phenotypes in the particular subsample. By repeating this process a large number of times we can estimate the fraction of times one a population subsample of the specified size achieves alpha-level significance, i.e., the power.

This form of power calculation is initiated using `grspwr` command. Suppose that we represent the genetic risk score as a data frame, which is layed out like this:

| SNP | WEIGHT | EAF | INFO |
| :--- | ---: | ---: | ---: |
| rs13387838 | 0.139 | 0.04 | 0.77 | 
| rs7550711	| 0.105 | 0.04 | 1.00 |
| rs4854349 | 0.090 | 0.83 | 0.99 |
| rs543874 | 0.077 | 0.20 | 1.00 |
| rs12429545 | 0.076 | 0.13 | 1.00 |
| rs11676272 | 0.068 | 0.46 | 0.96 |
| rs13130484 | 0.067 | 0.44 | 1.00 |
| rs7132908 | 0.066 | 0.39 | 1.00 |
| rs987237 | 0.062 | 0.19 | 1.00 |
| rs1421085 | 0.059 | 0.41 | 1.00 |
| rs6567160 | 0.050 | 0.23 | 0.99 |
| rs12041852 | 0.046 | 0.46 | 0.99 |
| rs8092503 | 0.045 | 0.27 | 0.99 |
| rs13253111 | 0.042 | 0.57 | 0.99 |
| rs3829849 | 0.041 | 0.36 | 0.97 |

Assume that this dataframe is loading into the variable `snps`. The SNP column uniquely identifies each SNP in the score. The WEIGHT indicates the weight of the SNP, which is usually beta/OR from GWAS summary statitics for the target trait. EAF is the frequency of the effect allele in _your_ data. 
The INFO column contains imputation quality esimates (estimate of R<sup>2</sup> between imputed genotypes and true genotypes) in _your_ data. 

To calculate power for this risk score we also need to state how many individuals we expect to have, what level of significance we are happy with. Suppose we have 500 individuals and would happy if the risk associates with out target trait a p<0.05 significance. To calculate power bases on 100 simulations we would run: 

```
n <- 500
alpha <- 0.05
num_simulations <- 100
result<-grspwr(snps, n, alpha, num_simulations)
```

This command returns a list that contain a vector of pvalues observed in the simulation (`result$pvalues`) and a field of the proportion of p-values less than 0.05 (`results$power`).

Note that the number simulations effectively limits the precision of `results$power`. If you use use a lower `alpha` value you also increase the number of simulations to accomodate this.
Another aspect that limits the precision is the assumptions about the population size, which can be specified to `grspwr` through the `popsize` argument. This does not need to actually be the size of the population, but merely needs to be large enough to approximate the variation in the population. By default `grspwr` assumes the `popsize` to be n<sup>2</sup>`. For  for very large values of n, this setting may become impractical. Depending on the heterogeneity of your population and the number of SNPs in your risk score you way wish to tune this parameter. By plotting the parameter against the obtained power estimates, it is possible to determine a point at which increase the popsize argument is subject to dimishing returns. 

### Weight adjustment

Documentation will follow soon

# grspwr

grspwr is a power calculator for Genetic Risk Scores (GRS) constructed from imputed genotypes.
It calculates power given a set of SNPs, their effect allele frequency, their presumed effect and the  $R^2$ measure of imputation quality.

To install
```{r}
devtools::install_github("cth/grspwr")
```

## How it works
The power calculation is based on a simulation that can be broken into several steps. First, we sample a population of individuals with a normally distribution random phenotype, genotypes for each SNP in the wGRS such that the per-allele effect on the phenotype mimics the weight specified in  the wGRS. In addition, for each SNP we sample dosages with a specific R2 correlation to the generated genotypes, mimicking the uncertainty of the genotype imputation. 

Given this population we can subsample to a specific sample size and test whether an alpha-level significant correlation exists between GRS constructed from the dosages and the phenotypes in the particular subsample. By repeating this process a large number of times we can estimate the fraction of times one a population subsample of the specified size achieves alpha-level significance, i.e., the power. 

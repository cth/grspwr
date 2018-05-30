#include <Rcpp.h>
using namespace Rcpp;

// This is a more efficient implementation of the code to shuffle
// genotypes so that they have the desired effect on the phenotype.
// Its written in C because we wish to exploit inplace swaps etc

// A summary data structure to keep track of counts, sums and means for each genotype class
struct PopDat {
  int genotype_cnts[3];
  double cumsums[3];
  double means[3];
};

void update_means(struct PopDat* pd) {
   for(int i = 0; i < 3; i++) {
    pd->means[i] = pd->cumsums[i] / pd->genotype_cnts[i];
  }
}

/**
 * swap_samples and update counts / cumsums / means
 */
void swap_samples(struct PopDat* pd, int a, int b, IntegerVector genotypes, NumericVector phenotypes) {
  int tmp;

  pd->cumsums[(int)(genotypes[a])] -= phenotypes[a];
  pd->cumsums[(int)(genotypes[b])] -= phenotypes[b];
  pd->cumsums[(int)(genotypes[a])] += phenotypes[b];
  pd->cumsums[(int)(genotypes[b])] += phenotypes[a];

  tmp = genotypes[a];
  genotypes[a] = genotypes[b];
  genotypes[b] = tmp;

  update_means(pd);
}

double effect_difference(struct PopDat* pd, double expected_effect) {
  double effect1 =  pd->means[1] - pd->means[0];
  double effect2 =  pd->means[2] - pd->means[1];
  return pow(effect1-expected_effect,2.0) + pow(effect2-expected_effect,2.0);
}

// [[Rcpp::export]]
IntegerVector permute_genotypes_C(double beta, IntegerVector genotypes, NumericVector phenotypes) {
  struct PopDat pd;
  int n = phenotypes.size();

  NumericVector out(n);

  memset(&pd,0,sizeof(struct PopDat));

  // Initialize counts
  for(int i = 0; i < n; i++) {
    pd.genotype_cnts[(int)genotypes[i]]++;
    pd.cumsums[(int)genotypes[i]] += phenotypes[i];
  }
  update_means(&pd);

  for(int i = 0; i < n*2; i++) {
    int s1 = rand() % n;
    int s2 = rand() % n;
    double error_before = effect_difference(&pd,beta);
    swap_samples(&pd,s1,s2,genotypes,phenotypes);
    if ( effect_difference(&pd,beta) > error_before)
      swap_samples(&pd, s1,s2, genotypes, phenotypes); // swap back

    //if (i % 100 == 0)
    //  printf("%f %f %f\n", pd.means[0], pd.means[1], pd.means[2]);
  }

  /*
  memset(&pd,0,sizeof(struct PopDat));
  for(int i = 0; i < n; i++) {
    pd.genotype_cnts[(int)genotypes[i]]++;
    pd.cumsums[(int)genotypes[i]] += phenotypes[i];
  }
  update_means(&pd);
  printf("%f %f %f\n", pd.means[0], pd.means[1], pd.means[2]);
  */

  return clone(genotypes);
}

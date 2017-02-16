# Permutation procedure from Cardoso-Moreira et al. (2016) Genome Research

The permutation procedure described in the methods of that paper claims that it controls the genome-wide FDR.  

Here, we use simple simulations to show that the FDR is **not** properly-controlled, and that the procedure has a
similar FDR to a per-gene test. 

## The model

The parameters are:

* $n$ is the sample size
* $nsites$ are the number of sites
* $a,b$ are the coefficients of the Beta distribution determining the sample frequency of the "presence" genotype.
* The mean trait value ("expression level") of the "absence" genotype is 1, and the mean trait value of the "presence"
  genotype is given by the parameter $upgreg$.
* $tpr$ is the true positive rate, e.g., the probability that a "presence" genotype has effect size $upgreg$. **Note:**
  this is applied *per site* and not *per sample*.
* Trait values for each individual at each site are given by the mean values described above plus a Gaussian noise term
  with mean 0 and standard deviation $noise$.

## Source files

* sim.R defines the functions relevant to carry out the simulation
* test.R simulates the specific case of 1,000 "CNVs" in a sample size of 20.  The sample frequencies of each CNV are
  i.i.d and follow a Beta(1,5) distribution, crudely mimicking a population-genetic sample.  A real sample would have an
  SFS shifted towards even rarer variants, resulting in further loss of power and higher FDR/FPR when the sample sizes
  are small.
* Makefile runs the workflow

### Statistical tests

The sim.R source file conducts the following three tests:

* A regression of normalized "expression level" (z-score) onto genotype. We call these "single-marker" tests. A
  single-marker test is declared "significant" if the regression $p <= 0.05$.
* Asking whether or not the median difference in expression level between "presence" and "absence" genotypes is outside
  of the 2.5 and 97.5 quantiles of the permutation distribution.  The permutation distribution is obtained using the
  procedure outlined in Cardoso-Moreira et al.
* Asking whether the single-marker test is significant at a "q-value" (Storey 2002, J. R. Statist. Soc. B 64: 479) of
  0.05 or smaller.  The simulation procedure here conforms to the assumptions of this FDR procedure.

For all of the above, we know which loci are true positives because that is recorded in the simulation procedure.  This
allows us to ask which fraction of loci "discovered" are **not** in the set of true positives.  This proportion is, by
definition, the FDR.

# A proper FDR

A controlled FDR remains constant even when the parameters change (true positive rate, etc.).

# Results

See the pdf in the repository.  The take-home message is that the FDR for the single-marker test and the permutation
procedure are very similar.  Further, the FDR for these two procedures depend on the parameters.  Thus, they are not
proper FDR.  Rather, they are single-marker tests.  As expected, only the q-value procedure used here qualifies as a proper FDR.

The permutation seems to control the per-marker FDR slightly better than the regression.  A reasonable hypothesis is
that the effect is due to the stability of the median vs the mean (the latter is what the regression is testing).

## Implications

In reality, the FDR of an experiment is a complex function of the true positive rate and the power.  The power is
affected by the distribution of effect sizes.  All of the results here for the single-marker regression are totally
expected from the GWAS literature (which is essentially what our study and Cordosa-Moreira are doing): small sample
sizes lead to false positives and it is critical to control the genome-wide FDR.

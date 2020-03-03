
# mmcm package


## Modified Maximum Contrast Method

An implementation of modified maximum contrast methods and
the maximum contrast method: Functions 'mmcm.mvt' and 'mcm.mvt'
give P-value by using randomized quasi-Monte Carlo method with
'pmvt' function of package 'mvtnorm', and 'mmcm.resamp' gives
P-value by using a permutation method.

- Sato Y, Laird NM, Nagashima K, et al. A new statistical screening approach for finding pharmacokinetics-related genes in genome-wide studies. *The Pharmacogenomics Journal* 2009; **9**(2): 137--146. [doi:10.1038/tpj.2008.17](https://doi.org/10.1038/tpj.2008.17).
- Nagashima K, Sato Y, Hamada C. A modified maximum contrast method for unequal sample sizes in pharmacogenomic studies. *Statistical Applications in Genetics and Molecular Biology* 2011; **10**(1): Article 41. [doi:10.2202/1544-6115.1560](https://doi.org/10.2202/1544-6115.1560).
- Yoshimura I, Wakana A, Hamada C. A performance comparison of maximum contrast methods to detect dose dependency. *Drug Informatoion Journal* 1997; **31**(2): 423--432. [doi:10.1177/009286159703100213](https://doi.org/10.1177/009286159703100213).


## Installation

Get the released version from CRAN:

```R
install.packages("mmcm")
```

Or the development version from github:

```R
# install.packages("devtools")
devtools::install_github("nshi-stat/mmcm")
```


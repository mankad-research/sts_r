# sts_r
## sts: An R Package for the Structural Topic and Sentiment-Discourse Model

Authors: [Shawn Mankad](http://mankad-research.github.io) and [Li Chen](https://business.cornell.edu/faculty-research/faculty/lc785/)

Please email all comments/questions to smankad [AT] ncsu.edu

[CRAN Version](https://CRAN.R-project.org/package=sts)

### Summary

This repository will host the development version of the package.  It is also available on CRAN. It estimates the Structural Topic and Sentiment-Discourse Model that was developed in Chen and Mankad (2024) <doi:10.1287/mnsc.2022.00261>. 

The package currently includes functionality to:
* estimate Structural Topic and Sentiment-Discourse Models
* calculate covariate effects on latent topic prevalence and sentiment-discourse
* create the main plots and tables used in our paper

### Installation Instructions
Assuming you already have R installed (if not see http://www.r-project.org/),
to install the CRAN version, simply use:
```
install.packages("sts")
```

### Acknowledgements
We want to thank our research assistant Nala Peng for her help in translating our original R code to C++.

Our work, both the implementation and the theoretical development, builds on the foundational [STM](https://github.com/bstewart/stm) library.




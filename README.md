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

You can install the most recent development version using the devtools package.  First you have 
to install devtools using the following code.  Note that you only have to do this once
```  
if(!require(devtools)) install.packages("devtools")
```   
Then you can load the package and use the function `install_github`

```
library(devtools)
install_github("mankad-research/sts_r",dependencies=TRUE)
```

Note that this will install all the packages suggested and required to run our package.  It may take a few minutes the first time, but this only needs to be done on the first use.  In the future you can update to the most recent development version using the same code. 

### Acknowledgements
We want to thank our research assistant Nala Peng for her help in translating our original R code to C++.

Our work, both the implementation and the theoretical development, builds on the foundational [STM](https://github.com/bstewart/stm) library.




# InfHMM:
Infinit hidden Markov model
The code interface is written in R, and for the sake of speed, most parts are written in C++. However, no prerequisite knowledge for both languages is required to run the code.
An R file called runInfHMM.R sources all needed functions to compile and run the code.

# Example: 
```
# Clean the memory
rm(list = ls())

# For reproducible results
set.seed(1990)

# You working directory to have all needed source files; e.g.:
setwd("~/")

# needed libs
library(Rcpp)
sourceCpp("FFBScplus.cpp")
sourceCpp("eXpandCplus.cpp")
sourceCpp("countCplus.cpp")
source("InfHMM.R")

### To run the model you need to feed to the infHmm function the following parameters: 
# meas: you data set as vector (time series), 
# n.iter: number of iterations, I recommand 25000, 
# n.burn: number of burning iterations---ignored iterations--- I recommand 10000, 
# ns.init: The starting number of states, I recommande 40 for clear sky index data, 
# pri.gam: the prior gamma, I used c(1,1) in the paper, 
# pri.alf: the prior alpha, I used c(1,1) in the paper.

# load your data, preprocess it and store it in vector.
Your_Data <- runif(500) #loaded.and.preprocessed.data

results <- infHmm(meas = Your_Data, n.iter = 25000, n.burn = 10000, ns.init = 40, pri.gam = c(1,1), pri.alf = c(1,1))

```

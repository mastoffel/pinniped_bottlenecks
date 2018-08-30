# Analysis workflow for "Recent demographic histories and genetic diversity across pinnipeds are shaped by anthropogenic interactions and mediated by ecology and life-history", Stoffel et al.

## Overview

To run the analyses, please download the complete folder. The folder named R contains the complete
analysis divided into 14 scripts, which are named 01_ to 14_ . As some of the analyses are 
computationally intensive, many datasets which are produced along the way are already saved
in data/ , so that the analysis can be started at any point. The coalescent simulations and
ABC analyses are outsourced into a second repository (but the results are available in data/).

To sum up, the complete analysis workflow is split into two repositories:

(1) Coalescent simulations and ABC analysis
- repository: https://github.com/mastoffel/Pinniped_bottleneck_CoalSimABC
- the scripts in this repository should run on a server with sufficient computing power and memory

(2) All further analyses are in the current repo.
- repository: https://github.com/mastoffel/pinniped_bottlenecks
- these scripts can be run a standard desktop machine (except for script 1, the STRUCTURE analysis)
- the scripts to reproduce the main results including the figures have been placed in the R folder
and are sequentially named from 01 to 14.

In addition, we packaged some specific functions into two packages:  
(a) sealABC:  devtools::install_github("mastoffel/sealABC")  
(b) mcmcR2:   devtools::install_github("mastoffel/mcmcR2")
  
Both packages are highly specific to the current analysis and probably have to be modified to 
be of use in other projects.



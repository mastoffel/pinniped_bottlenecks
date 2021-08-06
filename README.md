
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Code and data for:

Stoffel, M. A., Humble, E., Paijmans, A. J., Acevedo-Whitehouse, K.,
Chilvers, B. L., Dickerson, B., … & Hoffman, J. I. (2018). Demographic
histories and genetic diversity across pinnipeds are shaped by human
exploitation, ecology and life-history. *Nature Communications*,
**9**(1), 1-12.
[![](https://img.shields.io/badge/doi-https://doi.org/10.1038/s41467--018--06695--z-green.svg)](https://doi.org/https://doi.org/10.1038/s41467-018-06695-z)
[![](https://img.shields.io/badge/Altmetric-150-Darkorange.svg)](https://www.altmetric.com/details/51271319)

<!-- badges: start -->
<!-- badges: end -->

## Project

We analysed genetic data from 11,000 individuals spanning 30 pinniped
species around the globe to shed light on the consequences of large
scale seal-hunting in the 18th and 19th century. Around a third of these
species showed distinct genetic signatures of past population
bottlenecks. Pinniped species with harem-based mating systems and
land-breeding species were particularly affected by population declines
and loss of genetic diversity due to overhunting.

![](./other_stuff/pics_github/elephant_seal_weaners.jpg)

## Reproducing the analysis

To run the analyses, please download the complete folder. The folder
named R contains the complete analysis divided into 14 scripts, which
are named 01\_ to 14\_ . As some of the analyses are computationally
intensive, many datasets which are produced along the way are already
saved in data/ , so that the analysis can be started at any point. The
coalescent simulations and ABC analyses are outsourced into a second
repository (but the results are available in data/).

To sum up, the complete analysis workflow is split into two
repositories:

1.  Coalescent simulations and ABC analysis

-   repository:
    <https://github.com/mastoffel/Pinniped_bottleneck_CoalSimABC>
-   the scripts in this repository should run on a server with
    sufficient computing power and memory

2.  All further analyses are in the current repo.

-   repository: <https://github.com/mastoffel/pinniped_bottlenecks>
-   these scripts can be run a standard desktop machine (except for
    script 1, the STRUCTURE analysis)
-   the scripts to reproduce the main results including the figures have
    been placed in the R folder and are sequentially named from 01
    to 14.

In addition, we packaged some specific functions into two packages:  
(a) sealABC: devtools::install\_github(“mastoffel/sealABC”)  
(b) mcmcR2: devtools::install\_github(“mastoffel/mcmcR2”)

Both packages are highly specific to the current analysis and probably
have to be modified to be of use in other projects.

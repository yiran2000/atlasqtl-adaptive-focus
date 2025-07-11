<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- First time: run usethis::use_readme_rmd() to create a pre-commit hook that 
prevents from committing if the README.Rmd has changed, but has not been 
re-knitted to generate an updated README.md -->

## atlasqtl - variable selection in sparse regression with hierarchically-related responses <img src="man/figures/atlasqtl_logo.png" align="right" height="150"/>

<!-- Run for the R CMD checks, run usethis::use_github_actions() to set up the pipeline, possibly modify the .yaml file and then: -->
<!-- [![R build status](https://github.com/hruffieux/atlasqtl/workflows/R-CMD-check/badge.svg)](https://github.com/hruffieux/atlasqtl/actions) # TODO. not enabled yet, needs pre-install of GSL lib for windows -->
<!-- [![](https://travis-ci.org/hruffieux/atlasqtl.svg?branch=master)](https://travis-ci.org/hruffieux/atlasqtl) -->

[![License: GPL
v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![](https://img.shields.io/badge/devel%20version-0.1.4-blue.svg)](https://github.com/hruffieux/atlasqtl)
[![](https://img.shields.io/github/languages/code-size/hruffieux/atlasqtl.svg)](https://github.com/hruffieux/atlasqtl)
[![](https://img.shields.io/badge/doi-10.1214/20--AOAS1332-blue.svg)](https://doi.org/10.1214/20-AOAS1332)

## Overview

**atlasqtl-adaptive-focus** is the atlasQTL (see "hruffieux/atlasqtl") method with 
a improved adaptive-focus CAVI implementation. This strategy is transferale to 
various sparse Bayesian hiearchical models


## Installation

**Important note:** the R package depends on the GSL library which needs
to be manually installed. For example on Ubuntu,

``` bash
$ sudo apt-get install libgsl-dev
```

or on mac,

``` bash
$ brew install gsl
```

after having installed Homebrew.

Then, to install the package in R, run the following command:

``` r
if (!require(remotes)) install.packages("remotes")
remotes::install_github("yiran2000/atlasqtl-adaptive-focus")
```

## License and authors

This software uses the GPL v3 license, see [LICENSE](LICENSE). Authors
and copyright are provided in [DESCRIPTION](DESCRIPTION).

## Issues

To report an issue, please use the [atlasqtl issue
tracker](https://github.com/hruffieux/atlasqtl/issues) at github.com.

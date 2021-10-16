# Generalized Correlates Analysis Processing

## Summary

This repository houses modular workflows for the processing of correlates of
risk / protection data, exploratory tabular and graphical analyses of correlates
immunogenicity, and the automated reporting of analytic results. It serves as
a generalized suite of tools, based on the analyses originally designed for the
USG Biostatistics Response Team's analysis of COVID-19 vaccine efficacy trials
(archived
[here](https://github.com/CoVPN/correlates_reporting_usgcove_archive/)). See
below for brief descriptions of each of the analysis modules. This repository is
designed as the first part of an analytic pipeline, with the [correlates
reporting](https://github.com/CoVPN/correlates_reporting2) module serving as
a downstream component.

[![Build Status](https://app.travis-ci.com/CoVPN/correlates_processing.svg?branch=master)](https://app.travis-ci.com/CoVPN/correlates_processing)
_Note:_ automated builds of the correlates of risk and protection analyses are
evaluated by the [Travis CI](https://travis-ci.org/) continuous integration
service and the PDF reports posted to this repository's [`gh-pages`
branch](https://github.com/CoVPN/correlates_processing/tree/gh-pages).

## Contents

* Immunogenicity Description
  * `immuno_tabular`: Tabular descriptions of marker immunogenicity.
  * `immuno_graphical`: Graphical descriptions of marker immunogenicity.
* Baseline Risk Score Evaluation
  * `riskscore_baseline`: Super Learner baseline risk score for exposure.

## Collaboration Guide

* [Code style guide](https://style.tidyverse.org/), with some modifications;
  this will largely be enforcd with [`styler`](https://styler.r-lib.org/).
* Project organization: _mostly_ independent subdirectories, each incorporating
  [`here`](https://here.r-lib.org/) for path resolution.
* Package version control and virtual environments using
  [`renv`](https://rstudio.github.io/renv/).
* Code review procedure: see our [contribution
   guidelines](https://github.com/CoVPN/correlates_processing/blob/master/CONTRIBUTING.md).

---

## License

The contents of this repository are distributed under the GPL-3 license. See
file [`LICENSE.md`](https://github.com/CoVPN/correlates_processing/blob/master/LICENSE.md)
for details.

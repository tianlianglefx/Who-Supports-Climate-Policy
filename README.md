# Who Supports Climate Policy?

This repository contains the R code used for the term paper:

**Who Supports Climate Policy? A Model-Based Boosting Analysis of Fossil-Fuel-Tax Support in the European Social Survey**

The analysis uses data from **European Social Survey (ESS) Round 8** and examines individual-level characteristics associated with support for increasing taxes on fossil fuels across 23 countries.

## Analysis

The analysis includes:

- weighted descriptive statistics using the ESS analysis weight;
- weighted logistic regression;
- comparison of core and full model specifications;
- component-wise model-based boosting;
- tuning of the boosting stopping iteration by subsampling;
- stability selection;
- a robustness analysis comparing:
  - support versus explicit opposition; and
  - support versus all other responses;
- diagnostics of structural missingness in the ESS climate module.

The main analysis focuses on climate attitudes, socio-demographic characteristics, political orientation, political trust, and country-level heterogeneity.

## Data

The analysis uses the integrated European Social Survey Round 8,
Edition 2.3 dataset.

The ESS data is too large and is therefore not included in this repository. The dataset can be
downloaded from the European Social Survey Data Portal:

https://www.europeansocialsurvey.org/

After downloading the data, place the following file in the repository
ESS8e02_3.csv

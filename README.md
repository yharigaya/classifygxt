# ClassifyGxT

## About

ClassifyGxT is a method for classifying gene-by-treatment (GxT) interactions using Bayesian model selection (BMS).
The method is primarily designed for molecular count phenotypes, such as gene expression and chromatin accessibility, although it can be used for other types of phenotypes.
It takes as input a list of feature-SNP pairs (or SNPs) with significant GxT interactions that have already been identified by a standard method and assigns posterior probability of different types of GxT interactions to each of the feature-SNP pairs (or SNPs).
See [Overview](https://yharigaya.github.io/classifygxt/#overview) below for the types of GxT interactions.

See `Get started` and `Reference` for usage instructions.

## Publication

Yuriko Harigaya, 
Nana Matoba,
Brandon D. Le,
Jordan M. Valone,
Jason L. Stein,
Michael I. Love\*,
William Valdar\*.
"Probabilistic classification of gene-by-treatment interactions on molecular count phenotypes."
doi: https://doi.org/10.1101/2024.08.03.605142
(\* These authors contributed equally to this work.)

## Overview

In the simplest form, our BMS framework can be described as follows. 
For $n$ samples indexed by $i$, we consider a linear model shown below. 
$y_i$ is the phenotype, $g_i$ is the genotype, and $t_i$ is the
indicator variable as to whether the sample $i$ was treated.
$\varepsilon_i$ is the residual error (noise) and is assumed to follow a Normal distributiion with a zero mean and variance $\sigma^2$.
$\beta_0$ is the intercept.
$\beta_g$, $\beta_t$, and $\beta_{g \times t}$ are the genotype, treatment, and GxT interaction effects, respectively.
We use the $m_g$, $m_t$, and $m_{g \times t}$ indicator variables to specify the inclusion and exclusion of the genotype, treatment, and GxT interaction terms, respectively.
The $\mathbf{m}_j$ vector refers to the $j$-th model ($j = 1, \dots, 8$).
The eight models represent different types of GxT interactions.

<center><img src="man/figures/bms.png" align="center" width="600"></center>

To account for the inherent nonlinear relationship between the genotype and the tranformed molecular count phenotype, we use nonlinear regression.
The eight models can be formulated in the same way using the $\mathbf{m}$ vectors.
See our manuscript ([Publication](https://yharigaya.github.io/classifygxt/#publication)) for prior specification, posterior computation, and other details.

Although we primarily use the $\mathbf{m}$ vector notations to specify the types of GxT interaction, it is also useful to consider the following nomenclatures. 

|$\mathbf{m}$ vector|Nomenclature|
|--|--|
|(0,0,0)|No effect|
|(1,0,0)|Genotype main effect only|
|(0,1,0)|Treatment main effect only|
|(1,1,0)|Genotype and treatment main effects only|
|(0,0,1)|Treatment-induced genotype effect|
|(1,0,1)|Treatment-altered genotype effect|
|(0,1,1)|Treatment main effect and treatment-induced genotype effect|
|(1,1,1)|Genotype and treatment main effects and interaction|

## Data input

ClassifyGxT requires individual genotype and phenotype data for a list of feature-SNP pairs (or SNPs) with significant GxT interactions that have already been identified by a standard method, such as [TensorQTL](https://github.com/broadinstitute/tensorqtl). See `Get started` for the format of input data. See [Using ClassifyGxT with TensorQTL](https://yharigaya.github.io/classifygxt/articles/tensorqtl.html) under `Articles` for an instruction on formatting input data for ClassifyGxT from typical input/output data for [TensorQTL](https://github.com/broadinstitute/tensorqtl). 

## Installation

The *classifygxt* R package can be installed using *devtools*.

```
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("yharigaya/classifygxt")
```

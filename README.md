# sensitivityq
This package provides the implementation of the proposed method in 'Sensitivity Analysis for Quantiles of Hidden Biases in Matched Observational Studies',
which extends Rosenbaum’s sensitivity analysis to infer all quantiles of the biases across all sets, extending previous works that essentially focusing on the maximum bias.

## To install the package
First, you need to install the devtools package. You can do this from CRAN. Invoke R and then type

`install.packages("devtools")`

Load the devtools package.

`library(devtools)`

Then do

`install_github("dongxiao-wu/sensitivityq")`

To import the package, simply do `library('sensitivityq')`.

## Two main functions in this package

The two main functions in this package are:

`senmk()`: calculate p-value, statistics, expectation and variance for the extended sensitivity analysis of a specific quantile of biases and a specific sensitivity parameter Gamma.

`Gamma_seq()`: calculates lower confidence limits for the biases at rank k's(for a sequance of k as needed) across all matched sets.

For more informations, see the help pages in the package.

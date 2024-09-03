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

### `senmk()`
Calculate p-value, statistics, expectation and variance for the extended sensitivity analysis of a specific quantile of biases and a specific sensitivity parameter Gamma.

Usage: 

`senmk(y, z, mset, k, gamma = 1, inner = 0, trim = 3, lambda = 1/2, tau = 0, alternative = "greater", TonT = FALSE, precise = FALSE)`

Example:

```
I <- 500

mset <- as.vector(rbind(1:I,1:I))

z <- as.vector(rbind(rep(1,I),rep(0,I)))

y <- rnorm(1000,sd=sqrt(0.5))+0.5*z

senmk(y, z, mset, k = I)
```


### `Gamma_seq()`
Calculates lower confidence limits for the biases at rank k's(for a sequance of k as needed) across all matched sets.

Usage: 

`Gamma_seq(y, z, mset, inner = 0, trim = 3, thres, Ks, tol = 1e-04, precise = FALSE)`

Example:

```
I <- 500

mset <- as.vector(rbind(1:I,1:I))

z <- as.vector(rbind(rep(1,I),rep(0,I)))

y <- rnorm(1000,sd=sqrt(0.5))+0.5*z

Gamma_seq(y,z,mset,Ks=500:450,thres=0.05)
```

For more detailed descriptions of the arguments and output values, see the help pages in the package.

#' Sensitivity Analysis for a Matched Comparison in an Observational Study.
#'
#' Each matched set contains one treated individual and one or more controls. It is based on the function senm() in the sensitivitymult package. Uses Huber's M-statistic as the basis for the test with the loosen precondition that only parts of the biases are bounded. Performs either a randomization test or an analysis of sensitivity to departures from random assignment.
#'
#' @param y A vector of responses with no missing data.
#' @param z Treatment indicator, z=1 for treated, z=0 for control with length(z)==length(y).
#' @param mset Matched set indicator, 1, 2, ..., sum(z) with length(mset)==length(y). Matched set indicators should be either integers or a factor.
#' @param k k smallest biases are bounded by gamma
#' @param gamma gamma is the sensitivity parameter Γ, where Γ ≥ 1. Setting Γ = 1 is equivalent to assuming ignorable treatment assignment given the matched sets, and it performs a within-set randomization test.
#' @param inner inner and trim together define the ψ-function for the M-statistic. The default values yield a version of Huber's ψ-function, while setting inner = 0 and trim = Inf uses the mean within each matched set. The ψ-function is an odd function, so ψ(w) = -ψ(-w). For w ≥ 0, the ψ-function is ψ(w)=0 for 0 ≤ w ≤ inner, is ψ(w)= trim for w ≥ trim, and rises linearly from 0 to trim for inner < w < trim.
#' An error will result unless 0 ≤ inner ≤ trim.
#' Taking trim < Inf limits the influence of outliers; see Huber (1981). Taking trim < Inf and inner = 0 uses Huber's psi function. Taking trim = Inf does no trimming and is similar to a weighted mean; see TonT. Taking inner > 0 often increases design sensitivity; see Rosenbaum (2013).
#' @param trim inner and trim together define the ψ-function for the M-statistic. See inner.
#' @param lambda Before applying the ψ-function to treated-minus-control differences, the differences are scaled by dividing by the lambda quantile of all within set absolute differences. Typically, lambda = 1/2 for the median. The value of lambda has no effect if trim=Inf and inner=0. See Maritz (1979) for the paired case and Rosenbaum (2007) for matched sets.
#' An error will result unless 0 < lambda < 1.
#' @param tau The null hypothesis asserts that the treatment has an additive effect, tau. By default, tau=0, so by default the null hypothesis is Fisher's sharp null hypothesis of no treatment effect.
#' @param alternative If alternative="greater", the null hypothesis of a treatment effect of tau is tested against the alternative of a treatment effect larger than tau. If alternative="less", the null hypothesis of a treatment effect of tau is tested against the alternative of a treatment effect smaller than tau. In particular, alternative="less" is equivalent to: (i) alternative="greater", (ii) y replaced by -y, and (iii) tau replaced by -tau. See the note for discussion of two-sided sensitivity analyses.
#' @param TonT TonT refers to the effect of the treatment on the treated; see Rosenbaum and Rubin (1985, equation 1.1.1) The default is TonT=FALSE. If TonT=FALSE, then the total score in matched set i is divided by the number ni of individuals in set i, as in expression (8) in Rosenbaum (2007). This division by ni has few consequences when every matched set has the same number of individuals, but when set sizes vary, dividing by ni is intended to increase efficiency by weighting inversely as the variance; see the discussion in section 4.2 of Rosenbaum (2007). If TonT=TRUE, then the division is by ni-1, not by ni, and there is a further division by the total number of matched sets to make it a type of mean. If TonT=TRUE and trim=Inf, then the statistic is the mean over matched sets of the treated minus mean-control response, so it is weighted to estimate the average effect of the treatment on the treated. See the examples.
#' @param precise use precise method or not in the case of matched pairs
#'
#' @return A list.
#' \itemize{
#'   \item pval - Approximate upper bound on the one-sided P-value.
#'   \item deviate - Deviate that is compared to the upper tail of the standard Normal distribution to obtain the P-value.
#'   \item statistic - Value of the test statistic.
#'   \item expectation - Maximum null expectation of the test statistic for the given value of gamma.
#'   \item variance	- Among null distributions that yield the maximum expectation, variance is the maximum possible variance for the given value of gamma.
#' }
#'
#' @examples
#' I <- 500
#' mset <- as.vector(rbind(1:I,1:I))
#' z <- as.vector(rbind(rep(1,I),rep(0,I)))
#' y <- rnorm(1000,sd=sqrt(0.5))+0.5*z
#' senmk(y, z, mset, k = I)
#'
#' @export
senmk <- function (y, z, mset, k, gamma = 1, inner = 0, trim = 3, lambda = 1/2,
                   tau = 0, alternative = "greater", TonT = FALSE, precise=FALSE)
{
  # This is based on the senm function from sensitivitymult
  stopifnot((alternative == "greater") | (alternative == "less"))
  stopifnot(gamma >= 1)
  stopifnot((inner >= 0) & (inner <= trim))
  stopifnot((lambda > 0) & (lambda < 1))
  stopifnot(is.vector(y) & is.vector(z) & is.vector(mset))
  stopifnot((length(z) == length(y)))
  stopifnot((length(z) == length(mset)))
  stopifnot(all(!is.na(y)))
  stopifnot(all((z == 0) | (z == 1)))
  tbcheck <- table(z, mset)
  ck <- all(tbcheck[2, ] == 1) & all(tbcheck[1, ] >= 1)
  if (!ck) {
    warning("Every matched set must contain one treated subject and at least one control.")
    stopifnot(ck)
  }
  mset <- as.integer(mset)
  o <- order(mset, 1 - z)
  y <- y[o]
  z <- z[o]
  mset <- mset[o]
  tb <- table(mset)
  nset <- length(tb)
  setsize <- max(tb)
  makeymat <- function(yj) {
    ymat <- matrix(NA, nset, setsize)
    m <- 0
    for (i in 1:nset) {
      ymat[i, 1:tb[i]] <- yj[(m + 1):(m + tb[i])]
      m <- m + tb[i]
    }
    ymat
  }
  ymat <- makeymat(y)
  if (alternative == "less") {
    ymat <- (-ymat)
    tau <- (-tau)
  }
  if (!(tau == 0))
    ymat[, 1] <- ymat[, 1] - tau
  ms <- mscorev(ymat, inner = inner, trim = trim, qu = lambda, TonT = TonT)
  separable1vk(ms, gamma = gamma, k, precise)
}

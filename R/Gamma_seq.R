#' The sequence of sensitivity values under the loosen precondition that only parts of the biases are bounded.
#'
#' This function calculates lower confidence limits for the biases at rank k(k=I,I-1,...,I-K+1) across all matched sets.
#'
#' @param y A vector of responses with no missing data.
#' @param z Treatment indicator, z=1 for treated, z=0 for control with length(z)==length(y).
#' @param mset Matched set indicator, 1, 2, ..., sum(z) with length(mset)==length(y). Matched set indicators should be either integers or a factor.
#' @param inner inner and trim together define the ψ-function for the M-statistic. The default values yield a version of Huber's ψ-function, while setting inner = 0 and trim = Inf uses the mean within each matched set. The ψ-function is an odd function, so ψ(w) = -ψ(-w). For w ≥ 0, the ψ-function is ψ(w)=0 for 0 ≤ w ≤ inner, is ψ(w)= trim for w ≥ trim, and rises linearly from 0 to trim for inner < w < trim.
#' An error will result unless 0 ≤ inner ≤ trim.
#' Taking trim < Inf limits the influence of outliers; see Huber (1981). Taking trim < Inf and inner = 0 uses Huber's psi function. Taking trim = Inf does no trimming and is similar to a weighted mean; see TonT. Taking inner > 0 often increases design sensitivity; see Rosenbaum (2013).
#' @param trim inner and trim together define the ψ-function for the M-statistic. See inner.
#' @param thres the significance level for the test
#' @param K the needed length of the sequence
#' @param tol the desired accuracy (convergence tolerance).
#' @param precise use precise method or not in the case of matched pairs
#'
#' @return A list of the lower confidence limits for the biases at rank k(k=I, I-1, ..., I-K+1)
#'
#' @examples
#' I <- 500
#' mset <- as.vector(rbind(1:I,1:I))
#' z <- as.vector(rbind(rep(1,I),rep(0,I)))
#' y <- rnorm(1000,sd=sqrt(0.5))+0.5*z
#' Gamma_seq(y,z,mset,K=200,thres=0.05)
#'
#' @export
Gamma_seq <- function(y, z, mset, inner = 0, trim = 3, thres, K, tol=0.0001, precise=FALSE, alternative = "greater")
{
  # Get Gamma_{[k]} (k=I, I-1, ..., I-K+1)
  I <- sum(z)
  Gamma <- rep(1,K)
  for (i in 1:K){
    if ((i > 1) && (Gamma[i-1]==1)){Gamma[i]=1}
    else{
      f <- function(gamma){senmk(y, z, mset, k=I-i+1, gamma=gamma, inner=inner, trim=trim, precise=precise, alternative = "greater")$pval-thres}
      if (f(1) < 0){
        gamma_sol <- uniroot(f, lower = 1+tol, upper = 100000, extendInt = "no",tol=tol)$root
        if (f(gamma_sol) > 0){
          while( f(gamma_sol) > 0 ){
            gamma_sol <- gamma_sol - tol
          }
          gamma_sol <- gamma_sol + tol
        }
        else{
          while( f(gamma_sol) <= 0 ){
            gamma_sol = gamma_sol + tol
          }
        }
      }
      else{gamma_sol = 1}
      Gamma[i] <- gamma_sol
    }
  }
  Gamma
}

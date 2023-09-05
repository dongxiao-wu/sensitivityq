#' Computes M-scores for M-tests and estimates.
#'
#' Of limited interest to most users, function mscorev() computes the M-scores used by functions senm(), senmCI(), comparison(), and principal() that perform Huber-Maritz M-tests. The function is also in the package sensitivitymv.
#'
#' @param ymat ymat is a matrix as described in the documentation for senmk().
#' @param inner inner is the parameter described in the documentation for senmk().
#' The inner parameter is discussed in Rosenbaum (2013).
#' @param trim trim is the parameter described in the documentation for senmk().
#' Note that the default for mscorev is trim = 2.5, but other functions in this package reset the default to trim = 3.
#' This is for consistency with the sensitivitymv package which also contains the mscorev function.
#' @param qu qu is the lambda parameter described in the documentation for senmk().
#' @param TonT If TonT=FALSE, then the total score in set (row) i is divided by the number ni of individuals in row i, as in expression (8) in Rosenbaum (2007).
#' If TonT=TRUE, then the division is by ni-1, not by ni, and there is a further division by the total number of matched sets.
#' See the discussion of TonT in the documentation for senmk().
#'
#' @return Generally, a matrix with the same dimensions as ymat containing the M-scores.
#' Exception: if a matched set does not contain at least one treated subject and at least one control, then that set will not appear in the result, and the result will have fewer rows than ymat.
#' However, if a matched set has several controls but no treated subject, then these controls will contribute to the estimate of the scale parameter, typically the median absolute pair difference.
#'
#' @export
mscorev <-function (ymat, inner = 0, trim = 2.5, qu = 0.5, TonT = FALSE)
  {
    # This is the mscorev function from sensitivitymv version 1.3
    ymat <- as.matrix(ymat)
    n <- dim(ymat)[1]
    m <- dim(ymat)[2]
    out <- matrix(NA, n, m)
    one <- rep(1, m - 1)
    difs <- array(NA, c(n, m, m - 1))
    for (j in 1:m) {
      difs[, j, ] <- outer(as.vector(unlist(ymat[, j])), one,
                           "*") - ymat[, -j]
    }
    ms <- as.vector(difs)
    if ((trim < Inf) | (inner > 0)) {
      hqu <- as.numeric(quantile(abs(ms), qu, na.rm = TRUE))
      if (hqu > 0) {
        ms <- ms/hqu
        if ((trim < Inf) & (inner < trim)) {
          ab <- pmin(1, pmax(0, (abs(ms) - inner))/(trim -
                                                      inner))
        }
        else if ((trim < Inf) & (inner == trim)) {
          ab <- 1 * (abs(ms) > inner)
        }
        else {
          ab <- pmax(0, abs(ms) - inner)
        }
        ms <- sign(ms) * ab
      }
      else {
        warning("Error: Scale factor is zero.  Increase lambda.")
      }
    }
    ms <- array(ms, c(n, m, m - 1))
    ms <- apply(ms, c(1, 2), sum, na.rm = TRUE)
    ms[is.na(ymat)] <- NA
    colnames(ms) <- colnames(ymat)
    ni <- apply(!is.na(ymat), 1, sum)
    use <- (ni >= 2) & (!is.na(ms[, 1]))
    ms <- ms[use, ]
    ni <- ni[use]
    if (TonT) {
      ms <- (ms/outer(ni - 1, rep(1, m), "*"))/(dim(ms)[1])
    }
    else {
      ms <- ms/outer(ni, rep(1, m), "*")
    }
    ms
  }


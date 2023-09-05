#' Asymptotic separable calculations internal to other functions.
#'
#' Of limited interest to most users, this general purpose function is internal to other functions in the sensitivityq package. It is based on separable1v function in the sensitivitymv package, version 1.3. The function performs the asymptotic separable calculations with the loosen precondition that only parts of the biases are bounded.
#'
#' @param ymat ymat is a matrix whose rows are matched sets and whose columns are matched individuals. The first column describes treated individuals. Other columns describe controls. If matched sets contain variable numbers of controls, NAs fill in empty spaces in ymat; see the documentation for senmv. In senmk, the matrix ymat is created by mscorev. Instead, if there were no NAs and ranks within rows were used in ymat, then separable1v would perform a sensitivity analysis for the stratified Wilcoxon two-sample test. Applied directly to data, it performs a sensitivity analysis for the permutational t-test.
#' @param gamma gamma is the value of the sensitivity parameter; see the documentation for the senmk function. One should use a value of gamma >= 1.
#' @param k k smallest biases are bounded by gamma
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
#' @export
separable1vk <- function (ymat, gamma = 1, k, precise=FALSE)
{
  I <- dim(ymat)[1]
  n <- dim(ymat)[2]
  mu_inf <- rep(NA,I)
  mu_gamma <- rep(NA,I)
  var_inf <- rep(NA,I)
  var_gamma <- rep(NA,I)
  tstat <- sum(ymat[,1])
  if (precise==FALSE){
    for (i in 1:I) {
      if (!is.na(ymat[i, 1])) {
        rk <- sort(as.vector(unlist(ymat[i, ])))
        rk <- sort(rk[!is.na(rk)])
        ni <- length(rk)
        mu_gamma[i] <- (-Inf)
        var_gamma[i] <- (-Inf)
        mu_inf[i] <- rk[ni]
        var_inf[i] <- 0
        for (ai in 1:(ni - 1)) {
          if (gamma == Inf){
            mu_ia_gamma <- sum(rk[(ai + 1):ni])/(ni - ai)
            var_ia_gamma <- sum(rk[(ai + 1):ni]^2)/(ni - ai)-(mu_ia_gamma^2)}
          else{
            mu_ia_gamma <- (sum(rk[1:ai]) + gamma * sum(rk[(ai + 1):ni])) / (ai + gamma * (ni - ai))
            var_ia_gamma <- max(((sum(rk[1:ai]^2) + gamma * sum(rk[(ai + 1):ni]^2)) / (ai + gamma * (ni - ai))) - (mu_ia_gamma^2),0)
          }
          if (mu_ia_gamma > mu_gamma[i]) {
            mu_gamma[i] <- mu_ia_gamma
            var_gamma[i] <- var_ia_gamma
          }
          else if (mu_ia_gamma == mu_gamma[i])
            var_gamma[i] <- max(var_gamma[i], var_ia_gamma)
        }
      }
    }
    mu_diff <- mu_inf - mu_gamma
    v_sort <- sort(var_gamma,decreasing = T,index.return = T)
    diff_resort <- mu_diff[v_sort$ix]
    diff_sort_res <- sort(diff_resort,index.return=T)
    gamma_index <- v_sort$ix[diff_sort_res$ix[1:k]]
    mu_gamma_k <- sum(mu_gamma[gamma_index])+sum(mu_inf[-gamma_index])
    var_gamma_k <- sum(var_gamma[gamma_index])+sum(var_inf[-gamma_index])
    tstat <- as.vector(tstat)
    dev <- (tstat - mu_gamma_k)/sqrt(var_gamma_k)
    pval <- 1 - stats::pnorm(dev)
    if (var_gamma_k==0){pval <- 0}
    result <- list(pval = pval, deviate = dev, statistic = tstat, expectation = mu_gamma_k,
                   variance = var_gamma_k)
  }
  if (precise==TRUE){
    q_diff <- rep(0,I)
    Q <- matrix(nrow=2,ncol=I)
    for (i in 1:I) {
      rk <- sort(as.vector(unlist(ymat[i, ])))
      rk <- sort(rk[!is.na(rk)])
      Q[,i] <- rk
      q_diff[i] <- abs(rk[2]-rk[1])
    }
    Ik <- which(rank(q_diff,ties.method="first") <= k)
    monte <- matrix(nrow=I,ncol=10000)
    for (i in 1:I){
      if (i %in% Ik){
        monte[i,] <- sample(Q[,i],size=10000,replace=T,prob=c(1/(1+gamma),gamma/(1+gamma)))
      }
      else{
        monte[i,] <- rep(Q[2,i],10000)
      }
    }
    pval <- mean(colSums(monte)>=tstat)
    result <- list(pval = pval, statistic = tstat)
  }
  result
}

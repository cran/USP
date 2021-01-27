#' Test statistic for dependence in contingency table
#'
#' This function computes the value of the test statistic \eqn{T_n} measuring the strength of
#' dependence in a contingency table. See Section 3.1 of \insertCite{BKS2020}{USP}
#' for a definition.
#'
#' @param freq Two-way contingency table whose strength of dependence is to be measured.
#'
#' @return A list containing the value of the test statistic \eqn{T_n}, the table of expected
#' null counts, and the table of contributions to \eqn{T_n}.
#' @export
#'
#' @references \insertRef{BKS2020}{USP}
#'
#' @examples
#' freq=r2dtable(1,rep(10,5),rep(10,5))[[1]]; DiscStat(freq)
#'
#' freq=diag(1:5); DiscStat(freq)
#'
#' freq=r2dtable(1,rep(10,5),rep(10,5))[[1]] + 4*diag(rep(1,5))
#' DiscStat(freq)
DiscStat=function(freq){
  r=dim(freq)[1]; c=dim(freq)[2]; n=sum(freq)
  freqc=colSums(freq); freqr=rowSums(freq)
  exp=freqr%*%t(freqc)/n
  contributions=(freq-exp)^2/(n*(n-3))-4*freq*exp/(n*(n-2)*(n-3))
  return(list(sum(contributions),exp,contributions))
}

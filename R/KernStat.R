#' Test statistic calculated from two kernel matrices
#'
#' Calculate the U-statistic measure of dependence given two kernel matrices \eqn{J}
#' and \eqn{K}, as described in Section 7.1 of \insertCite{BKS2020}{USP}. For the featured
#' examples considered these matrices can be calculated using FourierKernel or
#' InfKern. Alternatively, if a different basis is to be used, then the kernels
#' can be entered separately.
#'
#' @param J \eqn{n \times n} kernel matrix corresponding to first sample.
#' @param K \eqn{n \times n} kernel matrix corresponding to second sample.
#'
#' @return Test statistic measure the strength of dependence between the two samples.
#' @export
#'
#' @references \insertRef{BKS2020}{USP}
#'
#' @examples
#' x=runif(100); y=runif(100); M=3
#' J=FourierKernel(x,M); K=FourierKernel(y,M)
#' KernStat(J,K)
KernStat=function(J,K){
  n=nrow(J)
  term1=sum(J*K)/(n*(n-3))

  J1=J%*%rep(1,n); K1=K%*%rep(1,n)
  term2=sum(J1*K1)/(n*(n-2)*(n-3))

  term3=sum(J)*sum(K)/(n*(n-1)*(n-2)*(n-3))
  return(term1-2*term2+term3)
}

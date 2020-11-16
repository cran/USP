#' Kernel matrix for Fourier basis
#'
#' Calculates the kernel matrix, described in \insertCite{BKS2020}{USP} for
#' univariate continuous data when using the Fourier basis.
#' This function is used in USPFourier.
#'
#' @param x A vector in \eqn{[0,1]^n} for some \eqn{n}, containing the observations.
#' @param M The maximum frequency of Fourier basis functions to compute.
#'
#' @return The kernel matrix \eqn{K}, to be used in independence testing.
#' @export
#'
#' @references \insertRef{BKS2020}{USP}
#'
#' @examples
#' n=10; x=runif(n)
#' FourierKernel(x,5)
FourierKernel=function(x,M){
  n=length(x)
  K=matrix(rep(0,n^2),ncol=n)
  for(m in 1:M){
    K=K+FourierBasis(0,m,x)%*%t(FourierBasis(0,m,x))
    K=K+FourierBasis(1,m,x)%*%t(FourierBasis(1,m,x))
  }
  diag(K)=0
  return(K)
}

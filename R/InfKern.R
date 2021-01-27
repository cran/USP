#' Kernel for infinite-dimensional example
#'
#' Function to produce the kernel matrices in the infinite dimensional example described
#' in Section 7.4 of \insertCite{BKS2020}{USP}. Here, a random function is converted to a sequence
#' of coefficients and we use the Fourier basis on these coefficients. This function is
#' an essential part of USPFunctional.
#'
#' @param X Matrix giving one of the samples to be tested. Each row corresponds to a
#' discretised function, with each column giving the values of the functions at the
#' corresponding grid point.
#' @param Ntrunc The total number of coefficients to look at in the basis expansion
#' of the functional data.
#' @param M The maximum frequency to look at in the Fourier basis.
#'
#' @return The kernel matrix for the sample \eqn{X}.
#' @export
#'
#' @references \insertRef{BKS2020}{USP}
#'
#' @examples
#' n=10  #number of observations
#' Ndisc=1000; t=1/Ndisc #functions represented at grid points 1/Ndisc, 2/Ndisc,...,1
#' X=matrix(rep(0,Ndisc*n),nrow=n)
#' for(i in 1:n){
#'  x=rnorm(Ndisc,mean=0,sd=1)
#'  X[i,]=cumsum(x*sqrt(t))
#' }
#' InfKern(X,2,2)
InfKern=function(X,Ntrunc,M){
  U=coeffs(X,Ntrunc)
  n=nrow(X)
  K=matrix(rep(0,n^2),ncol=n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      K[i,j]=sumbasis(Ntrunc,M,U[i,],U[j,])
      K[j,i]=sumbasis(Ntrunc,M,U[i,],U[j,])
    }
  }
  diag(K)=0
  return(K)
}

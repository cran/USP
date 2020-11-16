#' Kernel entries in infinite dimensional case
#'
#' Function to calculate each entry of the kernel matrix in the infinite dimensional
#' example described in Section 7.4 of \insertCite{BKS2020}{USP}. Here, a random function is converted
#' to a sequence of coefficients and we use the Fourier basis on these coefficients. This
#' function is only used in the function InfKern.
#'
#' @param Ntrunc The total number of coefficients to look at.
#' @param M The maximum frequency to look at in the Fourier basis.
#' @param x1 The coefficients of the first data point.
#' @param x2 The coefficients of the second data point.
#'
#' @return The entry of the kernel corresponding to the two data points.
#' @export
#'
#' @references \insertRef{BKS2020}{USP}
#'
#' @examples
#' x1=runif(5); x2=runif(5); sumbasis(5,2,x1,x2)
sumbasis=function(Ntrunc,M,x1,x2){
  cossum=rep(0,Ntrunc)
  for(j in 1:Ntrunc){
    x1j=x1[j]; x2j=x2[j]
    cossum[j]=sum(FourierBasis(0,1:M,x1j)*FourierBasis(0,1:M,x2j))
  }
  return(sum(cossum))
}

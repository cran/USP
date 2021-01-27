#' Fourier basis functions
#'
#' Computes the values of the one-dimensional Fourier basis functions at a vector of locations \eqn{x}
#' and with a vector of frequencies \eqn{m}. The scaling factor of \eqn{2\pi} is included, so that
#' the function returns, e.g., \eqn{\sqrt{2} \cos(2\pi m x)}.
#'
#' @param a Sine or cosine; \eqn{a=0} gives cosine and \eqn{a=1} gives sine.
#' @param m Vector of frequencies \eqn{m}.
#' @param x Vector of locations \eqn{x}.
#'
#' @return Returns the values of \eqn{\sqrt{2} \cos(2\pi m x)}.
#' @export
#'
#' @references \insertRef{BKS2020}{USP}
#' @examples
#' e=FourierBasis(1,1:100,0.01); plot(0.01*(1:100),e,type="l")
#' e=FourierBasis(0,1,0.01*(1:100)); plot(0.01*(1:100),e,type="l")
#' FourierBasis(1,1:3,0.1*(1:10))
#'
FourierBasis=function(a,m,x){
  matrixinput=matrix(x)%*%t(matrix(m))
  return(sqrt(2)*cos(2*pi*matrixinput-a*pi/2))
}

#' Calculate coefficients of a function's series expansion
#'
#' This function is used in InfKern to produce the kernel matrix from functional
#' data defined on the interval \eqn{[0,1]}. For further details see Section 7.4
#' of \insertCite{BKS2020}{USP}.
#'
#' @param X The discretised functions whose coefficients are required.
#' This should be a matrix with one row per function, and with \eqn{Ndisc} columns,
#' where \eqn{Ndisc} is the grid size of the discretisation.
#' @param Ntrunc The number of coefficients that are required.
#' The function returns coefficients 1,...,\eqn{Ntrunc}.
#'
#' @return The coefficients of \eqn{X} in its expansion in terms of sine functions.
#' See \insertCite{BKS2020}{USP} for more detail.
#' @export
#'
#' @importFrom stats pnorm
#'
#' @references \insertRef{BKS2020}{USP}
#'
#' @examples
#' t=seq(from=0,to=1,length.out=1000); X=t^2
#' U=coeffs(X,100)[1,]; L=5
#' plot(t,X,type="l")
#' approx=rep(0,1000)
#' for(l in 1:L){
#' approx=approx+qnorm(U[l])*sqrt(2)*sin((l-1/2)*pi*t)/((l-1/2)*pi)
#' lines(t,approx,col=l+1)
#' }
coeffs=function(X,Ntrunc){
  if(is.vector(X)) X=t(as.matrix(X))
  n=dim(X)[1]; Ndisc=dim(X)[2]; t=1/Ndisc
  x=t*(1:Ndisc)
  SineMat=sin(pi*x%*%t(1:Ntrunc-1/2))
  RescaledSine=sqrt(2)*pi*t(t(SineMat)*(1:Ntrunc-1/2))
  U=pnorm(t*X%*%RescaledSine)
  return(U)
}

#' Independence test for functional data
#'
#' We implement the permutation independence test described in \insertCite{BKS2020}{USP}
#' for functional data taking values in \eqn{L^2([0,1])}. The discretised functions are
#' expressed in a series expansion, and an independence test is carried out between the
#' coefficients of the functions, using a Fourier basis to define the test statistic.
#'
#' @param X A matrix of the discretised functional data from the first sample. There are \eqn{n} rows,
#' where \eqn{n} is the sample size, and Ndisc columns, where Ndisc is the grid size such
#' that the values of each function on 1/Ndisc, 2/Ndisc, ..., 1 are given.
#' @param Y A matrix of the discretised functional data from the second sample. The discretisation
#' grid may be different to the grid used for \eqn{X}, if required.
#' @param Ntrunc The number of coefficients to retain from the series expansions of \eqn{X} and \eqn{Y}.
#' @param M The maximum frequency to use in the Fourier basis when testing the independence of the
#' coefficients.
#' @param B The number of permutations used to calibrate the test.
#' @param ties.method If "standard" then calculate the p-value as in (5) of \insertCite{BKS2020}{USP},
#' which is slightly conservative. If "random" then break ties randomly. This preserves Type I error
#' control.
#'
#' @return A p-value for the test of the independence of \eqn{X} and \eqn{Y}.
#' @export
#'
#' @references \insertRef{BKS2020}{USP}
#'
#' @examples
#' n=50; r=0.6; Ndisc=1000; t=1/Ndisc
#' X=matrix(rep(0,Ndisc*n),nrow=n); Y=matrix(rep(0,Ndisc*n),nrow=n)
#' for(i in 1:n){
#'  x = rnorm(Ndisc, mean=0, sd= 1)
#'  se = sqrt(1 - r^2) #standard deviation of error
#'  e = rnorm(Ndisc, mean=0, sd=se)
#'  y = r*x + e
#'  X[i,] <- cumsum(x*sqrt(t))
#'  Y[i,] <- cumsum(y*sqrt(t))
#' }
#' USPFunctional(X,Y,2,1,999)
USPFunctional=function(X,Y,Ntrunc,M,B=999,ties.method="standard"){
  J=InfKern(X,Ntrunc,M); K=InfKern(Y,Ntrunc,M)
  return(USP(J,K,B,ties.method))
}

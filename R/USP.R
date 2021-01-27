#' Permutation test of independence.
#'
#' Carry out an independence test of the independence of two samples, give
#' two kernel matrices \eqn{J} and \eqn{K}, as described in Section 7.1 of \insertCite{BKS2020}{USP}.
#' We calculate the test statistic and null statistics using the function KernStat,
#' before comparing them to produce a p-value. For the featured
#' examples considered these matrices can be calculated using FourierKernel or
#' InfKern. Alternatively, if a different basis is to be used, then the kernels
#' can be entered separately.
#'
#' @param J \eqn{n \times n } kernel matrix corresponding to first sample.
#' @param K \eqn{n \times n } kernel matrix corresponding to second sample.
#' @param B The number of permutation used to calibrate the test.
#' @param nullstats If TRUE, returns a vector of the null statistic values.
#' @param ties.method If "standard" then calculate the p-value as in (5) of \insertCite{BKS2020}{USP},
#' which is slightly conservative. If "random" then break ties randomly. This preserves Type I error
#' control.
#'
#' @return Returns the p-value for this independence test and the value of the test statistic, \eqn{D_n},
#' as defined in \insertCite{BKS2020}{USP}. If nullstats=TRUE is used, then the function also
#' returns a vector of the null statistics.
#' @export
#'
#' @references \insertRef{BKS2020}{USP}
#'
#' @examples
#' x=runif(100); y=runif(100); M=3
#' J=FourierKernel(x,M); K=FourierKernel(y,M)
#' USP(J,K,999)
#'
#' n=50; r=0.6; Ndisc=1000; t=1/Ndisc
#' X=matrix(rep(0,Ndisc*n),nrow=n); Y=matrix(rep(0,Ndisc*n),nrow=n)
#' for(i in 1:n){
#'  x = rnorm(Ndisc, mean=0, sd= 1)
#'  se = sqrt(1 - r^2) #standard deviation of error
#'  e = rnorm(Ndisc, mean=0, sd=se)
#'  y = r*x + e
#'  X[i,] = cumsum(x*sqrt(t))
#'  Y[i,] = cumsum(y*sqrt(t))
#' }
#' J=InfKern(X,2,1); K=InfKern(Y,2,1)
#' USP(J,K,999)
USP=function(J,K,B=999,ties.method="standard",nullstats=FALSE){
  n=nrow(J)
  TestStat=KernStat(J,K)
  NullStats=rep(0,B)
  for(b in 1:B){
    perm=sample(n)
    NullStats[b]=KernStat(J,K[perm,perm])
  }
  if(ties.method=="random"){
    pval=(B+2-rank(c(TestStat,NullStats),ties.method = "random")[1])/(B+1)
  }else{
    pval=(1+sum(TestStat<=NullStats))/(B+1)
  }

  if(nullstats==FALSE){
    ans=list(pval,TestStat)
    names(ans)=c("p.value","TestStat")
  }else{
    ans=list(pval,TestStat,NullStats)
    names(ans)=c("p.value","TestStat","NullStats")
  }
  return(ans)
}

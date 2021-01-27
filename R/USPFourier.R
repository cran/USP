#' Independence test for continuous data
#'
#' Performs a permutation test of independence between two univariate continuous random
#' variables, using the Fourier basis to construct the test statistic, as described in
#' \insertCite{BKS2020}{USP}.
#'
#' @param x A vector containing the first sample, with each entry in \eqn{[0,1]}.
#' @param y A vector containing the second sample, with each entry in \eqn{[0,1]}.
#' @param M The maximum frequency to use in the Fourier basis.
#' @param B The number of permutation to use when calibrating the test.
#' @param nullstats If TRUE, returns a vector of the null statistic values.
#' @param ties.method If "standard" then calculate the p-value as in (5) of \insertCite{BKS2020}{USP},
#' which is slightly conservative. If "random" then break ties randomly. This preserves Type I error
#' control.
#'
#' @return Returns the p-value for this independence test and the value of the test statistic, \eqn{D_n},
#' as defined in \insertCite{BKS2020}{USP}. If nullstats=TRUE is used, then the function also
#' returns a vector of the null statistics.
#' @export
#' @references \insertRef{BKS2020}{USP}
#'
#' @examples
#' x=runif(10); y=x^2
#' USPFourier(x,y,1,999)
#'
#' n=100; w=2; x=integer(n); y=integer(n); m=300
#' unifdata=matrix(runif(2*m,min=0,max=1),ncol=2); x1=unifdata[,1]; y1=unifdata[,2]
#' unif=runif(m); prob=0.5*(1+sin(2*pi*w*x1)*sin(2*pi*w*y1)); accept=(unif<prob);
#' Data1=unifdata[accept,]; x=Data1[1:n,1]; y=Data1[1:n,2]
#' plot(x,y)
#' USPFourier(x,y,2,999)
#'
#' x=runif(100); y=runif(100)
#' test=USPFourier(x,y,3,999,nullstats=TRUE)
#' plot(density(test$NullStats,from=min(test$NullStats),to=max(max(test$NullStats),test$TestStat)),
#'          xlim=c(min(test$NullStats),max(max(test$NullStats),test$TestStat)),main="Test Statistics")
#' abline(v=test$TestStat,col=2); TestStats=c(test$TestStat,test$NullStats)
#' abline(v=quantile(TestStats,probs=0.95),lty=2)
USPFourier=function(x,y,M,B=999,ties.method="standard",nullstats=FALSE){
  if(length(x)!=length(y)) print("x and y should have the same length")
  J=FourierKernel(x,M); K=FourierKernel(y,M)
  return(USP(J,K,B,ties.method,nullstats))
}

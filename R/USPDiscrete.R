#' Independence test for discrete data
#'
#' Carry out a permutation independence test on a two-way contingency table.
#' The test statistic is \eqn{Tn}, as described in Sections 3.1 and 7.1 of \insertCite{BKS2020}{USP}.
#' The critical value is found by sampling null contingency tables,
#' with the same row and column totals as the input, via Patefield's algorithm, and recomputing
#' the test statistic.
#'
#' @param freq Contingency table whose independence is to be tested.
#' @param B The number of resampled null tables to be used to calibrate the test.
#' @param nullstats If TRUE, returns a vector of the null statistic values.
#'
#' @return Returns the p-value for this independence test and the value of the test statistic, \eqn{T_n},
#' as defined in \insertCite{BKS2020}{USP}. If nullstats=TRUE is used, then the function also
#' returns a vector of the null statistics.
#' @export
#'
#' @references \insertRef{BKS2020}{USP}
#'
#'
#'@importFrom Rdpack reprompt
#'@importFrom stats r2dtable
#'
#' @examples
#' freq=r2dtable(1,rep(10,5),rep(10,5))[[1]] + 4*diag(rep(1,5))
#' USPDiscrete(freq,99)
#'
#' freq=diag(1:5); USPDiscrete(freq,99)
#'
#' freq=r2dtable(1,rep(10,5),rep(10,5))[[1]];
#' test=USPDiscrete(freq,99,nullstats=TRUE)
#' plot(density(test$NullStats,from=0,
#' to=max(max(test$NullStats),test$TestStat)),
#'     xlim=c(min(test$NullStats),max(max(test$NullStats),test$TestStat)),
#'     main="Test Statistics")
#' abline(v=test$TestStat,col=2); TestStats=c(test$TestStat,test$NullStats)
#' abline(v=quantile(TestStats,probs=0.95),lty=2)
USPDiscrete=function(freq,B,nullstats=FALSE){
  freqc=colSums(freq); freqr=rowSums(freq)
  PermFreqs=r2dtable(B,freqr,freqc)
  TestStat=DiscStat(freq)
  NullStats=sapply(PermFreqs,DiscStat,simplify=TRUE)
  pval=(1+sum(NullStats>=TestStat))/(B+1)
  if(nullstats==FALSE){
    ans=list(pval,TestStat)
    names(ans)=c("pvalue","TestStat")
  }else{
    ans=list(pval,TestStat,NullStats)
    names(ans)=c("pvalue","TestStat","NullStats")
  }
  return(ans)
}

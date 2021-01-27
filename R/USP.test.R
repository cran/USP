#' Independence test for discrete data
#'
#' Carry out a permutation independence test on a two-way contingency table.
#' The test statistic is \eqn{Tn}, as described in Sections 3.1 and 7.1 of \insertCite{BKS2020}{USP}.
#' This also appears as \eqn{Un} in \insertCite{BS2021}{USP}.
#' The critical value is found by sampling null contingency tables,
#' with the same row and column totals as the input, via Patefield's algorithm, and recomputing
#' the test statistic.
#'
#' @param freq Two-way contingency table whose independence is to be tested.
#' @param B The number of resampled null tables to be used to calibrate the test.
#' @param nullstats If TRUE, returns a vector of the null statistic values.
#' @param ties.method If "standard" then calculate the p-value as in (5) of \insertCite{BKS2020}{USP},
#' which is slightly conservative. If "random" then break ties randomly. This preserves Type I error
#' control.
#'
#' @return Returns the p-value for this independence test and the value of the test statistic, \eqn{T_n},
#' as defined in \insertCite{BKS2020}{USP}. The third element of the list is the table of expected counts,
#' and the final element is the table of contributions to \eqn{T_n}. If nullstats=TRUE is used, then the function also
#' returns a vector of the null statistics.
#' @export
#'
#' @references \insertRef{BKS2020}{USP}
#' @references \insertRef{BS2021}{USP}
#'
#'
#'@importFrom Rdpack reprompt
#'@importFrom stats r2dtable
#'
#' @examples
#' freq=r2dtable(1,rep(10,5),rep(10,5))[[1]] + 4*diag(rep(1,5))
#' USP.test(freq,999)
#'
#' freq=diag(1:5); USP.test(freq,999)
#'
#' freq=r2dtable(1,rep(10,5),rep(10,5))[[1]];
#' test=USP.test(freq,999,nullstats=TRUE)
#' plot(density(test$NullStats,from=0,
#' to=max(max(test$NullStats),test$TestStat)),
#'     xlim=c(min(test$NullStats),max(max(test$NullStats),test$TestStat)),
#'     main="Test Statistics")
#' abline(v=test$TestStat,col=2); TestStats=c(test$TestStat,test$NullStats)
#' abline(v=quantile(TestStats,probs=0.95),lty=2)
USP.test=function(freq,B=999,ties.method="standard",nullstats=FALSE){
  if (is.data.frame(freq))
    freq <- as.matrix(freq)
  if(!is.matrix(freq))
    stop("'freq' must be a matrix")
  if (any(freq < 0) || anyNA(freq))
    stop("all entries of 'freq' must be nonnegative and finite")
  if ((sum(freq)) == 0)
    stop("at least one entry of 'freq' must be positive")

  freqc=colSums(freq); freqr=rowSums(freq)
  PermFreqs=r2dtable(B,freqr,freqc)
  Statistics=DiscStat(freq)
  TestStat=Statistics[[1]]
  exp=Statistics[[2]]
  contributions=Statistics[[3]]
  NullStats=unlist(sapply(PermFreqs,DiscStat,simplify=TRUE)[1,])
  if(ties.method=="random"){
    pval=(B+2-rank(c(TestStat,NullStats),ties.method = "random")[1])/(B+1)
  }else{
    pval=(1+sum(TestStat<=NullStats))/(B+1)
  }
  if(nullstats==FALSE){
    ans=list(pval,TestStat,exp,contributions)
    names(ans)=c("p.value","TestStat","expected","contributions")
  }else{
    ans=list(pval,TestStat,NullStats,exp,contributions)
    names(ans)=c("p.value","TestStat","NullStats","expected","contributions")
  }
  return(ans)
}

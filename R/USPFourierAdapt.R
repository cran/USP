#' Adaptive permutation test of independence for continuous data.
#'
#' We implement the adaptive version of the independence test for univariate continuous data
#' using the Fourier basis, as described in Section 4 of \insertCite{BKS2020}{USP}. This applies
#' USPFourier with a range of values of \eqn{M}, and a properly corrected significance
#' level.
#'
#' @param x The vector of data points from the first sample, each entry belonging to \eqn{[0,1]}.
#' @param y The vector of data points from the second sample, each entry belonging to \eqn{[0,1]}.
#' @param alpha The desired significance level of the test.
#' @param B Controls the number of permutations to be used. With a sample size of \eqn{n}  each
#' test uses \eqn{B \log_2 n} permutations. If \eqn{B+1 < 1/\alpha} then it is not possible to reject
#' the null hypothesis.
#'
#' @return Returns an indicator with value 1 if the null hypothesis of independence is rejected and
#'  0 otherwise. If the null hypothesis is rejected, the function also outputs the value of \eqn{M}
#'  at the which the null was rejected and the value of the test statistic.
#' @export
#'
#' @references \insertRef{BKS2020}{USP}
#'
#' @examples
#' n=100; w=2; x=integer(n); y=integer(n); m=300
#' unifdata=matrix(runif(2*m,min=0,max=1),ncol=2); x1=unifdata[,1]; y1=unifdata[,2]
#' unif=runif(m); prob=0.5*(1+sin(2*pi*w*x1)*sin(2*pi*w*y1)); accept=(unif<prob);
#' Data1=unifdata[accept,]; x=Data1[1:n,1]; y=Data1[1:n,2]
#' plot(x,y)
#' USPFourierAdapt(x,y,0.05,99)
USPFourierAdapt=function(x,y,alpha,B){
  a=0
  if(B+1<1/alpha){
    print("B is too small to reject the null")
    a=1
  }
  n=length(x)
  Gam=ceiling(log2(n))      # Largest model considered should have O(n^2) basis functions
  for(j in 1:Gam){
    rej=1
    M=2^(j-1)
    test=USPFourier(x,y,M,Gam*B)
    pval=test$pvalue
    TestStat=test$TestStat
    if(pval<=alpha/Gam) break
    rej=0
  }
  if(a==0){
    if(rej==1){
      print("Reject the null hypothesis")
      ans=list(rej,M,TestStat)
      names(ans)=c("Indicator","M","TestStat")
      return(ans)
    }else{
      print("Do not reject the null hypothesis")
      ans=list(rej)
      names(ans)="Indicator"
      return(ans)
    }
  }
}

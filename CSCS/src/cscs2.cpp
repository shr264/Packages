#include <cmath>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
double soft_thresh(double xx, double lambda) {
  double temp1;
  double temp2;
  temp1 = fabs(xx)-lambda;
  if(temp1<0){
    temp1 = 0.0;
  }
  if(xx>0){
    temp2 = 1.0;
  } else {
    temp2 = -1.0;
  }
  return(temp2*temp1);
}

double Tj(int j, NumericMatrix A, NumericMatrix x,double lambda, int k){
  double sumterm, out;
  sumterm = 0;
  for(int i = 1; i<=k; i++){
    if(i!=j)
      sumterm += A(i-1,j-1)*x(k-1,i-1);
      }
  out = soft_thresh(-2*sumterm, lambda)/(2*A(j-1,j-1));
  return(out);
}

double Tk(int k, NumericMatrix A, NumericMatrix x){
  double sumterm, out;
  sumterm = 0;
  for(int j = 1; j<=k; j++){
    if(j!=k)
      sumterm += A(j-1,k-1)*x(k-1,j-1);
      }
  out = (-sumterm + sqrt(sumterm*sumterm + 4*A(k-1,k-1)))/(2*A(k-1,k-1));
  return(out);
}

// [[Rcpp::export()]]
NumericMatrix CSCS(NumericMatrix S, NumericMatrix newL, double lambda, int maxitr, double tol){
  int p = S.ncol();
  NumericMatrix oldL;
  oldL= clone(newL);
  oldL(0,0) = 1/sqrt(S(1,1));
  newL(0,0) = 1/sqrt(S(1,1));
  for(int i = 1; i<=p; i++){
    int r = 1;
    int converged = 0;
    while((r<maxitr)&&(converged==0)){
      for(int j=1;j<i;j++){
        newL(i-1,j-1) = Tj(j,S,newL,lambda,i);
      }
      newL(i-1,i-1) = Tk(i,S,newL);
      if (Rcpp::max(abs(newL - oldL)) < tol)
        {
          converged = 1;
        } else {
        r += 1;
      }
      oldL = clone(newL);

    }
  }
  return(newL);
}




#include <RcppArmadillo.h> 
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;

// declare used Functions
int subCountsCpp(const double& x,const double& a,const double& b);
//*********************************************************************

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int countsCpp(const NumericVector& y,
              const double& z,
              const double& c){
  int s = y.size();
  NumericVector r(s);
  for(int j=0; j<s; j++) r[j] = subCountsCpp(y[j], z, c);
  return sum(r);
}

//-----------------------------------------------------
int subCountsCpp(const double& x,
                 const double& a,
                 const double& b){
  NumericVector res(x);
  if(x == 0){ 
    res[0] = 0;
  } else {
    for(int i=0; i<x; i++){ 
      res[i] = a*b/(i + a*b); // i - 1 + 1, C++ starts from 0
      res[i] = R::rbinom(1, res[i]);
    }
  }
  return sum(res);
}

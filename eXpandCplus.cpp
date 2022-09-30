
#include <RcppArmadillo.h> 
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;


// declare used Functions
NumericVector dirichletDistcpp(const NumericVector& x);
NumericVector resizingVec(const NumericVector& x, const double& a);
NumericMatrix resizingMat(const NumericMatrix& x, const NumericVector& a, const NumericVector& b, const double& c);


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List eXpand(const NumericMatrix& matTran,
            const NumericVector& auxVariables,
            const double& Gam,
            const double& alf,
            const NumericVector& distBase) {
  
  int supZ = 0;
  int Z = matTran.cols();
  double h(0);
  NumericVector tmpVec = distBase;
  NumericMatrix tmpMat = matTran;
  NumericVector cltmp;
  NumericVector newLtmp;
  NumericVector tvec;
  List output(4);
  
  if(max(tmpMat(_,Z-1)) > min(auxVariables)){
    do {
      supZ += 1;
      Z += 1;
      h = R::rbeta(1, Gam);
      tmpVec = resizingVec(tmpVec, h);
      h = R::rbeta(alf*tmpVec[Z-2], alf*tmpVec[Z-1]);
      cltmp = tmpMat(_,Z-2);
      tvec = alf*tmpVec;
      newLtmp = dirichletDistcpp(tvec);
      tmpMat = resizingMat(tmpMat, cltmp, newLtmp, h);
    } while(max(tmpMat(_,Z-1)) > min(auxVariables));
  }
  
  output[0] = tmpMat;
  output[1] = tmpVec;
  output[2] = Z-1;
  output[3] = supZ;
  return output;
}

// [[Rcpp::export]]
NumericVector dirichletDistcpp(const NumericVector& x){
  int n = x.size();
  NumericVector y(n);
  for(int i=0; i<n; i++) y[i] = R::rgamma(x[i], 1);
  y = y/sum(y);
  return y;
}

//****************
NumericVector resizingVec(const NumericVector& x,
                          const double& a){
  int n = x.size() ;
  NumericVector y(n+1);
  for( int i=0; i<n-1; i++) y[i] = x[i];
  y[n-1] = a*x[n-1];
  y[n] = (1-a)*x[n-1];
  return y;
}
//****************
NumericMatrix resizingMat(const NumericMatrix& x,
                          const NumericVector& a,
                          const NumericVector& b,
                          const double& c){
  int n = x.rows();
  int m = x.cols();
  NumericMatrix y(n+1,m+1);
  for(int i=0; i<n; i++){
    for(int j=0; j<m-1; j++){
      y(i,j) = x(i,j);
    }
  }
  y(_,m-1) = c * a;
  y(_,m) = (1 - c) * a;
  y(n,_) = b;
  
  return y;
}


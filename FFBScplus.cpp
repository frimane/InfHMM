
#include <RcppArmadillo.h> 
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

//*********************************************************************

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec FF_BS(const mat& obsL, 
          const mat& matTran,
          const vec& iniP, 
          const vec& auxVariables) {
  
  //   obsL: the likelihood of data
  //   matTran: the transition matrix
  //   iniP: the initial distribution of states
  //   auxVariables: the vector of auxiliary variables
  
  
  //   S: the number of states
  //   N: the length of data
  int N = obsL.n_cols;
  int S = obsL.n_rows;
  
  // Forward-filtering pass with normalization
  mat alphas = zeros(S,N);
  alphas.col(0) = obsL.col(0)%iniP;
  alphas.col(0) /= accu(alphas.col(0));
  
  mat tmp_idx = zeros(S,S);
  for (int n = 1; n < N; ++n) {
    
    for(int i = 0; i < S; ++i){ 
      for(int j = 0; j < S; ++j){
        if (matTran(i,j) > auxVariables[n]) {
          tmp_idx(i,j) = 1;
        } else {
          tmp_idx(i,j) = 0;
        }
      }
    }
    
    alphas.col(n) = obsL.col(n) % (tmp_idx.t()*alphas.col(n-1));
    alphas.col(n) /= accu(alphas.col(n));
  }
  
  
  // Backwards-sampling pass with running normalization
  vec hidden_states = ones(N);
  uvec ID = regspace<uvec>(1,1,S);
  NumericVector idx = wrap(ID);
  NumericVector proba = wrap(alphas.col(N-1));
  // // Sample last hidden state
  hidden_states[N-1] = sample(idx,1,FALSE,proba)[0];
  
  // // Sample nth hidden state conditional on (n+1)st hidden state
  vec tmp_idx_vec = zeros(S);
  for (int n = N-2; n >= 0; --n) {
    for (int i = 0; i < S; ++i) {
      if(matTran(i,hidden_states[n+1]-1) > auxVariables[n+1]) {
        tmp_idx_vec[i] = 1;
      } else {
        tmp_idx_vec[i] = 0;
      }
    }
    proba = wrap(alphas.col(n) % tmp_idx_vec);
    hidden_states[n] = sample(idx, 1, FALSE, proba)[0];
  }
  
  return hidden_states;
}


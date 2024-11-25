#include <RcppArmadillo.h>
using namespace arma;
using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
vec distance(const mat &pixel, const mat &ref){
  int n = pixel.n_rows;
  int m = ref.n_rows;
  vec output = zeros<vec>(n);
  for(int i = 0; i < n; i++){
    vec dist_i = zeros<vec>(m);
    for(int j = 0; j < m; j++){
      dist_i(j) = sqrt(sum(pow(pixel.row(i).t() - ref.row(j).t(), 2)));
    }
    output(i) = min(dist_i);
  }
  return output;
}

#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
long double dnbinomCpp(long double x, long double mu, long double size){
  long double prob = size / (size + mu);
  long double dres = tgamma(x + size)/(tgamma(size) * tgamma(x + 1.0)) * pow(prob, size) * pow(1.0-prob, x);
	return dres;
}

// [[Rcpp::export]]
long double productCpp(NumericVector arr, int n) { 
    long double result = 1.0; 
    for (int i = 0; i < n; i++) {
        result = result * arr[i]; 
    }
    return result; 
}

// [[Rcpp::export]]
long double mcintCpp(NumericVector x, NumericVector mu, long double size, int n){
  NumericVector result_arr(n);
  long double result;
  for (int i = 0; i < n; i++) {
  	result_arr[i] = dnbinomCpp(x[i], mu[i], size);
  }
  result = productCpp(result_arr, n);
  return result;
}


// #include <Rcpp.h>
// #include <RcppEigen.h>
// #include <RcppGSL.h>
// #include <gsl/gsl_randist.h>
// #include <gsl/gsl_cdf.h>
// #include <cmath>


#include <RcppArmadillo.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include<omp.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]


// [[Rcpp::export]]
int posel(int po, mat &X){
  vec tes = zeros(2*po);
  int inc = 0;
  for(int i = 0; i < (2*po); i++){
    tes(i) = accu(pow(cov(pow(atan(X),(i+1)/10))-cov(X),2)); //
    if(i > 0){
      if(tes(i) > tes(i-1)){
        inc = inc + 1;
      }
      if(tes(i) < tes(i-1)){
        inc = 0;
      }
    }
    if(inc > 10){ //the distance is increasing for last 10 iterations so stop.
      break;
    }
  }
  
  double ret = (index_min( tes.elem( find(tes > 0) ) )+1)/10;
  
  return ret;
}
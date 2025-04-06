// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppDist)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <RcppDist.h>

using namespace Rcpp;
using namespace arma;

// ll: Poisson log-likelihood with atan interaction
// [[Rcpp::export]]
double ll_cpp(const arma::vec& ind, const arma::umat& index, double lambda, double beta) {
  double sum_poisson = arma::accu(log(dpois(ind, lambda)) + lambda);
  int idx1 = index(0, 0) - 1;  // 1-based indexing in R
  int idx2 = index(1, 0) - 1;
  double interaction = beta * std::atan(ind(idx1)) * std::atan(ind(idx2));
  return sum_poisson + interaction;
}

// numerlike: summing exp(ll) over all rows of matrix
// [[Rcpp::export]]
double numerlike_cpp(const arma::mat& mat, const arma::umat& index, double lambda, double beta) {
  int rows = mat.n_rows;
  double total = 0.0;
  for (int i = 0; i < rows; ++i) {
    arma::vec row = mat.row(i).t();
    total += std::exp(ll_cpp(row, index, lambda, beta));
  }
  return total;
}

// numerlike1: numerical approx. of integral using poisson weights + atan interaction
// [[Rcpp::export]]
double numerlike1_cpp(int i, int j, const arma::mat& X, const arma::vec& lambda,
                      const arma::vec& beta, int po, int apnum = 100) {
  int c = X.n_cols;
  arma::mat Beta = arma::zeros(c, c);
  int bidx = 0;
  for (int row = 1; row < c; ++row) {
    for (int col = 0; col < row; ++col) {
      Beta(row, col) = beta(bidx);
      Beta(col, row) = beta(bidx);
      bidx++;
    }
  }
  
  arma::vec ks = arma::regspace(0, apnum);
  arma::vec poisd = dpois(ks, lambda(j)) * std::exp(lambda(j));
  double total = 0.0;
  
  for (int k = 0; k <= apnum; ++k) {
    double inner = 0.0;
    for (int col = 0; col < c; ++col) {
      if (col == j) continue;
      inner += -Beta(j, col) * std::pow(std::atan(X(i, col)), po);
    }
    total += poisd(k) * std::exp(inner * std::pow(std::atan(k), po));
  }
  
  return total;
}

// llhoodl: log-likelihood of lambda
// [[Rcpp::export]]
double llhoodl_cpp(int j, const arma::mat& X, const arma::vec& lambda, const arma::vec& beta,
                   int po, int apnum = 100) {
  int c = X.n_cols;
  int Ti = X.n_rows;
  arma::mat Beta = arma::zeros(c, c);
  int bidx = 0;
  for (int row = 1; row < c; ++row) {
    for (int col = 0; col < row; ++col) {
      Beta(row, col) = beta(bidx);
      Beta(col, row) = beta(bidx);
      bidx++;
    }
  }
  
  double sum_ret = 0.0;
  for (int k = 0; k < Ti; ++k) {
    double term = -lambda(j) + X(k, j) * std::log(lambda(j)) + lambda(j);
    term -= std::log(numerlike1_cpp(k, j, X, lambda, beta, po, apnum));
    sum_ret += term;
  }
  
  return sum_ret;
}

// llhoodb: log-likelihood of beta
// [[Rcpp::export]]
double llhoodb_cpp(int j, const arma::mat& X, const arma::vec& lambda, const arma::vec& beta,
                   int po, int apnum = 100) {
  int c = X.n_cols;
  int Ti = X.n_rows;
  arma::mat Beta = arma::zeros(c, c);
  int bidx = 0;
  for (int row = 1; row < c; ++row) {
    for (int col = 0; col < row; ++col) {
      Beta(row, col) = beta(bidx);
      Beta(col, row) = beta(bidx);
      bidx++;
    }
  }
  
  double total = 0.0;
  for (int i = 0; i < Ti; ++i) {
    double term = 0.0;
    for (int k = 0; k < c; ++k) {
      if (k == j) continue;
      term += -Beta(j, k) * std::pow(std::atan(X(i, k)), po);
    }
    term *= std::pow(std::atan(X(i, j)), po);
    total += term - std::log(numerlike1_cpp(i, j, X, lambda, beta, po, apnum));
  }
  
  return total;
}

// atanmean: Expected value of atan^po under Poisson(theta)
// [[Rcpp::export]]
double atanmean_cpp(double theta, int po, int apnum = 100) {
  arma::vec ks = arma::regspace(0, apnum);
  arma::vec atan_k = arma::pow(arma::atan(ks), po);
  arma::vec pois = dpois(ks, theta);
  return arma::dot(atan_k, pois);
}


// [[Rcpp::export]]
List update_beta_mh_cpp(arma::mat Beta, arma::vec beta, const arma::mat& Z,
                        double s1, double s0, const arma::mat& pdx,
                        const arma::vec& pdxid, const arma::mat& X,
                        const arma::vec& lambda, const arma::uvec& index,
                        int Ti, int po, double lambdashrk) {
  
  int c = X.n_cols;
  int betalen = beta.n_elem;
  arma::mat bsigma = arma::zeros(c, c);
  arma::vec s1vec = arma::ones(betalen) * s1;
  arma::vec s0vec = arma::ones(betalen) * s0;
  arma::vec sigmas = arma::vectorise(Z) % s1vec + (1 - arma::vectorise(Z)) % s0vec;
  
  int idx = 0;
  for (int row = 1; row < c; ++row) {
    for (int col = 0; col < row; ++col) {
      bsigma(row, col) = sigmas(idx);
      idx++;
    }
  }
  bsigma = bsigma + bsigma.t();
  
  arma::mat Betac = Beta;
  int acbeta = 0;
  
  for (int i = 0; i < c; ++i) {
    arma::vec mean = -pdx.row(i).t();
    mean.shed_row(i); // drop i-th entry
    
    arma::vec varc = bsigma.row(i).t();
    varc.shed_row(i);
    
    arma::mat varctemp = Beta;
    varctemp.shed_row(i);
    varctemp.shed_col(i);
    varctemp.diag() = pdxid;
    varctemp.diag().shed_row(i);
    
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, varctemp);
    
    arma::vec varcei = (arma::var(arma::pow(arma::atan(X.col(i)), po)) * Ti + lambdashrk) /
      arma::abs(eigval) + 1.0 / varc;
    
    arma::mat varceiU = eigvec.t();
    for (size_t r = 0; r < varceiU.n_rows; ++r) {
      varceiU.row(r) /= std::sqrt(std::abs(varcei(r)));
    }
    
    arma::mat varC = varceiU.t() * varceiU;
    
    arma::vec betac = rmvnorm(1, varC * mean, varC).t();
    
    // Replace any NaNs
    betac.elem(find_nonfinite(betac)).zeros();
    
    Betac.row(i).cols(0, i - 1) = betac.head(i).t();
    Betac.row(i).cols(i + 1, c - 1) = betac.tail(c - i - 1).t();
    Betac.col(i).rows(0, i - 1) = betac.head(i);
    Betac.col(i).rows(i + 1, c - 1) = betac.tail(c - i - 1);
    
    arma::vec betavec_old = Beta.row(i).t();
    betavec_old.shed_row(i);
    arma::vec betavec_new = Betac.row(i).t();
    betavec_new.shed_row(i);
    
    double R = llhoodb_cpp(i, X, lambda, Betac, po) - llhoodb_cpp(i, X, lambda, beta, po);
    
    R += arma::accu(dnorm(betavec_new, 0.0, varc, true) -
      dnorm(betavec_old, 0.0, varc, true));
    
    double Q = dmvnorm(betavec_old, varC * mean, varC, true) -
      dmvnorm(betavec_new, varC * mean, varC, true);
    
    R += Q;
    
    if (std::isnan(R) || std::isinf(R)) R = 1.0;
    
    if (std::log(R::runif(0, 1)) < R) {
      beta = vectorise(Betac(index));
      Beta = Betac;
      acbeta++;
    }
  }
  
  return List::create(
    Named("Beta") = Beta,
    Named("beta") = beta,
    Named("acbeta") = acbeta
  );
}

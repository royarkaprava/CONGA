// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppDist)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <RcppDist.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double ll_cpp(const arma::vec& ind, const arma::umat& index, double lambda, double beta) {
  double sum_poisson = 0.0;
  for (int i = 0; i < ind.n_elem; ++i) {
    sum_poisson += R::dpois(ind[i], lambda, true);  // log = true
  }
  sum_poisson += lambda * ind.n_elem;  // same as sum(lambda)
  
  int idx1 = index(0, 0) - 1;  // R is 1-based
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


// [[Rcpp::export]]
double numerlike1_cpp(int i, int j, const arma::mat& X, const arma::vec& lambda,
                      const arma::mat& Beta, int po, int apnum = 100) {
  
  int c = X.n_cols;
  
  // Create the k values and Poisson part (Poisson probabilities with lambda)
  arma::vec ks = arma::regspace(0, apnum);
  arma::vec poisd(apnum + 1);
  for (int k = 0; k <= apnum; ++k) {
    poisd[k] = R::dpois(ks[k], lambda(j), true);  // log Poisson probability
  }
  poisd *= std::exp(lambda(j));  // Multiply by exp(lambda)
  
  double total = 0.0;
  
  // Loop over k (Poisson counts)
  for (int k = 0; k <= apnum; ++k) {
    double inner = 0.0;
    
    // Precompute atan(X(i, col)) for efficiency
    arma::vec atans = arma::atan(X.row(i).t());
    
    for (int col = 0; col < c; ++col) {
      inner += -Beta(j, col) * std::pow(atans(col), po);
    }
    
    total += poisd[k] * std::exp(inner * std::pow(std::atan(k), po));
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
double llhoodb_cpp(int j, const arma::mat& X, const arma::vec& lambda, const arma::mat& Beta,
                   int po, int apnum = 100) {
  int c = X.n_cols;
  int Ti = X.n_rows;
  double total = 0.0;
  for (int i = 0; i < Ti; ++i) {
    double term = 0.0;
    for (int k = 0; k < c; ++k) {
      term += -Beta(j, k) * std::pow(std::atan(X(i, k)), po);
    }
    term *= std::pow(std::atan(X(i, j)), po);
    total += term - std::log(numerlike1_cpp(i, j, X, lambda, Beta, po, apnum));
  }
  
  return total;
}

// atanmean: Expected value of atan^po under Poisson(theta)
// [[Rcpp::export]]
double atanmean_cpp(double theta, int po, int apnum = 100) {
  arma::vec ks = arma::regspace(0, apnum);  // 0 to apnum
  arma::vec atan_k = arma::pow(arma::atan(ks), po);  // atan(k) ^ po
  arma::vec pois(apnum + 1);
  
  // Loop to compute Poisson probabilities for each k
  for (int k = 0; k <= apnum; ++k) {
    pois[k] = R::dpois(ks[k], theta, false);  // log = false (density, not log)
  }
  
  return arma::dot(atan_k, pois);  // Dot product of atan(k)^po and Poisson
}

// Function to calculate the log of the multivariate normal density
double dmvnorm_cpp(const arma::vec& x, const arma::vec& mean, const arma::mat& cov) {
  int d = x.n_elem;
  arma::vec diff = x - mean;
  
  // Calculate the log determinant of the covariance matrix and its inverse
  double log_det_cov;
  arma::mat inv_cov;
  
  // Calculate the Mahalanobis distance
  double mahalanobis = as_scalar(diff.t() * inv_cov * diff);
  
  // Log of the multivariate normal density
  double log_density = -0.5 * (mahalanobis);
  
  return log_density;
}

// Function to compute the difference of log-densities (Q)
double compute_Q(const arma::vec& betavec_old, const arma::vec& betavec_new,
                 const arma::vec& mean, const arma::mat& varC) {
  // Compute the log-density for betavec_old and betavec_new
  double log_density_old = dmvnorm_cpp(betavec_old, mean, varC);
  double log_density_new = dmvnorm_cpp(betavec_new, mean, varC);
  
  // Return the difference
  return log_density_old - log_density_new;
}

arma::vec sample_canonical_mvnorm(const arma::mat& Q, const arma::mat& Qhalf, const arma::vec& eta) {
  //arma::mat L = arma::chol(Q); // Q = L * L^T
  arma::vec mu = arma::solve(Q, eta); // mu = Q^{-1} * eta
  arma::vec z = arma::randn(Q.n_rows);
  arma::vec y = arma::solve(Qhalf, z); // L^T y = z
  return mu + y;
}

// [[Rcpp::export]]
List update_beta_mh_cpp(arma::mat Beta,
                        double s1, double s0, const arma::mat& pdx,
                        const arma::vec& pdxid, const arma::mat& X,
                        const arma::vec& lambda, const arma::uvec& index,
                        int Ti, double po, double lambdashrk) {
  
  int c = X.n_cols;
  int betalen = c*(c-1)/2;
  arma::mat bsigma = arma::zeros(c, c);
  
  // Create the sigmas vector for the bsigma matrix
  arma::vec s1vec = arma::ones(betalen) * s1;
  arma::vec s0vec = arma::ones(betalen) * s0;
  arma::vec sigmas = s1vec; //+ (1 - arma::vectorise(Z)) % s0vec;
  
  // Populate bsigma matrix symmetrically
  int idx = 0;
  for (int row = 1; row < c; ++row) {
    for (int col = 0; col < row; ++col) {
      bsigma(row, col) = sigmas(idx);
      bsigma(col, row) = sigmas(idx);
      idx++;
    }
  }
  bsigma = bsigma + bsigma.t();
  
  arma::mat Betac = Beta;
  int acbeta = 0;
  
  // Metropolis-Hastings updates
  for (int i = 0; i < c; ++i) {
    // Prepare mean and variance for the proposal distribution
    arma::vec mean = -pdx.row(i).t();
    mean.shed_row(i);  // Exclude the i-th element for the mean
    
    arma::vec varc = bsigma.row(i).t();
    varc.shed_row(i);  // Exclude the i-th element for the variance
    
    arma::mat varctemp = Beta;
    varctemp.diag() = pdxid;
    varctemp.shed_row(i);
    varctemp.shed_col(i);
    
    // Eigen decomposition of varctemp
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, varctemp);
    
    // Calculate inverse covariance
    arma::vec varcei = (arma::var(arma::pow(arma::atan(X.col(i)), po)) * Ti + lambdashrk) /
      arma::abs(eigval) + 1.0 / varc;
    
    arma::mat varceiU = eigvec.t();
    for (size_t r = 0; r < varceiU.n_rows; ++r) {
      varceiU.row(r) *= std::sqrt(std::abs(varcei(r)));
    }
    
    arma::mat varCinv = varceiU.t() * varceiU;
    
    // Sample from the proposal distribution
    arma::vec betac = sample_canonical_mvnorm(varCinv, varceiU, mean); //rmvnorm(1, varC * mean, varC).t();
    
    // Replace any NaNs in betac
    betac.elem(find_nonfinite(betac)).zeros();
    
    arma::vec betacc(c, arma::fill::zeros);
    if (i > 0) {
      betacc.subvec(0, i - 1) = betac.subvec(0, i - 1);
    }
    if (i < c - 1) {
      betacc.subvec(i + 1, c - 1) = betac.subvec(i, c - 2);
    }
    
    Betac.row(i) = betacc.t(); // betac has to be updated to match the size
    // Update the first column, excluding the diagonal
    Betac.col(i) = betacc;
    // Update the Beta matrix
   
    arma::vec betavec_old = Beta.row(i).t();
    betavec_old.shed_row(i);
    arma::vec betavec_new = Betac.row(i).t();
    betavec_new.shed_row(i);
    
    //Compute the acceptance ratio
    double R = llhoodb_cpp(i, X, lambda, Betac, po) - llhoodb_cpp(i, X, lambda, Beta, po);
    //double R = 0;
    // Correcting dnorm application for each element of betavec_new and betavec_old
    arma::vec norm_new(betavec_new.n_elem);
    arma::vec norm_old(betavec_old.n_elem);

    // Loop through each element of betavec_new and betavec_old and compute the normal densities
    // Loop through each element of betavec_new and betavec_old and compute the normal densities
    // for (size_t i = 0; i < betavec_new.n_elem; ++i) {
    //   norm_new(i) = R::dnorm(betavec_new(i), 0.0, varc(i), true);  // Use varc(i) as the std deviation
    //   norm_old(i) = R::dnorm(betavec_old(i), 0.0, varc(i), true);  // Use varc(i) as the std deviation
    // }
    // Sum the elements to get the scalar result for the acceptance ratio
    R += arma::accu(-pow(betavec_new,2)/varc + pow(betavec_old,2)/varc)/2;
    
    arma::vec Q = betavec_old%(-varCinv * betavec_old/2+mean)-betavec_new%(-varCinv * betac/2 +mean);
      
    //double Q = dmvnorm_cpp(betavec_old, varC * mean, varC) -
    //  dmvnorm_cpp(betavec_new, varC * mean, varC);

    R += arma::accu(Q);


    if (std::isnan(R) || std::isinf(R)) R = 1.0;

    // Accept or reject the proposal
    if (R::runif(0, 1) < R) {
      Beta = Betac;
      acbeta++;
    }
  }
  
  return List::create(
    Named("Beta") = Beta,
    Named("acbeta") = acbeta
  );
}

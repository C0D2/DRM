#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// helper: maximum over an IntegerMatrix ignoring NAs
int max_data(const IntegerMatrix& data) {
  int maxval = INT_MIN;
  for (int n = 0; n < data.nrow(); n++) {
    for (int i = 0; i < data.ncol(); i++) {
      int v = data(n, i);
      if (v != NA_INTEGER && v > maxval) {
        maxval = v;
      }
    }
  }
  if (maxval == INT_MIN) return -1; // indicate no non-NA values
  return maxval;
}
////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
NumericVector P_2PL_cpp(const NumericMatrix& theta,
                        const NumericVector& a,
                        double b) {
  int M = theta.nrow();
  int d = theta.ncol();
  NumericVector prob(M);

  for (int m = 0; m < M; m++) {
    double linpred = b;
    for (int j = 0; j < d; j++) {
      linpred += theta(m, j) * a[j];
    }
    prob[m] = 1.0 / (1.0 + std::exp(-linpred));
  }
  return prob;
}

// [[Rcpp::export]]
List P_DRM_cpp(const NumericMatrix& theta,
                        const NumericVector& a,
                        double b,
                        double nu,
                        const NumericVector& cut_score,
                        bool return_mu = false) {
  int M = theta.nrow();
  int d = theta.ncol();
  int ncats = cut_score.size() + 1;

  NumericMatrix probs(M, ncats);
  NumericVector mu(M);
  NumericVector one_minus_mu(M);

  const double LINPRED_MAX = 700.0;
  const double LINPRED_MIN = -700.0;
  const double P_EPS = 1e-12;

  // Compute 2PL probability (mu)
  for (int m = 0; m < M; m++) {
    double linpred = b;
    for (int j = 0; j < d; j++) linpred += theta(m, j) * a[j];
    linpred = std::max(std::min(linpred, LINPRED_MAX), LINPRED_MIN);
    mu[m] = 1.0 / (1.0 + std::exp(-linpred));
    one_minus_mu[m] = 1.0 - mu[m];
  }

  // Compute cumulative probabilities
  NumericMatrix cumprobs(M, ncats - 1);
  for (int i = 0; i < cut_score.size(); i++) {
    double q = cut_score[i];
    for (int m = 0; m < M; m++) {
      cumprobs(m, i) = R::pbeta(q, mu[m] * nu, one_minus_mu[m] * nu, 1, 0);
      if (cumprobs(m, i) < P_EPS) cumprobs(m, i) = P_EPS;
      if (cumprobs(m, i) > 1.0 - P_EPS) cumprobs(m, i) = 1.0 - P_EPS;
    }
  }

  // Convert cumulative probs to category probabilities
  for (int m = 0; m < M; m++) {
    probs(m, 0) = cumprobs(m, 0);
    for (int i = 1; i < ncats - 1; i++) {
      probs(m, i) = cumprobs(m, i) - cumprobs(m, i - 1);
    }
    probs(m, ncats - 1) = 1.0 - cumprobs(m, ncats - 2);
  }

  if (return_mu) {
    return List::create(
      _["probs"] = probs,
      _["mu"]    = mu
    );
  } else {
    return List::create(
      _["probs"] = probs
    );
  }
}

// [[Rcpp::export]]
NumericMatrix llik_cpp(const IntegerVector& data,        // N x 1
                       const NumericMatrix& theta,       // M x d
                       const NumericVector& item,        // d + 2
                       const NumericVector& cut_score) { // cut points
  int N = data.size();
  int M = theta.nrow();

  NumericVector a(theta.ncol());
  for (int j = 0; j < theta.ncol(); j++) a[j] = item[j];
  double b = item[theta.ncol()];
  double nu = item[theta.ncol() + 1];

  List tmp = P_DRM_cpp(theta, a, b, nu, cut_score, false);
  NumericMatrix p_matrix = tmp["probs"];
  int K = p_matrix.ncol();

  NumericMatrix loglik(N, M);

  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      int y = data[n];
      double prob;
      if (y == NA_INTEGER) {
        // treat missing response: set prob to 1 so it doesn't affect likelihood
        prob = 1.0;
      } else if (y >= 0 && y < K) {
        prob = p_matrix(m, y);
      } else {
        prob = 1.0;
      }
      loglik(n, m) = std::log(std::max(prob, 1e-12));
    }
  }


  return loglik;
}

// [[Rcpp::export]]
NumericVector llik_ratio_mhrm(const NumericMatrix& data,        // N x I
                              const NumericMatrix& theta_c,       // N x d - current theta
                              const NumericMatrix& theta_p,       // N x d - proposed theta
                              const NumericMatrix& item,        // I x (d + 2)
                              const List& cut_score,     // cut points
                              const NumericMatrix& Sigma_inv){
  int N = data.nrow();
  int I = data.ncol();
  int d = theta_c.ncol();


  NumericVector loglik_ratio(N);

  // initialize with log prior ratio
  for (int n = 0; n < N; n++) {
    double quad_p = 0.0, quad_c = 0.0;
    for (int j = 0; j < d; j++) {
      for (int k = 0; k < d; k++) {
        quad_p += theta_p(n,j) * Sigma_inv(j,k) * theta_p(n,k);
        quad_c += theta_c(n,j) * Sigma_inv(j,k) * theta_c(n,k);
      }
    }
    loglik_ratio[n] = -0.5 * (quad_p - quad_c);
  }

  for (int i = 0; i < I; i++) {
    NumericVector a(d);
    for (int j = 0; j < d; j++) a[j] = item(i,j);
    double b = item(i, d);
    double nu = item(i, d + 1);
    NumericVector cut_i = Rcpp::as<NumericVector>(cut_score[i]);

    List tmp_c = P_DRM_cpp(theta_c, a, b, nu, cut_i, false);
    List tmp_p = P_DRM_cpp(theta_p, a, b, nu, cut_i, false);
    const NumericMatrix& p_matrix_c = tmp_c["probs"];
    const NumericMatrix& p_matrix_p = tmp_p["probs"];
    int K = p_matrix_c.ncol();

    for (int n = 0; n < N; n++) {
      int y = data(n,i);
      if (IntegerVector::is_na(y) || y < 0 || y >= K) continue;

      double pc = std::max(p_matrix_c(n, y), 1e-12);
      double pp = std::max(p_matrix_p(n, y), 1e-12);

      loglik_ratio[n] += std::log(pp) - std::log(pc);
    }
  }

  return loglik_ratio;
}

// [[Rcpp::export]]
double  llik_mhrm(const NumericMatrix& data,        // N x I
                  const NumericMatrix& theta_c,       // N x d - current theta
                  const NumericMatrix& item,        // I x (d + 2)
                  const List& cut_score,     // cut points
                  const NumericMatrix& Sigma_inv){
  int N = data.nrow();
  int I = data.ncol();
  int d = theta_c.ncol();


  double loglik = 0.0;

  // initialize with log prior ratio
  for (int n = 0; n < N; n++) {
    double quad_c = 0.0;
    for (int j = 0; j < d; j++) {
      for (int k = 0; k < d; k++) {
        quad_c += theta_c(n,j) * Sigma_inv(j,k) * theta_c(n,k);
      }
    }
    loglik += -0.5 * quad_c;
  }

  for (int i = 0; i < I; i++) {
    NumericVector a(d);
    for (int j = 0; j < d; j++) a[j] = item(i,j);
    double b = item(i, d);
    double nu = item(i, d + 1);
    NumericVector cut_i = Rcpp::as<NumericVector>(cut_score[i]);

    List tmp_c = P_DRM_cpp(theta_c, a, b, nu, cut_i, false);
    const NumericMatrix& p_matrix_c = tmp_c["probs"];
    int K = p_matrix_c.ncol();

    for (int n = 0; n < N; n++) {
      int y = data(n,i);
      if (IntegerVector::is_na(y) || y < 0 || y >= K) continue;

      double pc = std::max(p_matrix_c(n, y), 1e-12);

      loglik += std::log(pc);
    }
  }

  return loglik;
}

// [[Rcpp::export]]
List Estep_cpp(IntegerMatrix data,
               NumericMatrix item,
               NumericMatrix grid,
               NumericVector prior,
               List cut_score) {
  int N = data.nrow();
  int I = data.ncol();
  int M = grid.nrow();

  // exponentiate the last column at once
  NumericMatrix item_exp = clone(item);
  item_exp(_, item_exp.ncol() - 1) = exp(item_exp(_, item_exp.ncol() - 1));

  // posterior accumulator
  NumericMatrix posterior(N, M);

  // accumulate over items
  for (int i = 0; i < I; i++) {
    IntegerVector data_col = data(_, i);
    NumericVector item_row = item_exp(i, _);

    // ith element of cut_score
    NumericVector cut_i = cut_score[i];

    // pass cut_i to llik_cpp
    NumericMatrix ll = llik_cpp(data_col, grid, item_row, cut_i); // N×M

    for (int n = 0; n < N; n++) {
      for (int m = 0; m < M; m++) {
        posterior(n, m) += ll(n, m);
      }
    }
  }

  // save logL
  NumericMatrix logL = clone(posterior);

  // posterior = exp(posterior + log(prior))
  for (int n = 0; n < N; n++) {
    for (int m = 0; m < M; m++) {
      posterior(n, m) = std::exp(posterior(n, m) + std::log(prior[m]));
    }
  }

  // normalize rows
  for (int n = 0; n < N; n++) {
    double rowsum = 0.0;
    for (int m = 0; m < M; m++) rowsum += posterior(n, m);
    if (rowsum > 0) {
      for (int m = 0; m < M; m++) posterior(n, m) /= rowsum;
    }
  }

  // expected response frequencies
  int categ = max_data(data) + 1;
  NumericVector e_response(I * M * categ);
  //#pragma omp parallel for
  for (int c = 0; c < categ; c++) {
    for (int i = 0; i < I; i++) {
      for (int m = 0; m < M; m++) {
        double sumval = 0.0;
        for (int n = 0; n < N; n++) {
          int val = data(n, i);
          if (val != NA_INTEGER && val == c) {
            sumval += posterior(n, m);
          }
        }
        e_response[i + I * m + I * M * c] = sumval;
      }
    }
  }

  // freq = column sums
  NumericVector freq(M);
  for (int m = 0; m < M; m++) {
    double s = 0.0;
    for (int n = 0; n < N; n++) s += posterior(n, m);
    freq[m] = s;
  }

  // logL
  double total_logL = 0.0;
  for (int n = 0; n < N; n++) {
    for (int m = 0; m < M; m++) {
      total_logL += logL(n, m) * posterior(n, m);
    }
  }
  for (int m = 0; m < M; m++) {
    total_logL += freq[m] * std::log(prior[m]);
  }
  for (int n = 0; n < N; n++) {
    for (int m = 0; m < M; m++) {
      double val = posterior(n, m);
      if (val > 0) total_logL -= val * std::log(val);
    }
  }

  // Ak
  double total_freq = std::accumulate(freq.begin(), freq.end(), 0.0);
  NumericVector Ak(M);
  if (total_freq > 0) {
    for (int m = 0; m < M; m++) Ak[m] = freq[m] / total_freq;
  }

  return List::create(
    _["posterior"]  = posterior,
    _["freq"]       = freq,
    _["e.response"] = e_response,
    _["grid"]       = grid,
    _["prior"]      = prior,
    _["logL"]       = total_logL,
    _["Ak"]         = Ak
  );
}






// [[Rcpp::export]]
NumericMatrix compute_pmat(const NumericVector& beta_grid,
                           const NumericVector& mu,
                           double nu) {
  int L = beta_grid.size();   // length(beta_grid)
  int M = mu.size();          // length(mu) == nrow(grid)
  NumericMatrix pmat(M, L);   // result: M x L (matches t(outer(...)) in R)

  for (int m = 0; m < M; ++m) {
    double a = nu * mu[m];
    double b = nu * (1.0 - mu[m]);

    // Beta(a,b) = Gamma(a)*Gamma(b)/Gamma(a+b)
    double beta_ab = R::gammafn(a) * R::gammafn(b) / R::gammafn(a + b);

    for (int l = 0; l < L; ++l) {
      // pow is safe here because beta_grid is assumed in (0,1)
      double p1 = std::pow(beta_grid[l], a - 1.0);
      double p2 = std::pow(1.0 - beta_grid[l], b - 1.0);
      pmat(m, l) = (p1 * p2) / beta_ab;
    }
  }

  return pmat; // M x L, same orientation as the R t(outer(...)) result
}

/*** R
# quick test in R
# beta_grid <- seq(1/2000, 1-1/2000, length.out = 10)
# mu <- c(0.2, 0.5, 0.8)
# nu <- 5
# pm <- compute_pmat(beta_grid, mu, nu)
# pm      # should be 3 x 10
# t(outer(beta_grid, nu*mu-1, FUN = "^")*outer(1-beta_grid, nu*(1-mu)-1, FUN = "^"))/beta(nu*mu,nu*(1-mu))
*/


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List compute_l1_cpp(const arma::mat& grid,
                          double nu,
                          const arma::vec& mu,
                          const arma::vec& beta_grid,
                          int ncats,
                          const arma::ivec& ind_cat) {

  int M = grid.n_rows;
  int d = grid.n_cols;
  int K = beta_grid.n_elem;

  // --- precompute digamma ---
  arma::vec dig_nu_mu(M), dig_nu_1mu(M);
  for (int m = 0; m < M; m++) {
    dig_nu_mu[m]  = R::digamma(nu * mu[m]);
    dig_nu_1mu[m] = R::digamma(nu * (1.0 - mu[m]));
  }
  double dig_nu = R::digamma(nu);

  // --- log transforms of beta_grid ---
  arma::vec logb = arma::log(beta_grid);
  arma::vec log1b = arma::log(1.0 - beta_grid);
  arma::vec logit = logb - log1b;

  const double p_floor = 1e-300;
  const double LOG_MIN = -700.0;

  // --- initialize matrices ---
  arma::mat pmat(M, K, fill::zeros);
  arma::mat l1mu(M, K, fill::zeros);
  arma::mat l1xi(M, K, fill::zeros);

  // --- compute pmat, l1mu, l1xi in parallel ---
#pragma omp parallel for
  for (int m = 0; m < M; m++) {
    double a = nu * mu[m];
    double b = nu * (1.0 - mu[m]);
    double a1 = a - 1.0;
    double b1 = b - 1.0;
    double log_beta_ab = R::lbeta(a, b);
    double diff_dig = dig_nu_mu[m] - dig_nu_1mu[m];

    for (int k = 0; k < K; k++) {
      double log_p = a1 * logb[k] + b1 * log1b[k] - log_beta_ab;
      if (!std::isfinite(log_p) || log_p < LOG_MIN) log_p = LOG_MIN;
      if (log_p < LOG_MIN) log_p = LOG_MIN;

      double p_val = std::exp(log_p);
      if (p_val < p_floor) p_val = p_floor;
      pmat(m, k) = p_val;

      l1mu(m, k) = nu * p_val * (logit[k] - diff_dig);
      double term = dig_nu + mu[m]*(logb[k]-dig_nu_mu[m]) + (1.0-mu[m])*(log1b[k]-dig_nu_1mu[m]);
      l1xi(m, k) = nu * term * p_val;
    }
  }

  // --- collapse to categories in a thread-safe way ---
  arma::mat l1m(M, ncats, fill::zeros);
  arma::mat l1x(M, ncats, fill::zeros);

#pragma omp parallel
{
  arma::mat l1m_private(M, ncats, fill::zeros);
  arma::mat l1x_private(M, ncats, fill::zeros);

#pragma omp for nowait
  for (int k = 0; k < K; k++) {
    int ct = ind_cat[k] - 1;
    for (int m = 0; m < M; m++) {
      l1m_private(m, ct) += l1mu(m, k);
      l1x_private(m, ct) += l1xi(m, k);
    }
  }

#pragma omp critical
{
  l1m += l1m_private;
  l1x += l1x_private;
}
}

// --- normalize by K and scale l1m by mu*(1-mu) ---
#pragma omp parallel for
for (int m = 0; m < M; m++) {
  double scale = mu[m]*(1.0-mu[m]);
  for (int ct = 0; ct < ncats; ct++) {
    l1m(m, ct) = (l1m(m, ct)/K) * scale;
    l1x(m, ct) /= K;
  }
}

// --- build l1s list ---
List l1s(d + 2);
arma::mat grid1(M, d+1, fill::ones);
grid1.cols(0, d-1) = grid;

for (int j = 0; j < d+1; j++) {
  arma::mat tmp = l1m.each_col() % grid1.col(j); // element-wise multiply with broadcasting
  l1s[j] = tmp;
}

l1s[d+1] = l1x;

return l1s;
}

/*** R
# quick test of compute_l1_cpp
# set.seed(1)
# M <- 4
# d <- 2
# ncats <- 3
# K <- 8
# cut_score <- c(.4,.6)
#
# grid <- matrix(runif(M*d), nrow=M, ncol=d)
# mu <- runif(M, 0.2, 0.8)
# nu <- 5
# beta_grid <- seq(0.05, 0.95, length.out=K)
# ind_cat <- as.numeric(cut(beta_grid,breaks = c(0,cut_score,1),labels = 1:ncats))
#
# res_cpp <- compute_l1_cpp(grid, nu, mu, beta_grid, ncats, ind_cat)
#
# # --- R reference implementation ---
# pmat <- t(outer(beta_grid, nu*mu-1, FUN="^") *
#             outer(1-beta_grid, nu*(1-mu)-1, FUN="^")) / beta(nu*mu, nu*(1-mu))
#
# l1mu <- nu*pmat * t(outer(log(beta_grid/(1-beta_grid)),
#                           digamma(nu*mu)-digamma(nu*(1-mu)), FUN="-"))
#
# l1xi <- nu*(digamma(nu) +
#               mu*t(outer(log(beta_grid), digamma(nu*mu), FUN="-")) +
#               (1-mu)*t(outer(log(1-beta_grid), digamma(nu*(1-mu)), FUN="-"))
# )*pmat
#
# l1m <- matrix(ncol=ncats, nrow=nrow(grid))
# l1x <- matrix(ncol=ncats, nrow=nrow(grid))
# for (ct in 1:ncats) {
#   l1m[,ct] <- rowSums(l1mu[,ind_cat==ct]) / length(beta_grid)
#   l1x[,ct] <- rowSums(l1xi[,ind_cat==ct]) / length(beta_grid)
# }
# l1m <- mu*(1-mu)*l1m
#
# l1s_ref <- list()
# for (i in 1:(d+1)) {
#   l1s_ref[[i]] <- sweep(l1m, 1, cbind(grid,1)[,i], FUN="*")
# }
# l1s_ref[[d+2]] <- l1x
#
# # --- compare results ---
# for (i in seq_along(res_cpp)) {
#   print(round(res_cpp[[i]] - l1s_ref[[i]], 15))
# }
*/


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List Grad_IM_cpp(List l1s,            // list of NumericMatrix, each M x ncats
                  NumericMatrix e_resp, // M x ncats
                  NumericMatrix p0,     // M x ncats
                  NumericVector f,      // length M
                  bool calculate = false) {

  int npar  = l1s.size();
  int M     = e_resp.nrow();
  int ncats = e_resp.ncol();
  int P     = M * ncats;   // flattened length
  const double eps = 1e-12;

  arma::mat A(npar, P, arma::fill::zeros);
  arma::vec Grad(npar, arma::fill::zeros);

  // --- precompute weights: sqrt(f[m] / p0(m,k)) ---
  arma::vec w(P);
  {
    int idx = 0;
    for (int k = 0; k < ncats; k++) {
      for (int m = 0; m < M; m++) {
        double p0val = p0(m, k);
        if (!std::isfinite(p0val) || p0val < eps) {
          p0val = eps; // clamp
        }
        w[idx++] = std::sqrt(f[m] / p0val);
      }
    }
  }

  // --- fill A (weighted l1s) and compute Grad ---
  for (int r = 0; r < npar; r++) {
    NumericMatrix lr = l1s[r];

    double g = 0.0;
    int idx = 0;
    for (int k = 0; k < ncats; k++) {
      for (int m = 0; m < M; m++) {
        double val = lr(m, k);
        double p0val = p0(m, k);
        if (!std::isfinite(p0val) || p0val < eps) {
          p0val = eps; // clamp
        }
        g += val * e_resp(m, k) / p0val;
        A(r, idx) = val * w[idx];
        idx++;
      }
    }
    Grad[r] = g;
  }

  // --- compute IM via BLAS: IM = A * A.t() ---
  arma::mat IM = A * A.t();

  if (calculate){
    arma::mat IMm(npar, npar, arma::fill::zeros);

    for (int r = 0; r < npar; r++) {
      for (int l = 0; l <= r; l++) {
        NumericMatrix lr1 = l1s[r];
        NumericMatrix lr2 = l1s[l];

        double g = 0.0;
        for (int m = 0; m < M; m++) {
          for (int k = 0; k < ncats; k++) {
            double p0val = p0(m, k);
            if (!std::isfinite(p0val) || p0val < eps) {
              continue; // clamp
            }
            g += lr2(m, k) * lr1(m, k) * e_resp(m, k) / p0val;
          }
        }
        IMm(r,l) = g;
        IMm(l,r) = g;
      }
    }
    return List::create(
      _["Grad"] = Grad,
      _["IM"]   = IM,
      _["IMm"]   = IMm
    );
  } else {
    return List::create(
      _["Grad"] = Grad,
      _["IM"]   = IM
    );
  }
}



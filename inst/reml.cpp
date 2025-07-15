#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// Helper: Convert Eigen::VectorXd to NumericVector for safe R return
NumericVector eigenVecToNumeric(const VectorXd& v) {
  return NumericVector(v.data(), v.data() + v.size());
}

// Helper: Convert Eigen::MatrixXd to NumericMatrix for safe R return
NumericMatrix eigenMatToNumeric(const MatrixXd& m) {
  NumericMatrix nm(m.rows(), m.cols());
  std::copy(m.data(), m.data() + m.size(), nm.begin());
  return nm;
}

// [[Rcpp::export]]
List reml_optimize_gamma_eigen(
    const Eigen::MappedSparseMatrix<double> &A,
    const Eigen::VectorXd &y,
    const Eigen::MappedSparseMatrix<double> &Q,
    const double sigma2 = 1.0,
    const double prior_alpha = 2.0,
    const double prior_beta = 0.1,
    const double log_gamma_init = 0.0,
    const double tol = 1e-6,
    const int maxit = 20,
    const bool verbose = true
) {
  const int n = A.cols();
  const int g = A.rows();

  // Precompute A'A and A'y (sparse)
  SparseMatrix<double> AtA = A.transpose() * A;
  VectorXd Aty = A.transpose() * y;

  double log_gamma = log_gamma_init;
  double gamma = std::exp(log_gamma);
  double final_nll = NA_REAL;

  for (int iter = 0; iter < maxit; ++iter) {
    if (verbose) Rcout << "\n[REML] Iter " << iter << " | log_gamma = " << log_gamma << ", gamma = " << gamma << "\n";

    // Build M = A'A + gamma * Q
    SparseMatrix<double> M = AtA + gamma * Q;

    if (verbose) Rcout << "[REML] Performing sparse Cholesky on M...\n";
    // Sparse Cholesky using Eigen's SimplicialLLT
    SimplicialLLT<SparseMatrix<double>> chol(M);
    if (chol.info() != Success) stop("Cholesky of M failed!");

    if (verbose) Rcout << "[REML] Solving M^{-1} * A^T...\n";
    // Dense solve for M^{-1} * A^T
    MatrixXd Minv_At = chol.solve(MatrixXd(A.transpose()));

    if (verbose) Rcout << "[REML] Computing Sigma_y = A * M^{-1} * A^T + sigma^2 I...\n";
    MatrixXd Sigma_y = A * Minv_At;
    Sigma_y.diagonal().array() += sigma2 + 1e-8;

    if (verbose) Rcout << "[REML] Cholesky decomposition of Sigma_y (dense)...\n";
    LLT<MatrixXd> Sigma_chol(Sigma_y.selfadjointView<Lower>());
    if (Sigma_chol.info() != Success) stop("Cholesky of Sigma_y failed!");

    VectorXd alpha = Sigma_chol.solve(y);
    double quad_form = alpha.squaredNorm();

    MatrixXd L = Sigma_chol.matrixL();
    double logdet = 2.0 * L.diagonal().array().log().sum();

    if (verbose) Rcout << "[REML] Computing gradient components...\n";
    VectorXd z_vec = A.transpose() * alpha;
    VectorXd m_inv_z = chol.solve(z_vec);
    VectorXd Q_m_inv_z = Q * m_inv_z;
    double quad_grad = -m_inv_z.dot(Q_m_inv_z);

    VectorXd qvec = Q_m_inv_z;
    VectorXd AQvec_tmp = chol.solve(qvec);
    VectorXd AQvec = AQvec_tmp;
    VectorXd S_inv_AQvec = Sigma_chol.solve(A * AQvec);
    double trace_term_approx = (A * AQvec).dot(S_inv_AQvec);

    // Add Gamma prior penalty and gradient (Gamma prior on gamma)
    double log_prior = (prior_alpha - 1.0) * log_gamma - prior_beta * gamma;
    double grad_prior = (prior_alpha - 1.0) - prior_beta * gamma;

    double grad_log_gamma = gamma * (trace_term_approx + quad_grad) - grad_prior;
    double nll = g * std::log(quad_form) + logdet - log_prior;
    final_nll = nll;

    if (verbose) {
      Rcout << "[REML] NLL = " << nll << ", grad_log_gamma = " << grad_log_gamma << "\n";
    }

    double step = -0.1 * grad_log_gamma;
    if (std::abs(step) < tol) {
      if (verbose) Rcout << "[REML] Converged.\n";
      break;
    }
    step = std::max(std::min(step, 1.0), -1.0);
    log_gamma += step;
    log_gamma = std::max(std::min(log_gamma, 20.0), -20.0);
    gamma = std::exp(log_gamma);
  }

  return List::create(
    Named("log_gamma") = log_gamma,
    Named("gamma") = gamma,
    Named("nll") = final_nll
  );
}

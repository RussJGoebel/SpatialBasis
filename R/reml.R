#' REML Estimation for General Prior
#'
#' Fits hyperparameters by maximizing the REML objective for any prior.
#'
#' @param X Design matrix (dense or sparse).
#' @param y Response vector.
#' @param prior_template A prior object (e.g., from make_sar_prior or make_matern_fourier_prior).
#' @param param_names Names of the hyperparameters to optimize.
#' @param param_start Named numeric vector of starting values on the natural scale.
#' @param control List of control options passed to optim().
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A list with components:
#' \describe{
#'   \item{opt}{The optim() result object.}
#'   \item{prior}{The prior object updated with estimated hyperparameters.}
#' }
#' @export
reml_estimate_general <- function(
    X,
    y,
    prior_template,
    param_names,
    param_start,
    control = list(),
    verbose = TRUE
) {
  stopifnot(length(param_names) == length(param_start))
  stopifnot(is.numeric(param_start))
  stopifnot(all(param_start > 0))
  stopifnot(all(param_names %in% names(prior_template$hyperparameters)))

  log_start <- log(param_start)

  if (verbose) message("Precomputing XtX and rhs...")
  XtX <- Matrix::crossprod(X)
  rhs <- Matrix::crossprod(X, y)

  # Objective function
  neg_reml <- function(log_params) {
    hyperparams <- as.list(exp(log_params))
    names(hyperparams) <- param_names

    if (verbose) {
      message("\n=== Evaluating objective ===")
      print(hyperparams)
    }

    prior <- prior_template
    prior$hyperparameters[names(hyperparams)] <- hyperparams

    Q <- compute_precision(prior)
    Q_post <- XtX + Q

    beta_hat <- solve_system(Q_post, rhs)

    r <- y - as.numeric(X %*% beta_hat)
    quad_form <- sum(r * r)

    chol_Q_post <- chol(Q_post)
    log_det <- 2 * sum(log(diag(chol_Q_post)))

    obj <- 0.5 * (log_det + quad_form)

    if (verbose) message("Objective value: ", signif(obj, 6))
    obj
  }

  # Gradient function
  grad_reml <- function(log_params) {
    hyperparams <- as.list(exp(log_params))
    names(hyperparams) <- param_names

    if (verbose) {
      message("\n=== Computing gradient ===")
      print(hyperparams)
    }

    prior <- prior_template
    prior$hyperparameters[names(hyperparams)] <- hyperparams

    Q <- compute_precision(prior)
    dQ_list <- compute_gradient(prior)
    Q_post <- XtX + Q

    beta_hat <- solve_system(Q_post, rhs)

    chol_Q_post <- chol(Q_post)
    inv_Q_post <- chol2inv(chol_Q_post)

    grad <- numeric(length(param_names))
    names(grad) <- param_names

    for (p in param_names) {
      if (verbose) message("Gradient wrt ", p, "...")
      dQ <- dQ_list[[p]]

      tr_term <- sum(diag(inv_Q_post %*% dQ))
      v_term <- crossprod(beta_hat, dQ %*% beta_hat)

      grad[p] <- 0.5 * (tr_term - v_term)
    }

    grad <- grad * exp(log_params)

    if (verbose) {
      message("Gradient vector:")
      print(grad)
    }

    grad
  }

  if (verbose) message("Starting optimization...")
  opt_result <- optim(
    par = log_start,
    fn = neg_reml,
    gr = grad_reml,
    method = "L-BFGS-B",
    control = control
  )

  est_params <- as.list(exp(opt_result$par))
  names(est_params) <- param_names

  prior_final <- prior_template
  prior_final$hyperparameters[names(est_params)] <- est_params

  list(
    opt = opt_result,
    prior = prior_final
  )
}

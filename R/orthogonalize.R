#' Orthogonalize one block of columns against the others (row-space projection)
#'
#' @param X [n x p] full design matrix
#' @param column_map named list: each entry is integer indices for a block
#' @param target_label which block to orthogonalize
#' @return list(X_new, column_map_new, A, Proj_other) where:
#'   A = (X_other' X_other)^{-1} X_other' X_target  [k x m], reuse for prediction
#'   Proj_other = (X_other' X_other)^{-1} X_other'  [k x n] (left pseudoinverse)
#' @export
orthogonalize_matrices <- function(X, column_map, target_label) {
  stopifnot(is.matrix(X) || inherits(X, "Matrix"))
  stopifnot(target_label %in% names(column_map))

  target_cols <- column_map[[target_label]]
  other_cols  <- unlist(column_map[names(column_map) != target_label])

  X_target <- X[, target_cols, drop = FALSE]

  if (length(other_cols) == 0L) {
    return(list(
      X_new = X,
      column_map_new = column_map,
      A = Matrix::Matrix(0, nrow = 0, ncol = ncol(X_target),
                         dimnames = list(character(0), colnames(X_target))),
      other_colnames = character(0)
    ))
  }

  X_other <- X[, other_cols, drop = FALSE]

  # Compute A without forming projector
  XtX_other <- Matrix::crossprod(X_other)
  chol_XtX <- Matrix::Cholesky(XtX_other, LDL = FALSE)  # factorization object
  RHS       <- Matrix::crossprod(X_other, X_target)

  # Sparse-safe Cholesky solve
  A <- Matrix::solve(chol_XtX, RHS, system = "A")
  #A <- Matrix::solve(chol_XtX, Matrix::solve(t(chol_XtX), RHS))

  # Residualize
  X_target_proj <- X_target - X_other %*% A

  # Rebuild without touching untouched cols
  all_idx <- seq_len(ncol(X))
  before_idx <- all_idx[all_idx < min(target_cols)]
  after_idx  <- all_idx[all_idx > max(target_cols)]

  X_new <- do.call(cbind, list(
    X[, before_idx, drop = FALSE],
    X_target_proj,
    X[, after_idx, drop = FALSE]
  ))

  list(
    X_new = X_new,
    column_map_new = column_map,
    A = A,
    other_colnames = colnames(X_other)
  )
}

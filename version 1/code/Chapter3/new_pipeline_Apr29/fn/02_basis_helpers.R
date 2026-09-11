build_trapezoid_weights <- function(grid) {
  n_grid <- length(grid)
  if (n_grid < 2) {
    stop("Grid must contain at least 2 points.")
  }

  step_size <- diff(grid)
  if (!all(abs(step_size - step_size[1]) < 1e-12)) {
    stop("Current implementation expects an equally spaced grid.")
  }

  delta <- step_size[1]
  weights <- rep(delta, n_grid)
  weights[c(1, n_grid)] <- delta / 2
  weights
}

build_raw_cubic_spline_basis <- function(grid, basis_k) {
  if (basis_k < 3) {
    stop("basis_k must be at least 3 for a cubic spline basis.")
  }

  if (basis_k == 3) {
    raw_basis <- splines::ns(
      x = grid,
      df = basis_k,
      intercept = TRUE,
      Boundary.knots = range(grid)
    )
    raw_basis <- as.matrix(raw_basis)
    if (ncol(raw_basis) != basis_k) {
      stop("Constructed natural spline basis has ", ncol(raw_basis), " columns; expected ", basis_k, ".")
    }
    return(raw_basis)
  }

  internal_knot_n <- basis_k - 4
  if (internal_knot_n > 0) {
    knot_probs <- seq_len(internal_knot_n) / (internal_knot_n + 1)
    knots <- as.numeric(stats::quantile(grid, probs = knot_probs))
  } else {
    knots <- NULL
  }

  raw_basis <- splines::bs(
    x = grid,
    knots = knots,
    degree = 3,
    intercept = TRUE,
    Boundary.knots = range(grid)
  )
  raw_basis <- as.matrix(raw_basis)
  if (ncol(raw_basis) != basis_k) {
    stop("Constructed spline basis has ", ncol(raw_basis), " columns; expected ", basis_k, ".")
  }
  raw_basis
}

orthonormalize_basis <- function(raw_basis, weights) {
  weighted_basis <- raw_basis * sqrt(weights)
  qr_fit <- qr(weighted_basis)
  q_weighted <- qr.Q(qr_fit)
  q_weighted / sqrt(weights)
}

build_spline_basis <- function(grid, basis_k = 5) {
  weights <- build_trapezoid_weights(grid)
  raw_basis <- build_raw_cubic_spline_basis(grid, basis_k = basis_k)
  basis_matrix <- orthonormalize_basis(raw_basis, weights = weights)
  basis_names <- paste0("Basis", seq_len(ncol(basis_matrix)))

  colnames(raw_basis) <- paste0("RawBasis", seq_len(ncol(raw_basis)))
  colnames(basis_matrix) <- basis_names

  list(
    grid = grid,
    weights = weights,
    raw_basis = raw_basis,
    basis_matrix = basis_matrix,
    basis_names = basis_names
  )
}

project_curves_to_basis <- function(curve_matrix, basis_obj) {
  weighted_basis <- sweep(basis_obj$basis_matrix, 1L, basis_obj$weights, "*")
  basis_scores <- curve_matrix %*% weighted_basis
  basis_scores <- as.matrix(basis_scores)
  colnames(basis_scores) <- basis_obj$basis_names
  basis_scores
}

verify_basis_orthonormality <- function(basis_obj) {
  weighted_basis <- sweep(basis_obj$basis_matrix, 1L, sqrt(basis_obj$weights), "*")
  cross_mat <- crossprod(weighted_basis)
  deviation <- max(abs(cross_mat - diag(ncol(cross_mat))))
  list(cross_mat = cross_mat, max_abs_deviation = deviation)
}

save_basis_plot <- function(basis_obj, out_file) {
  grDevices::png(filename = out_file, width = 1600, height = 900, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::matplot(
    basis_obj$grid,
    basis_obj$basis_matrix,
    type = "l",
    lty = 1,
    lwd = 2,
    xlab = "s",
    ylab = "Basis value",
    main = "Orthonormalized cubic spline basis"
  )
  graphics::legend(
    "topright",
    legend = basis_obj$basis_names,
    col = seq_along(basis_obj$basis_names),
    lty = 1,
    lwd = 2,
    bty = "n"
  )
}

make_basis_grid_table <- function(basis_obj) {
  data.frame(
    grid = basis_obj$grid,
    weight = basis_obj$weights,
    basis_obj$raw_basis,
    basis_obj$basis_matrix,
    check.names = FALSE
  )
}

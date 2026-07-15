spectral_distance <- function(
  Omega_hat_k0,
  epsilon = 1e-6,
  laplacian_target_type = c("identity", "diag_Omega_hat"),
  adj_threshold = 1e-4,
  laplacian_norm_type = c("symmetric", "unsymmetric")
) {
  laplacian_target_type <- match.arg(laplacian_target_type)
  laplacian_norm_type <- match.arg(laplacian_norm_type)
  Omega_hat_k0 <- as.matrix(Omega_hat_k0)
  if (!is.numeric(Omega_hat_k0) || !nrow(Omega_hat_k0) ||
      nrow(Omega_hat_k0) != ncol(Omega_hat_k0) ||
      any(!is.finite(Omega_hat_k0))) {
    stop("Omega_hat_k0 must be a finite, nonempty square matrix.",
         call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L ||
      !is.finite(epsilon) || epsilon <= 0) {
    stop("epsilon must be finite and strictly positive.", call. = FALSE)
  }
  if (!is.numeric(adj_threshold) || length(adj_threshold) != 1L ||
      !is.finite(adj_threshold) || adj_threshold < 0) {
    stop("adj_threshold must be finite and non-negative.", call. = FALSE)
  }

  # Threshold the precision matrix and form the requested graph Laplacian.
  get_laplacian <- function(matrix, threshold, normalization) {
    adjacency <- abs(matrix) > threshold
    diag(adjacency) <- FALSE
    graph <- igraph::graph_from_adjacency_matrix(
      adjacency,
      mode = "undirected"
    )
    igraph::laplacian_matrix(
      graph,
      normalization = if (identical(normalization, "symmetric")) {
        "symmetric"
      } else {
        "unnormalized"
      },
      sparse = FALSE
    )
  }

  laplacian <- get_laplacian(
    Omega_hat_k0,
    threshold = adj_threshold,
    normalization = laplacian_norm_type
  )
  spectrum <- sort(eigen(
    laplacian,
    symmetric = TRUE,
    only.values = TRUE
  )$values)

  # Both public targets represent an empty adjacency graph: diagonal entries
  # do not create edges in either construction.
  target <- diag(0, nrow(Omega_hat_k0))
  target_laplacian <- get_laplacian(
    target,
    threshold = adj_threshold,
    normalization = laplacian_norm_type
  )
  target_spectrum <- sort(eigen(
    target_laplacian,
    symmetric = TRUE,
    only.values = TRUE
  )$values)
  if (length(spectrum) != length(target_spectrum)) {
    stop("The fitted and target Laplacian spectra have different dimensions.",
         call. = FALSE)
  }

  distance <- sqrt(sum((spectrum - target_spectrum)^2))
  weight <- 1 / max(distance, epsilon)
  matrix(weight, nrow(Omega_hat_k0), ncol(Omega_hat_k0))
}

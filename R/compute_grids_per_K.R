compute_grids_per_K <- function(
  X, nbcluster, P_method,
  distance_method     = "Euclidean",
  eps_w0              = sqrt(.Machine$double.eps),
  L                   = 5,
  frac_min            = 0.05,
  init_method         = c(
    "kmeans", "mclust_hc", "mclust_em", "random", "user",
    "previous_fit", "deterministic_multistart", "hc"
  ),
  n.start             = 250,
  lambda_omega_0      = 50,
  penalize_diag       = FALSE,
  laplacian_target_type = c("identity", "diag_Omega_hat"),
  adj_threshold       = 1e-4,
  laplacian_norm_type = c("symmetric", "unsymmetric"),
  initializers        = NULL
) {
  stopifnot(L >= 2, frac_min > 0, frac_min < 1)
  if (length(eps_w0) != 1L || !is.numeric(eps_w0) ||
      !is.finite(eps_w0) || eps_w0 <= 0) {
    stop("eps_w0 must be a finite positive number.", call. = FALSE)
  }
  if (!is.logical(penalize_diag) || length(penalize_diag) != 1L ||
      is.na(penalize_diag)) {
    stop("penalize_diag must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(penalize_diag)) {
    stop("penalize_diag=TRUE is unsupported for ranking grids.", call. = FALSE)
  }
  X <- as.matrix(X)
  if (!is.numeric(X) || !length(X) || any(!is.finite(X))) {
    stop("X must be a non-empty finite numeric matrix.", call. = FALSE)
  }
  nbcluster <- as.integer(nbcluster)
  if (!length(nbcluster) || anyNA(nbcluster) || any(nbcluster < 1L) ||
      any(nbcluster > nrow(X))) {
    stop("nbcluster must contain positive integers no larger than nrow(X).", call. = FALSE)
  }
  if (!is.null(initializers) &&
      (!is.list(initializers) || length(initializers) != length(nbcluster))) {
    stop(
      "initializers must be NULL or a list aligned one-to-one with nbcluster.",
      call. = FALSE
    )
  }
  P_method <- match.arg(
    P_method,
    c("common", "weighted_by_W0", "weighted_by_dist_to_I",
      "weighted_by_dist_to_diag_W0", "laplacian_spectral")
  )
  init_method <- match.arg(init_method)
  laplacian_target_type <- match.arg(laplacian_target_type)
  laplacian_norm_type  <- match.arg(laplacian_norm_type)

  covariance_distance <- function(first, second, method) {
    if (!is.character(method) || length(method) != 1L || is.na(method)) {
      stop("distance_method must be one non-missing character string.",
           call. = FALSE)
    }
    if (identical(method, "Euclidean")) {
      # Euclidean covariance distance equals the Frobenius norm, so the
      # default does not require an optional dependency.
      return(sqrt(sum((first - second)^2)))
    }
    if (!requireNamespace("shapes", quietly = TRUE)) {
      stop(
        sprintf(
          paste0(
            "distance_method='%s' requires the optional package 'shapes'; ",
            "install it with install.packages('shapes')."
          ),
          method
        ),
        call. = FALSE
      )
    }
    shapes::distcov(first, second, method = method)
  }

  build_Pk <- function(Omega_cube) {
    dims <- dim(Omega_cube)
    p <- dims[1]; K <- dims[3]

    Pk <- switch(P_method,
      common = array(1, dim = dims),

      weighted_by_W0 = 1 / (eps_w0 + abs(Omega_cube)),

      weighted_by_dist_to_I = {
        d <- apply(Omega_cube, 3, function(Om)
          covariance_distance(Om, diag(p), method = distance_method)
        )
        d[d < eps_w0] <- eps_w0
        array(rep(1/d, each = p*p), dim = dims)
      },

      weighted_by_dist_to_diag_W0 = {
        d <- apply(Omega_cube, 3, function(Om) {
          D <- diag(diag(Om))
          covariance_distance(Om, D, method = distance_method)
        })
        d[d < eps_w0] <- eps_w0
        array(rep(1/d, each = p*p), dim = dims)
      },

      laplacian_spectral = {
        Pk_out <- array(0, dim = dims)
        for (k in seq_len(dims[3])) {
          Pk_out[,,k] <- spectral_distance(
            Omega_hat_k0         = Omega_cube[,,k],
            epsilon              = eps_w0,
            laplacian_target_type = laplacian_target_type,
            adj_threshold        = adj_threshold,
            laplacian_norm_type  = laplacian_norm_type
          )
        }
        Pk_out
      },

      stop("unknown P_method: ", P_method)
    )

    for (k in seq_len(dims[3])) diag(Pk[,,k]) <- 0
    if (any(!is.finite(Pk)) || any(Pk < 0)) {
      stop(
        "Computed penalty weights must be finite and non-negative.",
        call. = FALSE
      )
    }
    Pk
  }

  out <- vector("list", length(nbcluster))
  names(out) <- paste0("K", nbcluster)

  for (i in seq_along(nbcluster)) {
    K <- nbcluster[i]
    init <- if (is.null(initializers)) {
      InitParameter(
        data = X,
        nbClust = K,
        init = init_method,
        n.start = n.start,
        lambda_omega_0 = lambda_omega_0
      )
    } else {
      initializers[[i]]
    }

    expected_names <- c("prop", "Mu", "SigmaCube", "OmegaCube", "Z")
    if (!is.list(init) || !all(expected_names %in% names(init)) ||
        length(init$prop) != K || !identical(dim(init$Mu), c(ncol(X), K)) ||
        !identical(dim(init$SigmaCube), c(ncol(X), ncol(X), K)) ||
        !identical(dim(init$OmegaCube), c(ncol(X), ncol(X), K)) ||
        !identical(dim(init$Z), c(nrow(X), K)) ||
        any(!is.finite(unlist(init[expected_names], use.names = FALSE)))) {
      stop(sprintf("InitParameter returned an invalid initializer for K=%d.", K),
           call. = FALSE)
    }
    initializer_report <- init[setdiff(names(init), expected_names)]
    init <- init[expected_names]

    Z0     <- init$Z
    nk0    <- colSums(Z0)
    Sigma0 <- init$SigmaCube
    Omega0 <- init$OmegaCube

    Pk_cube <- build_Pk(Omega0)

    weighted_sums <- t(Z0) %*% X
    lambda_component_bounds <- vapply(seq_len(K), function(k) {
      max(abs(Omega0[, , k] %*% weighted_sums[k, ]))
    }, numeric(1))
    lam_max <- max(lambda_component_bounds)
    if (!is.finite(lam_max)) {
      stop(sprintf("Non-finite lambda path bound for K=%d.", K), call. = FALSE)
    }
    lambda_path_anchor <- if (lam_max > 0) lam_max else .Machine$double.eps

    rho_max <- 0
    off_diagonal <- row(Pk_cube[, , 1L]) != col(Pk_cube[, , 1L])
    for (k in seq_len(K)) {
      Sk  <- Sigma0[,,k]
      weights <- Pk_cube[, , k][off_diagonal]
      if (any(weights <= 0)) {
        stop(
          sprintf("Off-diagonal penalty weights must be positive for K=%d.", K),
          call. = FALSE
        )
      }
      tmp <- nk0[k] * abs(Sk[off_diagonal]) / weights
      if (any(!is.finite(tmp))) {
        stop(sprintf("Non-finite rho path bound for K=%d.", K), call. = FALSE)
      }
      rho_max <- max(rho_max, tmp)
    }
    if (!is.finite(rho_max)) {
      stop(sprintf("Non-finite rho path bound for K=%d.", K), call. = FALSE)
    }
    if (rho_max == 0) rho_max <- 1

    geo <- function(maxv) {
      minv <- maxv * frac_min
      maxv * (minv / maxv) ^ (seq(0, L - 1) / (L - 1))
    }

    out[[i]] <- list(
      K         = K,
      lambda_mu = geo(lambda_path_anchor),
      lambda_max_kkt = lam_max,
      lambda_component_bounds = lambda_component_bounds,
      rho       = geo(rho_max),
      initializer = init,
      initializer_report = initializer_report,
      Pk_cube   = Pk_cube,
      grid_order = "rho_outer_lambda_inner"
    )
  }

  attr(out, "grid_order") <- "rho_outer_lambda_inner"
  attr(out, "state_metadata") <- "stored_initializer_and_weights"
  out
}

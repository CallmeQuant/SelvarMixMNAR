SortvarClust <- function(x,
                         nbcluster,
                         type="lasso", 
                         lambda=seq(20, 100, by = 10),
                         rho=seq(0.1, 1, length=5),
                         group_shrinkage_method = "common", 
                         distance_method = "Euclidean",
                         lambda_omega_0 = 50,
                         epsilon_weighted_by_W0 = sqrt(.Machine$double.eps),
                         laplacian_target_type = c("identity", "diag_Omega_hat"),
                         adj_threshold = 1e-4,
                         laplacian_norm_type = c("symmetric", "unsymmetric"),
                          penalize_diag = FALSE,
                          initialize = "hc",
                          nbcores=min(2, parallel::detectCores(logical = FALSE)), 
                          n.start = 250,
                          scale = FALSE,
                           penalty_grids = NULL,
                           warm_start = c("none", "inner", "outer", "both"),
                           verbose = FALSE,
                           min_scorable_fraction = 0.5,
                           ...,
                           init_control = list())
{
  # Validate the complete-data ranking problem and numerical controls.
  if(missing(x)){ stop("x is required.", call. = FALSE) }
  if(!is.matrix(x) && !is.data.frame(x)) stop(sQuote("x"), " must be a matrix or data frame")
  x <- data.matrix(x) 
  if(any(!is.finite(x))) stop("Input data 'x' contains non-finite values (NA, NaN, Inf). Consider imputation or filtering.")
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("type must be one non-missing character token.", call. = FALSE)
  }
  if (identical(tolower(type), "likelihood")) {
    stop(structure(
      list(
        message = paste0(
          "type='likelihood' is unavailable because no likelihood-based ",
          "variable-ranking estimator is implemented. Use type='lasso'."
        ),
        call = NULL
      ),
      class = c("selvarmix_ranking_method_unavailable", "error", "condition")
    ))
  }
  if (!identical(tolower(type), "lasso")) {
    stop("type must be exactly 'lasso'.", call. = FALSE)
  }

  if(missing(nbcluster)){ stop("nbcluster is required.", call. = FALSE) }
  if(any(!is.wholenumber(nbcluster)) || any(nbcluster < 1)) stop("nbcluster must contain positive integers.", call. = FALSE)

  if (is.null(penalty_grids)) {
    if(!is.vector(lambda) || length(lambda) < 1) stop(sQuote("lambda"), " must be a vector with length >= 1")
    if (any(!is.finite(lambda)) || any(lambda < 0)) stop("lambda must be finite and non-negative.", call. = FALSE)

    if(!is.vector(rho) || length(rho) < 1) stop(sQuote("rho"), " must be a vector with length >= 1")
    if(any(!is.finite(rho)) || any(rho < 0)) stop("rho must be finite and non-negative.", call. = FALSE)
  } else if (!is.list(penalty_grids) || length(penalty_grids) != length(nbcluster)) {
    stop("penalty_grids must be a list aligned one-to-one with nbcluster.")
  }
  if (!is.null(penalty_grids) && length(init_control)) {
    stop(
      "Supply either penalty_grids with pinned initializers or init_control, not both.",
      call. = FALSE
    )
  }

  if(!is.wholenumber(nbcores) || nbcores < 1) stop(sQuote("nbcores"), " must be an integer > 0")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }

  group_shrinkage_method <- tryCatch(
      match.arg(group_shrinkage_method, c("common", "weighted_by_W0", "weighted_by_dist_to_I", "weighted_by_dist_to_diag_W0", "laplacian_spectral")),
      error = function(e) stop("Invalid 'group_shrinkage_method'. Choose from 'common', 'weighted_by_W0', 'weighted_by_dist_to_I', 'weighted_by_dist_to_diag_W0', 'laplacian_spectral'.")
  )

  laplacian_target_type <- match.arg(laplacian_target_type)
  laplacian_norm_type  <- match.arg(laplacian_norm_type)
  warm_start <- match.arg(warm_start)
  # Validate initialization.  The richer user/previous-fit modes enter through
  # a precomputed per-K initializer because they require an explicit state.
  initialize_methods <- c(
    "mclust_hc", "mclust_em", "kmeans", "random", "user",
    "previous_fit", "deterministic_multistart", "hc"
  )
  initialize <- tryCatch(
      match.arg(initialize, initialize_methods),
      error = function(e) stop(
        "Invalid 'initialize'. Choose one of: ",
        paste(initialize_methods, collapse = ", "), ".",
        call. = FALSE
      )
  )
  # Direct calls may request standardization; the main fitter supplies its
  # already transformed ranking matrix.
  if (!is.logical(scale) || length(scale) != 1L || is.na(scale)) {
    stop("scale must be TRUE or FALSE.", call. = FALSE)
  }
  column_sd <- apply(x, 2L, stats::sd)
  magnitude <- pmax(1, apply(abs(x), 2L, max))
  degenerate <- !is.finite(column_sd) |
    column_sd <= .Machine$double.eps * magnitude
  if (any(degenerate)) {
    stop(
      "Input data contain non-varying column(s): ",
      paste(which(degenerate), collapse = ", "),
      "; filter them before penalized ranking.",
      call. = FALSE
    )
  }
  if (scale){
    x <- scale(x, center = TRUE, scale = TRUE)
  }
  if (any(!is.finite(x))) {
      stop("Scaling resulted in non-finite values. Check for constant columns in the input data.")
  }
  p <- as.integer(ncol(x))

  # Aggregate component-mean activity over each validated penalty grid.
  OrderVariable <- matrix(NA, nrow = length(nbcluster), ncol = p)
  rownames(OrderVariable) <- paste("K", nbcluster, sep="=")
  colnames(OrderVariable) <- colnames(x) # Assign variable names if available

  if(tolower(type) == "lasso")
  {
    initializer_adapters <- NULL
    if (is.null(penalty_grids)) {
      initializer_adapters <- .build_initializer_registry(
        data = x,
        nbcluster = nbcluster,
        init_control = init_control,
        legacy_method = initialize,
        legacy_n_start = n.start,
        legacy_lambda_omega_0 = lambda_omega_0
      )
      penalty_grids <- lapply(seq_along(nbcluster), function(index) {
        list(
          K = as.integer(nbcluster[index]),
          lambda_mu = as.numeric(lambda),
          rho = as.numeric(rho),
          initializer = .initializer_state_projection(
            initializer_adapters[[index]]
          ),
          Pk_cube = NULL,
          grid_order = "rho_outer_lambda_inner"
        )
      })
      names(penalty_grids) <- names(initializer_adapters)
      attr(penalty_grids, "initializer_adapters") <- initializer_adapters
    } else {
      grid_initializers <- lapply(seq_along(nbcluster), function(index) {
        entry <- penalty_grids[[index]]
        if (!is.list(entry) || is.null(entry$initializer)) {
          stop(
            sprintf("penalty_grids[[%d]] lacks a pinned initializer.", index),
            call. = FALSE
          )
        }
        .initializer_state_projection(entry$initializer)
      })
      names(grid_initializers) <- paste0("K", nbcluster)
      supplied_adapters <- attr(
        penalty_grids, "initializer_adapters", exact = TRUE
      )
      if (is.null(supplied_adapters)) {
        per_k <- lapply(grid_initializers, function(initializer) {
          list(method = "user", user_init = initializer)
        })
        initializer_adapters <- .build_initializer_registry(
          data = x,
          nbcluster = nbcluster,
          init_control = list(method = "user", per_k = per_k),
          legacy_method = "user",
          legacy_n_start = n.start,
          legacy_lambda_omega_0 = lambda_omega_0
        )
      } else {
        valid_adapter_registry <- is.list(supplied_adapters) &&
          length(supplied_adapters) == length(nbcluster) &&
          all(vapply(
            supplied_adapters,
            inherits,
            logical(1),
            what = "selvarmix_numeric_initialization"
          ))
        if (!valid_adapter_registry) {
          stop(
            paste0(
              "The supplied-grid initializer_adapters attribute must contain ",
              "one validated numeric adapter per candidate K."
            ),
            call. = FALSE
          )
        }
        for (index in seq_along(nbcluster)) {
          adapter_state <- .initializer_state_projection(
            supplied_adapters[[index]]
          )
          if (!identical(adapter_state, grid_initializers[[index]])) {
            stop(
              sprintf(
                paste0(
                  "The reported initializer adapter for K=%d does not match ",
                  "the state pinned in penalty_grids."
                ),
                nbcluster[index]
              ),
              call. = FALSE
            )
          }
        }
        initializer_adapters <- supplied_adapters
      }
    }

    VarRole <- ClusteringEMGlassoWeighted(
        data = x,
        nbcluster = nbcluster,
        lambda = lambda,
        rho = rho,
        group_shrinkage_method = group_shrinkage_method, 
        distance_method = distance_method,           
        lambda_omega_0 = lambda_omega_0,           
        epsilon_weighted_by_W0 = epsilon_weighted_by_W0, 
        laplacian_target_type = laplacian_target_type,
        adj_threshold = adj_threshold,
        laplacian_norm_type = laplacian_norm_type,
        penalize_diag = penalize_diag,    
        initialize = initialize,
        nbcores = nbcores,
        n.start = n.start,
        penalty_grids = penalty_grids,
        warm_start = warm_start,
        verbose = verbose,
        min_scorable_fraction = min_scorable_fraction
        )

    used_grids <- attr(VarRole, "penalty_grids", exact = TRUE)
    expected_dim <- c(nrow(used_grids[[1L]]), p, length(nbcluster))
    if (!is.array(VarRole) || !identical(dim(VarRole), expected_dim)) {
        stop("ClusteringEMGlassoWeighted did not return the expected 3D activity array.")
    }
    if (!is.integer(VarRole)) {
        stop("ClusteringEMGlassoWeighted returned non-integer activity; no ranking was produced.")
    }

    MatrixScores <- matrix(0, nrow = length(nbcluster), ncol = p)
    scorable_grid_rows <- vector("list", length(nbcluster))
    names(scorable_grid_rows) <- paste0("K", nbcluster)
    for (k_idx in 1:length(nbcluster)) {
      # Failed fits remain all-NA rows and never contribute zero activity.
      activity_k <- matrix(
        VarRole[, , k_idx, drop = FALSE],
        nrow = dim(VarRole)[1L],
        ncol = p
      )
      row_observed <- rowSums(!is.na(activity_k))
      if (any(!(row_observed %in% c(0L, p)))) {
        stop("A failed grid fit returned partially missing activity; no ranking was produced.")
      }
      valid_rows <- which(row_observed == p)
      if (!length(valid_rows) ||
          any(!is.finite(activity_k[valid_rows, , drop = FALSE])) ||
          any(!(activity_k[valid_rows, , drop = FALSE] %in% 0:1))) {
        stop("No complete binary activity row is available for one candidate K; no ranking was produced.")
      }
      fit_status <- attr(VarRole, "fit_status", exact = TRUE)[[k_idx]]
      status_rows <- which(vapply(fit_status, function(status) {
        isTRUE(status$converged) && isTRUE(status$scorable)
      }, logical(1L)))
      if (!identical(valid_rows, status_rows)) {
        stop("Activity rows and explicit native fit statuses disagree; no ranking was produced.")
      }
      scorable_grid_rows[[k_idx]] <- valid_rows
      MatrixScores[k_idx, ] <- colSums(activity_k[valid_rows, , drop = FALSE])
      OrderVariable[k_idx, ] <- order(MatrixScores[k_idx, ], decreasing = TRUE)
    }

    for (attribute_name in c(
      "penalty_grids", "initializers", "fit_status", "fit_states",
      "objective_traces", "glasso_diagnostics", "path_metadata",
      "grid_failures",
      "fit_metadata"
    )) {
      attr(OrderVariable, attribute_name) <- attr(
        VarRole, attribute_name, exact = TRUE
      )
    }
    attr(OrderVariable, "initializer_adapters") <- initializer_adapters
    attr(OrderVariable, "scorable_grid_rows") <- scorable_grid_rows
    attr(OrderVariable, "ranking_scores") <- MatrixScores

  } else {
      stop("Internal ranking dispatch error.", call. = FALSE)
  }

  return(OrderVariable) # Return matrix of ranked variable indices
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

ClusteringEMGlassoWeighted <- function(data,
                                       nbcluster,
                                       lambda,
                                       rho,
                                       group_shrinkage_method = c(
                                         "common",
                                         "weighted_by_W0",
                                         "weighted_by_dist_to_I",
                                         "weighted_by_dist_to_diag_W0",
                                         "laplacian_spectral"
                                       ),
                                       distance_method = "Euclidean",
                                       lambda_omega_0 = 50,
                                       epsilon_weighted_by_W0 = sqrt(.Machine$double.eps),
                                       penalize_diag = FALSE,
                                       laplacian_target_type = c("identity", "diag_Omega_hat"),
                                       adj_threshold = 1e-4,
                                       laplacian_norm_type = c("symmetric", "unsymmetric"),
                                       initialize = c(
                                         "kmeans", "mclust_hc", "mclust_em",
                                         "random", "user", "previous_fit",
                                         "deterministic_multistart", "hc"
                                       ),
                                       nbcores = 1,
                                       n.start = 250,
                                       penalty_grids = NULL,
                                       warm_start = c("none", "inner", "outer", "both"),
                                       verbose = FALSE,
                                       min_scorable_fraction = 0.5)
{
  empty_grid <- function() {
    matrix(
      numeric(), nrow = 0L, ncol = 2L,
      dimnames = list(NULL, c("lambda", "rho"))
    )
  }

  stop_grid_fit <- function(message, K, failed_grid = NULL, n_grid = NULL,
                            failures = list()) {
    if (is.null(failed_grid)) failed_grid <- empty_grid()
    condition <- structure(
      list(
        message = message,
        call = NULL,
        nbcluster = K,
        failed_grid = failed_grid,
        n_failed = nrow(failed_grid),
        n_grid = n_grid,
        failures = failures
      ),
      class = c("selvarmix_grid_fit_error", "error", "condition")
    )
    stop(condition)
  }

  validate_initializer <- function(init, K, source) {
    expected_names <- c("prop", "Mu", "SigmaCube", "OmegaCube", "Z")
    if (!is.list(init) || length(init) != 5L) {
      stop(sprintf("%s must be a five-element initializer for K=%d.", source, K),
           call. = FALSE)
    }
    if (is.null(names(init))) names(init) <- expected_names
    if (!all(expected_names %in% names(init))) {
      stop(sprintf("%s has missing initializer fields for K=%d.", source, K),
           call. = FALSE)
    }
    init <- init[expected_names]
    valid <- is.numeric(init$prop) && length(init$prop) == K &&
      is.matrix(init$Mu) && identical(dim(init$Mu), c(p, K)) &&
      is.array(init$SigmaCube) &&
      identical(dim(init$SigmaCube), c(p, p, K)) &&
      is.array(init$OmegaCube) &&
      identical(dim(init$OmegaCube), c(p, p, K)) &&
      is.matrix(init$Z) && identical(dim(init$Z), c(n, K))
    if (!valid || any(!is.finite(unlist(init, use.names = FALSE))) ||
        any(init$prop <= 0) || abs(sum(init$prop) - 1) > 1e-8) {
      stop(sprintf("%s contains invalid initializer values for K=%d.", source, K),
           call. = FALSE)
    }
    init
  }

  validate_Pk <- function(Pk, K, source) {
    if (!is.array(Pk) || !identical(dim(Pk), c(p, p, K)) ||
        any(!is.finite(Pk)) || any(Pk < 0)) {
      stop(sprintf("%s contains invalid penalty weights for K=%d.", source, K),
           call. = FALSE)
    }
    diagonal_weights <- unlist(lapply(
      seq_len(K), function(component) diag(Pk[, , component])
    ), use.names = FALSE)
    if (any(diagonal_weights != 0)) {
      stop(
        sprintf(
          "%s has nonzero diagonal penalty weights for K=%d; diagonal precision penalization is unsupported.",
          source, K
        ),
        call. = FALSE
      )
    }
    symmetric <- vapply(seq_len(K), function(component) {
      value <- Pk[, , component]
      norm(value - t(value), type = "I") <=
        1e-12 * max(1, norm(value, type = "I"))
    }, logical(1))
    if (!all(symmetric)) {
      stop(sprintf("%s must be symmetric for every component of K=%d.",
                   source, K), call. = FALSE)
    }
    Pk
  }

  valid_covariance_precision_pair <- function(covariance, precision) {
    matrix_valid <- function(value) {
      scale_value <- max(1, norm(value, type = "I"))
      symmetric <- norm(value - t(value), type = "I") <= 1e-8 * scale_value
      positive_definite <- tryCatch({
        chol(0.5 * (value + t(value)))
        TRUE
      }, error = function(condition) FALSE)
      symmetric && positive_definite
    }
    if (!matrix_valid(covariance) || !matrix_valid(precision)) return(FALSE)
    identity <- diag(nrow(covariance))
    residual <- max(
      norm(covariance %*% precision - identity, type = "I"),
      norm(precision %*% covariance - identity, type = "I")
    )
    normalized_residual <- residual / (
      1 + norm(covariance, type = "I") * norm(precision, type = "I")
    )
    tolerance <- max(1e-10, 25 * sqrt(1 + nrow(covariance)) * 1e-7)
    is.finite(normalized_residual) && normalized_residual <= tolerance
  }

  valid_ranking_state <- function(state, K) {
    if (any(state$prop <= 0) || abs(sum(state$prop) - 1) > 1e-8 ||
        any(state$Z < 0) ||
        any(abs(rowSums(state$Z) - 1) > 1e-8) ||
        any(colSums(state$Z) <= 0)) {
      return(FALSE)
    }
    all(vapply(seq_len(K), function(component) {
      valid_covariance_precision_pair(
        state$SigmaCube[, , component],
        state$OmegaCube[, , component]
      )
    }, logical(1)))
  }

  make_grid <- function(lambda_values, rho_values, source) {
    if (!is.numeric(lambda_values) || !length(lambda_values) ||
        any(!is.finite(lambda_values)) || any(lambda_values < 0) ||
        !is.numeric(rho_values) || !length(rho_values) ||
        any(!is.finite(rho_values)) || any(rho_values < 0)) {
      stop(sprintf("%s must contain finite, non-negative lambda_mu and rho paths.",
                   source), call. = FALSE)
    }
    grid <- expand.grid(
      lambda = as.double(lambda_values),
      rho = as.double(rho_values),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    as.matrix(grid)
  }

  compute_Pk <- function(omega0_cube,
                         method = "common",
                         distance_method = "Euclidean",
                         laplacian_target_type_pk = "diag_Omega_hat",
                         adj_threshold_pk = 1e-4,
                         laplacian_norm_type_pk = "symmetric",
                         eps_w0 = sqrt(.Machine$double.eps),
                         penalize_diag_pk = FALSE) {
    dims <- dim(omega0_cube)
    p_local <- dims[1L]
    K_local <- dims[3L]
    covariance_distance <- function(first, second, requested_method) {
      if (!is.character(requested_method) ||
          length(requested_method) != 1L || is.na(requested_method)) {
        stop("distance_method must be one non-missing character string.",
             call. = FALSE)
      }
      if (identical(requested_method, "Euclidean")) {
        return(sqrt(sum((first - second)^2)))
      }
      if (!requireNamespace("shapes", quietly = TRUE)) {
        stop(
          sprintf(
            paste0(
              "distance_method='%s' requires the optional package 'shapes'; ",
              "install it with install.packages('shapes')."
            ),
            requested_method
          ),
          call. = FALSE
        )
      }
      shapes::distcov(first, second, method = requested_method)
    }

    out <- switch(
      method,
      common = array(1, dim = dims),
      weighted_by_W0 = 1 / (eps_w0 + abs(omega0_cube)),
      weighted_by_dist_to_I = {
        d <- apply(omega0_cube, 3L, function(omega) {
          covariance_distance(
            omega, diag(p_local), requested_method = distance_method
          )
        })
        d[d < eps_w0] <- eps_w0
        array(rep(1 / d, each = p_local * p_local), dim = dims)
      },
      weighted_by_dist_to_diag_W0 = {
        d <- apply(omega0_cube, 3L, function(omega) {
          covariance_distance(
            omega, diag(diag(omega)), requested_method = distance_method
          )
        })
        d[d < eps_w0] <- eps_w0
        array(rep(1 / d, each = p_local * p_local), dim = dims)
      },
      laplacian_spectral = {
        result <- array(0, dim = dims)
        for (component in seq_len(K_local)) {
          result[, , component] <- spectral_distance(
            omega0_cube[, , component],
            epsilon = eps_w0,
            laplacian_target_type = laplacian_target_type_pk,
            adj_threshold = adj_threshold_pk,
            laplacian_norm_type = laplacian_norm_type_pk
          )
        }
        result
      },
      stop("Unknown group_shrinkage_method.", call. = FALSE)
    )

    if (!penalize_diag_pk) {
      for (component in seq_len(K_local)) diag(out[, , component]) <- 0
    }
    if (any(!is.finite(out)) || any(out < 0)) {
      stop("Penalty weights must be finite and non-negative.", call. = FALSE)
    }
    out
  }

  normalize_native_fit <- function(raw_fit, job, K) {
    prefix <- sprintf(
      "K=%d, lambda=%.4g, rho=%.4g", K, job$lambda, job$rho
    )
    if (!isTRUE(raw_fit$ok)) {
      return(list(
        ok = FALSE,
        message = sprintf("%s: %s", prefix, raw_fit$message),
        status = NULL
      ))
    }

    result <- raw_fit$result
    status <- attr(result, "fit_status", exact = TRUE)
    state <- attr(result, "fit_state", exact = TRUE)
    trace <- attr(result, "objective_trace", exact = TRUE)
    fit_metadata <- attr(result, "fit_metadata", exact = TRUE)
    glasso_diagnostics <- attr(result, "glasso_diagnostics", exact = TRUE)
    mean_diagnostics <- attr(result, "mean_diagnostics", exact = TRUE)
    component_diagnostics <- attr(
      result, "component_diagnostics", exact = TRUE
    )
    expected_state <- c("prop", "Mu", "SigmaCube", "OmegaCube", "Z")

    structurally_valid <- is.integer(result) && length(result) == p &&
      is.list(status) && is.list(state) &&
      all(expected_state %in% names(state)) && is.numeric(trace) &&
      length(trace) == as.integer(status$iterations) + 1L &&
      all(is.finite(trace)) &&
      length(state$prop) == K &&
      identical(dim(state$Mu), c(p, K)) &&
      identical(dim(state$SigmaCube), c(p, p, K)) &&
      identical(dim(state$OmegaCube), c(p, p, K)) &&
      identical(dim(state$Z), c(n, K)) &&
      all(is.finite(unlist(state[expected_state], use.names = FALSE))) &&
      is.list(fit_metadata) &&
      identical(fit_metadata$state_ownership, "deep_copy")

    if (!isTRUE(structurally_valid)) {
      return(list(
        ok = FALSE,
        message = sprintf("%s: native fit returned incomplete fit metadata", prefix),
        status = status
      ))
    }
    trace_changes <- diff(trace)
    trace_tolerance <- if (length(trace_changes)) {
      1e-8 * (1 + abs(utils::head(trace, -1L)))
    } else {
      numeric()
    }
    status_valid <- isTRUE(status$state_valid) &&
      isTRUE(status$glasso_validated) &&
      !identical(status$mean_kkt_validated, FALSE) &&
      !identical(status$component_support_validated, FALSE) &&
      isTRUE(status$no_material_decrease) &&
      all(trace_changes >= -trace_tolerance)
    if (!status_valid || !valid_ranking_state(state, K)) {
      return(list(
        ok = FALSE,
        message = sprintf(
          "%s: native fit returned an invalid numerical state or objective trace",
          prefix
        ),
        status = status
      ))
    }
    if (!isTRUE(status$converged) || !isTRUE(status$scorable)) {
      return(list(
        ok = FALSE,
        message = sprintf(
          "%s: native fit is not scorable (reason=%s, iterations=%s)",
          prefix, as.character(status$reason), as.character(status$iterations)
        ),
        status = status
      ))
    }
    if (anyNA(result) || any(!(result %in% 0:1))) {
      return(list(
        ok = FALSE,
        message = sprintf("%s: native activity is incomplete or non-binary", prefix),
        status = status
      ))
    }

    list(
      ok = TRUE,
      activity = as.integer(result),
      status = status,
      state = state,
      objective_trace = trace,
      fit_metadata = fit_metadata,
      glasso_diagnostics = glasso_diagnostics,
      mean_diagnostics = mean_diagnostics,
      component_diagnostics = component_diagnostics
    )
  }

  run_native_job <- function(job, input_state, Pk, native_function = NULL,
                             inner_warm_start = FALSE,
                             predecessor_fitted = FALSE) {
    tryCatch(
      {
        if (is.null(native_function)) {
          native_function <- utils::getFromNamespace(
            "rcppClusteringEMGlassoWeighted", "SelvarMixMNAR"
          )
        }
        list(
          ok = TRUE,
          result = native_function(
            input_state, job$lambda, job$rho, Pk,
            inner_warm_start = inner_warm_start,
            predecessor_fitted = predecessor_fitted
          )
        )
      },
      error = function(condition) {
        list(ok = FALSE, message = conditionMessage(condition))
      }
    )
  }

  make_continuation_plan <- function(grid) {
    pair_key <- function(lambda_value, rho_value) {
      paste(
        formatC(lambda_value, digits = 17L, format = "g"),
        formatC(rho_value, digits = 17L, format = "g"),
        sep = "|"
      )
    }
    input_keys <- pair_key(grid[, "lambda"], grid[, "rho"])
    if (anyDuplicated(input_keys)) {
      stop(
        "Warm-start penalty paths require unique (lambda, rho) pairs.",
        call. = FALSE
      )
    }

    lambda_values <- sort(unique(grid[, "lambda"]), decreasing = TRUE)
    rho_values <- sort(unique(grid[, "rho"]), decreasing = TRUE)
    canonical_grid <- as.matrix(expand.grid(
      lambda = lambda_values,
      rho = rho_values,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ))
    if (nrow(canonical_grid) != nrow(grid)) {
      stop(
        "Warm-start penalty paths must be a complete lambda/rho cross-product.",
        call. = FALSE
      )
    }
    canonical_keys <- pair_key(
      canonical_grid[, "lambda"], canonical_grid[, "rho"]
    )
    execution_to_grid_index <- match(canonical_keys, input_keys)
    if (anyNA(execution_to_grid_index)) {
      stop(
        "Warm-start canonical path could not be mapped to the input grid.",
        call. = FALSE
      )
    }

    n_lambda <- length(lambda_values)
    parent_execution_index <- integer(nrow(canonical_grid))
    for (rho_index in seq_along(rho_values)) {
      for (lambda_index in seq_along(lambda_values)) {
        execution_index <- (rho_index - 1L) * n_lambda + lambda_index
        if (execution_index == 1L) {
          parent_execution_index[execution_index] <- 0L
        } else if (lambda_index == 1L) {
          parent_execution_index[execution_index] <-
            (rho_index - 2L) * n_lambda + 1L
        } else {
          parent_execution_index[execution_index] <- execution_index - 1L
        }
      }
    }
    parent_grid_index <- rep.int(NA_integer_, nrow(canonical_grid))
    has_parent <- parent_execution_index > 0L
    parent_grid_index[has_parent] <- execution_to_grid_index[
      parent_execution_index[has_parent]
    ]

    list(
      order = "lambda_descending_spokes_on_rho_descending_backbone",
      execution_grid = canonical_grid,
      execution_to_grid_index = execution_to_grid_index,
      parent_execution_index = parent_execution_index,
      parent_grid_index = parent_grid_index
    )
  }

  state_to_input <- function(data_matrix, state) {
    list(
      data_matrix,
      state$prop,
      state$Mu,
      state$SigmaCube,
      state$OmegaCube,
      state$Z
    )
  }

  run_continuation_context <- function(context, native_function = NULL) {
    plan <- context$plan
    outer_enabled <- warm_start %in% c("outer", "both")
    inner_enabled <- warm_start %in% c("inner", "both")
    fits_execution <- vector("list", nrow(plan$execution_grid))

    for (execution_index in seq_len(nrow(plan$execution_grid))) {
      parent_index <- if (outer_enabled) {
        plan$parent_execution_index[execution_index]
      } else {
        0L
      }
      grid_index <- plan$execution_to_grid_index[execution_index]
      job <- list(
        index = grid_index,
        lambda = plan$execution_grid[execution_index, "lambda"],
        rho = plan$execution_grid[execution_index, "rho"]
      )

      if (parent_index > 0L && !isTRUE(fits_execution[[parent_index]]$ok)) {
        fits_execution[[execution_index]] <- list(
          ok = FALSE,
          message = sprintf(
            paste0(
              "K=%d, lambda=%.4g, rho=%.4g: not executed because ",
              "continuation predecessor %d failed"
            ),
            context$K, job$lambda, job$rho, parent_index
          ),
          status = list(
            converged = FALSE,
            scorable = FALSE,
            reason = "predecessor_failure",
            predecessor_execution_index = parent_index
          ),
          state = NULL,
          objective_trace = numeric(),
          fit_metadata = NULL,
          glasso_diagnostics = NULL,
          mean_diagnostics = NULL,
          component_diagnostics = NULL,
          aborted = TRUE
        )
        next
      }

      predecessor_fitted <- parent_index > 0L
      input_state <- if (predecessor_fitted) {
        state_to_input(data, fits_execution[[parent_index]]$state)
      } else {
        context$base_input_state
      }
      raw_fit <- run_native_job(
        job = job,
        input_state = input_state,
        Pk = context$Pk_cube,
        native_function = native_function,
        inner_warm_start = inner_enabled,
        predecessor_fitted = predecessor_fitted
      )
      fit <- normalize_native_fit(raw_fit, job, context$K)
      initial_state_source <- if (predecessor_fitted) {
        "path_predecessor"
      } else {
        "stored_initializer"
      }
      if (is.list(fit$status)) {
        fit$status$initial_state_source <- initial_state_source
        fit$status$path_execution_index <- execution_index
        fit$status$path_grid_index <- grid_index
        fit$status$parent_execution_index <- if (parent_index) {
          parent_index
        } else {
          NA_integer_
        }
        fit$status$parent_grid_index <- plan$parent_grid_index[execution_index]
      }
      if (is.list(fit$fit_metadata)) {
        fit$fit_metadata$warm_start_mode <- warm_start
        fit$fit_metadata$initial_state_source <- initial_state_source
      }
      if (is.list(fit$glasso_diagnostics)) {
        fit$glasso_diagnostics$initial_state_source <- initial_state_source
        fit$glasso_diagnostics$parent_execution_index <- if (parent_index) {
          parent_index
        } else {
          NA_integer_
        }
      }
      fit$aborted <- FALSE
      fits_execution[[execution_index]] <- fit
    }

    fits_output <- vector("list", length(fits_execution))
    fits_output[plan$execution_to_grid_index] <- fits_execution
    list(
      fits = fits_output,
      path = c(
        plan,
        list(
          warm_start_mode = warm_start,
          outer_continuation = outer_enabled,
          inner_glasso_warm_start = inner_enabled
        )
      )
    )
  }

  data <- as.matrix(data)
  if (!is.numeric(data) || !length(data) || any(!is.finite(data))) {
    stop("data must be a non-empty finite numeric matrix.", call. = FALSE)
  }
  n <- nrow(data)
  p <- ncol(data)
  nbcluster <- as.integer(nbcluster)
  if (!length(nbcluster) || anyNA(nbcluster) || any(nbcluster < 1L) ||
      any(nbcluster > n)) {
    stop("nbcluster must contain positive integers no larger than nrow(data).",
         call. = FALSE)
  }
  group_shrinkage_method <- match.arg(group_shrinkage_method)
  laplacian_target_type <- match.arg(laplacian_target_type)
  laplacian_norm_type <- match.arg(laplacian_norm_type)
  initialize <- match.arg(initialize)
  warm_start <- match.arg(warm_start)
  n.start <- as.integer(n.start)
  lambda_omega_0 <- as.double(lambda_omega_0)
  nbcores <- max(1L, as.integer(nbcores))
  if (!is.logical(penalize_diag) || length(penalize_diag) != 1L ||
      is.na(penalize_diag)) {
    stop("penalize_diag must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(penalize_diag)) {
    stop(
      paste0(
        "penalize_diag=TRUE is unsupported: the native ranking ",
        "objective penalizes off-diagonal precision entries only."
      ),
      call. = FALSE
    )
  }

  if (is.null(penalty_grids)) {
    shared_grid <- make_grid(lambda, rho, "lambda/rho")
    grid_specs <- lapply(seq_along(nbcluster), function(index) {
      list(
        K = nbcluster[index],
        grid = shared_grid,
        initializer = NULL,
        Pk_cube = NULL,
        source = "shared_user_grid"
      )
    })
  } else {
    if (!is.list(penalty_grids) || length(penalty_grids) != length(nbcluster)) {
      stop("penalty_grids must be a list aligned one-to-one with nbcluster.",
           call. = FALSE)
    }
    grid_specs <- lapply(seq_along(nbcluster), function(index) {
      K <- nbcluster[index]
      entry <- penalty_grids[[index]]
      source <- sprintf("penalty_grids[[%d]]", index)
      if (!is.list(entry) ||
          (!is.null(entry$K) && !identical(as.integer(entry$K), K))) {
        stop(sprintf("%s is not aligned with K=%d.", source, K), call. = FALSE)
      }
      list(
        K = K,
        grid = make_grid(entry$lambda_mu, entry$rho, source),
        initializer = validate_initializer(entry$initializer, K, source),
        Pk_cube = if (is.null(entry$Pk_cube)) {
          NULL
        } else {
          validate_Pk(entry$Pk_cube, K, source)
        },
        source = "stored_per_K_grid"
      )
    })
    grid_sizes <- vapply(grid_specs, function(spec) nrow(spec$grid), integer(1))
    if (length(unique(grid_sizes)) != 1L) {
      stop(
        "All per-K penalty grids must have equal cardinality for the activity array.",
        call. = FALSE
      )
    }
  }

  names(grid_specs) <- paste0("K", nbcluster)
  n_grid <- nrow(grid_specs[[1L]]$grid)
  parallel_axis <- if (identical(warm_start, "none")) {
    "independent_grid_points_within_K"
  } else {
    "independent_K_paths"
  }
  nbcores <- if (identical(warm_start, "none")) {
    min(nbcores, n_grid)
  } else {
    min(nbcores, length(nbcluster))
  }
  is_windows <- identical(unname(Sys.info()["sysname"]), "Windows")
  cl <- NULL
  on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)

  if (nbcores > 1L && is_windows) {
    cl <- parallel::makePSOCKcluster(nbcores)
    parallel::clusterCall(cl, function(library_paths) {
      .libPaths(library_paths)
      loadNamespace("SelvarMixMNAR")
      NULL
    }, .libPaths())
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(min_scorable_fraction) ||
      length(min_scorable_fraction) != 1L ||
      !is.finite(min_scorable_fraction) ||
      min_scorable_fraction <= 0 || min_scorable_fraction > 1) {
    stop("min_scorable_fraction must be in (0, 1].", call. = FALSE)
  }
  if (verbose) message("Initializing parameters...")
  initializers <- lapply(seq_along(nbcluster), function(index) {
    spec <- grid_specs[[index]]
    if (!is.null(spec$initializer)) return(spec$initializer)
    init <- InitParameter(
      data = data,
      nbClust = spec$K,
      init = initialize,
      n.start = n.start,
      lambda_omega_0 = lambda_omega_0
    )
    validate_initializer(init, spec$K, "InitParameter")
  })
  names(initializers) <- paste0("K", nbcluster)
  if (verbose) message("Initialization complete.")

  row_labels <- if (all(vapply(
    grid_specs,
    function(spec) identical(spec$grid, grid_specs[[1L]]$grid),
    logical(1)
  ))) {
    paste(
      "L", grid_specs[[1L]]$grid[, "lambda"],
      "R", grid_specs[[1L]]$grid[, "rho"]
    )
  } else {
    paste0("grid_", seq_len(n_grid))
  }
  VarRole <- array(
    NA_integer_,
    dim = c(n_grid, p, length(nbcluster)),
    dimnames = list(row_labels, colnames(data), paste("K", nbcluster))
  )
  all_status <- all_states <- all_traces <- all_glasso <- all_means <-
    all_components <- all_paths <- all_failures <-
    vector("list", length(nbcluster))
  result_names <- paste0("K", nbcluster)
  names(all_status) <- names(all_states) <- names(all_traces) <-
    names(all_glasso) <- names(all_means) <- names(all_components) <-
    names(all_paths) <- names(all_failures) <- result_names

  contexts <- lapply(seq_along(nbcluster), function(k_idx) {
    K_current <- nbcluster[k_idx]
    spec <- grid_specs[[k_idx]]
    P_init <- initializers[[k_idx]]
    Pk_cube <- if (!is.null(spec$Pk_cube)) {
      spec$Pk_cube
    } else {
      compute_Pk(
        omega0_cube = P_init$OmegaCube,
        method = group_shrinkage_method,
        distance_method = distance_method,
        laplacian_target_type_pk = laplacian_target_type,
        adj_threshold_pk = as.double(adj_threshold),
        laplacian_norm_type_pk = laplacian_norm_type,
        eps_w0 = epsilon_weighted_by_W0,
        penalize_diag_pk = penalize_diag
      )
    }
    Pk_cube <- validate_Pk(Pk_cube, K_current, "ranking penalty weights")
    list(
      K = K_current,
      pen_grid = spec$grid,
      Pk_cube = Pk_cube,
      base_input_state = list(
        data,
        P_init$prop,
        P_init$Mu,
        P_init$SigmaCube,
        P_init$OmegaCube,
        P_init$Z
      ),
      plan = if (identical(warm_start, "none")) {
        NULL
      } else {
        make_continuation_plan(spec$grid)
      }
    )
  })
  names(contexts) <- result_names

  if (identical(warm_start, "none")) {
    path_results <- lapply(contexts, function(context) {
      jobs <- lapply(seq_len(nrow(context$pen_grid)), function(grid_index) {
        list(
          index = grid_index,
          lambda = context$pen_grid[grid_index, "lambda"],
          rho = context$pen_grid[grid_index, "rho"]
        )
      })
      if (verbose) message(sprintf(
        "Processing K = %d with %d independent cold-start grid fits...",
        context$K, length(jobs)
      ))
      if (nbcores > 1L && is_windows) {
        raw_fits <- parallel::parLapply(
          cl, jobs, run_native_job,
          input_state = context$base_input_state,
          Pk = context$Pk_cube,
          native_function = NULL,
          inner_warm_start = FALSE,
          predecessor_fitted = FALSE
        )
      } else {
        raw_fits <- parallel::mclapply(
          jobs, run_native_job,
          input_state = context$base_input_state,
          Pk = context$Pk_cube,
          native_function = rcppClusteringEMGlassoWeighted,
          inner_warm_start = FALSE,
          predecessor_fitted = FALSE,
          mc.cores = nbcores,
          mc.preschedule = TRUE,
          mc.cleanup = TRUE
        )
      }
      fits <- Map(
        function(raw_fit, job) {
          normalize_native_fit(raw_fit, job, context$K)
        },
        raw_fits, jobs
      )
      list(
        fits = fits,
        path = list(
          order = "input_grid_index",
          execution_grid = context$pen_grid,
          execution_to_grid_index = seq_len(nrow(context$pen_grid)),
          parent_execution_index = rep.int(0L, nrow(context$pen_grid)),
          parent_grid_index = rep.int(NA_integer_, nrow(context$pen_grid)),
          warm_start_mode = "none",
          outer_continuation = FALSE,
          inner_glasso_warm_start = FALSE
        )
      )
    })
  } else {
    if (verbose) message(sprintf(
      "Processing %d independent K continuation path%s with warm_start='%s'...",
      length(contexts), if (length(contexts) == 1L) "" else "s", warm_start
    ))
    if (nbcores > 1L && is_windows) {
      path_results <- parallel::parLapply(
        cl, contexts, run_continuation_context, native_function = NULL
      )
    } else {
      path_results <- parallel::mclapply(
        contexts, run_continuation_context,
        native_function = NULL,
        mc.cores = nbcores,
        mc.preschedule = TRUE,
        mc.cleanup = TRUE
      )
    }
  }
  names(path_results) <- result_names

  for (k_idx in seq_along(nbcluster)) {
    K_current <- nbcluster[k_idx]
    pen_grid <- contexts[[k_idx]]$pen_grid
    fits <- path_results[[k_idx]]$fits
    failed_rows <- which(!vapply(fits, `[[`, logical(1), "ok"))
    if (length(failed_rows)) {
      for (failed_index in failed_rows) {
        warning(fits[[failed_index]]$message, call. = FALSE)
      }
      failed_grid <- pen_grid[failed_rows, , drop = FALSE]
      pair_text <- paste(
        sprintf(
          "(lambda=%.4g, rho=%.4g)",
          failed_grid[, "lambda"], failed_grid[, "rho"]
        ),
        collapse = ", "
      )
      if (length(failed_rows) == nrow(pen_grid)) {
        stop_grid_fit(
          sprintf(
            paste0(
              "Penalized EM grid failed for K=%d at %d of %d parameter pairs; ",
              "failed fits were not scored and nonconverged fits were not scored: %s"
            ),
            K_current, length(failed_rows), nrow(pen_grid), pair_text
          ),
          K = K_current,
          failed_grid = failed_grid,
          n_grid = nrow(pen_grid),
          failures = fits[failed_rows]
        )
      }
      scorable_fraction <- 1 - length(failed_rows) / nrow(pen_grid)
      if (scorable_fraction < min_scorable_fraction) {
        stop_grid_fit(
          sprintf(
            paste0(
              "Penalized EM grid for K=%d retained only %.3f of parameter ",
              "pairs, below min_scorable_fraction=%.3f; failed fits were ",
              "not scored: %s"
            ),
            K_current, scorable_fraction, min_scorable_fraction, pair_text
          ),
          K = K_current,
          failed_grid = failed_grid,
          n_grid = nrow(pen_grid),
          failures = fits[failed_rows]
        )
      }
    }

    successful_rows <- setdiff(seq_len(nrow(pen_grid)), failed_rows)
    var_role_k <- matrix(NA_integer_, nrow = nrow(pen_grid), ncol = p)
    for (successful_index in successful_rows) {
      var_role_k[successful_index, ] <- fits[[successful_index]]$activity
    }
    expected_dim <- c(nrow(pen_grid), p)
    if (!identical(dim(var_role_k), expected_dim) ||
        anyNA(var_role_k[successful_rows, , drop = FALSE]) ||
        any(!(var_role_k[successful_rows, , drop = FALSE] %in% 0:1)) ||
        (length(failed_rows) &&
         any(!is.na(var_role_k[failed_rows, , drop = FALSE])))) {
      stop_grid_fit(
        sprintf(
          "Penalized EM grid for K=%d returned an invalid activity matrix; no ranking was produced.",
          K_current
        ),
        K = K_current,
        failed_grid = pen_grid,
        n_grid = nrow(pen_grid),
        failures = fits
      )
    }

    VarRole[, , k_idx] <- var_role_k
    all_status[[k_idx]] <- lapply(fits, `[[`, "status")
    all_states[[k_idx]] <- lapply(fits, `[[`, "state")
    all_traces[[k_idx]] <- lapply(fits, `[[`, "objective_trace")
    all_glasso[[k_idx]] <- lapply(fits, `[[`, "glasso_diagnostics")
    all_means[[k_idx]] <- lapply(fits, `[[`, "mean_diagnostics")
    all_components[[k_idx]] <- lapply(
      fits, `[[`, "component_diagnostics"
    )
    all_paths[[k_idx]] <- path_results[[k_idx]]$path
    all_failures[[k_idx]] <- if (length(failed_rows)) {
      list(
        rows = failed_rows,
        grid = pen_grid[failed_rows, , drop = FALSE],
        messages = vapply(fits[failed_rows], `[[`, character(1), "message")
      )
    } else {
      list(
        rows = integer(),
        grid = pen_grid[FALSE, , drop = FALSE],
        messages = character()
      )
    }
  }

  attr(VarRole, "penalty_grids") <- lapply(grid_specs, `[[`, "grid")
  attr(VarRole, "initializers") <- initializers
  attr(VarRole, "fit_status") <- all_status
  attr(VarRole, "fit_states") <- all_states
  attr(VarRole, "objective_traces") <- all_traces
  attr(VarRole, "glasso_diagnostics") <- all_glasso
  attr(VarRole, "mean_diagnostics") <- all_means
  attr(VarRole, "component_diagnostics") <- all_components
  attr(VarRole, "path_metadata") <- all_paths
  attr(VarRole, "grid_failures") <- all_failures
  inner_enabled <- warm_start %in% c("inner", "both")
  outer_enabled <- warm_start %in% c("outer", "both")
  attr(VarRole, "fit_metadata") <- list(
    state_ownership = "deep_copy_per_grid_fit",
    initialization_order = "nbcluster_input_order",
    grid_order = "rho_outer_lambda_inner",
    grid_result_order = "input_grid_index",
    grid_fits = if (outer_enabled) {
      "canonical_continuation_tree"
    } else if (inner_enabled) {
      "independent_inner_warm_starts"
    } else {
      "independent_cold_starts"
    },
    path_warm_start = outer_enabled,
    inner_glasso_warm_start = inner_enabled,
    warm_start_mode = warm_start,
    continuation_order = if (outer_enabled) {
      "lambda_descending_spokes_on_rho_descending_backbone"
    } else {
      NA_character_
    },
    parallel_axis = parallel_axis,
    effective_workers = nbcores,
    glasso_start = if (identical(warm_start, "none")) {
      "cold"
    } else {
      "cold_then_previous_em_warm"
    },
    mean_solver = "gauss_seidel_coordinate_descent_to_subgradient_kkt",
    glasso_solver_tolerance = 1e-7,
    glasso_validation =
      "spd_plus_normwise_inverse_backward_error_plus_kkt",
    component_support_rule =
      "p_effective_unpenalized_or_positive_variance_fully_penalized",
    partial_grid_rule = paste0(
      "retain grid alignment with all-NA activity rows for explicit failed ",
      "fits; rank only scorable rows; enforce min_scorable_fraction"
    ),
    min_scorable_fraction = min_scorable_fraction
  )
  if (verbose) message("Variable ranking complete.")
  VarRole
}

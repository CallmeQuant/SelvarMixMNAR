.numeric_initialization_schema_version <- "1.0.0"

.numeric_initialization_state_names <- c(
  "prop", "Mu", "SigmaCube", "OmegaCube", "Z"
)

.numeric_init_integer <- function(value, name, minimum = 1L) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
      value != round(value) || value < minimum ||
      value > .Machine$integer.max) {
    stop(sprintf("%s must be one finite integer >= %d.", name, minimum),
         call. = FALSE)
  }
  as.integer(value)
}

.numeric_init_seed <- function(seed, name = "seed") {
  if (is.null(seed)) return(NULL)
  .numeric_init_integer(seed, name, minimum = 0L)
}

.numeric_init_clone_state <- function(state) {
  list(
    prop = structure(as.double(state$prop), names = names(state$prop)),
    Mu = matrix(
      as.double(state$Mu), nrow = nrow(state$Mu), ncol = ncol(state$Mu),
      dimnames = dimnames(state$Mu)
    ),
    SigmaCube = array(
      as.double(state$SigmaCube), dim = dim(state$SigmaCube),
      dimnames = dimnames(state$SigmaCube)
    ),
    OmegaCube = array(
      as.double(state$OmegaCube), dim = dim(state$OmegaCube),
      dimnames = dimnames(state$OmegaCube)
    ),
    Z = matrix(
      as.double(state$Z), nrow = nrow(state$Z), ncol = ncol(state$Z),
      dimnames = dimnames(state$Z)
    )
  )
}

.numeric_init_abort <- function(message, method, seed = NA_integer_,
                                code = "initialization_failed",
                                candidate_provenance = NULL,
                                parent = NULL) {
  condition <- structure(
    list(
      message = as.character(message),
      call = NULL,
      method = as.character(method),
      seed = seed,
      failure = list(
        code = as.character(code),
        message = as.character(message),
        parent_class = if (is.null(parent)) NULL else class(parent)
      ),
      candidate_provenance = candidate_provenance,
      parent = parent
    ),
    class = c("selvarmix_initialization_error", "error", "condition")
  )
  stop(condition)
}

.numeric_init_with_rng <- function(seed, code) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (!is.null(seed)) set.seed(seed)
  force(code)
}

.numeric_init_canonical_method <- function(method) {
  choices <- c(
    "mclust_hc", "mclust_em", "kmeans", "random", "user",
    "previous_fit", "deterministic_multistart", "hc"
  )
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("init must be one initialization method token.", call. = FALSE)
  }
  requested <- match.arg(method, choices)
  canonical <- if (identical(requested, "hc")) "mclust_hc" else requested
  list(
    requested = requested,
    canonical = canonical,
    legacy_alias = identical(requested, "hc")
  )
}

.numeric_init_validate_data <- function(data, nbClust) {
  data <- as.matrix(data)
  if (!is.numeric(data) || length(dim(data)) != 2L ||
      any(dim(data) < 1L) || any(!is.finite(data))) {
    stop("data must be a non-empty finite numeric matrix.", call. = FALSE)
  }
  K <- .numeric_init_integer(nbClust, "nbClust")
  if (K > nrow(data)) {
    stop("nbClust cannot exceed the number of observations.", call. = FALSE)
  }
  centered <- sweep(data, 2L, colMeans(data), FUN = "-")
  column_scale <- sqrt(colSums(centered^2) / nrow(data))
  magnitude <- pmax(1, apply(abs(data), 2L, max))
  degenerate <- !is.finite(column_scale) |
    column_scale <= .Machine$double.eps * magnitude
  if (any(degenerate)) {
    stop(
      sprintf(
        "data contain non-varying column(s): %s; filter them before initialization.",
        paste(which(degenerate), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  list(data = data, K = K, n = nrow(data), p = ncol(data))
}

.numeric_init_extract_supplied <- function(object, source) {
  if (is.null(object)) {
    stop(sprintf("%s must be supplied.", source), call. = FALSE)
  }

  # A fitted object may supply its reusable numeric state directly when
  # `previous_fit` is requested.
  if (is.list(object) && is.list(object$initialization) &&
      !is.null(object$initialization$selected)) {
    selected <- .numeric_init_extract_supplied(
      object$initialization$selected,
      paste0(source, "$initialization$selected")
    )
    selected$source_class <- class(object)
    selected$source_path <- "initialization$selected"
    return(selected)
  }

  state <- NULL
  if (is.list(object) &&
      all(.numeric_initialization_state_names %in% names(object))) {
    state <- object[.numeric_initialization_state_names]
  } else if (is.list(object) && is.list(object$state) &&
             all(.numeric_initialization_state_names %in% names(object$state))) {
    state <- object$state[.numeric_initialization_state_names]
  } else {
    fitted_state <- attr(object, "fit_state", exact = TRUE)
    if (is.list(fitted_state) &&
        all(.numeric_initialization_state_names %in% names(fitted_state))) {
      state <- fitted_state[.numeric_initialization_state_names]
    }
  }
  if (is.null(state)) {
    stop(
      sprintf(
        "%s does not contain prop, Mu, SigmaCube, OmegaCube, and Z.",
        source
      ),
      call. = FALSE
    )
  }

  report <- attr(object, "initialization_report", exact = TRUE)
  partition <- if (is.list(object) && !is.null(object$partition)) {
    object$partition
  } else if (is.list(report) && !is.null(report$partition)) {
    report$partition
  } else {
    NULL
  }

  list(
    state = state,
    partition = partition,
    source_class = class(object),
    source_path = "direct_state",
    source_schema_version = if (is.list(object)) {
      object$schema_version
    } else {
      NULL
    }
  )
}

.numeric_init_validate_state <- function(state, data, K, partition = NULL,
                                         source = "initializer",
                                         tolerance = 1e-8) {
  n <- nrow(data)
  p <- ncol(data)
  if (!is.list(state) ||
      !all(.numeric_initialization_state_names %in% names(state))) {
    stop(sprintf("%s has an incomplete numeric state.", source), call. = FALSE)
  }
  state <- state[.numeric_initialization_state_names]

  shape_ok <- is.numeric(state$prop) && length(state$prop) == K &&
    is.matrix(state$Mu) && identical(dim(state$Mu), c(p, K)) &&
    is.array(state$SigmaCube) &&
    identical(dim(state$SigmaCube), c(p, p, K)) &&
    is.array(state$OmegaCube) &&
    identical(dim(state$OmegaCube), c(p, p, K)) &&
    is.matrix(state$Z) && identical(dim(state$Z), c(n, K))
  if (!shape_ok) {
    stop(sprintf("%s has dimensions inconsistent with data and nbClust.", source),
         call. = FALSE)
  }
  if (any(!is.finite(unlist(state, use.names = FALSE)))) {
    stop(sprintf("%s contains non-finite numeric values.", source),
         call. = FALSE)
  }
  if (any(state$prop <= 0) ||
      abs(sum(state$prop) - 1) > tolerance * (1 + K)) {
    stop(sprintf("%s has invalid mixture proportions.", source),
         call. = FALSE)
  }
  if (any(state$Z < 0) ||
      any(abs(rowSums(state$Z) - 1) > tolerance * (1 + K)) ||
      any(colSums(state$Z) <= 0)) {
    stop(sprintf("%s has invalid responsibilities.", source), call. = FALSE)
  }
  if (max(abs(colMeans(state$Z) - state$prop)) >
      1e-6 * (1 + max(state$prop))) {
    stop(
      sprintf("%s has proportions inconsistent with its responsibilities.", source),
      call. = FALSE
    )
  }

  inverse_residual <- numeric(K)
  minimum_eigenvalue <- numeric(K)
  for (component in seq_len(K)) {
    covariance <- matrix(
      state$SigmaCube[, , component], nrow = p, ncol = p
    )
    precision <- matrix(
      state$OmegaCube[, , component], nrow = p, ncol = p
    )
    covariance_scale <- max(norm(covariance, type = "I"),
                            .Machine$double.xmin)
    precision_scale <- max(norm(precision, type = "I"),
                           .Machine$double.xmin)
    if (norm(covariance - t(covariance), type = "I") >
        tolerance * covariance_scale ||
        norm(precision - t(precision), type = "I") >
        tolerance * precision_scale) {
      stop(sprintf("%s has a non-symmetric covariance/precision pair.", source),
           call. = FALSE)
    }
    covariance <- 0.5 * (covariance + t(covariance))
    precision <- 0.5 * (precision + t(precision))
    chol_covariance <- tryCatch(chol(covariance), error = function(e) NULL)
    chol_precision <- tryCatch(chol(precision), error = function(e) NULL)
    if (is.null(chol_covariance) || is.null(chol_precision)) {
      stop(sprintf("%s has a non-positive-definite covariance/precision pair.",
                   source), call. = FALSE)
    }
    minimum_eigenvalue[component] <- min(eigen(
      covariance, symmetric = TRUE, only.values = TRUE
    )$values)
    identity <- diag(p)
    inverse_residual[component] <- max(
      norm(covariance %*% precision - identity, type = "I"),
      norm(precision %*% covariance - identity, type = "I")
    )
    if (!is.finite(inverse_residual[component]) ||
        inverse_residual[component] > 1e-6 * (1 + p)) {
      stop(sprintf("%s has covariance and precision matrices that are not inverses.",
                   source), call. = FALSE)
    }
  }

  implied_partition <- max.col(state$Z, ties.method = "first")
  if (is.null(partition)) {
    partition <- implied_partition
  } else {
    if (!is.numeric(partition) || length(partition) != n ||
        any(!is.finite(partition)) || any(partition != round(partition)) ||
        any(!(partition %in% seq_len(K)))) {
      stop(sprintf("%s has an invalid partition.", source), call. = FALSE)
    }
    partition <- as.integer(partition)
    if (!identical(partition, as.integer(implied_partition))) {
      stop(sprintf("%s has a partition inconsistent with max.col(Z).", source),
           call. = FALSE)
    }
  }

  list(
    state = .numeric_init_clone_state(state),
    partition = as.integer(partition),
    validation = list(
      inverse_residual = inverse_residual,
      minimum_covariance_eigenvalue = minimum_eigenvalue,
      tolerance = tolerance
    )
  )
}

.numeric_init_loglik <- function(data, state) {
  n <- nrow(data)
  p <- ncol(data)
  K <- length(state$prop)
  component_log_density <- matrix(NA_real_, nrow = n, ncol = K)
  constant <- p * log(2 * pi)

  for (component in seq_len(K)) {
    covariance <- matrix(
      state$SigmaCube[, , component], nrow = p, ncol = p
    )
    upper <- chol(covariance)
    centered <- sweep(data, 2L, state$Mu[, component], FUN = "-")
    standardized <- backsolve(upper, t(centered), transpose = TRUE)
    quadratic <- colSums(standardized * standardized)
    log_determinant <- 2 * sum(log(diag(upper)))
    component_log_density[, component] <-
      log(state$prop[component]) -
      0.5 * (constant + log_determinant + quadratic)
  }

  row_maximum <- apply(component_log_density, 1L, max)
  value <- sum(row_maximum + log(rowSums(exp(
    component_log_density - row_maximum
  ))))
  if (!is.finite(value)) {
    stop("The Gaussian-mixture initialization objective is non-finite.",
         call. = FALSE)
  }
  as.double(value)
}

.numeric_init_hc_pairs <- function(data, model_name) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("mclust is required for mclust_hc and mclust_em initialization.",
         call. = FALSE)
  }
  hc_object <- do.call(
    mclust::hc,
    list(data = data, modelName = model_name, use = "SVD"),
    envir = asNamespace("mclust")
  )
  hc_pairs <- unclass(hc_object)
  storage.mode(hc_pairs) <- "double"
  if (!is.numeric(hc_pairs) || !is.matrix(hc_pairs) ||
      !identical(dim(hc_pairs), c(2L, nrow(data) - 1L)) ||
      any(!is.finite(hc_pairs)) || inherits(hc_pairs, "hc")) {
    stop("mclust::hc did not return valid numeric hcPairs.", call. = FALSE)
  }
  hc_pairs
}

.numeric_init_mclust_em <- function(data, K, model_name, hc_pairs) {
  if (!is.numeric(hc_pairs) || !is.matrix(hc_pairs) ||
      inherits(hc_pairs, "hc")) {
    stop("Mclust initialization requires numeric hcPairs.", call. = FALSE)
  }
  warnings <- character()
  fit <- withCallingHandlers(
    do.call(
      mclust::Mclust,
      list(
        data = data,
        G = K,
        modelNames = model_name,
        initialization = list(hcPairs = hc_pairs),
        warn = TRUE,
        verbose = FALSE
      ),
      envir = asNamespace("mclust")
    ),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  if (is.null(fit) || !inherits(fit, "Mclust")) {
    stop("mclust::Mclust did not return a fitted Mclust object.",
         call. = FALSE)
  }
  list(fit = fit, warnings = warnings)
}

.numeric_init_covariances_from_partition <- function(data, partition, K,
                                                     lambda_omega_0,
                                                     epsilon_pd) {
  n <- nrow(data)
  p <- ncol(data)
  Z <- matrix(0, nrow = n, ncol = K)
  Z[cbind(seq_len(n), partition)] <- 1
  component_mass <- colSums(Z)
  if (any(component_mass <= 0)) {
    stop("The proposed partition contains an empty component.", call. = FALSE)
  }
  prop <- component_mass / n
  Mu <- matrix(NA_real_, nrow = p, ncol = K)
  SigmaCube <- array(NA_real_, dim = c(p, p, K))
  OmegaCube <- array(NA_real_, dim = c(p, p, K))
  covariance_diagnostics <- vector("list", K)
  global_centered <- sweep(data, 2L, colMeans(data), FUN = "-")
  global_covariance <- crossprod(global_centered) / n
  global_covariance <- 0.5 * (global_covariance + t(global_covariance))
  global_scale <- norm(global_covariance, type = "I")
  if (!is.finite(global_scale) || global_scale <= 0) {
    stop(
      "The data have zero covariance scale; a Gaussian initializer is not identifiable.",
      call. = FALSE
    )
  }

  for (component in seq_len(K)) {
    rows <- which(partition == component)
    Mu[, component] <- colMeans(data[rows, , drop = FALSE])
    centered <- sweep(
      data[rows, , drop = FALSE], 2L, Mu[, component], FUN = "-"
    )
    empirical <- crossprod(centered) / length(rows)
    empirical <- matrix(empirical, nrow = p, ncol = p)
    empirical <- 0.5 * (empirical + t(empirical))
    empirical_scale <- norm(empirical, type = "I")
    relative_scale <- if (is.finite(empirical_scale) && empirical_scale > 0) {
      empirical_scale
    } else {
      global_scale
    }
    empirical_minimum <- min(eigen(
      empirical, symmetric = TRUE, only.values = TRUE
    )$values)

    if (empirical_minimum > epsilon_pd * relative_scale) {
      covariance <- empirical
      precision <- chol2inv(chol(covariance))
      covariance_diagnostics[[component]] <- list(
        component = component,
        solver = "empirical_inverse",
        regularized = FALSE,
        penalty = 0,
        iterations = 0L,
        errflag = 0L,
        empirical_minimum_eigenvalue = empirical_minimum,
        covariance_scale = relative_scale
      )
    } else {
      # Treat lambda_omega_0 as a dimensionless multiplier. Scaling the
      # regularizer by the empirical covariance scale makes initialization
      # equivariant under a global change of measurement units.
      penalty <- max(
        2 * lambda_omega_0 * relative_scale / component_mass[component],
        epsilon_pd * relative_scale
      )
      if (!is.finite(penalty) || penalty < 0) {
        stop("The singular-covariance regularization penalty is invalid.",
             call. = FALSE)
      }
      off_diagonal <- row(empirical) != col(empirical)
      if (!length(empirical[off_diagonal]) ||
          all(empirical[off_diagonal] == 0)) {
        diagonal <- diag(empirical) + penalty
        if (any(!is.finite(diagonal)) || any(diagonal <= 0)) {
          stop("Analytic diagonal regularization was not positive definite.",
               call. = FALSE)
        }
        covariance <- diag(diagonal, nrow = p)
        precision <- diag(1 / diagonal, nrow = p)
        covariance_diagnostics[[component]] <- list(
          component = component,
          solver = "analytic_diagonal_glasso",
          regularized = TRUE,
          penalty = penalty,
          iterations = 0L,
          errflag = 0L,
          empirical_minimum_eigenvalue = empirical_minimum,
          covariance_scale = relative_scale
        )
      } else {
        if (!requireNamespace("glassoFast", quietly = TRUE)) {
          stop("glassoFast is required to regularize a singular initializer.",
               call. = FALSE)
        }
        maximum_iterations <- 1000L
        glasso_fit <- glassoFast::glassoFast(
          S = empirical,
          rho = penalty,
          start = "cold",
          thr = 1e-4,
          maxIt = maximum_iterations
        )
        valid_fit <- is.list(glasso_fit) && is.matrix(glasso_fit$w) &&
          is.matrix(glasso_fit$wi) &&
          identical(dim(glasso_fit$w), c(p, p)) &&
          identical(dim(glasso_fit$wi), c(p, p)) &&
          all(is.finite(glasso_fit$w)) && all(is.finite(glasso_fit$wi)) &&
          length(glasso_fit$errflag) == 1L && glasso_fit$errflag == 0L &&
          length(glasso_fit$niter) == 1L && is.finite(glasso_fit$niter) &&
          glasso_fit$niter < maximum_iterations
        if (!valid_fit) {
          stop("glassoFast failed or reached its cap during initialization.",
               call. = FALSE)
        }
        covariance <- glasso_fit$w
        precision <- glasso_fit$wi
        covariance_diagnostics[[component]] <- list(
          component = component,
          solver = "glassoFast",
          regularized = TRUE,
          penalty = penalty,
          iterations = as.integer(glasso_fit$niter),
          errflag = as.integer(glasso_fit$errflag),
          empirical_minimum_eigenvalue = empirical_minimum,
          covariance_scale = relative_scale
        )
      }
    }
    SigmaCube[, , component] <- covariance
    OmegaCube[, , component] <- precision
  }

  list(
    state = list(
      prop = prop,
      Mu = Mu,
      SigmaCube = SigmaCube,
      OmegaCube = OmegaCube,
      Z = Z
    ),
    partition = as.integer(partition),
    details = list(
      component_mass = component_mass,
      covariance = covariance_diagnostics
    )
  )
}

.numeric_init_state_from_mclust <- function(fit, data, K) {
  parameters <- fit$parameters
  if (!is.list(parameters) || is.null(parameters$variance)) {
    stop("The Mclust object has no usable Gaussian parameters.", call. = FALSE)
  }
  prop <- as.double(parameters$pro)
  Mu <- parameters$mean
  if (is.null(dim(Mu))) Mu <- matrix(Mu, ncol = K)
  Mu <- matrix(as.double(Mu), nrow = ncol(data), ncol = K,
               dimnames = dimnames(Mu))

  SigmaCube <- parameters$variance$sigma
  if (is.null(SigmaCube)) SigmaCube <- parameters$variance$Sigma
  if (is.null(SigmaCube)) {
    stop("The Mclust object has no covariance array.", call. = FALSE)
  }
  SigmaCube <- array(
    as.double(SigmaCube), dim = c(ncol(data), ncol(data), K),
    dimnames = if (length(dim(SigmaCube)) == 3L) dimnames(SigmaCube) else NULL
  )
  OmegaCube <- array(NA_real_, dim = dim(SigmaCube))
  for (component in seq_len(K)) {
    covariance <- matrix(
      SigmaCube[, , component], nrow = ncol(data), ncol = ncol(data)
    )
    OmegaCube[, , component] <- chol2inv(chol(covariance))
  }
  Z <- fit$z
  if (is.null(dim(Z))) Z <- matrix(Z, ncol = K)
  Z <- matrix(as.double(Z), nrow = nrow(data), ncol = K,
              dimnames = dimnames(Z))
  partition <- as.integer(fit$classification)

  list(
    state = list(
      prop = prop,
      Mu = Mu,
      SigmaCube = SigmaCube,
      OmegaCube = OmegaCube,
      Z = Z
    ),
    partition = partition,
    details = list(
      backend_object_class = class(fit),
      requested_model = NULL,
      fitted_model = fit$modelName,
      backend_loglik = as.double(fit$loglik),
      backend_bic = as.double(fit$bic)
    )
  )
}

.numeric_init_run_candidate <- function(method, data, K, n.start,
                                        lambda_omega_0, epsilon_pd,
                                        user_init, previous_fit,
                                        mclust_model) {
  result <- switch(
    method,
    kmeans = {
      fit <- stats::kmeans(
        data, centers = K, nstart = n.start, iter.max = 1000L
      )
      result <- .numeric_init_covariances_from_partition(
        data, as.integer(fit$cluster), K, lambda_omega_0, epsilon_pd
      )
      result$details$partition_method <- "stats::kmeans"
      result$details$withinss <- as.double(fit$tot.withinss)
      result
    },
    random = {
      labels <- rep(seq_len(K), length.out = nrow(data))
      labels <- sample(labels, length(labels), replace = FALSE)
      result <- .numeric_init_covariances_from_partition(
        data, as.integer(labels), K, lambda_omega_0, epsilon_pd
      )
      result$details$partition_method <- "balanced_random_partition"
      result
    },
    mclust_hc = {
      hc_pairs <- .numeric_init_hc_pairs(data, mclust_model)
      labels <- as.integer(mclust::hclass(hc_pairs, K))
      result <- .numeric_init_covariances_from_partition(
        data, labels, K, lambda_omega_0, epsilon_pd
      )
      result$details$partition_method <- "mclust::hc/hclass"
      result$details$hc_pairs <- list(
        type = typeof(hc_pairs),
        numeric = is.numeric(hc_pairs),
        inherits_hc = inherits(hc_pairs, "hc"),
        dim = dim(hc_pairs)
      )
      result
    },
    mclust_em = {
      hc_pairs <- .numeric_init_hc_pairs(data, mclust_model)
      fitted <- .numeric_init_mclust_em(
        data, K, mclust_model, hc_pairs
      )
      result <- .numeric_init_state_from_mclust(fitted$fit, data, K)
      result$details$requested_model <- mclust_model
      result$details$mclust_warnings <- fitted$warnings
      result$details$mclust_initialization <-
        "list(hcPairs = numeric_matrix)"
      result$details$hc_pairs <- list(
        type = typeof(hc_pairs),
        numeric = is.numeric(hc_pairs),
        inherits_hc = inherits(hc_pairs, "hc"),
        dim = dim(hc_pairs)
      )
      result
    },
    user = {
      supplied <- .numeric_init_extract_supplied(user_init, "user_init")
      list(
        state = supplied$state,
        partition = supplied$partition,
        details = list(
          partition_method = "user",
          source_class = supplied$source_class,
          source_path = supplied$source_path,
          source_schema_version = supplied$source_schema_version
        )
      )
    },
    previous_fit = {
      supplied <- .numeric_init_extract_supplied(previous_fit, "previous_fit")
      list(
        state = supplied$state,
        partition = supplied$partition,
        details = list(
          partition_method = "previous_fit",
          source_class = supplied$source_class,
          source_path = supplied$source_path,
          source_schema_version = supplied$source_schema_version
        )
      )
    },
    stop("Unsupported single initialization method.", call. = FALSE)
  )

  validated <- .numeric_init_validate_state(
    result$state, data, K, result$partition,
    source = sprintf("%s initialization", method)
  )
  objective <- .numeric_init_loglik(data, validated$state)
  backend_loglik <- result$details$backend_loglik
  if (!is.null(backend_loglik) &&
      abs(objective - backend_loglik) > 1e-6 * (1 + abs(backend_loglik))) {
    stop(
      "The independently recomputed objective disagrees with Mclust loglik.",
      call. = FALSE
    )
  }

  list(
    state = validated$state,
    partition = validated$partition,
    objective = objective,
    converged = TRUE,
    convergence_reason = if (identical(method, "mclust_em")) {
      "mclust_returned_valid_fit"
    } else {
      "non_iterative_state_constructed"
    },
    details = c(result$details, list(validation = validated$validation))
  )
}

.numeric_init_candidate_row <- function(candidate_id, method, seed,
                                        success, converged, objective,
                                        failure_message = NA_character_) {
  data.frame(
    candidate_id = as.integer(candidate_id),
    method = as.character(method),
    seed = as.integer(seed),
    success = isTRUE(success),
    converged = isTRUE(converged),
    objective = as.double(objective),
    failure_message = as.character(failure_message),
    stringsAsFactors = FALSE
  )
}

.numeric_init_build_adapter <- function(candidate, method, requested_method,
                                        seed, seed_source,
                                        candidate_provenance,
                                        multistart_provenance) {
  state <- .numeric_init_clone_state(candidate$state)
  out <- c(
    state,
    list(
      partition = as.integer(candidate$partition),
      seed = if (is.null(seed)) NA_integer_ else as.integer(seed),
      seed_source = as.character(seed_source),
      method = as.character(method),
      requested_method = as.character(requested_method),
      objective = as.double(candidate$objective),
      objective_name = "gaussian_mixture_observed_loglik",
      objective_direction = "maximize",
      converged = isTRUE(candidate$converged),
      convergence_reason = as.character(candidate$convergence_reason),
      failure = NULL,
      candidate_provenance = candidate_provenance,
      multistart_provenance = multistart_provenance,
      selected_candidate_details = candidate$details,
      schema_version = .numeric_initialization_schema_version
    )
  )
  structure(out, class = c("selvarmix_numeric_initialization", "list"))
}

.numeric_init_derived_seeds <- function(base_seed, count) {
  modulus <- as.double(.Machine$integer.max) + 1
  as.integer((as.double(base_seed) + seq_len(count) - 1) %% modulus)
}

.numeric_initialization_adapter <- function(
    data,
    nbClust,
    init,
    n.start,
    lambda_omega_0,
    epsilon_pd,
    seed,
    user_init,
    previous_fit,
    multistart_methods,
    multistart_replicates,
    multistart_seeds,
    mclust_model) {
  checked <- .numeric_init_validate_data(data, nbClust)
  data <- checked$data
  K <- checked$K
  method_info <- .numeric_init_canonical_method(init)
  method <- method_info$canonical
  seed <- .numeric_init_seed(seed)
  n.start <- .numeric_init_integer(n.start, "n.start")
  if (length(lambda_omega_0) != 1L || !is.numeric(lambda_omega_0) ||
      !is.finite(lambda_omega_0) || lambda_omega_0 < 0) {
    stop("lambda_omega_0 must be one finite non-negative number.",
         call. = FALSE)
  }
  if (length(epsilon_pd) != 1L || !is.numeric(epsilon_pd) ||
      !is.finite(epsilon_pd) || epsilon_pd <= 0) {
    stop("epsilon_pd must be one finite positive number.", call. = FALSE)
  }
  if (!is.character(mclust_model) || length(mclust_model) != 1L ||
      is.na(mclust_model) || !nzchar(mclust_model)) {
    stop("mclust_model must be one non-empty model name.", call. = FALSE)
  }

  if (!identical(method, "deterministic_multistart")) {
    candidate_seed <- seed
    seed_source <- if (is.null(seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        "caller_rng_state"
      } else {
        "system_initialized_rng"
      }
    } else {
      "explicit"
    }
    candidate <- tryCatch(
      .numeric_init_with_rng(
        candidate_seed,
        .numeric_init_run_candidate(
          method, data, K, n.start, lambda_omega_0, epsilon_pd,
          user_init, previous_fit, mclust_model
        )
      ),
      error = identity
    )
    if (inherits(candidate, "error")) {
      failed_row <- .numeric_init_candidate_row(
        1L, method,
        if (is.null(candidate_seed)) NA_integer_ else candidate_seed,
        FALSE, FALSE, NA_real_, conditionMessage(candidate)
      )
      .numeric_init_abort(
        sprintf("%s initialization failed: %s", method,
                conditionMessage(candidate)),
        method = method,
        seed = if (is.null(candidate_seed)) NA_integer_ else candidate_seed,
        code = "candidate_failed",
        candidate_provenance = failed_row,
        parent = candidate
      )
    }
    candidate_row <- .numeric_init_candidate_row(
      1L, method,
      if (is.null(candidate_seed)) NA_integer_ else candidate_seed,
      TRUE, candidate$converged, candidate$objective
    )
    return(.numeric_init_build_adapter(
      candidate = candidate,
      method = method,
      requested_method = method_info$requested,
      seed = candidate_seed,
      seed_source = seed_source,
      candidate_provenance = candidate_row,
      multistart_provenance = list(
        enabled = FALSE,
        selected_candidate_id = 1L,
        selected_method = method,
        requested_methods = method,
        replicates = 1L,
        tie_break = "first_candidate_in_declared_order"
      )
    ))
  }

  multistart_replicates <- .numeric_init_integer(
    multistart_replicates, "multistart_replicates"
  )
  if (!is.character(multistart_methods) || !length(multistart_methods) ||
      anyNA(multistart_methods)) {
    stop("multistart_methods must contain method tokens.", call. = FALSE)
  }
  canonical_methods <- vapply(multistart_methods, function(value) {
    .numeric_init_canonical_method(value)$canonical
  }, character(1))
  unsupported <- setdiff(
    unique(canonical_methods),
    c("mclust_hc", "mclust_em", "kmeans", "random")
  )
  if (length(unsupported)) {
    stop(
      paste0(
        "deterministic_multistart supports only mclust_hc, mclust_em, ",
        "kmeans, and random candidates."
      ),
      call. = FALSE
    )
  }
  plan <- expand.grid(
    replicate = seq_len(multistart_replicates),
    method_index = seq_along(canonical_methods),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  plan <- plan[order(plan$method_index, plan$replicate), , drop = FALSE]
  plan$method <- canonical_methods[plan$method_index]
  candidate_count <- nrow(plan)
  if (is.null(multistart_seeds)) {
    base_seed <- if (is.null(seed)) 1L else seed
    candidate_seeds <- .numeric_init_derived_seeds(base_seed, candidate_count)
    seed_source <- if (is.null(seed)) "deterministic_method_default" else "explicit"
  } else {
    if (!is.numeric(multistart_seeds) ||
        length(multistart_seeds) != candidate_count ||
        any(!is.finite(multistart_seeds)) ||
        any(multistart_seeds != round(multistart_seeds)) ||
        any(multistart_seeds < 0) ||
        any(multistart_seeds > .Machine$integer.max)) {
      stop(
        "multistart_seeds must contain one valid integer per candidate.",
        call. = FALSE
      )
    }
    candidate_seeds <- as.integer(multistart_seeds)
    base_seed <- if (is.null(seed)) candidate_seeds[1L] else seed
    seed_source <- "explicit_candidate_seeds"
  }

  candidates <- vector("list", candidate_count)
  candidate_details <- vector("list", candidate_count)
  candidate_rows <- vector("list", candidate_count)
  for (index in seq_len(candidate_count)) {
    candidate <- tryCatch(
      .numeric_init_with_rng(
        candidate_seeds[index],
        .numeric_init_run_candidate(
          plan$method[index], data, K, n.start, lambda_omega_0,
          epsilon_pd, NULL, NULL, mclust_model
        )
      ),
      error = identity
    )
    if (inherits(candidate, "error")) {
      candidate_rows[[index]] <- .numeric_init_candidate_row(
        index, plan$method[index], candidate_seeds[index],
        FALSE, FALSE, NA_real_, conditionMessage(candidate)
      )
      candidate_details[[index]] <- list(
        failure = list(
          message = conditionMessage(candidate),
          class = class(candidate)
        )
      )
    } else {
      candidates[[index]] <- candidate
      candidate_rows[[index]] <- .numeric_init_candidate_row(
        index, plan$method[index], candidate_seeds[index],
        TRUE, candidate$converged, candidate$objective
      )
      candidate_details[[index]] <- candidate$details
    }
  }
  candidate_provenance <- do.call(rbind, candidate_rows)
  valid <- candidate_provenance$success & candidate_provenance$converged &
    is.finite(candidate_provenance$objective)
  if (!any(valid)) {
    .numeric_init_abort(
      "deterministic_multistart produced no valid converged candidate.",
      method = method,
      seed = base_seed,
      code = "all_candidates_failed",
      candidate_provenance = candidate_provenance
    )
  }
  valid_indices <- which(valid)
  selected <- valid_indices[which.max(
    candidate_provenance$objective[valid_indices]
  )]
  selected_candidate <- candidates[[selected]]

  .numeric_init_build_adapter(
    candidate = selected_candidate,
    method = method,
    requested_method = method_info$requested,
    seed = base_seed,
    seed_source = seed_source,
    candidate_provenance = candidate_provenance,
    multistart_provenance = list(
      enabled = TRUE,
      selected_candidate_id = as.integer(selected),
      selected_method = plan$method[selected],
      requested_methods = multistart_methods,
      canonical_methods = canonical_methods,
      replicates = multistart_replicates,
      candidate_seeds = candidate_seeds,
      tie_break = "first_candidate_in_declared_order",
      candidate_details = candidate_details
    )
  )
}

.initializer_state_projection <- function(adapter) {
  extracted <- .numeric_init_extract_supplied(adapter, "adapter")
  K <- length(extracted$state$prop)
  n <- nrow(extracted$state$Z)
  p <- nrow(extracted$state$Mu)
  if (!is.finite(K) || K < 1L || !is.finite(n) || n < 1L ||
      !is.finite(p) || p < 1L) {
    stop("adapter does not contain inferable state dimensions.", call. = FALSE)
  }
  dummy_data <- matrix(0, nrow = n, ncol = p)
  validated <- .numeric_init_validate_state(
    extracted$state, dummy_data, K, extracted$partition,
    source = "adapter"
  )
  .numeric_init_clone_state(validated$state)
}

.build_initializer_registry <- function(
    data,
    nbcluster,
    init_control = NULL,
    legacy_method = "hc",
    legacy_n_start = 250,
    legacy_lambda_omega_0 = 50) {
  if (is.null(init_control)) init_control <- list()
  if (!is.list(init_control) || is.null(names(init_control)) && length(init_control)) {
    stop("init_control must be NULL or a named list.", call. = FALSE)
  }
  if (anyDuplicated(names(init_control))) {
    stop("init_control names must be unique.", call. = FALSE)
  }
  allowed <- c(
    "method", "seed", "n_start", "lambda_omega_0", "epsilon_pd",
    "user_init", "previous_fit", "multistart_methods",
    "multistart_replicates", "multistart_seeds", "mclust_model", "per_k"
  )
  unknown <- setdiff(names(init_control), allowed)
  if (length(unknown)) {
    stop(sprintf("Unknown init_control field(s): %s.",
                 paste(unknown, collapse = ", ")), call. = FALSE)
  }

  checked_data <- .numeric_init_validate_data(data, 1L)$data
  if (!is.numeric(nbcluster) || !length(nbcluster) ||
      any(!is.finite(nbcluster)) || any(nbcluster != round(nbcluster)) ||
      any(nbcluster < 1) || any(nbcluster > nrow(checked_data))) {
    stop("nbcluster must contain valid component counts.", call. = FALSE)
  }
  nbcluster <- as.integer(nbcluster)
  if (anyDuplicated(nbcluster)) {
    stop("nbcluster values must be unique in an initializer registry.",
         call. = FALSE)
  }

  per_k <- init_control$per_k
  init_control$per_k <- NULL
  per_k_map <- list()
  if (!is.null(per_k)) {
    if (!is.list(per_k) || !length(per_k) || is.null(names(per_k)) ||
        any(!nzchar(names(per_k))) || anyDuplicated(names(per_k))) {
      stop("init_control$per_k must be a uniquely named non-empty list.",
           call. = FALSE)
    }
    parsed <- suppressWarnings(as.integer(sub("^[Kk]", "", names(per_k))))
    if (anyNA(parsed) || any(!(parsed %in% nbcluster)) || anyDuplicated(parsed)) {
      stop("init_control$per_k keys must uniquely match requested K values.",
           call. = FALSE)
    }
    for (index in seq_along(per_k)) {
      override <- per_k[[index]]
      if (!is.list(override) || is.null(names(override)) && length(override) ||
          anyDuplicated(names(override))) {
        stop("Each per-K initializer override must be a uniquely named list.",
             call. = FALSE)
      }
      unknown_override <- setdiff(names(override), setdiff(allowed, "per_k"))
      if (length(unknown_override)) {
        stop(sprintf("Unknown per-K init_control field(s): %s.",
                     paste(unknown_override, collapse = ", ")), call. = FALSE)
      }
      per_k_map[[paste0("K", parsed[index])]] <- override
    }
  }

  defaults <- list(
    method = legacy_method,
    n_start = legacy_n_start,
    lambda_omega_0 = legacy_lambda_omega_0,
    epsilon_pd = sqrt(.Machine$double.eps),
    seed = NULL,
    user_init = NULL,
    previous_fit = NULL,
    multistart_methods = c("mclust_hc", "mclust_em", "kmeans", "random"),
    multistart_replicates = 1L,
    multistart_seeds = NULL,
    mclust_model = "VVV"
  )
  global <- utils::modifyList(defaults, init_control, keep.null = TRUE)
  registry <- vector("list", length(nbcluster))
  names(registry) <- paste0("K", nbcluster)
  for (index in seq_along(nbcluster)) {
    key <- names(registry)[index]
    controls <- utils::modifyList(
      global,
      if (is.null(per_k_map[[key]])) list() else per_k_map[[key]],
      keep.null = TRUE
    )
    registry[[index]] <- InitParameter(
      data = checked_data,
      nbClust = nbcluster[index],
      init = controls$method,
      n.start = controls$n_start,
      lambda_omega_0 = controls$lambda_omega_0,
      epsilon_pd = controls$epsilon_pd,
      seed = controls$seed,
      user_init = controls$user_init,
      previous_fit = controls$previous_fit,
      multistart_methods = controls$multistart_methods,
      multistart_replicates = controls$multistart_replicates,
      multistart_seeds = controls$multistart_seeds,
      mclust_model = controls$mclust_model,
      return_adapter = TRUE
    )
  }
  structure(
    registry,
    class = c("selvarmix_initializer_registry", "list"),
    schema_version = .numeric_initialization_schema_version,
    K = nbcluster
  )
}

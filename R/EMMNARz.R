.mnar_abort <- function(message, subclass = "selvarmix_mnar_fit_error") {
  condition <- structure(
    list(message = message, call = NULL),
    class = unique(c(subclass, "selvarmix_mnar_fit_error", "error", "condition"))
  )
  stop(condition)
}

.mnar_warn <- function(message, subclass) {
  condition <- structure(
    list(message = message, call = NULL),
    class = c(subclass, "warning", "condition")
  )
  warning(condition)
}

.mnar_validate_gaussian_design <- function(x, diag) {
  observed_count <- colSums(!is.na(x))
  observed_sd <- vapply(seq_len(ncol(x)), function(j) {
    stats::sd(x[, j], na.rm = TRUE)
  }, numeric(1))
  observed_scale <- vapply(seq_len(ncol(x)), function(j) {
    values <- x[!is.na(x[, j]), j]
    if (!length(values)) .Machine$double.xmin else
      max(.Machine$double.xmin, abs(values))
  }, numeric(1))
  invalid <- observed_count < 2L | !is.finite(observed_sd) |
    observed_sd <= sqrt(.Machine$double.eps) * observed_scale
  if (any(invalid)) {
    .mnar_abort(
      paste0(
        "each Gaussian coordinate requires at least two distinct observed ",
        "values; invalid column(s): ", paste(which(invalid), collapse = ", ")
      ),
      "selvarmix_mnar_input_error"
    )
  }

  if (!isTRUE(diag) && ncol(x) > 1L) {
    joint_count <- crossprod(!is.na(x))
    unidentified <- which(
      upper.tri(joint_count) & joint_count < 2L,
      arr.ind = TRUE
    )
    if (nrow(unidentified)) {
      labels <- apply(unidentified, 1L, function(pair) {
        paste0("(", pair[[1L]], ",", pair[[2L]], ")")
      })
      .mnar_abort(
        paste0(
          "a full Gaussian covariance requires at least two jointly observed ",
          "rows for every variable pair; unidentified pair(s): ",
          paste(labels, collapse = ", ")
        ),
        "selvarmix_mnar_input_error"
      )
    }
  }
  invisible(TRUE)
}

.mnar_parameter_count <- function(K, d, diag = TRUE,
                                  mecha = c("MNARz", "mixed"),
                                  is_mnar = NULL) {
  mecha <- match.arg(mecha)
  if (length(K) != 1L || !is.finite(K) || K < 1 || K != as.integer(K)) {
    stop("K must be a positive integer")
  }
  if (length(d) != 1L || !is.finite(d) || d < 1 || d != as.integer(d)) {
    stop("d must be a positive integer")
  }
  if (length(diag) != 1L || is.na(diag) || !is.logical(diag)) {
    stop("diag must be TRUE or FALSE")
  }

  covariance_count <- if (diag) K * d else K * d * (d + 1) / 2
  missingness_count <- if (mecha == "MNARz") {
    # Pure MNARz: one alpha_k per component, repeated across variables in
    # storage. The distinct MNARzj model would contribute K*d here.
    K
  } else {
    if (is.null(is_mnar) || length(is_mnar) != d || anyNA(is_mnar)) {
      stop("is_mnar must be a non-missing logical vector of length d")
    }
    is_mnar <- as.logical(is_mnar)
    # Ignorable mask factors are not parameterized; all selected MNARz
    # coordinates share one component-specific probability.
    if (any(is_mnar)) K else 0
  }

  as.numeric((K - 1) + K * d + covariance_count + missingness_count)
}

.mnar_initial_rho <- function(x, K, responsibilities = NULL,
                              coordinates = seq_len(ncol(x))) {
  mask <- is.na(x[, coordinates, drop = FALSE])
  overall <- mean(mask)
  if (!is.finite(overall) || overall <= 0 || overall >= 1) {
    .mnar_abort(
      "the modeled MNARz mask is on the all-observed or all-missing boundary",
      "selvarmix_mnar_boundary_error"
    )
  }
  if (is.null(responsibilities)) return(rep(overall, K))
  responsibilities <- as.matrix(responsibilities)
  if (!identical(dim(responsibilities), c(nrow(x), as.integer(K))) ||
      any(!is.finite(responsibilities)) || any(responsibilities < 0) ||
      any(abs(rowSums(responsibilities) - 1) > 1e-8) ||
      any(colSums(responsibilities) <= 0)) {
    .mnar_abort(
      "initializer responsibilities are invalid for class-specific MNARz starts",
      "selvarmix_mnar_input_error"
    )
  }
  missing_per_row <- rowSums(mask)
  trials <- colSums(responsibilities) * ncol(mask)
  successes <- as.numeric(crossprod(responsibilities, missing_per_row))
  # Jeffreys smoothing is used only for the initial state.  It keeps a class
  # with no initial missing entries finite without clipping any EM update.
  as.numeric((successes + 0.5) / (trials + 1))
}

.mnar_default_gaussian_init <- function(x, K, diag) {
  observed_count <- colSums(!is.na(x))
  observed_sd <- vapply(seq_len(ncol(x)), function(j) {
    stats::sd(x[, j], na.rm = TRUE)
  }, numeric(1))
  invalid <- observed_count < 2L | !is.finite(observed_sd) | observed_sd <= 0
  if (any(invalid)) {
    .mnar_abort(
      paste0(
        "automatic MNARz initialization requires at least two distinct ",
        "observations per variable; invalid column(s): ",
        paste(which(invalid), collapse = ", ")
      ),
      "selvarmix_mnar_input_error"
    )
  }
  completed <- x
  for (j in seq_len(ncol(x))) {
    completed[is.na(completed[, j]), j] <- mean(x[, j], na.rm = TRUE)
  }
  standardized <- sweep(completed, 2L, colMeans(completed), FUN = "-")
  completed_sd <- apply(standardized, 2L, stats::sd)
  completed_sd[!is.finite(completed_sd) | completed_sd <= 0] <- observed_sd[
    !is.finite(completed_sd) | completed_sd <= 0
  ]
  standardized <- sweep(standardized, 2L, completed_sd, FUN = "/")

  minimum_mass <- if (isTRUE(diag)) 1L else max(1L, ncol(x))
  if (nrow(x) < K * (minimum_mass + 1L)) {
    .mnar_abort(
      paste0(
        "automatic initialization cannot form K components larger than the ",
        "dimension-aware component floor (", minimum_mass, ")"
      ),
      "selvarmix_mnar_component_error"
    )
  }
  if (K == 1L) {
    labels <- rep(1L, nrow(x))
  } else {
    labels <- stats::cutree(
      stats::hclust(stats::dist(standardized), method = "ward.D2"),
      k = K
    )
    if (any(tabulate(labels, nbins = K) <= minimum_mass)) {
      principal_score <- tryCatch(
        as.numeric(standardized %*% svd(standardized, nu = 0L, nv = 1L)$v[, 1L]),
        error = function(e) rowSums(standardized)
      )
      order_index <- order(principal_score, seq_along(principal_score))
      sizes <- rep(nrow(x) %/% K, K)
      sizes[seq_len(nrow(x) %% K)] <- sizes[seq_len(nrow(x) %% K)] + 1L
      labels <- integer(nrow(x))
      labels[order_index] <- rep(seq_len(K), times = sizes)
    }
  }
  initialized <- .numeric_init_covariances_from_partition(
    completed, as.integer(labels), as.integer(K),
    lambda_omega_0 = 1, epsilon_pd = sqrt(.Machine$double.eps)
  )
  initialized$state
}

.mnar_adapt_gaussian_init <- function(init, K, d, x, mecha, is_mnar = NULL) {
  if (is.null(init)) return(NULL)
  if (!is.list(init)) return(init)

  responsibilities <- init$Z
  if (is.null(responsibilities) && !is.null(init$tik)) {
    responsibilities <- init$tik
  }
  if (is.null(responsibilities) && !is.null(init$partition)) {
    partition <- as.integer(init$partition)
    if (length(partition) == nrow(x) && !anyNA(partition) &&
        all(partition >= 1L & partition <= K)) {
      responsibilities <- matrix(0, nrow(x), K)
      responsibilities[cbind(seq_len(nrow(x)), partition)] <- 1
    }
  }

  legacy_names <- c("pik", "mu", "sigma")
  if (all(legacy_names %in% names(init))) {
    adapted <- init
  } else {
    if (!all(c("prop", "Mu") %in% names(init))) return(init)
    sigma_cube <- init$SigmaCube
    if (is.null(sigma_cube)) sigma_cube <- init$Sigma
    if (is.null(sigma_cube)) return(init)

    prop <- as.numeric(init$prop)
    Mu <- as.matrix(init$Mu)
    sigma_dim <- dim(sigma_cube)
    if (length(prop) != K || !identical(dim(Mu), c(d, as.integer(K))) ||
        length(sigma_dim) != 3L ||
        !identical(as.integer(sigma_dim), c(d, d, as.integer(K)))) {
      .mnar_abort(
        paste0(
          "initializer fields prop, Mu, and SigmaCube/Sigma must have dimensions ",
          "K, d by K, and d by d by K"
        ),
        "selvarmix_mnar_input_error"
      )
    }

    adapted <- list(
      pik = prop,
      mu = lapply(seq_len(K), function(k) as.numeric(Mu[, k])),
      sigma = lapply(seq_len(K), function(k) {
        matrix(sigma_cube[, , k], nrow = d, ncol = d)
      })
    )
    if (!is.null(init$rho)) adapted$rho <- as.numeric(init$rho)
    if (!is.null(init$alpha)) adapted$alpha <- as.matrix(init$alpha)
    if (!is.null(init$beta)) adapted$beta <- as.matrix(init$beta)
  }

  if (identical(mecha, "MNARz") && is.null(adapted$alpha)) {
    if (!is.null(adapted$rho)) {
      if (length(adapted$rho) != K || any(!is.finite(adapted$rho)) ||
          any(adapted$rho <= 0 | adapted$rho >= 1)) {
        .mnar_abort(
          "initializer rho must have length K with values strictly between zero and one",
          "selvarmix_mnar_input_error"
        )
      }
      adapted$alpha <- matrix(stats::qnorm(adapted$rho), nrow = K, ncol = d)
    } else {
      probability <- .mnar_initial_rho(
        x, K, responsibilities = responsibilities
      )
      adapted$rho <- probability
      adapted$alpha <- matrix(stats::qnorm(probability), nrow = K, ncol = d)
    }
  }
  if (identical(mecha, "mixed") &&
      is.null(adapted$alpha) && is.null(adapted$rho)) {
    if (any(is_mnar)) {
      adapted$rho <- .mnar_initial_rho(
        x, K, responsibilities = responsibilities,
        coordinates = which(is_mnar)
      )
    } else {
      adapted$alpha <- matrix(NA_real_, nrow = K, ncol = d)
      adapted$rho <- rep(NA_real_, K)
    }
  }
  adapted
}

.mnar_native_error <- function(error) {
  message <- conditionMessage(error)
  subclass <- if (grepl("Effective component size", message, fixed = TRUE)) {
    "selvarmix_mnar_component_error"
  } else if (grepl("log-likelihood decreased materially", message, fixed = TRUE)) {
    "selvarmix_mnar_monotonicity_error"
  } else if ((grepl("no finite", message, fixed = TRUE) &&
              grepl("intercept", message, fixed = TRUE)) ||
             grepl("boundary", message, fixed = TRUE)) {
    "selvarmix_mnar_boundary_error"
  } else if (grepl("Initial ", message, fixed = TRUE) ||
             grepl("must be", message, fixed = TRUE) ||
             grepl("Unsupported missingness", message, fixed = TRUE)) {
    "selvarmix_mnar_input_error"
  } else {
    "selvarmix_mnar_numerical_error"
  }
  .mnar_abort(paste0("EMGaussian failed: ", message), subclass)
}

EMClustMNARz <- function(x,
                         K,
                         mecha     = "MNARz",
                         criterion = "BIC",
                         diag      = TRUE,
                         rmax      = 100,
                         init      = NULL,
                         tol       = 1e-4,
                         is_mnar   = NULL) {
  x <- as.matrix(x)
  if (!is.numeric(x) || !length(x) || nrow(x) < 1L || ncol(x) < 1L) {
    .mnar_abort("x must be a non-empty numeric matrix", "selvarmix_mnar_input_error")
  }
  if (any(is.nan(x)) || any(is.infinite(x))) {
    .mnar_abort(
      "x may contain NA values but not NaN or infinite values",
      "selvarmix_mnar_input_error"
    )
  }
  if (length(K) != 1L || !is.numeric(K) || !is.finite(K) ||
      K < 1 || K != as.integer(K)) {
    .mnar_abort("K must be a positive integer", "selvarmix_mnar_input_error")
  }
  if (K > nrow(x)) {
    .mnar_abort(
      "K cannot exceed the number of observations",
      "selvarmix_mnar_input_error"
    )
  }
  all_missing <- which(colSums(!is.na(x)) == 0L)
  if (length(all_missing)) {
    .mnar_abort(
      paste0(
        "all-missing column(s) are unidentified under the Gaussian MNARz model: ",
        paste(all_missing, collapse = ", ")
      ),
      "selvarmix_mnar_input_error"
    )
  }
  if (length(diag) != 1L || !is.logical(diag) || is.na(diag)) {
    .mnar_abort("diag must be TRUE or FALSE", "selvarmix_mnar_input_error")
  }
  .mnar_validate_gaussian_design(x, diag)
  if (length(rmax) != 1L || !is.numeric(rmax) || !is.finite(rmax) ||
      rmax < 1 || rmax != as.integer(rmax)) {
    .mnar_abort("rmax must be a positive integer", "selvarmix_mnar_input_error")
  }
  if (length(tol) != 1L || !is.numeric(tol) || !is.finite(tol) || tol < 0) {
    .mnar_abort("tol must be finite and non-negative", "selvarmix_mnar_input_error")
  }
  if (!is.null(init) && !is.list(init)) {
    .mnar_abort("init must be NULL or a list", "selvarmix_mnar_input_error")
  }
  if (length(criterion) != 1L || !is.character(criterion) || is.na(criterion)) {
    .mnar_abort("criterion must be one non-missing string", "selvarmix_mnar_input_error")
  }
  if (!criterion %in% c("BIC", "ICL")) {
    .mnar_abort(
      "criterion must be exactly one of 'BIC' or 'ICL'",
      "selvarmix_mnar_input_error"
    )
  }

  if (length(mecha) != 1L || !is.character(mecha) || is.na(mecha) ||
      !is.null(attributes(mecha))) {
    .mnar_abort(
      "mecha must be one non-missing character string without attributes",
      "selvarmix_mnar_input_error"
    )
  }
  if (!mecha %in% c("MNARz", "mixed")) {
    .mnar_abort(
      paste0("Unknown or unsupported mechanism '", mecha, "'"),
      "selvarmix_mnar_input_error"
    )
  }

  n <- nrow(x)
  d <- ncol(x)
  if (identical(mecha, "mixed")) {
    if (is.null(is_mnar) || !is.logical(is_mnar) || length(is_mnar) != d ||
        anyNA(is_mnar)) {
      .mnar_abort(
        "is_mnar must be a non-missing logical vector of length ncol(x) for mecha = 'mixed'",
        "selvarmix_mnar_input_error"
      )
    }
    is_mnar <- as.logical(is_mnar)
  }
  automatic_initialization <- is.null(init)
  if (automatic_initialization) {
    init <- .mnar_default_gaussian_init(x, as.integer(K), diag)
  }
  init <- .mnar_adapt_gaussian_init(
    init, K = as.integer(K), d = d, x = x,
    mecha = mecha, is_mnar = is_mnar
  )

  result <- tryCatch(
    if (identical(mecha, "mixed")) {
      EMGaussianMixed(x, K, mecha, is_mnar, diag, rmax, init, tol)
    } else {
      EMGaussian(x, K, mecha, diag, rmax, init, tol)
    },
    error = .mnar_native_error
  )

  required_status <- c(
    "iterations", "converged", "termination_reason",
    "final_loglik_improvement", "loglik_monotone",
    "no_material_loglik_decrease",
    "min_effective_component_size", "component_floor",
    "covariance_adjustments"
  )
  if (!all(required_status %in% names(result))) {
    .mnar_abort(
      "EMGaussian did not return the required convergence diagnostics",
      "selvarmix_mnar_result_error"
    )
  }

  loglikelihood <- as.numeric(result$loglik_vec)
  if (!length(loglikelihood) || any(!is.finite(loglikelihood))) {
    .mnar_abort(
      "EMGaussian returned a non-finite observed log-likelihood",
      "selvarmix_mnar_numerical_error"
    )
  }
  loglik_final <- loglikelihood[[length(loglikelihood)]]

  tik <- as.matrix(result$tik)
  if (!identical(dim(tik), c(n, as.integer(K))) || any(!is.finite(tik)) ||
      any(tik < 0) || any(abs(rowSums(tik) - 1) > 1e-8)) {
    .mnar_abort(
      "EMGaussian returned invalid posterior probabilities",
      "selvarmix_mnar_result_error"
    )
  }

  parameter_count <- .mnar_parameter_count(
    K, d, diag, mecha = mecha, is_mnar = is_mnar
  )
  missingness_parameter_count <- if (identical(mecha, "mixed")) {
    if (any(is_mnar)) as.numeric(K) else 0
  } else {
    as.numeric(K)
  }
  converged <- isTRUE(result$converged)
  criterion_available <- converged && isTRUE(result$no_material_loglik_decrease)

  entropy <- if (any(tik > 0)) {
    -sum(tik[tik > 0] * log(tik[tik > 0]))
  } else {
    NA_real_
  }
  if (criterion_available) {
    bic <- -2 * loglik_final + parameter_count * log(n)
    icl <- bic + 2 * entropy
  } else {
    bic <- NA_real_
    icl <- NA_real_
    .mnar_warn(
      paste0(
        if (identical(mecha, "mixed")) "Mixed EM" else "MNARz EM",
        " stopped with termination_reason = '",
        result$termination_reason,
        "'; BIC and ICL are unavailable for this non-converged fit."
      ),
      "selvarmix_mnar_nonconvergence_warning"
    )
  }

  crit_list <- list(BIC = bic, ICL = icl)

  diagnostics <- list(
    iterations = as.integer(result$iterations),
    converged = converged,
    termination_reason = as.character(result$termination_reason),
    loglik_monotone = isTRUE(result$loglik_monotone),
    no_material_loglik_decrease = isTRUE(result$no_material_loglik_decrease),
    final_loglik_improvement = as.numeric(result$final_loglik_improvement),
    min_effective_component_size = as.numeric(result$min_effective_component_size),
    component_floor = as.numeric(result$component_floor),
    covariance_adjustments = as.integer(result$covariance_adjustments),
    parameter_count = parameter_count,
    missingness_parameter_count = missingness_parameter_count,
    criterion_available = criterion_available,
    initialization = if (automatic_initialization) {
      "deterministic_data_and_mask_aware"
    } else {
      "supplied"
    },
    criterion_scope = paste0(
      if (identical(mecha, "mixed")) {
        paste0(
          "Schwarz-type criterion for the modeled observed likelihood: Gaussian ",
          "Y_obs plus class-only MNARz masks on is_mnar = TRUE coordinates; ",
          "the unrestricted ignorable-MAR mask factor is omitted and contributes ",
          "no counted parameter. "
        )
      } else {
        "Schwarz-type criterion for the joint observed likelihood of (Y_obs, C); "
      },
      "the parameter count is the nominal dimension of the returned local fit.",
      if (result$covariance_adjustments > 0L) {
        paste0(
          " This fit used a machine-scale covariance eigenvalue floor; ",
          "the criterion evaluates the returned covariance-floored fit."
        )
      } else {
        ""
      }
    )
  )
  optional_native_diagnostics <- c(
    "absolute_tolerance", "relative_tolerance", "decrease_tolerance",
    "component_mass_rule", "covariance_model_dimension",
    "missing_pattern_count", "sufficient_statistic_storage",
    "conditional_moment_storage", "final_relative_loglik_improvement",
    "convergence_threshold"
  )
  for (name in optional_native_diagnostics) {
    if (!is.null(result[[name]])) diagnostics[[name]] <- result[[name]]
  }

  parameters <- list(
    pik = result$pik,
    mu = result$mu,
    sigma = result$sigma,
    alpha = result$alpha
  )
  if (identical(mecha, "mixed")) {
    parameters$rho <- result$rho
    parameters$beta <- result$beta
    parameters$is_mnar <- is_mnar
    parameters$mechanism_specification <- result$mechanism_specification
    diagnostics$mnar_coordinate_count <- as.integer(sum(is_mnar))
    diagnostics$ignored_mar_coordinate_count <- as.integer(sum(!is_mnar))
    diagnostics$alpha_scope <- result$alpha_scope
    diagnostics$beta_deprecated <- isTRUE(result$beta_deprecated)
  }

  list(
    loglik_obs = loglik_final,
    partition = max.col(tik, ties.method = "first"),
    imputedData = result$imputedData,
    criterionValue = crit_list,
    parameters = parameters,
    proba = tik,
    diagnostics = diagnostics
  )
}

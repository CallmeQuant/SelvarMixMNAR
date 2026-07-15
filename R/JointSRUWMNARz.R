# Fixed-role SRUW--MNARz estimation combines the Gaussian SRUW factorization
# with one class-specific missingness probability shared across coordinates.
# Variable ranking is external to this likelihood.

.joint_mnar_abort <- function(message,
                              subclass = "selvarmix_joint_mnar_fit_error",
                              data = list()) {
  condition <- structure(
    c(list(message = message, call = NULL), data),
    class = unique(c(
      subclass,
      "selvarmix_joint_mnar_fit_error",
      "error",
      "condition"
    ))
  )
  stop(condition)
}

.joint_deep_clone <- function(object) {
  unserialize(serialize(object, NULL))
}

.joint_check_scalar_integer <- function(x, name, minimum = 1L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x != as.integer(x) || x < minimum) {
    .joint_mnar_abort(
      paste0(name, " must be an integer greater than or equal to ", minimum),
      "selvarmix_joint_mnar_input_error"
    )
  }
  as.integer(x)
}

.joint_validate_roles <- function(roles, d) {
  if (!is.list(roles) || !all(c("S", "R", "U", "W") %in% names(roles))) {
    .joint_mnar_abort(
      "roles must be a list with S, R, U, and W entries",
      "selvarmix_joint_mnar_input_error"
    )
  }

  clean_role <- function(value, name, allow_empty = TRUE) {
    if (!length(value) && allow_empty) return(integer())
    if (!is.numeric(value) || anyNA(value) || any(!is.finite(value)) ||
        any(value != as.integer(value)) || any(value < 1L) || any(value > d) ||
        anyDuplicated(value)) {
      .joint_mnar_abort(
        paste0("roles$", name, " must contain distinct column indices in 1:d"),
        "selvarmix_joint_mnar_input_error"
      )
    }
    sort(as.integer(value))
  }

  S <- clean_role(roles$S, "S", allow_empty = FALSE)
  R <- clean_role(roles$R, "R")
  U <- clean_role(roles$U, "U")
  W <- clean_role(roles$W, "W")
  if (!length(S)) {
    .joint_mnar_abort(
      "roles$S must contain at least one clustering variable",
      "selvarmix_joint_mnar_input_error"
    )
  }
  if (!all(R %in% S)) {
    .joint_mnar_abort(
      "roles$R must be a subset of roles$S",
      "selvarmix_joint_mnar_input_error"
    )
  }
  partition <- c(S, U, W)
  if (length(partition) != d || anyDuplicated(partition) ||
      !identical(sort(partition), seq_len(d))) {
    .joint_mnar_abort(
      "roles$S, roles$U, and roles$W must partition all columns of x",
      "selvarmix_joint_mnar_input_error"
    )
  }

  list(S = S, R = R, U = U, W = W)
}

.joint_validate_observation_design <- function(x, roles, diagonal) {
  observed_count <- colSums(!is.na(x))
  observed_sd <- vapply(seq_len(ncol(x)), function(j) {
    stats::sd(x[, j], na.rm = TRUE)
  }, numeric(1))
  observed_scale <- vapply(seq_len(ncol(x)), function(j) {
    values <- x[!is.na(x[, j]), j]
    if (!length(values)) {
      .Machine$double.xmin
    } else {
      max(.Machine$double.xmin, abs(values))
    }
  }, numeric(1))
  invalid <- observed_count < 2L | !is.finite(observed_sd) |
    observed_sd <= sqrt(.Machine$double.eps) * observed_scale
  if (any(invalid)) {
    .joint_mnar_abort(
      paste0(
        "each modeled Gaussian coordinate requires at least two distinct ",
        "observed values; invalid column(s): ",
        paste(which(invalid), collapse = ", ")
      ),
      "selvarmix_joint_mnar_input_error"
    )
  }

  joint_count <- crossprod(!is.na(x))
  require_pairs <- function(left, right = left, within = TRUE, label) {
    if (!length(left) || !length(right)) return(invisible(NULL))
    candidates <- expand.grid(left = left, right = right)
    if (within) candidates <- candidates[candidates$left < candidates$right, ]
    if (!nrow(candidates)) return(invisible(NULL))
    count <- joint_count[cbind(candidates$left, candidates$right)]
    bad <- candidates[count < 2L, , drop = FALSE]
    if (nrow(bad)) {
      pairs <- paste0("(", bad$left, ",", bad$right, ")")
      .joint_mnar_abort(
        paste0(
          label, " requires at least two jointly observed rows for pair(s): ",
          paste(pairs, collapse = ", ")
        ),
        "selvarmix_joint_mnar_input_error"
      )
    }
    invisible(NULL)
  }

  if (!isTRUE(diagonal)) {
    require_pairs(roles$S, label = "the full S covariance")
    require_pairs(roles$U, label = "the full regression-residual covariance")
    require_pairs(roles$W, label = "the full W covariance")
  }
  require_pairs(
    roles$R, roles$U, within = FALSE,
    label = "the U-on-R regression"
  )
  invisible(TRUE)
}

.joint_stabilize_covariance <- function(sigma,
                                        diagonal,
                                        covariance_floor,
                                        label,
                                        subclass = "selvarmix_joint_mnar_numerical_error") {
  sigma <- as.matrix(sigma)
  p <- nrow(sigma)
  if (!identical(ncol(sigma), p) || any(!is.finite(sigma))) {
    .joint_mnar_abort(
      paste0(label, " must be a finite square matrix"),
      subclass
    )
  }
  if (!p) return(list(value = sigma, adjustments = 0L))

  sigma <- (sigma + t(sigma)) / 2
  scale <- max(.Machine$double.xmin, max(abs(sigma)))
  floor_value <- covariance_floor * scale
  negative_tolerance <- 64 * .Machine$double.eps * scale

  if (diagonal) {
    values <- base::diag(sigma)
    if (any(values < -negative_tolerance)) {
      .joint_mnar_abort(
        paste0(label, " has a materially negative diagonal entry"),
        subclass
      )
    }
    adjusted <- values < floor_value
    values[adjusted] <- floor_value
    return(list(
      value = base::diag(values, nrow = p, ncol = p),
      adjustments = as.integer(sum(adjusted))
    ))
  }

  decomposition <- eigen(sigma, symmetric = TRUE)
  if (any(decomposition$values < -negative_tolerance)) {
    .joint_mnar_abort(
      paste0(label, " is materially indefinite"),
      subclass
    )
  }
  adjusted <- decomposition$values < floor_value
  values <- pmax(decomposition$values, floor_value)
  stabilized <- decomposition$vectors %*%
    (values * t(decomposition$vectors))
  stabilized <- (stabilized + t(stabilized)) / 2
  list(value = stabilized, adjustments = as.integer(sum(adjusted)))
}

.joint_chol_solve <- function(a,
                              b,
                              label,
                              rank_check = FALSE,
                              covariance_floor = sqrt(.Machine$double.eps),
                              subclass = "selvarmix_joint_mnar_numerical_error") {
  a <- as.matrix(a)
  b_matrix <- as.matrix(b)
  if (!nrow(a)) return(b_matrix)
  if (!identical(dim(a), c(nrow(a), nrow(a))) || any(!is.finite(a)) ||
      nrow(b_matrix) != nrow(a) || any(!is.finite(b_matrix))) {
    .joint_mnar_abort(paste0(label, " has invalid dimensions or values"), subclass)
  }
  a <- (a + t(a)) / 2
  if (rank_check) {
    values <- eigen(a, symmetric = TRUE, only.values = TRUE)$values
    scale <- max(.Machine$double.xmin, max(abs(values)))
    if (min(values) <= covariance_floor * scale) {
      .joint_mnar_abort(
        paste0(label, " is rank deficient under the declared numerical tolerance"),
        "selvarmix_joint_mnar_regression_error",
        list(minimum_eigenvalue = min(values), numerical_scale = scale)
      )
    }
  }
  factor <- tryCatch(chol(a), error = identity)
  if (inherits(factor, "error")) {
    .joint_mnar_abort(
      paste0(label, " is not positive definite: ", conditionMessage(factor)),
      subclass
    )
  }
  solution <- backsolve(factor, forwardsolve(t(factor), b_matrix))
  if (is.null(dim(b))) as.numeric(solution) else solution
}

.joint_log_sum_exp <- function(x) {
  maximum <- max(x)
  if (!is.finite(maximum)) return(maximum)
  maximum + log(sum(exp(x - maximum)))
}

.joint_log_mvn_observed <- function(x, mean, sigma, label) {
  p <- length(x)
  if (!p) return(0)
  factor <- tryCatch(chol(sigma), error = identity)
  if (inherits(factor, "error")) {
    .joint_mnar_abort(
      paste0(label, " is not positive definite: ", conditionMessage(factor)),
      "selvarmix_joint_mnar_numerical_error"
    )
  }
  residual <- x - mean
  standardized <- forwardsolve(t(factor), residual)
  -0.5 * (p * log(2 * pi) + 2 * sum(log(base::diag(factor))) +
            sum(standardized^2))
}

.joint_global_parameters <- function(state, roles, d) {
  K <- length(state$pik)
  s <- length(roles$S)
  r <- length(roles$R)
  u <- length(roles$U)
  r_in_s <- match(roles$R, roles$S)
  means <- vector("list", K)
  covariances <- vector("list", K)

  for (k in seq_len(K)) {
    mean_k <- numeric(d)
    covariance_k <- matrix(0, nrow = d, ncol = d)
    mean_k[roles$S] <- state$mu_S[[k]]
    covariance_k[roles$S, roles$S] <- state$sigma_S[[k]]

    if (u) {
      mean_u <- state$a
      covariance_su <- matrix(0, nrow = s, ncol = u)
      covariance_uu <- state$Omega
      if (r) {
        mean_u <- mean_u + as.numeric(t(state$beta) %*%
                                        state$mu_S[[k]][r_in_s])
        covariance_su <- state$sigma_S[[k]][, r_in_s, drop = FALSE] %*%
          state$beta
        covariance_uu <- covariance_uu + t(state$beta) %*%
          state$sigma_S[[k]][r_in_s, r_in_s, drop = FALSE] %*%
          state$beta
      }
      mean_k[roles$U] <- mean_u
      covariance_k[roles$S, roles$U] <- covariance_su
      covariance_k[roles$U, roles$S] <- t(covariance_su)
      covariance_k[roles$U, roles$U] <- covariance_uu
    }
    if (length(roles$W)) {
      mean_k[roles$W] <- state$gamma
      covariance_k[roles$W, roles$W] <- state$Gamma
    }
    means[[k]] <- mean_k
    covariances[[k]] <- covariance_k
  }
  list(mean = means, covariance = covariances)
}

.joint_observed_state_rowwise <- function(x, mechanism_mask, state, roles) {
  n <- nrow(x)
  d <- ncol(x)
  K <- length(state$pik)
  mask_counts <- rowSums(mechanism_mask)
  mask_dimension <- ncol(mechanism_mask)
  global <- .joint_global_parameters(state, roles, d)
  log_components <- matrix(NA_real_, nrow = n, ncol = K)

  for (i in seq_len(n)) {
    observed <- which(!is.na(x[i, ]))
    for (k in seq_len(K)) {
      gaussian_term <- .joint_log_mvn_observed(
        x[i, observed],
        global$mean[[k]][observed],
        global$covariance[[k]][observed, observed, drop = FALSE],
        paste0("observed covariance for row ", i, ", component ", k)
      )
      missingness_term <- mask_counts[i] * log(state$rho[k]) +
        (mask_dimension - mask_counts[i]) * log1p(-state$rho[k])
      log_components[i, k] <- log(state$pik[k]) + gaussian_term +
        missingness_term
    }
  }

  row_loglik <- apply(log_components, 1L, .joint_log_sum_exp)
  if (any(!is.finite(row_loglik))) {
    .joint_mnar_abort(
      "the joint observed likelihood has a non-finite row contribution",
      "selvarmix_joint_mnar_numerical_error"
    )
  }
  tik <- exp(log_components - row_loglik)
  if (any(!is.finite(tik)) || any(tik < 0) ||
      any(abs(rowSums(tik) - 1) > 1e-8)) {
    .joint_mnar_abort(
      "the joint E-step produced invalid responsibilities",
      "selvarmix_joint_mnar_numerical_error"
    )
  }
  list(
    loglik = sum(row_loglik),
    row_loglik = row_loglik,
    log_components = log_components,
    tik = tik,
    global = global
  )
}

.joint_missing_patterns <- function(x) {
  keys <- apply(is.na(x), 1L, function(mask) {
    paste0(as.integer(mask), collapse = "")
  })
  unique_keys <- unique(keys)
  groups <- lapply(unique_keys, function(key) which(keys == key))
  patterns <- lapply(seq_along(groups), function(pattern_index) {
    rows <- groups[[pattern_index]]
    missing <- which(is.na(x[rows[1L], ]))
    list(
      index = as.integer(pattern_index),
      key = unique_keys[pattern_index],
      rows = as.integer(rows),
      observed = setdiff(seq_len(ncol(x)), missing),
      missing = missing
    )
  })
  list(
    count = as.integer(length(patterns)),
    row_pattern = match(keys, unique_keys),
    patterns = patterns
  )
}

.joint_pattern_component_cache <- function(mean,
                                           sigma,
                                           pattern,
                                           label) {
  observed <- pattern$observed
  missing <- pattern$missing
  factor <- NULL
  log_normalizer <- 0
  gain <- matrix(0, nrow = length(missing), ncol = length(observed))
  conditional_missing_covariance <-
    matrix(0, nrow = length(missing), ncol = length(missing))
  adjustments <- 0L

  if (!length(observed)) {
    conditional_missing_covariance <-
      sigma[missing, missing, drop = FALSE]
  } else {
    sigma_oo <- sigma[observed, observed, drop = FALSE]
    factor <- tryCatch(chol(sigma_oo), error = identity)
    if (inherits(factor, "error")) {
      .joint_mnar_abort(
        paste0(label, " is not positive definite: ", conditionMessage(factor)),
        "selvarmix_joint_mnar_numerical_error"
      )
    }
    log_normalizer <- -0.5 * (
      length(observed) * log(2 * pi) +
        2 * sum(log(base::diag(factor)))
    )

    if (length(missing)) {
      solved_cross <- backsolve(
        factor,
        forwardsolve(
          t(factor),
          sigma[observed, missing, drop = FALSE]
        )
      )
      gain <- t(solved_cross)
      conditional_missing_covariance <-
        sigma[missing, missing, drop = FALSE] -
        sigma[missing, observed, drop = FALSE] %*% solved_cross
      conditional_missing_covariance <-
        (conditional_missing_covariance +
           t(conditional_missing_covariance)) / 2
      decomposition <- eigen(
        conditional_missing_covariance,
        symmetric = TRUE
      )
      numerical_scale <- max(
        .Machine$double.xmin,
        max(abs(conditional_missing_covariance))
      )
      negative_tolerance <- 64 * .Machine$double.eps * numerical_scale
      if (any(decomposition$values < -negative_tolerance)) {
        .joint_mnar_abort(
          paste0(
            label,
            " conditional missing covariance is materially indefinite"
          ),
          "selvarmix_joint_mnar_numerical_error"
        )
      }
      adjusted <- decomposition$values < 0
      conditional_missing_covariance <- decomposition$vectors %*%
        (pmax(decomposition$values, 0) * t(decomposition$vectors))
      conditional_missing_covariance <-
        (conditional_missing_covariance +
           t(conditional_missing_covariance)) / 2
      adjustments <- as.integer(sum(adjusted))
    }
  }

  list(
    factor = factor,
    log_normalizer = log_normalizer,
    gain = gain,
    conditional_missing_covariance = conditional_missing_covariance,
    adjustments = adjustments,
    factorization_count = as.integer(length(observed) > 0L)
  )
}

.joint_pattern_conditional_means <- function(x,
                                             rows,
                                             mean,
                                             pattern,
                                             component_cache) {
  row_count <- length(rows)
  conditional_means <- matrix(
    rep(mean, each = row_count),
    nrow = row_count,
    ncol = length(mean)
  )
  observed <- pattern$observed
  missing <- pattern$missing
  if (length(observed)) {
    observed_values <- x[rows, observed, drop = FALSE]
    conditional_means[, observed] <- observed_values
    if (length(missing)) {
      residual <- sweep(observed_values, 2L, mean[observed], "-")
      predicted_missing <- residual %*% t(component_cache$gain)
      conditional_means[, missing] <- sweep(
        predicted_missing, 2L, mean[missing], "+"
      )
    }
  }
  conditional_means
}

.joint_observed_state <- function(x,
                                  mechanism_mask,
                                  state,
                                  roles,
                                  patterns = NULL) {
  n <- nrow(x)
  d <- ncol(x)
  K <- length(state$pik)
  mask_counts <- rowSums(mechanism_mask)
  mask_dimension <- ncol(mechanism_mask)
  global <- .joint_global_parameters(state, roles, d)
  log_components <- matrix(NA_real_, nrow = n, ncol = K)
  if (is.null(patterns)) patterns <- .joint_missing_patterns(x)
  pattern_cache <- vector("list", patterns$count)
  factorization_count <- 0L

  for (pattern_index in seq_len(patterns$count)) {
    pattern <- patterns$patterns[[pattern_index]]
    rows <- pattern$rows
    observed <- pattern$observed
    components <- vector("list", K)
    for (k in seq_len(K)) {
      component_cache <- .joint_pattern_component_cache(
        global$mean[[k]],
        global$covariance[[k]],
        pattern,
        paste0(
          "observed covariance for pattern ", pattern_index,
          ", component ", k
        )
      )
      factorization_count <- factorization_count +
        component_cache$factorization_count
      components[[k]] <- component_cache

      if (!length(observed)) {
        gaussian_term <- rep(0, length(rows))
      } else {
        residual <- sweep(
          x[rows, observed, drop = FALSE],
          2L,
          global$mean[[k]][observed],
          "-"
        )
        standardized <- forwardsolve(
          t(component_cache$factor),
          t(residual)
        )
        gaussian_term <- component_cache$log_normalizer -
          0.5 * colSums(standardized^2)
      }
      missingness_term <- mask_counts[rows] * log(state$rho[k]) +
        (mask_dimension - mask_counts[rows]) * log1p(-state$rho[k])
      log_components[rows, k] <- log(state$pik[k]) + gaussian_term +
        missingness_term
    }
    pattern_cache[[pattern_index]] <- c(
      pattern,
      list(components = components)
    )
  }

  row_loglik <- apply(log_components, 1L, .joint_log_sum_exp)
  if (any(!is.finite(row_loglik))) {
    .joint_mnar_abort(
      "the joint observed likelihood has a non-finite row contribution",
      "selvarmix_joint_mnar_numerical_error"
    )
  }
  tik <- exp(log_components - row_loglik)
  if (any(!is.finite(tik)) || any(tik < 0) ||
      any(abs(rowSums(tik) - 1) > 1e-8)) {
    .joint_mnar_abort(
      "the joint E-step produced invalid responsibilities",
      "selvarmix_joint_mnar_numerical_error"
    )
  }
  list(
    loglik = sum(row_loglik),
    row_loglik = row_loglik,
    log_components = log_components,
    tik = tik,
    global = global,
    patterns = patterns,
    pattern_cache = pattern_cache,
    factorization_count = as.integer(factorization_count)
  )
}

.joint_conditional_moments <- function(x_row,
                                       mean,
                                       sigma,
                                       covariance_floor,
                                       label) {
  d <- length(x_row)
  observed <- which(!is.na(x_row))
  missing <- which(is.na(x_row))
  conditional_mean <- mean
  conditional_covariance <- matrix(0, nrow = d, ncol = d)
  adjustments <- 0L

  if (!length(observed)) {
    conditional_covariance <- sigma
  } else {
    conditional_mean[observed] <- x_row[observed]
    if (length(missing)) {
      sigma_oo <- sigma[observed, observed, drop = FALSE]
      sigma_mo <- sigma[missing, observed, drop = FALSE]
      residual <- x_row[observed] - mean[observed]
      solved_residual <- .joint_chol_solve(
        sigma_oo, residual, paste0(label, " observed covariance")
      )
      solved_cross <- .joint_chol_solve(
        sigma_oo,
        sigma[observed, missing, drop = FALSE],
        paste0(label, " observed covariance")
      )
      conditional_mean[missing] <- mean[missing] +
        as.numeric(sigma_mo %*% solved_residual)
      conditional_missing_covariance <- sigma[missing, missing, drop = FALSE] -
        sigma_mo %*% solved_cross
      conditional_missing_covariance <-
        (conditional_missing_covariance + t(conditional_missing_covariance)) / 2
      decomposition <- eigen(conditional_missing_covariance, symmetric = TRUE)
      numerical_scale <- max(
        .Machine$double.xmin,
        max(abs(conditional_missing_covariance))
      )
      negative_tolerance <- 64 * .Machine$double.eps * numerical_scale
      if (any(decomposition$values < -negative_tolerance)) {
        .joint_mnar_abort(
          paste0(label, " conditional missing covariance is materially indefinite"),
          "selvarmix_joint_mnar_numerical_error"
        )
      }
      adjusted <- decomposition$values < 0
      conditional_missing_covariance <- decomposition$vectors %*%
        (pmax(decomposition$values, 0) * t(decomposition$vectors))
      conditional_covariance[missing, missing] <-
        (conditional_missing_covariance + t(conditional_missing_covariance)) / 2
      adjustments <- as.integer(sum(adjusted))
    }
  }

  list(
    mean = conditional_mean,
    covariance = conditional_covariance,
    second = conditional_covariance + tcrossprod(conditional_mean),
    adjustments = adjustments
  )
}

.joint_expected_moments_rowwise <- function(x, observed_state, covariance_floor) {
  n <- nrow(x)
  K <- ncol(observed_state$tik)
  values <- vector("list", K)
  adjustments <- 0L
  for (k in seq_len(K)) {
    component <- vector("list", n)
    for (i in seq_len(n)) {
      component[[i]] <- .joint_conditional_moments(
        x[i, ],
        observed_state$global$mean[[k]],
        observed_state$global$covariance[[k]],
        covariance_floor,
        paste0("row ", i, ", component ", k)
      )
      adjustments <- adjustments + component[[i]]$adjustments
    }
    values[[k]] <- component
  }
  list(values = values, adjustments = as.integer(adjustments))
}

# Rowwise reference implementation for checking the streamed sufficient
# statistics used by the fitted model.
.joint_expected_moments <- function(x, observed_state, covariance_floor) {
  .joint_expected_moments_rowwise(x, observed_state, covariance_floor)
}

.joint_weighted_block_moments_rowwise <- function(expected, tik, columns) {
  K <- ncol(tik)
  n <- nrow(tik)
  p <- length(columns)
  first <- numeric(p)
  second <- matrix(0, nrow = p, ncol = p)
  total <- 0
  for (k in seq_len(K)) {
    for (i in seq_len(n)) {
      weight <- tik[i, k]
      first <- first + weight * expected[[k]][[i]]$mean[columns]
      second <- second + weight *
        expected[[k]][[i]]$second[columns, columns, drop = FALSE]
      total <- total + weight
    }
  }
  list(first = first, second = second, total = total)
}

.joint_cached_covariance_block <- function(pattern,
                                           component_cache,
                                           rows,
                                           columns = rows) {
  answer <- matrix(0, nrow = length(rows), ncol = length(columns))
  if (!length(rows) || !length(columns) || !length(pattern$missing)) {
    return(answer)
  }
  row_positions <- match(rows, pattern$missing)
  column_positions <- match(columns, pattern$missing)
  valid_rows <- which(!is.na(row_positions))
  valid_columns <- which(!is.na(column_positions))
  if (length(valid_rows) && length(valid_columns)) {
    answer[valid_rows, valid_columns] <-
      component_cache$conditional_missing_covariance[
        row_positions[valid_rows],
        column_positions[valid_columns],
        drop = FALSE
      ]
  }
  answer
}

.joint_stream_block_moments <- function(conditional_means,
                                        weights,
                                        columns,
                                        pattern,
                                        component_cache) {
  p <- length(columns)
  total <- sum(weights)
  if (!p) {
    return(list(
      first = numeric(),
      second = matrix(numeric(), nrow = 0L, ncol = 0L),
      total = total
    ))
  }
  block <- conditional_means[, columns, drop = FALSE]
  weighted_block <- sweep(block, 1L, weights, "*")
  conditional_covariance <- .joint_cached_covariance_block(
    pattern, component_cache, columns
  )
  list(
    first = colSums(weighted_block),
    second = crossprod(block, weighted_block) +
      total * conditional_covariance,
    total = total
  )
}

.joint_stream_cross_moment <- function(conditional_means,
                                       weights,
                                       left,
                                       right,
                                       pattern,
                                       component_cache) {
  if (!length(left) || !length(right)) {
    return(matrix(0, nrow = length(left), ncol = length(right)))
  }
  left_block <- conditional_means[, left, drop = FALSE]
  weighted_right <- sweep(
    conditional_means[, right, drop = FALSE],
    1L,
    weights,
    "*"
  )
  crossprod(left_block, weighted_right) +
    sum(weights) * .joint_cached_covariance_block(
      pattern, component_cache, left, right
    )
}

.joint_expected_statistics <- function(x, observed_state, roles) {
  K <- ncol(observed_state$tik)
  s <- length(roles$S)
  r <- length(roles$R)
  u <- length(roles$U)
  w <- length(roles$W)
  q <- r + 1L
  component_sizes <- colSums(observed_state$tik)
  s_first <- lapply(seq_len(K), function(k) numeric(s))
  s_second <- lapply(
    seq_len(K),
    function(k) matrix(0, nrow = s, ncol = s)
  )
  a_moment <- matrix(0, nrow = q, ncol = q)
  b_moment <- matrix(0, nrow = q, ncol = u)
  d_moment <- matrix(0, nrow = u, ncol = u)
  w_first <- numeric(w)
  w_second <- matrix(0, nrow = w, ncol = w)
  w_total <- 0
  adjustments <- 0L

  for (pattern_index in seq_len(observed_state$patterns$count)) {
    pattern <- observed_state$pattern_cache[[pattern_index]]
    rows <- pattern$rows
    for (k in seq_len(K)) {
      component_cache <- pattern$components[[k]]
      weights <- observed_state$tik[rows, k]
      conditional_means <- .joint_pattern_conditional_means(
        x,
        rows,
        observed_state$global$mean[[k]],
        pattern,
        component_cache
      )
      adjustments <- adjustments +
        length(rows) * component_cache$adjustments

      s_moments <- .joint_stream_block_moments(
        conditional_means,
        weights,
        roles$S,
        pattern,
        component_cache
      )
      s_first[[k]] <- s_first[[k]] + s_moments$first
      s_second[[k]] <- s_second[[k]] + s_moments$second

      if (u) {
        r_moments <- .joint_stream_block_moments(
          conditional_means,
          weights,
          roles$R,
          pattern,
          component_cache
        )
        u_moments <- .joint_stream_block_moments(
          conditional_means,
          weights,
          roles$U,
          pattern,
          component_cache
        )
        total <- u_moments$total
        a_moment[1L, 1L] <- a_moment[1L, 1L] + total
        b_moment[1L, ] <- b_moment[1L, ] + u_moments$first
        if (r) {
          a_moment[1L, -1L] <- a_moment[1L, -1L] + r_moments$first
          a_moment[-1L, 1L] <- a_moment[-1L, 1L] + r_moments$first
          a_moment[-1L, -1L] <-
            a_moment[-1L, -1L, drop = FALSE] + r_moments$second
          b_moment[-1L, ] <- b_moment[-1L, , drop = FALSE] +
            .joint_stream_cross_moment(
              conditional_means,
              weights,
              roles$R,
              roles$U,
              pattern,
              component_cache
            )
        }
        d_moment <- d_moment + u_moments$second
      }

      if (w) {
        w_moments <- .joint_stream_block_moments(
          conditional_means,
          weights,
          roles$W,
          pattern,
          component_cache
        )
        w_first <- w_first + w_moments$first
        w_second <- w_second + w_moments$second
        w_total <- w_total + w_moments$total
      }
    }
  }

  list(
    component_sizes = component_sizes,
    S = list(first = s_first, second = s_second),
    regression = list(A = a_moment, B = b_moment, D = d_moment),
    W = list(first = w_first, second = w_second, total = w_total),
    adjustments = as.integer(adjustments),
    storage = paste0(
      "streamed by missingness pattern and component; retains only ",
      "role-block sufficient statistics and one pattern-row conditional-mean ",
      "matrix, not K*n full conditional covariance matrices"
    )
  )
}

.joint_m_step_rowwise <- function(x,
                                  mechanism_mask,
                                  tik,
                                  expected,
                                  roles,
                                  diagonal,
                                  component_floor,
                                  covariance_floor) {
  n <- nrow(x)
  K <- ncol(tik)
  s <- length(roles$S)
  r <- length(roles$R)
  u <- length(roles$U)
  w <- length(roles$W)
  component_sizes <- colSums(tik)
  if (any(!is.finite(component_sizes)) || any(component_sizes <= component_floor)) {
    failed <- which(component_sizes <= component_floor | !is.finite(component_sizes))[1L]
    .joint_mnar_abort(
      paste0(
        "effective component size for component ", failed,
        " is not greater than component_floor"
      ),
      "selvarmix_joint_mnar_component_error",
      list(
        component = failed,
        effective_component_size = component_sizes[failed],
        component_floor = component_floor
      )
    )
  }

  covariance_adjustments <- 0L
  mu_S <- vector("list", K)
  sigma_S <- vector("list", K)
  for (k in seq_len(K)) {
    first <- numeric(s)
    second <- matrix(0, nrow = s, ncol = s)
    for (i in seq_len(n)) {
      weight <- tik[i, k]
      first <- first + weight * expected[[k]][[i]]$mean[roles$S]
      second <- second + weight *
        expected[[k]][[i]]$second[roles$S, roles$S, drop = FALSE]
    }
    updated_mean <- first / component_sizes[k]
    updated_covariance <- second / component_sizes[k] -
      tcrossprod(updated_mean)
    stabilized <- .joint_stabilize_covariance(
      updated_covariance,
      diagonal,
      covariance_floor,
      paste0("updated S covariance for component ", k)
    )
    mu_S[[k]] <- updated_mean
    sigma_S[[k]] <- stabilized$value
    covariance_adjustments <- covariance_adjustments + stabilized$adjustments
  }

  a <- numeric(u)
  beta <- matrix(numeric(), nrow = r, ncol = u)
  Omega <- matrix(numeric(), nrow = u, ncol = u)
  if (u) {
    q <- r + 1L
    a_moment <- matrix(0, nrow = q, ncol = q)
    b_moment <- matrix(0, nrow = q, ncol = u)
    d_moment <- matrix(0, nrow = u, ncol = u)
    for (k in seq_len(K)) {
      for (i in seq_len(n)) {
        weight <- tik[i, k]
        mean_i <- expected[[k]][[i]]$mean
        second_i <- expected[[k]][[i]]$second
        predictor_mean <- c(1, mean_i[roles$R])
        predictor_second <- matrix(1, nrow = q, ncol = q)
        if (r) {
          predictor_second[1L, -1L] <- mean_i[roles$R]
          predictor_second[-1L, 1L] <- mean_i[roles$R]
          predictor_second[-1L, -1L] <-
            second_i[roles$R, roles$R, drop = FALSE]
        }
        cross_moment <- matrix(mean_i[roles$U], nrow = q, ncol = u,
                               byrow = TRUE)
        if (r) {
          cross_moment[-1L, ] <-
            second_i[roles$R, roles$U, drop = FALSE]
        }
        a_moment <- a_moment + weight * predictor_second
        b_moment <- b_moment + weight * cross_moment
        d_moment <- d_moment + weight *
          second_i[roles$U, roles$U, drop = FALSE]
      }
    }
    coefficients <- .joint_chol_solve(
      a_moment,
      b_moment,
      "expected regression design cross-product",
      rank_check = TRUE,
      covariance_floor = covariance_floor,
      subclass = "selvarmix_joint_mnar_regression_error"
    )
    a <- as.numeric(coefficients[1L, ])
    if (r) beta <- coefficients[-1L, , drop = FALSE]
    updated_omega <- (d_moment - t(coefficients) %*% b_moment -
                        t(b_moment) %*% coefficients +
                        t(coefficients) %*% a_moment %*% coefficients) / n
    stabilized <- .joint_stabilize_covariance(
      updated_omega,
      diagonal,
      covariance_floor,
      "updated regression covariance"
    )
    Omega <- stabilized$value
    covariance_adjustments <- covariance_adjustments + stabilized$adjustments
  }

  gamma <- numeric(w)
  Gamma <- matrix(numeric(), nrow = w, ncol = w)
  if (w) {
    moments_w <- .joint_weighted_block_moments_rowwise(
      expected, tik, roles$W
    )
    gamma <- moments_w$first / moments_w$total
    updated_gamma <- moments_w$second / moments_w$total - tcrossprod(gamma)
    stabilized <- .joint_stabilize_covariance(
      updated_gamma,
      diagonal,
      covariance_floor,
      "updated W covariance"
    )
    Gamma <- stabilized$value
    covariance_adjustments <- covariance_adjustments + stabilized$adjustments
  }

  mask_counts <- rowSums(mechanism_mask)
  mask_dimension <- ncol(mechanism_mask)
  rho <- as.numeric(crossprod(tik, mask_counts) /
                      (component_sizes * mask_dimension))
  if (any(!is.finite(rho)) || any(rho <= 0) || any(rho >= 1)) {
    failed <- which(!is.finite(rho) | rho <= 0 | rho >= 1)[1L]
    .joint_mnar_abort(
      paste0(
        "component ", failed,
        " has no finite class-only MNARz intercept because its updated ",
        "missingness probability is on the boundary"
      ),
      "selvarmix_joint_mnar_boundary_error",
      list(component = failed, missingness_probability = rho[failed])
    )
  }

  list(
    state = list(
      pik = component_sizes / n,
      mu_S = mu_S,
      sigma_S = sigma_S,
      a = a,
      beta = beta,
      Omega = Omega,
      gamma = gamma,
      Gamma = Gamma,
      rho = rho
    ),
    component_sizes = component_sizes,
    covariance_adjustments = as.integer(covariance_adjustments)
  )
}

.joint_m_step <- function(x,
                          mechanism_mask,
                          tik,
                          expected,
                          roles,
                          diagonal,
                          component_floor,
                          covariance_floor) {
  # Analytic fixtures may supply rowwise moments; fitted models use the
  # equivalent streamed sufficient statistics.
  if (is.null(expected$component_sizes)) {
    return(.joint_m_step_rowwise(
      x,
      mechanism_mask,
      tik,
      expected,
      roles,
      diagonal,
      component_floor,
      covariance_floor
    ))
  }
  n <- nrow(x)
  K <- ncol(tik)
  r <- length(roles$R)
  u <- length(roles$U)
  w <- length(roles$W)
  component_sizes <- as.numeric(expected$component_sizes)
  direct_component_sizes <- colSums(tik)
  size_tolerance <- 64 * .Machine$double.eps *
    max(.Machine$double.xmin, max(abs(direct_component_sizes)))
  if (length(component_sizes) != K || any(!is.finite(component_sizes)) ||
      any(abs(component_sizes - direct_component_sizes) > size_tolerance)) {
    .joint_mnar_abort(
      "streamed component sizes are inconsistent with the responsibilities",
      "selvarmix_joint_mnar_numerical_error"
    )
  }
  if (any(component_sizes <= component_floor)) {
    failed <- which(component_sizes <= component_floor)[1L]
    .joint_mnar_abort(
      paste0(
        "effective component size for component ", failed,
        " is not greater than component_floor"
      ),
      "selvarmix_joint_mnar_component_error",
      list(
        component = failed,
        effective_component_size = component_sizes[failed],
        component_floor = component_floor
      )
    )
  }

  covariance_adjustments <- 0L
  mu_S <- vector("list", K)
  sigma_S <- vector("list", K)
  for (k in seq_len(K)) {
    updated_mean <- expected$S$first[[k]] / component_sizes[k]
    updated_covariance <- expected$S$second[[k]] / component_sizes[k] -
      tcrossprod(updated_mean)
    stabilized <- .joint_stabilize_covariance(
      updated_covariance,
      diagonal,
      covariance_floor,
      paste0("updated S covariance for component ", k)
    )
    mu_S[[k]] <- updated_mean
    sigma_S[[k]] <- stabilized$value
    covariance_adjustments <- covariance_adjustments + stabilized$adjustments
  }

  a <- numeric(u)
  beta <- matrix(numeric(), nrow = r, ncol = u)
  Omega <- matrix(numeric(), nrow = u, ncol = u)
  if (u) {
    a_moment <- expected$regression$A
    b_moment <- expected$regression$B
    d_moment <- expected$regression$D
    coefficients <- .joint_chol_solve(
      a_moment,
      b_moment,
      "expected regression design cross-product",
      rank_check = TRUE,
      covariance_floor = covariance_floor,
      subclass = "selvarmix_joint_mnar_regression_error"
    )
    a <- as.numeric(coefficients[1L, ])
    if (r) beta <- coefficients[-1L, , drop = FALSE]
    updated_omega <- (d_moment - t(coefficients) %*% b_moment -
                        t(b_moment) %*% coefficients +
                        t(coefficients) %*% a_moment %*% coefficients) / n
    stabilized <- .joint_stabilize_covariance(
      updated_omega,
      diagonal,
      covariance_floor,
      "updated regression covariance"
    )
    Omega <- stabilized$value
    covariance_adjustments <- covariance_adjustments + stabilized$adjustments
  }

  gamma <- numeric(w)
  Gamma <- matrix(numeric(), nrow = w, ncol = w)
  if (w) {
    if (!is.finite(expected$W$total) || expected$W$total <= 0) {
      .joint_mnar_abort(
        "streamed W sufficient statistics have non-positive total weight",
        "selvarmix_joint_mnar_numerical_error"
      )
    }
    gamma <- expected$W$first / expected$W$total
    updated_gamma <- expected$W$second / expected$W$total -
      tcrossprod(gamma)
    stabilized <- .joint_stabilize_covariance(
      updated_gamma,
      diagonal,
      covariance_floor,
      "updated W covariance"
    )
    Gamma <- stabilized$value
    covariance_adjustments <- covariance_adjustments + stabilized$adjustments
  }

  mask_counts <- rowSums(mechanism_mask)
  mask_dimension <- ncol(mechanism_mask)
  rho <- as.numeric(crossprod(tik, mask_counts) /
                      (component_sizes * mask_dimension))
  if (any(!is.finite(rho)) || any(rho <= 0) || any(rho >= 1)) {
    failed <- which(!is.finite(rho) | rho <= 0 | rho >= 1)[1L]
    .joint_mnar_abort(
      paste0(
        "component ", failed,
        " has no finite class-only MNARz intercept because its updated ",
        "missingness probability is on the boundary"
      ),
      "selvarmix_joint_mnar_boundary_error",
      list(component = failed, missingness_probability = rho[failed])
    )
  }

  list(
    state = list(
      pik = component_sizes / n,
      mu_S = mu_S,
      sigma_S = sigma_S,
      a = a,
      beta = beta,
      Omega = Omega,
      gamma = gamma,
      Gamma = Gamma,
      rho = rho
    ),
    component_sizes = component_sizes,
    covariance_adjustments = as.integer(covariance_adjustments)
  )
}

.joint_parameter_count <- function(K, roles, diagonal) {
  s <- length(roles$S)
  r <- length(roles$R)
  u <- length(roles$U)
  w <- length(roles$W)
  sigma_count <- if (diagonal) K * s else K * s * (s + 1) / 2
  omega_count <- if (diagonal) u else u * (u + 1) / 2
  gamma_count <- if (diagonal) w else w * (w + 1) / 2
  as.numeric(
    (K - 1) + K * s + sigma_count +
      u * (r + 1) + omega_count +
      w + gamma_count + K
  )
}

.joint_mean_impute <- function(x) {
  completed <- x
  for (j in seq_len(ncol(x))) {
    observed <- x[, j][!is.na(x[, j])]
    if (!length(observed)) {
      .joint_mnar_abort(
        paste0("column ", j, " is entirely missing and cannot be initialized"),
        "selvarmix_joint_mnar_input_error"
      )
    }
    completed[is.na(completed[, j]), j] <- mean(observed)
  }
  completed
}

.joint_validate_initial_labels <- function(labels, n, K, component_floor) {
  if (length(labels) != n || !is.numeric(labels) || anyNA(labels) ||
      any(!is.finite(labels)) || any(labels != as.integer(labels)) ||
      any(labels < 1L) || any(labels > K)) {
    .joint_mnar_abort(
      "initial_labels must contain one integer in 1:K for every row",
      "selvarmix_joint_mnar_input_error"
    )
  }
  labels <- as.integer(labels)
  sizes <- tabulate(labels, nbins = K)
  if (any(sizes <= component_floor)) {
    failed <- which(sizes <= component_floor)[1L]
    .joint_mnar_abort(
      paste0("initial component ", failed, " is not larger than component_floor"),
      "selvarmix_joint_mnar_component_error",
      list(component = failed, effective_component_size = sizes[failed])
    )
  }
  labels
}

.joint_initial_state_from_labels <- function(x,
                                             mechanism_mask,
                                             K,
                                             roles,
                                             labels,
                                             diagonal,
                                             component_floor,
                                             covariance_floor) {
  completed <- .joint_mean_impute(x)
  labels <- .joint_validate_initial_labels(labels, nrow(x), K, component_floor)
  n <- nrow(x)
  s <- length(roles$S)
  r <- length(roles$R)
  u <- length(roles$U)
  w <- length(roles$W)
  sizes <- tabulate(labels, nbins = K)
  adjustments <- 0L
  mu_S <- vector("list", K)
  sigma_S <- vector("list", K)
  for (k in seq_len(K)) {
    block <- completed[labels == k, roles$S, drop = FALSE]
    mu_S[[k]] <- colMeans(block)
    centered <- sweep(block, 2L, mu_S[[k]], "-")
    covariance <- crossprod(centered) / nrow(block)
    stabilized <- .joint_stabilize_covariance(
      covariance,
      diagonal,
      covariance_floor,
      paste0("initial S covariance for component ", k)
    )
    sigma_S[[k]] <- stabilized$value
    adjustments <- adjustments + stabilized$adjustments
  }

  a <- numeric(u)
  beta <- matrix(numeric(), nrow = r, ncol = u)
  Omega <- matrix(numeric(), nrow = u, ncol = u)
  if (u) {
    design <- cbind(1, completed[, roles$R, drop = FALSE])
    response <- completed[, roles$U, drop = FALSE]
    coefficients <- .joint_chol_solve(
      crossprod(design),
      crossprod(design, response),
      "initial regression design cross-product",
      rank_check = TRUE,
      covariance_floor = covariance_floor,
      subclass = "selvarmix_joint_mnar_regression_error"
    )
    a <- as.numeric(coefficients[1L, ])
    if (r) beta <- coefficients[-1L, , drop = FALSE]
    residual <- response - design %*% coefficients
    covariance <- crossprod(residual) / n
    stabilized <- .joint_stabilize_covariance(
      covariance,
      diagonal,
      covariance_floor,
      "initial regression covariance"
    )
    Omega <- stabilized$value
    adjustments <- adjustments + stabilized$adjustments
  }

  gamma <- numeric(w)
  Gamma <- matrix(numeric(), nrow = w, ncol = w)
  if (w) {
    block <- completed[, roles$W, drop = FALSE]
    gamma <- colMeans(block)
    centered <- sweep(block, 2L, gamma, "-")
    covariance <- crossprod(centered) / n
    stabilized <- .joint_stabilize_covariance(
      covariance,
      diagonal,
      covariance_floor,
      "initial W covariance"
    )
    Gamma <- stabilized$value
    adjustments <- adjustments + stabilized$adjustments
  }

  overall_rho <- mean(mechanism_mask)
  if (!is.finite(overall_rho) || overall_rho <= 0 || overall_rho >= 1) {
    .joint_mnar_abort(
      paste0(
        "the mechanism mask is on the all-observed or all-missing boundary; ",
        "no finite class-only MNARz intercept exists"
      ),
      "selvarmix_joint_mnar_boundary_error"
    )
  }
  missing_per_row <- rowSums(mechanism_mask)
  class_missing <- vapply(seq_len(K), function(k) {
    sum(missing_per_row[labels == k])
  }, numeric(1))
  class_trials <- sizes * ncol(mechanism_mask)
  # Jeffreys smoothing applies only to the initial state, avoiding infinite
  # intercepts for an initially pure class without clipping any later M-step.
  initial_rho <- (class_missing + 0.5) / (class_trials + 1)

  list(
    state = list(
      pik = sizes / n,
      mu_S = mu_S,
      sigma_S = sigma_S,
      a = a,
      beta = beta,
      Omega = Omega,
      gamma = gamma,
      Gamma = Gamma,
      rho = as.numeric(initial_rho)
    ),
    covariance_adjustments = as.integer(adjustments)
  )
}

.joint_extract_state <- function(init) {
  candidate <- if (!is.null(init$parameters)) init$parameters else init
  if (is.null(candidate$pik) && !is.null(candidate$pi)) candidate$pik <- candidate$pi
  candidate
}

.joint_validate_supplied_state <- function(init,
                                           K,
                                           roles,
                                           diagonal,
                                           covariance_floor) {
  state <- .joint_deep_clone(.joint_extract_state(init))
  required <- c(
    "pik", "mu_S", "sigma_S", "a", "beta", "Omega",
    "gamma", "Gamma"
  )
  if (!all(required %in% names(state)) ||
      (!"rho" %in% names(state) && !"alpha" %in% names(state))) {
    .joint_mnar_abort(
      paste0("init does not contain the complete joint SRUW--MNARz state: ",
             paste(c(required, "rho or alpha"), collapse = ", ")),
      "selvarmix_joint_mnar_input_error"
    )
  }
  s <- length(roles$S)
  r <- length(roles$R)
  u <- length(roles$U)
  w <- length(roles$W)
  state$pik <- as.numeric(state$pik)
  if (length(state$pik) != K || any(!is.finite(state$pik)) ||
      any(state$pik <= 0) || abs(sum(state$pik) - 1) > 1e-8) {
    .joint_mnar_abort("init$pik must be positive and sum to one",
                      "selvarmix_joint_mnar_input_error")
  }
  if (!is.list(state$mu_S) || length(state$mu_S) != K ||
      !is.list(state$sigma_S) || length(state$sigma_S) != K) {
    .joint_mnar_abort("init must provide K S means and covariances",
                      "selvarmix_joint_mnar_input_error")
  }
  adjustments <- 0L
  for (k in seq_len(K)) {
    state$mu_S[[k]] <- as.numeric(state$mu_S[[k]])
    if (length(state$mu_S[[k]]) != s || any(!is.finite(state$mu_S[[k]]))) {
      .joint_mnar_abort("each init$mu_S entry must be finite with length |S|",
                        "selvarmix_joint_mnar_input_error")
    }
    if (!identical(dim(as.matrix(state$sigma_S[[k]])), c(s, s))) {
      .joint_mnar_abort("each init$sigma_S entry must be an |S| square matrix",
                        "selvarmix_joint_mnar_input_error")
    }
    stabilized <- .joint_stabilize_covariance(
      state$sigma_S[[k]], diagonal, covariance_floor,
      paste0("supplied S covariance for component ", k),
      "selvarmix_joint_mnar_input_error"
    )
    state$sigma_S[[k]] <- stabilized$value
    adjustments <- adjustments + stabilized$adjustments
  }
  state$a <- as.numeric(state$a)
  state$beta <- as.matrix(state$beta)
  state$Omega <- as.matrix(state$Omega)
  state$gamma <- as.numeric(state$gamma)
  state$Gamma <- as.matrix(state$Gamma)
  if (length(state$a) != u || any(!is.finite(state$a)) ||
      !identical(dim(state$beta), c(r, u)) || any(!is.finite(state$beta)) ||
      !identical(dim(state$Omega), c(u, u)) || any(!is.finite(state$Omega)) ||
      length(state$gamma) != w || any(!is.finite(state$gamma)) ||
      !identical(dim(state$Gamma), c(w, w)) || any(!is.finite(state$Gamma))) {
    .joint_mnar_abort("the supplied SRUW block dimensions do not match roles",
                      "selvarmix_joint_mnar_input_error")
  }
  if (u) {
    stabilized <- .joint_stabilize_covariance(
      state$Omega, diagonal, covariance_floor, "supplied regression covariance",
      "selvarmix_joint_mnar_input_error"
    )
    state$Omega <- stabilized$value
    adjustments <- adjustments + stabilized$adjustments
  }
  if (w) {
    stabilized <- .joint_stabilize_covariance(
      state$Gamma, diagonal, covariance_floor, "supplied W covariance",
      "selvarmix_joint_mnar_input_error"
    )
    state$Gamma <- stabilized$value
    adjustments <- adjustments + stabilized$adjustments
  }
  if (is.null(state$rho)) {
    alpha <- as.matrix(state$alpha)
    if (nrow(alpha) != K) {
      .joint_mnar_abort("init$alpha must have K rows",
                        "selvarmix_joint_mnar_input_error")
    }
    if (ncol(alpha) > 1L &&
        any(abs(alpha - alpha[, 1L, drop = FALSE]) > 1e-12)) {
      .joint_mnar_abort("each init$alpha row must contain one shared MNARz intercept",
                        "selvarmix_joint_mnar_input_error")
    }
    state$rho <- stats::pnorm(alpha[, 1L])
  }
  state$rho <- as.numeric(state$rho)
  if (length(state$rho) != K || any(!is.finite(state$rho)) ||
      any(state$rho <= 0) || any(state$rho >= 1)) {
    .joint_mnar_abort("init$rho must contain K probabilities strictly between zero and one",
                      "selvarmix_joint_mnar_boundary_error")
  }
  state$alpha <- NULL
  list(state = state, covariance_adjustments = as.integer(adjustments))
}

.joint_impute_rowwise <- function(x, observed_state, covariance_floor) {
  imputed <- x
  n <- nrow(x)
  K <- ncol(observed_state$tik)
  adjustments <- 0L
  for (i in seq_len(n)) {
    missing <- which(is.na(x[i, ]))
    if (!length(missing)) next
    conditional <- vector("list", K)
    for (k in seq_len(K)) {
      conditional[[k]] <- .joint_conditional_moments(
        x[i, ],
        observed_state$global$mean[[k]],
        observed_state$global$covariance[[k]],
        covariance_floor,
        paste0("final imputation row ", i, ", component ", k)
      )
      adjustments <- adjustments + conditional[[k]]$adjustments
    }
    for (j in missing) {
      imputed[i, j] <- sum(vapply(
        seq_len(K),
        function(k) observed_state$tik[i, k] * conditional[[k]]$mean[j],
        numeric(1)
      ))
    }
  }
  list(value = imputed, adjustments = as.integer(adjustments))
}

.joint_impute <- function(x, observed_state, covariance_floor) {
  imputed <- x
  K <- ncol(observed_state$tik)
  adjustments <- 0L
  for (pattern_index in seq_len(observed_state$patterns$count)) {
    pattern <- observed_state$pattern_cache[[pattern_index]]
    rows <- pattern$rows
    missing <- pattern$missing
    if (!length(missing)) next
    imputed_values <- matrix(
      0,
      nrow = length(rows),
      ncol = length(missing)
    )
    for (k in seq_len(K)) {
      component_cache <- pattern$components[[k]]
      conditional_means <- .joint_pattern_conditional_means(
        x,
        rows,
        observed_state$global$mean[[k]],
        pattern,
        component_cache
      )
      imputed_values <- imputed_values + sweep(
        conditional_means[, missing, drop = FALSE],
        1L,
        observed_state$tik[rows, k],
        "*"
      )
      adjustments <- adjustments +
        length(rows) * component_cache$adjustments
    }
    imputed[rows, missing] <- imputed_values
  }
  list(value = imputed, adjustments = as.integer(adjustments))
}

#' Joint fixed-role SRUW Gaussian mixture with class-only MNARz missingness
#'
#' This internal estimator fits a prescribed role partition. `roles$S`,
#' `roles$U`, and `roles$W` partition the columns of `x`; `roles$R` is a common
#' subset of `roles$S` used to regress U.
#' The missingness mask must equal `is.na(x)`;
#' pure joint MNARz models the mask of every variable in the fitted role model.
#'
#' @keywords internal
EMJointSRUWMNARz <- function(
    x,
    K,
    roles = list(
      S = seq_len(ncol(as.matrix(x))),
      R = integer(),
      U = integer(),
      W = integer()
    ),
    diag = TRUE,
    rmax = 100L,
    tol = 1e-4,
    init = NULL,
    initial_labels = NULL,
    mechanism_mask = NULL,
    component_floor = 1,
    covariance_floor = sqrt(.Machine$double.eps),
    absolute_tolerance = sqrt(.Machine$double.eps),
    monotonicity_tolerance = sqrt(.Machine$double.eps)) {
  x <- as.matrix(x)
  if (!is.numeric(x) || !length(x) || nrow(x) < 1L || ncol(x) < 1L) {
    .joint_mnar_abort("x must be a non-empty numeric matrix",
                      "selvarmix_joint_mnar_input_error")
  }
  if (any(is.nan(x)) || any(is.infinite(x))) {
    .joint_mnar_abort("x may contain NA values but not NaN or infinite values",
                      "selvarmix_joint_mnar_input_error")
  }
  K <- .joint_check_scalar_integer(K, "K")
  rmax <- .joint_check_scalar_integer(rmax, "rmax")
  if (K > nrow(x)) {
    .joint_mnar_abort("K cannot exceed the number of rows",
                      "selvarmix_joint_mnar_input_error")
  }
  if (length(diag) != 1L || !is.logical(diag) || is.na(diag)) {
    .joint_mnar_abort("diag must be TRUE or FALSE",
                      "selvarmix_joint_mnar_input_error")
  }
  if (length(tol) != 1L || !is.numeric(tol) || !is.finite(tol) || tol < 0) {
    .joint_mnar_abort("tol must be finite and non-negative",
                      "selvarmix_joint_mnar_input_error")
  }
  if (length(component_floor) != 1L || !is.numeric(component_floor) ||
      !is.finite(component_floor) || component_floor < 0) {
    .joint_mnar_abort("component_floor must be finite and non-negative",
                      "selvarmix_joint_mnar_input_error")
  }
  requested_component_floor <- as.numeric(component_floor)
  if (length(covariance_floor) != 1L || !is.numeric(covariance_floor) ||
      !is.finite(covariance_floor) || covariance_floor <= 0) {
    .joint_mnar_abort("covariance_floor must be finite and positive",
                      "selvarmix_joint_mnar_input_error")
  }
  if (length(absolute_tolerance) != 1L || !is.numeric(absolute_tolerance) ||
      !is.finite(absolute_tolerance) || absolute_tolerance < 0) {
    .joint_mnar_abort("absolute_tolerance must be finite and non-negative",
                      "selvarmix_joint_mnar_input_error")
  }
  if (length(monotonicity_tolerance) != 1L ||
      !is.numeric(monotonicity_tolerance) ||
      !is.finite(monotonicity_tolerance) || monotonicity_tolerance < 0) {
    .joint_mnar_abort("monotonicity_tolerance must be finite and non-negative",
                      "selvarmix_joint_mnar_input_error")
  }
  if (!is.null(init) && !is.list(init)) {
    .joint_mnar_abort("init must be NULL or a complete fitted-state list",
                      "selvarmix_joint_mnar_input_error")
  }
  if (!is.null(init) && !is.null(initial_labels)) {
    .joint_mnar_abort("supply at most one of init and initial_labels",
                      "selvarmix_joint_mnar_input_error")
  }

  roles <- .joint_validate_roles(roles, ncol(x))
  .joint_validate_observation_design(x, roles, diag)
  component_floor <- max(
    requested_component_floor,
    if (isTRUE(diag)) 1 else max(1L, length(roles$S))
  )
  if (is.null(mechanism_mask)) mechanism_mask <- is.na(x)
  mechanism_mask <- as.matrix(mechanism_mask)
  if (!identical(dim(mechanism_mask), dim(x)) ||
      anyNA(mechanism_mask) ||
      any(!(mechanism_mask %in% c(FALSE, TRUE, 0, 1)))) {
    .joint_mnar_abort(
      "mechanism_mask must be a complete binary matrix with dim(x)",
      "selvarmix_joint_mnar_input_error"
    )
  }
  mechanism_mask <- matrix(
    as.logical(mechanism_mask),
    nrow = nrow(mechanism_mask),
    ncol = ncol(mechanism_mask)
  )
  if (!all(mechanism_mask == is.na(x))) {
    .joint_mnar_abort(
      "mechanism_mask must equal is.na(x) for the pure joint MNARz likelihood",
      "selvarmix_joint_mnar_input_error"
    )
  }
  missing_patterns <- .joint_missing_patterns(x)

  initialization <- NULL
  if (!is.null(init)) {
    initialized <- .joint_validate_supplied_state(
      init, K, roles, diag, covariance_floor
    )
    initialization <- "supplied_complete_state"
  } else {
    completed <- .joint_mean_impute(x)
    if (is.null(initial_labels)) {
      if (K == 1L) {
        initial_labels <- rep(1L, nrow(x))
      } else {
        clustering <- stats::hclust(
          stats::dist(completed[, roles$S, drop = FALSE]),
          method = "ward.D2"
        )
        initial_labels <- stats::cutree(clustering, k = K)
      }
      initialization <- "deterministic_hierarchical"
    } else {
      initialization <- "pinned_labels"
    }
    initialized <- .joint_initial_state_from_labels(
      x,
      mechanism_mask,
      K,
      roles,
      initial_labels,
      diag,
      component_floor,
      covariance_floor
    )
  }

  state <- initialized$state
  covariance_adjustments <- initialized$covariance_adjustments
  conditional_adjustments <- 0L
  observed <- .joint_observed_state(
    x, mechanism_mask, state, roles, patterns = missing_patterns
  )
  observed_block_factorizations <- observed$factorization_count
  loglik_trace <- observed$loglik
  min_effective_component_size <- min(colSums(observed$tik))
  converged <- FALSE
  final_improvement <- NA_real_

  for (iteration in seq_len(rmax)) {
    expected <- .joint_expected_statistics(x, observed, roles)
    conditional_adjustments <- conditional_adjustments + expected$adjustments
    updated <- .joint_m_step(
      x,
      mechanism_mask,
      observed$tik,
      expected,
      roles,
      diag,
      component_floor,
      covariance_floor
    )
    candidate <- .joint_observed_state(
      x,
      mechanism_mask,
      updated$state,
      roles,
      patterns = missing_patterns
    )
    observed_block_factorizations <- observed_block_factorizations +
      candidate$factorization_count
    candidate_component_sizes <- colSums(candidate$tik)
    if (any(!is.finite(candidate_component_sizes)) ||
        any(candidate_component_sizes <= component_floor)) {
      failed <- which(
        !is.finite(candidate_component_sizes) |
          candidate_component_sizes <= component_floor
      )[1L]
      .joint_mnar_abort(
        paste0(
          "post-update effective component size for component ", failed,
          " is not greater than component_floor; the candidate state is unscorable"
        ),
        "selvarmix_joint_mnar_component_error",
        list(
          iteration = iteration,
          component = failed,
          effective_component_size = candidate_component_sizes[failed],
          component_floor = component_floor,
          phase = "candidate_e_step"
        )
      )
    }
    improvement <- candidate$loglik - observed$loglik
    decrease_tolerance <- monotonicity_tolerance *
      (1 + abs(observed$loglik))
    if (!is.finite(improvement) || improvement < -decrease_tolerance) {
      .joint_mnar_abort(
        paste0(
          "joint observed log-likelihood decreased materially at iteration ",
          iteration
        ),
        "selvarmix_joint_mnar_monotonicity_error",
        list(
          iteration = iteration,
          previous_loglik = observed$loglik,
          candidate_loglik = candidate$loglik,
          improvement = improvement,
          tolerance = decrease_tolerance
        )
      )
    }
    state <- updated$state
    observed <- candidate
    covariance_adjustments <- covariance_adjustments +
      updated$covariance_adjustments
    min_effective_component_size <- min(
      min_effective_component_size,
      updated$component_sizes,
      candidate_component_sizes
    )
    loglik_trace <- c(loglik_trace, observed$loglik)
    final_improvement <- improvement
    convergence_threshold <- max(
      absolute_tolerance,
      tol * (1 + abs(loglik_trace[length(loglik_trace) - 1L]))
    )
    if (abs(improvement) <= convergence_threshold) {
      converged <- TRUE
      break
    }
  }

  iterations <- length(loglik_trace) - 1L
  if (!converged) {
    .joint_mnar_abort(
      paste0(
        "joint SRUW--MNARz EM reached rmax = ", rmax,
        " without satisfying the convergence tolerance; no BIC is available"
      ),
      "selvarmix_joint_mnar_nonconvergence_error",
      list(
        iterations = iterations,
        termination_reason = "max_iterations",
        loglik_trace = loglik_trace,
        final_loglik_improvement = final_improvement
      )
    )
  }

  imputed <- .joint_impute(x, observed, covariance_floor)
  conditional_adjustments <- conditional_adjustments + imputed$adjustments
  parameter_count <- .joint_parameter_count(K, roles, diag)
  bic_min <- -2 * observed$loglik + parameter_count * log(nrow(x))
  score <- -bic_min
  alpha <- matrix(
    stats::qnorm(state$rho),
    nrow = K,
    ncol = ncol(mechanism_mask)
  )
  global <- .joint_global_parameters(state, roles, ncol(x))

  list(
    estimator_mode = "joint_mnarz",
    roles = roles,
    K = K,
    covariance = if (diag) "diagonal" else "full",
    loglik_obs = observed$loglik,
    bic_min = bic_min,
    score = score,
    criterionValue = list(BIC = bic_min, BIC_score = score),
    criterion_convention = list(
      BIC = "minimize -2 * logLik + q * log(n)",
      score = "maximize 2 * logLik - q * log(n)"
    ),
    parameter_count = parameter_count,
    parameters = list(
      pik = state$pik,
      mu_S = state$mu_S,
      sigma_S = state$sigma_S,
      a = state$a,
      beta = state$beta,
      Omega = state$Omega,
      gamma = state$gamma,
      Gamma = state$Gamma,
      rho = state$rho,
      alpha = alpha,
      global_mu = global$mean,
      global_sigma = global$covariance
    ),
    proba = observed$tik,
    partition = max.col(observed$tik, ties.method = "first"),
    imputedData = imputed$value,
    loglik_trace = loglik_trace,
    diagnostics = list(
      iterations = as.integer(iterations),
      converged = TRUE,
      termination_reason = "tolerance",
      final_loglik_improvement = final_improvement,
      loglik_monotone = all(diff(loglik_trace) >= 0),
      no_material_loglik_decrease = TRUE,
      min_effective_component_size = min_effective_component_size,
      component_floor = component_floor,
      requested_component_floor = requested_component_floor,
      component_mass_rule = if (isTRUE(diag)) {
        "max(user floor, 1)"
      } else {
        "max(user floor, dimension of component-specific S covariance)"
      },
      covariance_floor = covariance_floor,
      covariance_adjustments = as.integer(covariance_adjustments),
      conditional_covariance_adjustments = as.integer(conditional_adjustments),
      missingness_parameter_count = as.numeric(K),
      convergence_absolute_tolerance = absolute_tolerance,
      convergence_relative_tolerance = tol,
      monotonicity_relative_tolerance = monotonicity_tolerance,
      final_convergence_threshold = convergence_threshold,
      mechanism_dimension = ncol(mechanism_mask),
      unique_missingness_pattern_count = missing_patterns$count,
      observed_block_factorizations_per_e_step =
        as.integer(observed$factorization_count),
      total_observed_block_factorizations =
        as.integer(observed_block_factorizations),
      observed_factorization_bound_per_e_step =
        as.integer(K * missing_patterns$count),
      factorization_complexity = paste0(
        "at most one observed-block Cholesky factorization per component and ",
        "unique missingness pattern in each E-step; the same factor serves ",
        "the likelihood, conditional means, and conditional covariance"
      ),
      conditional_moment_storage = expected$storage,
      retained_full_conditional_moment_matrices = 0L,
      parameter_count = parameter_count,
      initialization = initialization,
      criterion_available = TRUE,
      criterion_scope = paste0(
        "Schwarz-type criterion for one joint observed likelihood of ",
        "(Y_obs, C) under fixed SRUW roles and class-only MNARz; parameter ",
        "dimension is the nominal dimension of the returned local fit."
      )
    ),
    paper_alignment = list(
      model = "fixed-role joint SRUW Gaussian factorization with pure class-only MNARz",
      missingness = paste0(
        "one rho_k per component, shared across all ",
        ncol(mechanism_mask), " mask coordinates"
      ),
      link = "rho is likelihood-primary; alpha = qnorm(rho) is stored for probit compatibility",
      criterion = "one total observed-data BIC; no additive complete-case SRUW decomposition",
      status = "fixed_role_joint_sruw_mnarz"
    )
  )
}

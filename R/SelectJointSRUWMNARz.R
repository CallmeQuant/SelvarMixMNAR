# Joint SRUW--MNARz role selector. Dependent candidates are evaluated
# sequentially; parallel execution is valid only across independent values of
# K or independent data sets.

.select_joint_initial_labels <- function(x, K, seed_variable) {
  seed <- .joint_mean_impute(x[, seed_variable, drop = FALSE])
  if (K == 1L) return(rep(1L, nrow(x)))
  labels <- stats::cutree(
    stats::hclust(stats::dist(seed), method = "ward.D2"),
    k = K
  )
  .joint_validate_initial_labels(labels, nrow(x), K, component_floor = 1)
}

.select_joint_role_key <- function(S, R, U, W) {
  encode <- function(x) paste(sort(as.integer(x)), collapse = ",")
  paste0(
    "S=", encode(S), "|R=", encode(R),
    "|U=", encode(U), "|W=", encode(W)
  )
}

.select_joint_sruw_mnarz_impl <- function(
    x,
    K,
    rank,
    hsize,
    diag,
    rmax,
    tol,
    r_search_cap,
    initial_labels = NULL,
    .fit_function = EMJointSRUWMNARz) {
  x <- as.matrix(x)
  if (!is.numeric(x) || !length(x) || nrow(x) < 1L || ncol(x) < 1L ||
      any(is.nan(x)) || any(is.infinite(x))) {
    .joint_mnar_abort(
      "x must be a non-empty numeric matrix containing only finite values or NA",
      "selvarmix_joint_mnar_selection_input_error"
    )
  }
  K <- .joint_check_scalar_integer(K, "K")
  hsize <- .joint_check_scalar_integer(hsize, "hsize")
  rmax <- .joint_check_scalar_integer(rmax, "rmax")
  r_search_cap <- .joint_check_scalar_integer(r_search_cap, "r_search_cap")
  if (K > nrow(x)) {
    .joint_mnar_abort(
      "K cannot exceed the number of rows",
      "selvarmix_joint_mnar_selection_input_error"
    )
  }
  if (length(diag) != 1L || !is.logical(diag) || is.na(diag)) {
    .joint_mnar_abort("diag must be TRUE or FALSE",
                      "selvarmix_joint_mnar_selection_input_error")
  }
  if (length(tol) != 1L || !is.numeric(tol) || !is.finite(tol) || tol < 0) {
    .joint_mnar_abort("tol must be finite and non-negative",
                      "selvarmix_joint_mnar_selection_input_error")
  }
  if (length(rank) != ncol(x) || !is.numeric(rank) || anyNA(rank) ||
      any(!is.finite(rank)) || any(rank != as.integer(rank)) ||
      !identical(sort(as.integer(rank)), seq_len(ncol(x)))) {
    .joint_mnar_abort(
      "rank must be an exact permutation of seq_len(ncol(x))",
      "selvarmix_joint_mnar_selection_input_error"
    )
  }
  if (!is.function(.fit_function)) {
    .joint_mnar_abort(".fit_function must be a function",
                      "selvarmix_joint_mnar_selection_input_error")
  }
  rank <- as.integer(rank)
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("V", seq_len(ncol(x)))
  }

  # Use one validated initial partition for every S-versus-U and U-versus-W
  # candidate. The default partition uses the first ranked variable.
  initialization <- if (is.null(initial_labels)) {
    initial_labels <- .select_joint_initial_labels(x, K, rank[1L])
    "deterministic seed-S hierarchical partition"
  } else {
    initial_labels <- .joint_validate_initial_labels(
      initial_labels, nrow(x), K, component_floor = 1
    )
    "validated adapter-supplied partition"
  }
  cache <- new.env(parent = emptyenv(), hash = TRUE)
  cache_hits <- 0L
  cache_misses <- 0L
  r_search_records <- list()

  fit_roles <- function(S, R = integer(), U = integer(), W = integer(),
                        context = "candidate") {
    S <- sort(as.integer(S))
    R <- sort(as.integer(R))
    U <- sort(as.integer(U))
    W <- sort(as.integer(W))
    if (!length(U)) R <- integer()
    key <- .select_joint_role_key(S, R, U, W)
    if (exists(key, envir = cache, inherits = FALSE)) {
      cache_hits <<- cache_hits + 1L
      return(get(key, envir = cache, inherits = FALSE))
    }
    cache_misses <<- cache_misses + 1L
    active <- sort(c(S, U, W))
    local_roles <- list(
      S = match(S, active),
      R = match(R, active),
      U = match(U, active),
      W = match(W, active)
    )
    fit <- tryCatch(
      .fit_function(
        x = x[, active, drop = FALSE],
        K = K,
        roles = local_roles,
        diag = diag,
        rmax = rmax,
        tol = tol,
        initial_labels = initial_labels
      ),
      error = identity
    )
    if (inherits(fit, "error")) {
      .joint_mnar_abort(
        paste0(
          "joint role selection failed in ", context, " for ", key,
          ": ", conditionMessage(fit)
        ),
        "selvarmix_joint_mnar_selection_error",
        list(parent = fit, context = context,
             role_specification = list(S = S, R = R, U = U, W = W))
      )
    }
    if (!is.list(fit) || length(fit$score) != 1L ||
        !is.finite(fit$score) || length(fit$bic_min) != 1L ||
        !is.finite(fit$bic_min) || !is.list(fit$diagnostics) ||
        !isTRUE(fit$diagnostics$converged) ||
        !isTRUE(fit$diagnostics$criterion_available)) {
      .joint_mnar_abort(
        paste0(
          "joint role selection received an unscorable fit in ", context,
          " for ", key
        ),
        "selvarmix_joint_mnar_selection_error",
        list(context = context,
             role_specification = list(S = S, R = R, U = U, W = W))
      )
    }
    sign_tolerance <- 64 * .Machine$double.eps *
      (1 + abs(fit$score) + abs(fit$bic_min))
    if (abs(fit$score + fit$bic_min) > sign_tolerance) {
      .joint_mnar_abort(
        paste0(
          "joint role selection received inconsistent score/BIC signs in ",
          context, " for ", key
        ),
        "selvarmix_joint_mnar_selection_error",
        list(
          context = context,
          role_specification = list(S = S, R = R, U = U, W = W),
          score = as.numeric(fit$score),
          bic_min = as.numeric(fit$bic_min),
          sign_tolerance = sign_tolerance
        )
      )
    }
    value <- list(
      fit = fit,
      score = as.numeric(fit$score),
      bic_min = as.numeric(fit$bic_min),
      roles = list(S = S, R = R, U = U, W = W),
      key = key
    )
    assign(key, value, envir = cache)
    value
  }

  best_R <- function(S, U, W, context) {
    S_ordered <- rank[rank %in% S]
    subset_count <- if (length(S_ordered) < 1024L) {
      2^length(S_ordered)
    } else {
      Inf
    }
    if (!length(U)) {
      selected <- fit_roles(S, integer(), U, W, context)
      record <- list(
        context = context,
        method = "not_applicable",
        evaluated = 1L,
        subset_count = 1,
        r_search_cap = r_search_cap,
        selected_R = integer(),
        selected_score = selected$score,
        candidates = list(list(R = integer(), score = selected$score))
      )
      r_search_records[[length(r_search_records) + 1L]] <<- record
      return(c(selected, list(
        R = integer(), method = "not_applicable", search_record = record
      )))
    }

    candidates <- list()
    evaluate_R <- function(R) {
      evaluated <- fit_roles(S, R, U, W, context)
      candidates[[length(candidates) + 1L]] <<- list(
        R = as.integer(R), score = evaluated$score, bic_min = evaluated$bic_min
      )
      evaluated
    }

    maximum_exhaustive_dimension <- floor(log(r_search_cap, base = 2))
    if (length(S_ordered) <= maximum_exhaustive_dimension) {
      subsets <- list(integer())
      for (size in seq_along(S_ordered)) {
        index_subsets <- utils::combn(
          length(S_ordered), size, simplify = FALSE
        )
        subsets <- c(
          subsets,
          lapply(index_subsets, function(index) S_ordered[index])
        )
      }
      selected <- NULL
      selected_R <- integer()
      for (R in subsets) {
        candidate <- evaluate_R(R)
        if (is.null(selected) || candidate$score > selected$score) {
          selected <- candidate
          selected_R <- as.integer(R)
        }
      }
      method <- "exhaustive"
    } else {
      selected_R <- integer()
      selected <- evaluate_R(selected_R)
      remaining <- S_ordered
      repeat {
        if (!length(remaining)) break
        step_best <- NULL
        step_variable <- NA_integer_
        for (variable in remaining) {
          candidate <- evaluate_R(c(selected_R, variable))
          if (is.null(step_best) || candidate$score > step_best$score) {
            step_best <- candidate
            step_variable <- variable
          }
        }
        if (is.null(step_best) || step_best$score <= selected$score) break
        selected <- step_best
        selected_R <- c(selected_R, step_variable)
        remaining <- remaining[remaining != step_variable]
      }
      method <- "deterministic_forward"
    }
    record <- list(
      context = context,
      method = method,
      evaluated = length(candidates),
      subset_count = subset_count,
      r_search_cap = r_search_cap,
      selected_R = selected_R,
      selected_score = selected$score,
      candidates = candidates,
      tie_rule = paste0(
        "strict score improvement only; ties retain the first subset in ",
        "size-then-rank order"
      )
    )
    r_search_records[[length(r_search_records) + 1L]] <<- record
    c(selected, list(R = selected_R, method = method, search_record = record))
  }

  S <- rank[1L]
  non_S <- integer()
  forward_trace <- list()
  forward_consecutive <- 0L
  forward_stopped <- FALSE
  if (length(rank) > 1L) {
    for (position in 2:length(rank)) {
      variable <- rank[position]
      remaining_W <- setdiff(rank, c(S, variable))
      s_candidate <- fit_roles(
        c(S, variable), integer(), integer(), remaining_W,
        context = paste0("forward S candidate for variable ", variable)
      )
      u_candidate <- best_R(
        S, variable, remaining_W,
        context = paste0("forward U candidate for variable ", variable)
      )
      increment <- s_candidate$score - u_candidate$score
      selected_as_S <- increment > 0
      if (selected_as_S) {
        S <- c(S, variable)
        forward_consecutive <- 0L
      } else {
        non_S <- c(non_S, variable)
        forward_consecutive <- forward_consecutive + 1L
      }
      forward_trace[[length(forward_trace) + 1L]] <- list(
        position = position,
        variable = variable,
        score_as_S = s_candidate$score,
        bic_as_S = s_candidate$bic_min,
        score_as_U = u_candidate$score,
        bic_as_U = u_candidate$bic_min,
        R_as_U = u_candidate$R,
        R_search_method = u_candidate$method,
        increment_S_minus_U = increment,
        decision = if (selected_as_S) "S" else "non_S",
        consecutive_nonpositive = forward_consecutive
      )
      if (forward_consecutive >= hsize) {
        if (position < length(rank)) {
          non_S <- c(non_S, rank[(position + 1L):length(rank)])
        }
        forward_stopped <- TRUE
        break
      }
    }
  }
  S <- rank[rank %in% S]
  non_S <- rank[rank %in% unique(non_S)]

  # Reverse pass begins with every non-S variable in U.  A candidate comparison
  # moves exactly one variable U -> W; variables not visited after rolling-c
  # stopping remain U by construction.
  U <- non_S
  W <- integer()
  reverse_trace <- list()
  reverse_consecutive <- 0L
  reverse_stopped <- FALSE
  if (length(non_S)) {
    reverse_order <- rev(non_S)
    for (position in seq_along(reverse_order)) {
      variable <- reverse_order[position]
      current_u <- best_R(
        S, U, W,
        context = paste0("reverse current-U model for variable ", variable)
      )
      moved_u <- U[U != variable]
      moved_w <- c(W, variable)
      moved_to_w <- best_R(
        S, moved_u, moved_w,
        context = paste0("reverse U-to-W candidate for variable ", variable)
      )
      increment <- moved_to_w$score - current_u$score
      selected_as_W <- increment > 0
      if (selected_as_W) {
        U <- moved_u
        W <- moved_w
        reverse_consecutive <- 0L
      } else {
        reverse_consecutive <- reverse_consecutive + 1L
      }
      reverse_trace[[length(reverse_trace) + 1L]] <- list(
        position = position,
        variable = variable,
        score_in_U = current_u$score,
        bic_in_U = current_u$bic_min,
        R_in_U = current_u$R,
        score_in_W = moved_to_w$score,
        bic_in_W = moved_to_w$bic_min,
        R_in_W = moved_to_w$R,
        increment_W_minus_U = increment,
        decision = if (selected_as_W) "W" else "U",
        consecutive_nonpositive = reverse_consecutive
      )
      if (reverse_consecutive >= hsize) {
        reverse_stopped <- TRUE
        break
      }
    }
  }
  U <- rank[rank %in% U]
  W <- rank[rank %in% W]
  final <- best_R(S, U, W, context = "final full-role fit")
  final_fit <- final$fit
  roles <- list(S = S, R = final$R, U = U, W = W)
  search_methods <- unique(vapply(
    r_search_records,
    function(record) record$method,
    character(1)
  ))

  list(
    estimator_mode = "joint_mnarz",
    S = S,
    R = final$R,
    U = U,
    W = W,
    roles = roles,
    K = K,
    rank = rank,
    hsize = hsize,
    covariance = if (diag) "diagonal" else "full",
    final_fit = final_fit,
    loglik_obs = final_fit$loglik_obs,
    bic_min = final_fit$bic_min,
    score = final_fit$score,
    criterionValue = final_fit$criterionValue,
    criterion_convention = final_fit$criterion_convention,
    parameter_count = final_fit$parameter_count,
    parameters = final_fit$parameters,
    proba = final_fit$proba,
    partition = final_fit$partition,
    imputedData = final_fit$imputedData,
    selection_trace = list(
      forward = forward_trace,
      reverse = reverse_trace,
      R_search = r_search_records
    ),
    diagnostics = list(
      converged = TRUE,
      criterion_available = TRUE,
      initial_labels = initial_labels,
      initialization = paste0(
        "one ", initialization, " reused by every candidate"
      ),
      cache_hits = cache_hits,
      cache_misses = cache_misses,
      cached_role_fits = length(ls(cache, all.names = TRUE)),
      r_search_cap = r_search_cap,
      r_search_cap_interpretation = "maximum number of exhaustive R subsets",
      maximum_exhaustive_R_dimension = floor(log(r_search_cap, base = 2)),
      r_search_methods = search_methods,
      r_tie_rule = "strict improvement; first subset in size-then-rank order wins ties",
      forward_stopped_by_rolling_c = forward_stopped,
      reverse_stopped_by_rolling_c = reverse_stopped,
      parallelization = paste0(
        "dependent role candidates are sequential; parallelize only across ",
        "independent K values, starts, or simulation replicates"
      )
    ),
    paper_alignment = list(
      status = "deterministic_total_observed_bic_path",
      forward = paste0(
        "candidate-as-S versus candidate-as-U under one total observed BIC; ",
        "all other variables are held in the same W role so every comparison ",
        "uses the full observed likelihood"
      ),
      reverse = paste0(
        "all non-S variables start in U and are tested in reverse rank order ",
        "for a U-to-W move with common R reselected"
      ),
      stopping = "literal rolling count of consecutive nonpositive increments",
      R_search = paste0(
        "all subsets including empty R when 2^|S| <= r_search_cap; otherwise ",
        "reported deterministic forward search; r_search_cap limits the ",
        "number of exhaustive subset fits"
      ),
      note = paste0(
        "forward and reverse role updates are compared using one total ",
        "observed-data BIC"
      )
    )
  )
}

# Internal entry point shared with the `joint_mnarz` dispatch.
select_joint_sruw_mnarz <- function(
    x,
    K,
    rank,
    hsize,
    diag,
    rmax,
    tol,
    r_search_cap,
    initial_labels = NULL) {
  .select_joint_sruw_mnarz_impl(
    x = x,
    K = K,
    rank = rank,
    hsize = hsize,
    diag = diag,
    rmax = rmax,
    tol = tol,
    r_search_cap = r_search_cap,
    initial_labels = initial_labels
  )
}

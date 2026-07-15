.selvarmix_schema_version <- "1.0"

.selvarmix_value_or <- function(value, fallback) {
  if (is.null(value) || !length(value)) fallback else value
}

.selvarmix_scalar_logical <- function(value) {
  if (is.logical(value) && length(value) == 1L && !is.na(value)) {
    value
  } else {
    NA
  }
}

.selvarmix_result_diagnostics <- function(model) {
  engine <- NULL
  if (is.list(model$clust_result) && is.list(model$clust_result$diagnostics)) {
    engine <- model$clust_result$diagnostics
  } else if (is.list(model$parametersMNARz) &&
             is.list(model$parametersMNARz$diagnostics)) {
    engine <- model$parametersMNARz$diagnostics
  }
  if (is.null(engine)) engine <- list()

  converged <- .selvarmix_scalar_logical(engine$converged)
  criterion_available <- .selvarmix_scalar_logical(
    engine$criterion_available
  )
  if (is.na(criterion_available)) {
    criterion_available <- is.numeric(model$criterionValue) &&
      length(model$criterionValue) == 1L &&
      is.finite(model$criterionValue)
  }
  status <- if (identical(converged, FALSE)) {
    "nonconverged"
  } else if (!isTRUE(criterion_available)) {
    "criterion_unavailable"
  } else if (identical(converged, TRUE)) {
    "converged"
  } else {
    "completed"
  }

  ranking_failures <- attr(model$ranking, "grid_failures", exact = TRUE)
  failed_grid_rows <- if (is.list(ranking_failures)) {
    failure_counts <- vapply(
      ranking_failures,
      function(record) {
        if (is.list(record) && "rows" %in% names(record)) {
          return(length(record$rows))
        }
        if (is.atomic(record) || is.null(record)) {
          return(length(record))
        }
        NA_integer_
      },
      integer(1L)
    )
    if (anyNA(failure_counts)) NA_integer_ else sum(failure_counts)
  } else {
    NA_integer_
  }

  list(
    status = status,
    converged = converged,
    iterations = if (length(engine$iterations) == 1L) {
      as.integer(engine$iterations)
    } else {
      NA_integer_
    },
    terminationReason = if (length(engine$termination_reason) == 1L) {
      as.character(engine$termination_reason)
    } else {
      NA_character_
    },
    criterionAvailable = criterion_available,
    workflow = .selvarmix_value_or(model$workflow, NA_character_),
    criterionConvention = .selvarmix_value_or(
      model$criterionConvention, NA_character_
    ),
    criterionScope = .selvarmix_value_or(
      model$criterionScope,
      .selvarmix_value_or(engine$criterion_scope, NA_character_)
    ),
    covarianceAdjustments = if (length(engine$covariance_adjustments) == 1L) {
      as.integer(engine$covariance_adjustments)
    } else {
      NA_integer_
    },
    minEffectiveComponentSize = if (
      length(engine$min_effective_component_size) == 1L
    ) {
      as.numeric(engine$min_effective_component_size)
    } else {
      NA_real_
    },
    selection = list(
      stoppingRule = .selvarmix_value_or(
        model$stoppingRule, NA_character_
      ),
      stoppingThreshold = .selvarmix_value_or(
        model$stoppingThreshold, NA_integer_
      ),
      evaluated = .selvarmix_value_or(model$nEvaluated, NA_integer_),
      stopReason = .selvarmix_value_or(model$stopReason, NA_character_),
      wStoppingRule = .selvarmix_value_or(
        model$wStoppingRule, NA_character_
      ),
      wStoppingThreshold = .selvarmix_value_or(
        model$wStoppingThreshold, NA_integer_
      ),
      wEvaluated = .selvarmix_value_or(model$wNEvaluated, NA_integer_),
      wStopReason = .selvarmix_value_or(model$wStopReason, NA_character_)
    ),
    ranking = list(
      available = !is.null(model$ranking),
      failedGridRows = as.integer(failed_grid_rows)
    )
  )
}

.selvarmix_role_vector_is_valid <- function(value) {
  if (is.null(value) || !length(value)) return(TRUE)
  if (is.character(value)) {
    return(!anyNA(value) && all(nzchar(value)) && !anyDuplicated(value))
  }
  is.numeric(value) && !anyNA(value) && all(is.finite(value)) &&
    all(value == round(value)) && all(value > 0) && !anyDuplicated(value)
}

.selvarmix_same_role_scale <- function(first, second) {
  if (!length(first) || !length(second)) return(TRUE)
  (is.character(first) && is.character(second)) ||
    (is.numeric(first) && is.numeric(second))
}

.validate_selvarmix_result <- function(object, strict = TRUE) {
  if (!is.list(object) || !inherits(object, "selvarmix")) {
    stop("object must inherit from class 'selvarmix'.", call. = FALSE)
  }
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("strict must be TRUE or FALSE.", call. = FALSE)
  }

  required <- c(
    "S", "R", "U", "W", "criterionValue", "criterion", "model",
    "parameters", "nbcluster", "partition", "proba", "regparameters",
    "imputedData", "workflow", "workflowRequested", "workflowEffective",
    "criterionConvention", "criterionScope"
  )
  missing_fields <- setdiff(required, names(object))
  if (strict && length(missing_fields)) {
    stop(
      "selvarmix result is missing required fields: ",
      paste(missing_fields, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!strict && length(missing_fields)) return(invisible(object))

  roles <- object[c("S", "R", "U", "W")]
  if (!all(vapply(roles, .selvarmix_role_vector_is_valid, logical(1L)))) {
    stop(
      "S, R, U, and W must contain unique positive indices or names.",
      call. = FALSE
    )
  }
  nonempty_roles <- roles[lengths(roles) > 0L]
  if (length(nonempty_roles) > 1L) {
    reference <- nonempty_roles[[1L]]
    if (!all(vapply(
      nonempty_roles,
      .selvarmix_same_role_scale,
      logical(1L),
      second = reference
    ))) {
      stop("Variable roles must use either indices or names consistently.",
           call. = FALSE)
    }
  }
  if (length(intersect(object$S, object$U)) ||
      length(intersect(object$S, object$W)) ||
      length(intersect(object$U, object$W))) {
    stop("S, U, and W must be pairwise disjoint.", call. = FALSE)
  }
  if (length(setdiff(object$R, object$S))) {
    stop("R must be a subset of S.", call. = FALSE)
  }

  nbcluster <- object$nbcluster
  if (!is.numeric(nbcluster) || length(nbcluster) != 1L ||
      !is.finite(nbcluster) || nbcluster != round(nbcluster) ||
      nbcluster < 1L) {
    stop("nbcluster must be one positive integer.", call. = FALSE)
  }
  if (!is.character(object$criterion) || length(object$criterion) != 1L ||
      is.na(object$criterion) || !nzchar(object$criterion)) {
    stop("criterion must be one non-empty character value.", call. = FALSE)
  }
  if (!is.character(object$model) || length(object$model) != 1L ||
      is.na(object$model) || !nzchar(object$model)) {
    stop("model must be one non-empty character value.", call. = FALSE)
  }
  if (!is.numeric(object$criterionValue) ||
      length(object$criterionValue) != 1L) {
    stop("criterionValue must be one numeric value.", call. = FALSE)
  }

  diagnostics <- object$diagnostics
  criterion_available <- if (is.list(diagnostics)) {
    .selvarmix_scalar_logical(diagnostics$criterionAvailable)
  } else {
    NA
  }
  if (!is.finite(object$criterionValue) &&
      !identical(criterion_available, FALSE)) {
    stop(
      "A non-finite criterionValue requires diagnostics$criterionAvailable = FALSE.",
      call. = FALSE
    )
  }

  partition <- object$partition
  if (!is.numeric(partition) || !length(partition) || anyNA(partition) ||
      any(!is.finite(partition)) || any(partition != round(partition)) ||
      any(partition < 1L) || any(partition > nbcluster)) {
    stop(
      "partition must contain component labels from 1 through nbcluster.",
      call. = FALSE
    )
  }
  probabilities <- object$proba
  if (!(is.matrix(probabilities) || is.data.frame(probabilities))) {
    stop("proba must be a numeric probability matrix.", call. = FALSE)
  }
  probabilities <- as.matrix(probabilities)
  if (!(is.numeric(probabilities) || is.logical(probabilities))) {
    stop("proba must be a numeric probability matrix.", call. = FALSE)
  }
  storage.mode(probabilities) <- "double"
  if (
      !identical(dim(probabilities), c(length(partition), as.integer(nbcluster))) ||
      anyNA(probabilities) || any(!is.finite(probabilities)) ||
      any(probabilities < -1e-12) || any(probabilities > 1 + 1e-12) ||
      any(abs(rowSums(probabilities) - 1) > 1e-6)) {
    stop(
      "proba must have one finite probability row per observation and one column per component.",
      call. = FALSE
    )
  }
  assigned_probability <- probabilities[
    cbind(seq_along(partition), as.integer(partition))
  ]
  row_maximum <- apply(probabilities, 1L, max)
  if (any(assigned_probability < row_maximum - 1e-10)) {
    stop(
      "partition must assign each observation to a posterior-probability maximum.",
      call. = FALSE
    )
  }

  completed <- as.matrix(object$imputedData)
  if (!is.numeric(completed) || nrow(completed) != length(partition) ||
      ncol(completed) < 1L || anyNA(completed) || any(!is.finite(completed))) {
    stop(
      "imputedData must be a finite complete matrix with one row per observation.",
      call. = FALSE
    )
  }
  primary_roles <- unique(c(object$S, object$U, object$W))
  if (is.numeric(primary_roles)) {
    if (!identical(sort(as.integer(primary_roles)), seq_len(ncol(completed)))) {
      stop("S, U, and W must partition the columns of imputedData.",
           call. = FALSE)
    }
  } else {
    completed_names <- colnames(completed)
    if (is.null(completed_names) || anyNA(completed_names) ||
        any(!nzchar(completed_names)) || anyDuplicated(completed_names) ||
        !identical(sort(primary_roles), sort(completed_names))) {
      stop(
        "Character-valued S, U, and W must partition colnames(imputedData).",
        call. = FALSE
      )
    }
  }

  if (strict) {
    if (!identical(object$schemaVersion, .selvarmix_schema_version)) {
      stop(
        "schemaVersion must equal '", .selvarmix_schema_version, "'.",
        call. = FALSE
      )
    }
    if (!is.list(object$diagnostics) ||
        !all(c(
          "status", "converged", "iterations", "terminationReason",
          "criterionAvailable", "workflow", "criterionConvention",
          "criterionScope", "selection", "ranking"
        ) %in% names(object$diagnostics))) {
      stop("diagnostics does not satisfy the selvarmix 1.0 schema.",
           call. = FALSE)
    }
  }
  invisible(object)
}

.new_selvarmix_result <- function(fields, validate = TRUE) {
  if (!is.list(fields)) {
    stop("fields must be a list.", call. = FALSE)
  }
  if (is.null(fields$schemaVersion)) {
    fields$schemaVersion <- .selvarmix_schema_version
  }
  if (is.null(fields$diagnostics)) {
    fields$diagnostics <- .selvarmix_result_diagnostics(fields)
  }
  class(fields) <- "selvarmix"
  if (isTRUE(validate)) .validate_selvarmix_result(fields, strict = TRUE)
  fields
}

.validate_selvarmix_collection <- function(object) {
  if (!is.list(object) || !inherits(object, "selvarmix_collection") ||
      !length(object)) {
    stop(
      "object must be a non-empty 'selvarmix_collection'.",
      call. = FALSE
    )
  }
  if (is.null(names(object)) || any(!nzchar(names(object))) ||
      anyDuplicated(names(object))) {
    stop("A selvarmix_collection requires unique criterion names.",
         call. = FALSE)
  }
  invisible(lapply(object, .validate_selvarmix_result, strict = TRUE))
}

.new_selvarmix_collection <- function(results) {
  if (!is.list(results) || !length(results)) {
    stop("results must be a non-empty list.", call. = FALSE)
  }
  class(results) <- c("selvarmix_collection", "list")
  .validate_selvarmix_collection(results)
  results
}

# Convert backend result lists to the validated public schema.
selvarmix <- function(bestModel) {
  if (inherits(bestModel, "selvarmix")) {
    .validate_selvarmix_result(bestModel, strict = TRUE)
    return(bestModel)
  }
  if (!is.list(bestModel) || !length(bestModel)) {
    stop("bestModel must be a non-empty model list.", call. = FALSE)
  }
  models <- if (!is.null(bestModel$S)) list(bestModel) else bestModel
  results <- lapply(models, function(model) {
    if (inherits(model, "selvarmix")) model else ProcessModelOutput(model)
  })
  if (length(results) == 1L) return(results[[1L]])
  if (is.null(names(results)) || any(!nzchar(names(results)))) {
    names(results) <- paste0("result", seq_along(results))
  }
  .new_selvarmix_collection(results)
}

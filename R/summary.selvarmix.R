.selvarmix_component_count <- function(object) {
  value <- .selvarmix_value_or(object$nbcluster, object$nbCluster)
  if (length(value) == 1L && is.numeric(value) && is.finite(value)) {
    as.integer(value)
  } else {
    NA_integer_
  }
}

.selvarmix_variable_count <- function(object) {
  if (!is.null(object$imputedData)) return(as.integer(ncol(object$imputedData)))
  primary <- unique(c(object$S, object$U, object$W))
  if (!length(primary)) return(0L)
  if (is.numeric(primary)) return(as.integer(max(primary)))
  as.integer(length(primary))
}

.selvarmix_format_number <- function(value, digits = 4L) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    return("NA")
  }
  trimws(formatC(value, digits = as.integer(digits), format = "g"))
}

.selvarmix_format_set <- function(value, max_items = 12L) {
  if (is.null(value) || !length(value)) return("none")
  max_items <- max(1L, as.integer(max_items))
  shown <- utils::head(value, max_items)
  text <- paste(shown, collapse = ", ")
  if (length(value) > max_items) {
    text <- paste0(text, ", ... (+", length(value) - max_items, ")")
  }
  text
}

.selvarmix_text_or_unknown <- function(value) {
  if (is.null(value) || !length(value) || is.na(value[[1L]]) ||
      !nzchar(as.character(value[[1L]]))) {
    "not recorded"
  } else {
    as.character(value[[1L]])
  }
}

#' Summarize a SelvarMixMNAR fit
#'
#' @param object A fitted object inheriting from class `selvarmix`.
#' @param ... Reserved for extensions.
#' @return A structured object of class `summary.selvarmix`.
#' @export
summary.selvarmix <- function(object, ...) {
  .validate_selvarmix_result(object, strict = FALSE)
  diagnostics <- if (is.list(object$diagnostics)) {
    object$diagnostics
  } else {
    .selvarmix_result_diagnostics(object)
  }
  roles <- list(
    S = .selvarmix_value_or(object$S, vector()),
    R = .selvarmix_value_or(object$R, vector()),
    U = .selvarmix_value_or(object$U, vector()),
    W = .selvarmix_value_or(object$W, vector())
  )
  model <- list(
    requested = .selvarmix_value_or(object$requestedModel, object$model),
    effective = .selvarmix_value_or(object$effectiveModel, object$model),
    framework = .selvarmix_value_or(object$framework, NA_character_),
    regression = .selvarmix_value_or(object$rmodel, NA_character_),
    independent = .selvarmix_value_or(object$imodel, NA_character_)
  )
  output <- list(
    schemaVersion = .selvarmix_value_or(object$schemaVersion, "legacy"),
    workflow = .selvarmix_value_or(object$workflow, NA_character_),
    criterion = .selvarmix_value_or(object$criterion, NA_character_),
    criterionValue = .selvarmix_value_or(object$criterionValue, NA_real_),
    criterionConvention = .selvarmix_value_or(
      object$criterionConvention,
      .selvarmix_value_or(diagnostics$criterionConvention, NA_character_)
    ),
    criterionScope = .selvarmix_value_or(
      object$criterionScope,
      .selvarmix_value_or(diagnostics$criterionScope, NA_character_)
    ),
    nbcluster = .selvarmix_component_count(object),
    nobs = as.integer(length(object$partition)),
    nvar = .selvarmix_variable_count(object),
    model = model,
    roles = roles,
    roleCounts = vapply(roles, length, integer(1L)),
    completedDataReturned = !is.null(object$imputedData),
    diagnostics = diagnostics,
    selectionCriterionValue = .selvarmix_value_or(
      object$selectionCriterionValue, NA_real_
    )
  )
  class(output) <- c("summary.selvarmix", "list")
  output
}

#' @export
print.summary.selvarmix <- function(x, ..., digits = 4L, max.roles = 12L) {
  if (!inherits(x, "summary.selvarmix")) {
    stop("x must inherit from class 'summary.selvarmix'.", call. = FALSE)
  }
  cat("SelvarMixMNAR fit summary\n")
  cat("  Workflow: ", .selvarmix_text_or_unknown(x$workflow), "\n", sep = "")
  cat(
    "  Criterion: ", .selvarmix_text_or_unknown(x$criterion), " = ",
    .selvarmix_format_number(x$criterionValue, digits),
    if (!is.na(x$criterionConvention)) {
      paste0(" (", x$criterionConvention, ")")
    } else {
      ""
    },
    "\n",
    sep = ""
  )
  cat("  Components: ", x$nbcluster, "\n", sep = "")
  cat("  Observations: ", x$nobs, "\n", sep = "")
  cat("  Variables: ", x$nvar, "\n", sep = "")
  cat(
    "  Model: ", .selvarmix_text_or_unknown(x$model$effective),
    if (!is.na(x$model$framework)) paste0(" [", x$model$framework, "]") else "",
    "\n",
    sep = ""
  )
  cat("\nVariable roles\n")
  role_labels <- c(
    S = "S (clustering)",
    R = "R (regression predictors within S)",
    U = "U (redundant)",
    W = "W (independent)"
  )
  for (role in names(role_labels)) {
    cat(
      "  ", role_labels[[role]], ": ",
      .selvarmix_format_set(x$roles[[role]], max.roles), "\n",
      sep = ""
    )
  }
  cat("\nDiagnostics\n")
  cat(
    "  Status: ", .selvarmix_text_or_unknown(x$diagnostics$status), "\n",
    sep = ""
  )
  if (!is.na(x$diagnostics$converged)) {
    cat(
      "  Converged: ", if (isTRUE(x$diagnostics$converged)) "yes" else "no",
      "\n",
      sep = ""
    )
  }
  if (!is.na(x$diagnostics$iterations)) {
    cat("  Iterations: ", x$diagnostics$iterations, "\n", sep = "")
  }
  if (!is.na(x$diagnostics$terminationReason)) {
    cat(
      "  Termination: ", x$diagnostics$terminationReason, "\n",
      sep = ""
    )
  }
  if (!isTRUE(x$diagnostics$criterionAvailable)) {
    cat("  Criterion available: no\n")
  }
  invisible(x)
}

#' @export
summary.selvarmix_collection <- function(object, ...) {
  .validate_selvarmix_collection(object)
  summaries <- lapply(object, summary, ...)
  class(summaries) <- c("summary.selvarmix_collection", "list")
  summaries
}

.selvarmix_collection_table <- function(object, digits = 4L) {
  data.frame(
    Result = names(object),
    Workflow = vapply(
      object,
      function(x) .selvarmix_text_or_unknown(x$workflow),
      character(1L)
    ),
    Criterion = vapply(
      object,
      function(x) .selvarmix_text_or_unknown(x$criterion),
      character(1L)
    ),
    Value = vapply(
      object,
      function(x) .selvarmix_format_number(x$criterionValue, digits),
      character(1L)
    ),
    K = vapply(object, function(x) as.integer(x$nbcluster), integer(1L)),
    S = vapply(object, function(x) length(x$roles$S), integer(1L)),
    U = vapply(object, function(x) length(x$roles$U), integer(1L)),
    W = vapply(object, function(x) length(x$roles$W), integer(1L)),
    Status = vapply(
      object,
      function(x) .selvarmix_text_or_unknown(x$diagnostics$status),
      character(1L)
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' @export
print.summary.selvarmix_collection <- function(x, ..., digits = 4L) {
  if (!inherits(x, "summary.selvarmix_collection")) {
    stop(
      "x must inherit from class 'summary.selvarmix_collection'.",
      call. = FALSE
    )
  }
  cat("SelvarMixMNAR model collection\n")
  print(.selvarmix_collection_table(x, digits), row.names = FALSE)
  invisible(x)
}

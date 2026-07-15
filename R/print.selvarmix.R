#' Print a SelvarMixMNAR fit
#'
#' @param x A fitted object inheriting from class `selvarmix`.
#' @param ... Reserved for extensions.
#' @param digits Number of significant digits used for the criterion.
#' @param max.roles Maximum number of indices printed for each variable role.
#' @return `x`, invisibly.
#' @export
print.selvarmix <- function(x, ..., digits = 4L, max.roles = 12L) {
  .validate_selvarmix_result(x, strict = FALSE)
  fit_summary <- summary.selvarmix(x)
  cat("SelvarMixMNAR fit\n")
  cat(
    "  ", .selvarmix_text_or_unknown(fit_summary$criterion), " = ",
    .selvarmix_format_number(fit_summary$criterionValue, digits),
    if (!is.na(fit_summary$criterionConvention)) {
      paste0(" (", fit_summary$criterionConvention, ")")
    } else {
      ""
    },
    "\n",
    sep = ""
  )
  cat(
    "  Workflow: ", .selvarmix_text_or_unknown(fit_summary$workflow), "\n",
    sep = ""
  )
  cat(
    "  Model: ", .selvarmix_text_or_unknown(fit_summary$model$effective),
    "; components: ", fit_summary$nbcluster, "\n",
    sep = ""
  )
  cat(
    "  Roles: S={", .selvarmix_format_set(fit_summary$roles$S, max.roles),
    "}; U={", .selvarmix_format_set(fit_summary$roles$U, max.roles),
    "}; W={", .selvarmix_format_set(fit_summary$roles$W, max.roles),
    "}\n",
    sep = ""
  )
  if (length(fit_summary$roles$R)) {
    cat(
      "  Regression predictors R within S: {",
      .selvarmix_format_set(fit_summary$roles$R, max.roles), "}\n",
      sep = ""
    )
  }
  cat(
    "  Status: ", .selvarmix_text_or_unknown(fit_summary$diagnostics$status),
    "\n",
    sep = ""
  )
  invisible(x)
}

#' @export
print.selvarmix_collection <- function(x, ..., digits = 4L) {
  .validate_selvarmix_collection(x)
  print(summary.selvarmix_collection(x), digits = digits)
  invisible(x)
}

CheckInputsC <- function(x, nbcluster, lambda, rho, type, hsize, criterion,
                         models, rmodel, imodel, nbcores) {
  if (missing(x) || (!is.matrix(x) && !is.data.frame(x))) {
    stop("x must be a numeric matrix or data frame.", call. = FALSE)
  }
  x_matrix <- as.matrix(x)
  if (!is.numeric(x_matrix) || !length(x_matrix) || nrow(x_matrix) < 2L ||
      ncol(x_matrix) < 1L) {
    stop("x must contain at least two observations and one numeric variable.",
         call. = FALSE)
  }
  if (any(is.nan(x_matrix)) || any(is.infinite(x_matrix))) {
    stop("x may contain NA values but not NaN or infinite values.",
         call. = FALSE)
  }

  if (missing(nbcluster) || !is.numeric(nbcluster) || !length(nbcluster) ||
      any(!is.finite(nbcluster)) || any(!is.wholenumber(nbcluster)) ||
      any(nbcluster < 1L) || any(nbcluster > nrow(x_matrix))) {
    stop("nbcluster must contain positive integers no larger than nrow(x).",
         call. = FALSE)
  }
  if (anyDuplicated(nbcluster)) {
    stop("nbcluster values must be unique.", call. = FALSE)
  }

  validate_penalty <- function(value, name) {
    if (is.null(value)) return(invisible(NULL))
    if (!is.numeric(value) || !length(value) || any(!is.finite(value)) ||
        any(value < 0)) {
      stop(name, " must be NULL or a finite non-negative numeric vector.",
           call. = FALSE)
    }
    invisible(NULL)
  }
  validate_penalty(lambda, "lambda")
  validate_penalty(rho, "rho")

  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("type must be one non-missing character token.", call. = FALSE)
  }
  if (identical(tolower(type), "likelihood")) {
    condition <- structure(
      list(
        message = paste0(
          "type='likelihood' is unavailable: SelvarMixMNAR has no ",
          "likelihood-based variable-ranking estimator. Use type='lasso'."
        ),
        call = NULL
      ),
      class = c(
        "selvarmix_ranking_method_unavailable", "error", "condition"
      )
    )
    stop(condition)
  }
  if (!identical(tolower(type), "lasso")) {
    stop("type must be exactly 'lasso'.", call. = FALSE)
  }

  if (length(hsize) != 1L || !is.numeric(hsize) || !is.finite(hsize) ||
      !is.wholenumber(hsize) || hsize < 1L || hsize > ncol(x_matrix)) {
    stop("hsize must be one positive integer no larger than ncol(x).",
         call. = FALSE)
  }
  if (!is.character(criterion) || !length(criterion) || anyNA(criterion) ||
      any(!criterion %in% c("BIC", "ICL"))) {
    stop("criterion must contain only 'BIC' and/or 'ICL'.", call. = FALSE)
  }

  for (requested_criterion in criterion) {
    .sruw_backend_specification(
      models, requested_criterion, supervised = FALSE
    )
  }
  if (!is.character(rmodel) || !length(rmodel) || anyNA(rmodel) ||
      any(!rmodel %in% c("LI", "LB", "LC"))) {
    stop("rmodel must contain only LI, LB, and/or LC.", call. = FALSE)
  }
  if (!is.character(imodel) || !length(imodel) || anyNA(imodel) ||
      any(!imodel %in% c("LI", "LB"))) {
    stop("imodel must contain only LI and/or LB.", call. = FALSE)
  }
  if (length(nbcores) != 1L || !is.numeric(nbcores) || !is.finite(nbcores) ||
      !is.wholenumber(nbcores) || nbcores < 1L) {
    stop("nbcores must be one positive integer.", call. = FALSE)
  }
  invisible(TRUE)
}

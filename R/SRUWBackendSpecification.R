.sruw_backend_condition <- function(message, subclass, ...) {
  structure(
    c(list(message = message, call = NULL), list(...)),
    class = c(
      subclass,
      "selvarmix_sruw_backend_error",
      "error",
      "condition"
    )
  )
}

.stop_sruw_backend <- function(message, subclass, ...) {
  stop(.sruw_backend_condition(message, subclass, ...))
}

.sruw_mclust_models <- function() {
  c(
    "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE",
    "EVE", "VVE", "EEV", "VEV", "EVV", "VVV"
  )
}

.sruw_mixall_models <- function() {
  c(
    "gaussian_pk_sjk", "gaussian_pk_sj", "gaussian_pk_sk",
    "gaussian_pk_s", "gaussian_p_sjk", "gaussian_p_sj",
    "gaussian_p_sk", "gaussian_p_s"
  )
}

.sruw_rmixmod_family <- function(model) {
  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    return(NA_character_)
  }

  compact <- gsub("[[:space:]]+", "", model)
  if (identical(compact, "mixmodGaussianModel()")) {
    return("general")
  }

  families <- c("general", "diagonal", "spherical", "all")
  double_quoted <- sprintf(
    "mixmodGaussianModel(family=\"%s\")", families
  )
  single_quoted <- sprintf(
    "mixmodGaussianModel(family='%s')", families
  )
  match_index <- match(compact, c(double_quoted, single_quoted))
  if (is.na(match_index)) {
    return(NA_character_)
  }
  families[((match_index - 1L) %% length(families)) + 1L]
}

.sruw_resolve_model <- function(models) {
  if (!is.character(models) || length(models) != 1L || is.na(models) ||
      !nzchar(models)) {
    .stop_sruw_backend(
      "The SRUW model must be one non-empty character token.",
      "selvarmix_sruw_model_unavailable",
      model = models
    )
  }

  if (models %in% .sruw_mclust_models()) {
    return(list(
      framework = "Mclust", model = models, package = "mclust",
      rmixmod_family = NULL
    ))
  }
  if (models %in% .sruw_mixall_models()) {
    return(list(
      framework = "MixAll", model = models, package = "MixAll",
      rmixmod_family = NULL
    ))
  }

  family <- .sruw_rmixmod_family(models)
  if (!is.na(family)) {
    return(list(
      framework = "Rmixmod",
      model = sprintf("mixmodGaussianModel(family=\"%s\")", family),
      package = "Rmixmod",
      rmixmod_family = family
    ))
  }

  .stop_sruw_backend(
    paste0(
      "Unsupported SRUW model token '", models, "'. Use an Mclust model ",
      "name, a MixAll Gaussian model name, or ",
      "mixmodGaussianModel(family=\"<general|diagonal|spherical|all>\")."
    ),
    "selvarmix_sruw_model_unavailable",
    model = models
  )
}

.sruw_backend_specification <- function(models, criterion,
                                        supervised = FALSE) {
  resolved <- .sruw_resolve_model(models)
  if (!is.character(criterion) || length(criterion) != 1L ||
      is.na(criterion) || !nzchar(criterion)) {
    .stop_sruw_backend(
      "The SRUW criterion must be one non-empty character token.",
      "selvarmix_sruw_criterion_unavailable",
      framework = resolved$framework,
      criterion = criterion
    )
  }
  if (!is.logical(supervised) || length(supervised) != 1L ||
      is.na(supervised)) {
    .stop_sruw_backend(
      "The SRUW supervised flag must be one non-missing logical value.",
      "selvarmix_sruw_supervision_unavailable",
      framework = resolved$framework,
      supervised = supervised
    )
  }

  if (identical(resolved$framework, "MixAll")) {
    .stop_sruw_backend(
      paste0(
        "MixAll model '", resolved$model, "' is not executed in-process ",
        "because an STK failure may terminate the R session. Choose an ",
        "Mclust or Rmixmod model."
      ),
      "selvarmix_sruw_backend_uncontained",
      framework = resolved$framework,
      model = resolved$model,
      containment = "unavailable_in_process"
    )
  }

  if (identical(resolved$framework, "Mclust")) {
    if (isTRUE(supervised)) {
      .stop_sruw_backend(
        "Mclust SRUW implements unsupervised fitting only.",
        "selvarmix_sruw_supervision_unavailable",
        framework = resolved$framework,
        supervised = TRUE
      )
    }
    supported_criteria <- "BIC"
  } else if (identical(resolved$framework, "Rmixmod")) {
    supported_criteria <- if (isTRUE(supervised)) {
      c("BIC", "CV")
    } else {
      c("BIC", "ICL")
    }
  } else {
    supported_criteria <- c("BIC", "ICL")
  }

  if (!criterion %in% supported_criteria) {
    .stop_sruw_backend(
      paste0(
        resolved$framework, " SRUW does not implement criterion '",
        criterion, "'. Available choices: ",
        paste(supported_criteria, collapse = ", "), "."
      ),
      "selvarmix_sruw_criterion_unavailable",
      framework = resolved$framework,
      criterion = criterion,
      supported_criteria = supported_criteria
    )
  }

  resolved$criterion <- criterion
  resolved$supervised <- supervised
  resolved
}

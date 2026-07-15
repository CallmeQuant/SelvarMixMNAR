.sruw_refit_empty_r <- function(fit, selection_fit, data, imodel) {
  if (!is.list(fit) || !is.null(fit$error)) return(fit)
  if (length(fit$R) || !length(fit$U)) return(fit)

  target_w <- sort(unique(c(as.integer(fit$U), as.integer(fit$W))))
  if (!length(target_w)) return(fit)
  clustering_criterion <- as.numeric(selection_fit$criterionValue)[1L]
  if (!is.finite(clustering_criterion)) {
    stop("The clustering contribution is unavailable for the empty-R refit.",
         call. = FALSE)
  }

  candidates <- lapply(imodel, function(model) {
    tryCatch({
      independent <- rcppRegressionBIC(
        data, response = target_w, predictors = integer(), model = model
      )
      value <- clustering_criterion + as.numeric(independent$bicvalue)[1L]
      if (!is.finite(value)) stop("non-finite independent criterion")
      list(
        valid = TRUE,
        model = model,
        criterionValue = value,
        independentCriterionValue = as.numeric(independent$bicvalue)[1L]
      )
    }, error = function(e) {
      list(valid = FALSE, model = model, error = conditionMessage(e))
    })
  })
  valid <- which(vapply(candidates, `[[`, logical(1), "valid"))
  if (!length(valid)) {
    stop(
      "No independent covariance model could refit U union W after R became empty: ",
      paste(vapply(candidates, function(x) x$error, character(1)),
            collapse = " | "),
      call. = FALSE
    )
  }
  values <- vapply(candidates[valid], `[[`, numeric(1), "criterionValue")
  selected <- candidates[[valid[which.max(values)]]]

  fit$emptyRRefit <- list(
    applied = TRUE,
    reason = paste0(
      "U conditional on an empty R is an independent block; U and W were ",
      "jointly refitted before comparing K candidates"
    ),
    previousCriterionValue = fit$criterionValue,
    clusteringCriterionValue = clustering_criterion,
    independentCriterionValue = selected$independentCriterionValue,
    selectedIndependentModel = selected$model,
    candidateFits = candidates
  )
  fit$R <- integer()
  fit$U <- integer()
  fit$W <- target_w
  fit$rmodel <- ""
  fit$imodel <- selected$model
  fit$regparameters <- matrix(numeric(), nrow = 0L, ncol = 0L)
  fit$criterionValue <- selected$criterionValue
  fit
}

ModelSelectionClust <- function(VariableSelectRes,
                                data,
                                rmodel,
                                imodel,
                                nbcores,
                                rng_seed = 1L) {
  if (!is.list(VariableSelectRes) || !length(VariableSelectRes)) {
    stop("VariableSelectRes must contain at least one selection candidate.",
         call. = FALSE)
  }
  data <- as.matrix(data)
  if (!is.null(rng_seed) &&
      (length(rng_seed) != 1L || !is.numeric(rng_seed) ||
       !is.finite(rng_seed) || rng_seed != as.integer(rng_seed) ||
       rng_seed < 0 || rng_seed > .Machine$integer.max)) {
    stop("rng_seed must be NULL or one non-negative integer.", call. = FALSE)
  }
  if (!is.null(rng_seed)) rng_seed <- as.integer(rng_seed)
  candidate_count <- length(VariableSelectRes)
  nbcores <- min(max(1L, as.integer(nbcores)), candidate_count)
  task_seeds <- lapply(seq_len(candidate_count), function(index) {
    .sruw_task_seed(rng_seed, index)
  })

  evaluate_candidate <- function(index) {
    .sruw_with_seed(task_seeds[[index]], {
      selection_fit <- VariableSelectRes[[index]]
      fit <- tryCatch(
        rcppCrit(data, selection_fit, rmodel, imodel),
        error = function(e) list(error = conditionMessage(e))
      )
      if (!is.list(fit) || !is.null(fit$error)) {
        message <- if (is.list(fit) && !is.null(fit$error)) {
          as.character(fit$error)[1L]
        } else {
          "rcppCrit returned an invalid object"
        }
        return(list(error = message, candidate = index))
      }
      fit <- tryCatch(
        .sruw_refit_empty_r(fit, selection_fit, data, imodel),
        error = function(e) list(error = conditionMessage(e), candidate = index)
      )
      if (is.null(fit$error)) fit$rngSeed <- task_seeds[[index]]
      fit
    })
  }

  evaluate_sequentially <- function() {
    lapply(seq_len(candidate_count), evaluate_candidate)
  }
  results <- NULL
  if (nbcores > 1L && candidate_count > 1L) {
    if (.Platform$OS.type == "windows") {
      cl <- tryCatch(parallel::makeCluster(nbcores), error = identity)
      if (!inherits(cl, "error")) {
        cluster_active <- TRUE
        on.exit({
          if (cluster_active) try(parallel::stopCluster(cl), silent = TRUE)
        }, add = TRUE)
        parallel_result <- tryCatch({
          parallel::clusterEvalQ(cl, {
            if (!requireNamespace("SelvarMixMNAR", quietly = TRUE)) {
              stop("SelvarMixMNAR namespace is unavailable on a model worker.")
            }
            TRUE
          })
          parallel::clusterExport(
            cl,
            c(
              "VariableSelectRes", "data", "rmodel", "imodel", "task_seeds",
              "rcppCrit", "rcppRegressionBIC", ".sruw_with_seed",
              ".sruw_refit_empty_r", "evaluate_candidate"
            ),
            envir = environment()
          )
          parallel::parLapply(cl, seq_len(candidate_count), evaluate_candidate)
        }, error = identity)
        try(parallel::stopCluster(cl), silent = TRUE)
        cluster_active <- FALSE
        if (!inherits(parallel_result, "error")) {
          results <- parallel_result
        } else {
          warning(
            "Parallel model evaluation failed; rerunning the complete candidate set sequentially: ",
            conditionMessage(parallel_result),
            call. = FALSE
          )
        }
      } else {
        warning(
          "Parallel cluster creation failed; evaluating models sequentially: ",
          conditionMessage(cl),
          call. = FALSE
        )
      }
    } else {
      parallel_result <- tryCatch(
        parallel::mclapply(
          seq_len(candidate_count), evaluate_candidate,
          mc.cores = nbcores, mc.silent = TRUE, mc.preschedule = TRUE,
          mc.cleanup = TRUE, mc.set.seed = FALSE
        ),
        error = identity
      )
      if (!inherits(parallel_result, "error")) {
        results <- parallel_result
      } else {
        warning(
          "Parallel model evaluation failed; rerunning sequentially: ",
          conditionMessage(parallel_result),
          call. = FALSE
        )
      }
    }
  }
  if (is.null(results)) results <- evaluate_sequentially()

  valid <- which(vapply(results, function(result) {
    is.list(result) && is.null(result$error) &&
      length(result$criterionValue) == 1L &&
      is.finite(result$criterionValue)
  }, logical(1)))
  if (!length(valid)) {
    messages <- vapply(results, function(result) {
      if (is.list(result) && !is.null(result$error)) {
        as.character(result$error)[1L]
      } else {
        "invalid or non-finite criterion"
      }
    }, character(1))
    stop(
      "No valid models found. Candidate failures: ",
      paste(messages, collapse = " | "),
      call. = FALSE
    )
  }
  criterion_values <- vapply(
    results[valid], function(result) as.numeric(result$criterionValue), numeric(1)
  )
  bestModel <- results[[valid[which.max(criterion_values)]]]
  if (!length(bestModel$R)) bestModel$R <- NULL
  if (!length(bestModel$U)) bestModel$U <- NULL
  if (!length(bestModel$W)) bestModel$W <- NULL
  bestModel
}

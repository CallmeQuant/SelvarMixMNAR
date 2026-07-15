.make_sruw_cluster <- function(nbcores) {
  parallel::makeCluster(nbcores)
}

.stop_sruw_cluster <- function(cl) {
  parallel::stopCluster(cl)
}

.sruw_task_seed <- function(base_seed, index) {
  if (is.null(base_seed)) return(NULL)
  modulus <- as.double(.Machine$integer.max) + 1
  as.integer((as.double(base_seed) + as.double(index) - 1) %% modulus)
}

.sruw_with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.run_windows_sruw_parallel <- function(cl, input_list, export_env,
                                       backend_package) {
  parallel::clusterExport(
    cl,
    varlist = c(
      "data", "framework", "model_name", "hsize", "criterion",
      "supervised", "z", "stopping", "rcppSelectS", "rcppSelectW",
      "wrapper.selectVar", "backend_package", ".sruw_with_seed"
    ),
    envir = export_env
  )

  parallel::clusterEvalQ(cl, {
    if (!requireNamespace("SelvarMixMNAR", quietly = TRUE)) {
      stop("SelvarMixMNAR namespace is unavailable on an SRUW worker.")
    }
    if (!requireNamespace(backend_package, quietly = TRUE)) {
      stop("Required SRUW backend package is unavailable: ", backend_package)
    }
    TRUE
  })

  parallel::parLapply(cl, input_list, function(x) {
    worker <- get("wrapper.selectVar", envir = .GlobalEnv, inherits = FALSE)
    tryCatch(
      worker(x$nbcluster, x$ordervar, x$task_seed),
      error = function(e) {
        list(
          error = paste("Parallel SRUW selection failed:", e$message),
          nbcluster = x$nbcluster
        )
      }
    )
  })
}

VariableSelection <- function(data,
                              nbcluster,
                              models,
                              criterion,
                              OrderVariable,
                              hsize,
                              supervised,
                              z,
                              nbcores,
                              stopping = c("consecutive", "legacy_block"),
                              rng_seed = 1L,
                              verbose = FALSE) {
  data <- as.matrix(data)
  nbcluster <- as.integer(nbcluster)
  criterion <- as.character(criterion)
  stopping <- match.arg(stopping)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(rng_seed) &&
      (length(rng_seed) != 1L || !is.numeric(rng_seed) ||
       !is.finite(rng_seed) || rng_seed != as.integer(rng_seed) ||
       rng_seed < 0 || rng_seed > .Machine$integer.max)) {
    stop("rng_seed must be NULL or one non-negative integer.", call. = FALSE)
  }
  if (!is.null(rng_seed)) rng_seed <- as.integer(rng_seed)
  p <- ncol(data)

  backend_spec <- .sruw_backend_specification(
    models, criterion, supervised
  )
  framework <- backend_spec$framework
  model_name <- backend_spec$model
  backend_package <- backend_spec$package
  if (!requireNamespace(backend_package, quietly = TRUE)) {
    .stop_sruw_backend(
      paste0("Required SRUW backend package is unavailable: ", backend_package),
      "selvarmix_sruw_backend_unavailable",
      framework = framework,
      package = backend_package
    )
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    .stop_sruw_backend(
      "The required parallel package is unavailable.",
      "selvarmix_sruw_backend_unavailable",
      package = "parallel"
    )
  }

  stop_sruw_fit <- function(results, input_list) {
    failed <- which(vapply(
      results,
      function(result) !is.list(result) || !is.null(result$error),
      logical(1)
    ))
    if (!length(failed)) return(invisible(NULL))
    messages <- vapply(failed, function(i) {
      result <- results[[i]]
      if (is.list(result) && !is.null(result$error)) {
        as.character(result$error)[1L]
      } else {
        "invalid SRUW result"
      }
    }, character(1))
    failed_k <- vapply(
      failed,
      function(i) as.integer(input_list[[i]]$nbcluster),
      integer(1)
    )
    backend_error_classes <- vapply(
      failed,
      function(i) {
        result <- results[[i]]
        if (!is.list(result)) return(NA_character_)
        value <- result$errorClass
        if (is.null(value)) NA_character_ else as.character(value)[1L]
      },
      character(1)
    )
    condition_classes <- unique(backend_error_classes[!is.na(
      backend_error_classes
    )])
    condition <- structure(
      list(
        message = sprintf(
          "SRUW selection failed for K=%s; incomplete candidate sets were not scored: %s",
          paste(failed_k, collapse = ","),
          paste(messages, collapse = " | ")
        ),
        call = NULL,
        failed_k = failed_k,
        errors = messages,
        stopping = stopping,
        backend_error_classes = backend_error_classes
      ),
      class = c(
        "selvarmix_sruw_fit_error", condition_classes, "error", "condition"
      )
    )
    stop(condition)
  }
  
  # Each candidate K uses a deterministic substream and its corresponding
  # ranked-variable order, so candidates remain independent across K.
  wrapper.selectVar <- function(current_nbcluster, current_ordervar,
                                task_seed = NULL) {
    .sruw_with_seed(task_seed, {
    result <- tryCatch({
      rcppSelectS(data,
                  current_ordervar,
                  current_nbcluster,
                  framework, model_name,
                  hsize, criterion,
                  as.integer(z), supervised,
                  stopping) 
    }, error = function(e) {
       list(error = paste("Error in rcppSelectS call:", e$message),
            S = integer(0),
            W = integer(0),
            U = integer(0)) 
    })

    if (!is.null(result$error) || is.null(result$S)) {
       return(list(
         S = integer(0), W = integer(0), U = 1:p,
         error = result$error,
         errorClass = result$errorClass,
         nbcluster = current_nbcluster
       ))
    }
    
    # Conditional on S, scan the complement in reverse order to identify W.
    OrderAux <- setdiff(current_ordervar, result$S)
    w_result <- tryCatch(
      list(value = rcppSelectW(data, OrderAux, result$S, hsize, stopping)),
      error = function(e) list(error = paste("Error in rcppSelectW call:", e$message))
    )
    if (!is.null(w_result$error)) {
      return(list(
        S = integer(0), W = integer(0), U = integer(0),
        error = w_result$error, nbcluster = current_nbcluster
      ))
    }
    w_value <- w_result$value
    result$wStoppingRule <- attr(w_value, "stoppingRule", exact = TRUE)
    result$wStoppingThreshold <- attr(
      w_value, "stoppingThreshold", exact = TRUE
    )
    result$wNEvaluated <- attr(w_value, "nEvaluated", exact = TRUE)
    result$wStopReason <- attr(w_value, "stopReason", exact = TRUE)
    result$W <- as.integer(w_value)
    
    # Variables in neither S nor W form the redundant block U.
    result$U <- setdiff(seq_len(ncol(data)), union(result$S, result$W))

    result$nbcluster_run <- current_nbcluster
    result$rngSeed <- task_seed
    return(result)
    })
  }

  max_cores <- parallel::detectCores(logical = FALSE)
  if (length(max_cores) != 1L || !is.finite(max_cores) || max_cores < 1L) {
    max_cores <- 1L
  }
  nbcores <- min(max(1L, as.integer(nbcores)), max_cores)
  # Bind each candidate K to its ranking and reproducible task seed.
  input_list <- lapply(seq_along(nbcluster), function(i) {
    list(
      nbcluster = nbcluster[i],
      ordervar = if (is.matrix(OrderVariable) && nrow(OrderVariable) >= length(nbcluster)) OrderVariable[i, ] else if (is.matrix(OrderVariable)) OrderVariable[1, ] else OrderVariable,
      task_seed = .sruw_task_seed(rng_seed, i)
    )
  })
  results <- list()

  if (nbcores > 1) {
    if (.Platform$OS.type == "windows") {
      cl <- tryCatch({
        .make_sruw_cluster(nbcores)
      }, error = function(e) {
        warning("Parallel cluster creation failed; rerunning sequentially: ", e$message, call. = FALSE)
        NULL
      })

      if (!is.null(cl)) {
        on.exit({
          if (!is.null(cl)) parallel::stopCluster(cl)
        }, add = TRUE)
        parallel_attempt <- tryCatch(
          list(
            ok = TRUE,
            results = .run_windows_sruw_parallel(
              cl,
              input_list,
              environment(),
              backend_package
            )
          ),
          error = function(e) {
            list(ok = FALSE, error = conditionMessage(e))
          }
        )
        if (isTRUE(parallel_attempt$ok)) {
          results <- parallel_attempt$results
        } else {
          warning(
            "Parallel SRUW selection failed; rerunning sequentially: ",
            parallel_attempt$error,
            call. = FALSE
          )
          try(.stop_sruw_cluster(cl), silent = TRUE)
          cl <- NULL
          nbcores <- 1L
        }
      }
      if (is.null(cl)) nbcores <- 1L

    } else {
      results <- tryCatch({
        parallel::mclapply(input_list,
                         FUN = function(x) {
                            tryCatch({
                               wrapper.selectVar(
                                 x$nbcluster, x$ordervar, x$task_seed
                               )
                            }, error = function(e) {
                               list(error = paste("Parallel SRUW selection failed:", e$message), nbcluster=x$nbcluster)
                            })
                         },
                         mc.cores = nbcores,
                         mc.silent = TRUE,
                         mc.set.seed = FALSE
                         )
      }, error = function(e) {
         warning("Parallel SRUW selection failed; rerunning sequentially: ", e$message, call. = FALSE)
         list()
      })
       if(length(results) == 0) nbcores <- 1 
    }
  } else {
     nbcores <- 1
  }


  if (nbcores <= 1) {
    if (verbose) cat("Running sequentially...\n")
    results <- lapply(input_list, function(x) {
       tryCatch({
          wrapper.selectVar(x$nbcluster, x$ordervar, x$task_seed)
       }, error = function(e) {
          list(error = paste("Sequential SRUW selection failed:", e$message), nbcluster=x$nbcluster)
       })
    })
  }

  # Failed candidates remain invalid. Surviving fits retain backend, stopping,
  # and random-number provenance for downstream model selection.
  stop_sruw_fit(results, input_list)
  VariableSelectRes <- list()
  valid_results_count <- 0

  for (ll in seq_along(results)) {
    current_result <- results[[ll]]
    current_nbcluster_val <- if(!is.null(current_result$nbcluster_run)) current_result$nbcluster_run else input_list[[ll]]$nbcluster

    if (!is.null(current_result$error)) {
      if (verbose) cat("Run for", current_nbcluster_val, "clusters failed:", current_result$error, "\n")
      next
    }
    if (is.null(current_result$S) || !is.numeric(current_result$S)) {
       if (verbose) cat("Run for", current_nbcluster_val, "clusters produced invalid 'S' component. Skipping.\n")
       next
    }

    VariableSelectRes[[length(VariableSelectRes) + 1]] <- list(
      S = current_result$S,
      W = if (!is.null(current_result$W)) current_result$W else integer(0), 
      U = if (!is.null(current_result$U)) current_result$U else integer(0),
      criterionValue = if (!is.null(current_result$criterionValue)) current_result$criterionValue else NA_real_,
      criterion = if (!is.null(current_result$criterion)) current_result$criterion else criterion, 
      model = if (!is.null(current_result$model)) current_result$model else model_name,
      framework = if (!is.null(current_result$framework)) current_result$framework else framework,
      requestedModel = if (!is.null(current_result$requestedModel)) current_result$requestedModel else model_name,
      effectiveModel = if (!is.null(current_result$effectiveModel)) current_result$effectiveModel else current_result$model,
      nbcluster = current_nbcluster_val,
      parameters = current_result$parameters,
      partition = current_result$partition, 
      proba = current_result$proba,      
      missingValues = current_result$missingValues,
      stoppingRule = if (!is.null(current_result$stoppingRule)) current_result$stoppingRule else stopping,
      stoppingThreshold = if (!is.null(current_result$stoppingThreshold)) current_result$stoppingThreshold else hsize,
      nEvaluated = if (!is.null(current_result$nEvaluated)) current_result$nEvaluated else NA_integer_,
      stopReason = if (!is.null(current_result$stopReason)) current_result$stopReason else NA_character_,
      wStoppingRule = if (!is.null(current_result$wStoppingRule)) current_result$wStoppingRule else stopping,
      wStoppingThreshold = if (!is.null(current_result$wStoppingThreshold)) current_result$wStoppingThreshold else hsize,
      wNEvaluated = if (!is.null(current_result$wNEvaluated)) current_result$wNEvaluated else NA_integer_,
      wStopReason = if (!is.null(current_result$wStopReason)) current_result$wStopReason else NA_character_,
      rngSeed = current_result$rngSeed
    )
    valid_results_count <- valid_results_count + 1
  }

  if (valid_results_count == 0) {
    stop("All model fitting runs during variable selection failed. Please check data and parameters.", call. = FALSE)
  }

  if (length(VariableSelectRes) > 1) {
     cluster_order <- order(sapply(VariableSelectRes, `[[`, "nbcluster"))
     VariableSelectRes <- VariableSelectRes[cluster_order]
  }


  return(VariableSelectRes)
}





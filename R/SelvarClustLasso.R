SelvarClustLasso <- function(
  x,
  nbcluster,
  strategy = NULL,
  lambda                    = NULL,
  rho                       = NULL,
  num_vals_penalty = 5, 
  type                      = "lasso",
  hsize = 3,             
  criterion = "BIC",
  models = "VVI",
  rmodel = c("LI", "LB"), 
  imodel = c("LI", "LB"),       
  nbcores = min(2, detectCores(all.tests = FALSE, logical = FALSE)),
  impute_missing = TRUE,
  use_copula = TRUE,           
  scale_data = "always",
  scale_check_method = "pairwise.complete.obs",
  use_missing_pattern = FALSE,
  use_diag = TRUE,
  true_labels = NULL,
  sd_ratio_threshold    = 10,
  cond_number_threshold = 30,
  rank           = NULL,
  rank_control   = list(),
  mnarz_control  = list(),
  selection_control = list(),
  workflow = c(
    "auto", "decoupled_mnarz", "imputed_sruw", "joint_mnarz", "decoupled"
  ),
  joint_control = list(),
  init_control = list(),
  verbose = FALSE
) {
  workflow_was_missing <- missing(workflow)
  legacy_missingness_supplied <- !missing(use_missing_pattern)

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(true_labels)) {
    warning(
      paste0(
        "true_labels is ignored during fitting. Unsupervised partitions and ",
        "model choices must be computed without outcome information; calculate ",
        "ARI or other truth-based metrics after fitting."
      ),
      call. = FALSE
    )
  }

  # Validate the public fitting problem before resolving its workflow.
  CheckInputsC(x, nbcluster, lambda, rho, type, hsize, criterion, models, rmodel, imodel, nbcores)
  
  # Center and scale once; all fitted stages share this transformation.
  x <- as.matrix(x)
  n <- as.integer(nrow(x))
  p <- as.integer(ncol(x))
  nbcluster <- as.integer(nbcluster)
  workflow_choice <- .selvar_resolve_workflow(
    workflow = workflow,
    workflow_was_missing = workflow_was_missing,
    use_missing_pattern = use_missing_pattern,
    legacy_missingness_supplied = legacy_missingness_supplied,
    has_missing = anyNA(x)
  )
  workflow_requested <- workflow_choice$requested
  workflow <- workflow_choice$effective

  preprocessing <- .selvar_preprocess(
    x = x,
    scale_data = scale_data,
    sd_ratio_threshold = sd_ratio_threshold,
    cond_number_threshold = cond_number_threshold,
    use = scale_check_method
  )
  centers <- preprocessing$center
  sds <- preprocessing$scale
  x_scaled <- preprocessing$data
  
  # Variable ranking uses one completed matrix; the original mask enters the
  # subsequent MNARz likelihood.
  if (impute_missing && any(is.na(x_scaled))) {
    if (use_copula) {
      .selvar_require_namespace("gcimputeR", "copula ranking imputation")
      x_imp_scaled <- as.matrix(gcimputeR::impute_GC(
        as.data.frame(x_scaled), verbose = FALSE
      )$Ximp)
    } else {
      .selvar_require_namespace("missRanger", "random-forest ranking imputation")
      x_imp_scaled <- as.matrix(missRanger::missRanger(
        as.data.frame(x_scaled), verbose = 0
      ))
    }
  } else if (any(is.na(x_scaled))) {
    stop(
      paste0(
        "SRUW ranking and role selection require one complete working data set; ",
        "set impute_missing=TRUE for incomplete x."
      ),
      call. = FALSE
    )
  } else {
    x_imp_scaled <- x_scaled
  }

  x_imp_orig <- sweep(x_imp_scaled, 2, sds, "*")
  x_imp_orig <- sweep(x_imp_orig, 2, centers, "+")
  dimnames(x_imp_orig) <- dimnames(x)


  .rank_defaults <- list(
    type = type,
    lambda = lambda,
    rho = rho,
    group_shrinkage_method = "weighted_by_dist_to_diag_W0",
    distance_method        = "Euclidean",
    lambda_omega_0         = 50,
    epsilon_weighted_by_W0 = sqrt(.Machine$double.eps),
    penalize_diag          = FALSE,
    laplacian_target_type  = "diag_Omega_hat",
    adj_threshold          = 1e-4,
    laplacian_norm_type    = "symmetric",
    initialize             = "hc",
    nbcores                = nbcores,
    n.start                = 250,
    warm_start             = "none",
    min_scorable_fraction  = 0.5,
    verbose                = verbose
  )
  if (missing(rank_control) || is.null(rank_control)) rank_control <- list()               

  rank_control <- utils::modifyList(.rank_defaults, rank_control)

  rank_is_fitted <- missing(rank) || is.null(rank)
  # Ranking controls are irrelevant when the caller supplies the rank.  Do not
  # validate or conflict-check an unused grid in that route.
  supplied_penalty_grids <- rank_is_fitted &&
    !is.null(rank_control$penalty_grids)
  if (supplied_penalty_grids && length(init_control)) {
    stop(
      paste0(
        "Supply either rank_control$penalty_grids with pinned initializers ",
        "or init_control, not both."
      ),
      call. = FALSE
    )
  }
  registry_control <- init_control
  if (supplied_penalty_grids) {
    if (!is.list(rank_control$penalty_grids) ||
        length(rank_control$penalty_grids) != length(nbcluster)) {
      stop(
        "rank_control$penalty_grids must align one-to-one with nbcluster.",
        call. = FALSE
      )
    }
    per_k <- lapply(seq_along(nbcluster), function(index) {
      entry <- rank_control$penalty_grids[[index]]
      if (!is.list(entry) || is.null(entry$initializer)) {
        stop(
          sprintf(
            "rank_control$penalty_grids[[%d]] lacks a pinned initializer.",
            index
          ),
          call. = FALSE
        )
      }
      list(method = "user", user_init = entry$initializer)
    })
    names(per_k) <- paste0("K", nbcluster)
    registry_control <- list(method = "user", per_k = per_k)
  }

  complete_initializer_data <- !anyNA(x_imp_scaled)
  # A supplied rank has no ranking initializer to reuse. The separate MNAR fit
  # therefore uses its own initialization unless init_control supplies a
  # reusable Gaussian state.
  needs_initializer <- rank_is_fitted ||
    workflow %in% c("decoupled_mnarz", "joint_mnarz") ||
    length(init_control) > 0L
  initializer_registry <- initializer_states <- NULL
  if (needs_initializer) {
    if (!complete_initializer_data) {
      stop(
        paste0(
          "A shared initializer requires complete initialization data. ",
          "Set impute_missing=TRUE or omit init_control when a supplied rank ",
          "is followed only by the standalone missing-data fit."
        ),
        call. = FALSE
      )
    }
    initializer_registry <- .build_initializer_registry(
      data = x_imp_scaled,
      nbcluster = nbcluster,
      init_control = registry_control,
      legacy_method = rank_control$initialize,
      legacy_n_start = rank_control$n.start,
      legacy_lambda_omega_0 = rank_control$lambda_omega_0
    )
    initializer_states <- lapply(
      initializer_registry, .initializer_state_projection
    )
    names(initializer_states) <- names(initializer_registry)
  }

  .mnarz_defaults <- list(
    mecha     = "MNARz",
    diag      = use_diag,
    rmax      = 100,
    tol       = 1e-4,
    init      = NULL
  )
  mnarz_control <- utils::modifyList(.mnarz_defaults, mnarz_control)

  .selection_defaults <- list(stopping = "consecutive", seed = 1L)
  if (is.null(selection_control)) selection_control <- list()
  if (!is.list(selection_control) ||
      (length(selection_control) && is.null(names(selection_control)))) {
    stop("selection_control must be a named list.", call. = FALSE)
  }
  unknown_selection_controls <- setdiff(
    names(selection_control), names(.selection_defaults)
  )
  if (length(unknown_selection_controls)) {
    stop(
      "Unknown selection_control entries: ",
      paste(unknown_selection_controls, collapse = ", "),
      call. = FALSE
    )
  }
  selection_control <- utils::modifyList(
    .selection_defaults,
    selection_control
  )
  selection_control$stopping <- match.arg(
    selection_control$stopping,
    c("consecutive", "legacy_block")
  )
  selection_control$seed <- .numeric_init_seed(
    selection_control$seed, "selection_control$seed"
  )

  .joint_defaults <- list(
    diag = use_diag,
    rmax = 200L,
    tol = 1e-4,
    r_search_cap = 256L
  )
  if (is.null(joint_control)) joint_control <- list()
  if (!is.list(joint_control)) {
    stop("joint_control must be a list.", call. = FALSE)
  }
  unknown_joint_controls <- setdiff(names(joint_control), names(.joint_defaults))
  if (length(unknown_joint_controls)) {
    stop(
      "Unknown joint_control entries: ",
      paste(unknown_joint_controls, collapse = ", "),
      call. = FALSE
    )
  }
  joint_control <- utils::modifyList(.joint_defaults, joint_control)

  if (!is.null(lambda)) rank_control$lambda <- lambda
  if (!is.null(rho))    rank_control$rho    <- rho
  rank_control$type <- type                         

  rank_control$penalize_diag <- isTRUE(rank_control$penalize_diag)

  if (rank_is_fitted) {
    .selvar_require_namespace("glassoFast", "penalized variable ranking")
  }
  if (rank_is_fitted && !supplied_penalty_grids) {
    if (rank_control$group_shrinkage_method %in%
        c("weighted_by_dist_to_I", "weighted_by_dist_to_diag_W0") &&
        !identical(rank_control$distance_method, "Euclidean")) {
      .selvar_require_namespace(
        "shapes", "distance-weighted precision penalties"
      )
    }
    if (identical(rank_control$group_shrinkage_method, "laplacian_spectral")) {
      .selvar_require_namespace("igraph", "spectral precision penalties")
    }
    grid_list <- compute_grids_per_K(
        X               = as.matrix(x_imp_scaled),
        nbcluster       = nbcluster,
        P_method        = rank_control$group_shrinkage_method,
        distance_method = rank_control$distance_method,
        eps_w0          = rank_control$epsilon_weighted_by_W0,
        L               = num_vals_penalty,
        frac_min        = 0.05,
        init_method     = rank_control$initialize,
        n.start         = rank_control$n.start,
        lambda_omega_0  = rank_control$lambda_omega_0,
        initializers    = initializer_states,
        penalize_diag   = rank_control$penalize_diag,
        laplacian_target_type = rank_control$laplacian_target_type,
        adj_threshold   = rank_control$adj_threshold,
        laplacian_norm_type = rank_control$laplacian_norm_type)

    # Preserve the K-specific paths and the exact initialization used to
    # construct each path. A user-supplied dimension remains common across K.
    if (!is.null(rank_control$lambda) && length(rank_control$lambda)) {
      grid_list <- lapply(grid_list, function(grid) {
        grid$lambda_mu <- as.numeric(rank_control$lambda)
        grid
      })
    }
    if (!is.null(rank_control$rho) && length(rank_control$rho)) {
      grid_list <- lapply(grid_list, function(grid) {
        grid$rho <- as.numeric(rank_control$rho)
        grid
      })
    }
    attr(grid_list, "initializer_adapters") <- initializer_registry
    rank_control$penalty_grids <- grid_list
  } else if (rank_is_fitted) {
    attr(rank_control$penalty_grids, "initializer_adapters") <-
      initializer_registry
  }
  
  if (!is.null(strategy)) {
    warning(
      paste0(
        "strategy is deprecated and is not used by the SRUW backend ",
        "interface. Configure initialization through init_control instead."
      ),
      call. = FALSE
    )
  }
  
  OrderVariable <- matrix(NA, nrow = length(nbcluster), ncol = p)

  # Rank variables by component-mean activity over the penalized path.
  if (rank_is_fitted) {
    rank_ctrl <- rank_control
    rank_ctrl$x         <- x_imp_scaled      
    rank_ctrl$nbcluster <- nbcluster   
    if (verbose) cat("Performing variable ranking\n")
    OrderVariable <- do.call(SortvarClust, rank_ctrl)
  } else {
    for (r in seq_len(nrow(OrderVariable))) {
      OrderVariable[r, ] <- rank
    }
  }
  if (verbose) {
    cat("Variable Ranks: \n")
    print(OrderVariable)
  }

  if (identical(workflow, "joint_mnarz")) {
    if (length(criterion) != 1L || !identical(as.character(criterion), "BIC")) {
      stop(
        "workflow='joint_mnarz' accepts the total observed-data BIC only.",
        call. = FALSE
      )
    }
    if (!anyNA(x_scaled)) {
      stop(
        paste0(
          "workflow='joint_mnarz' requires a non-degenerate missingness mask; ",
          "all-observed data have no finite class-only missingness intercept."
        ),
        call. = FALSE
      )
    }
    if (anyNA(x_imp_scaled)) {
      stop(
        paste0(
          "Joint-workflow ranking requires one fixed complete imputation. ",
          "Set impute_missing=TRUE or supply complete ranking data upstream."
        ),
        call. = FALSE
      )
    }

    joint_fits <- lapply(seq_along(nbcluster), function(index) {
      select_joint_sruw_mnarz(
        x = x_scaled,
        K = nbcluster[index],
        rank = as.integer(OrderVariable[index, ]),
        hsize = hsize,
        diag = isTRUE(joint_control$diag),
        rmax = as.integer(joint_control$rmax),
        tol = as.double(joint_control$tol),
        r_search_cap = as.integer(joint_control$r_search_cap),
        initial_labels = initializer_registry[[index]]$partition
      )
    })
    joint_bic <- vapply(joint_fits, `[[`, numeric(1), "bic_min")
    if (length(joint_bic) != length(nbcluster) || any(!is.finite(joint_bic))) {
      stop(
        "At least one joint K candidate was unscorable; partial K sets are not selected.",
        call. = FALSE
      )
    }
    best_index <- which.min(joint_bic)
    best_joint <- joint_fits[[best_index]]
    fixed_fit <- best_joint$final_fit
    joint_imputed_original <- sweep(fixed_fit$imputedData, 2, sds, "*")
    joint_imputed_original <- sweep(
      joint_imputed_original, 2, centers, "+"
    )
    dimnames(joint_imputed_original) <- dimnames(x)
    joint_regression <- NULL
    if (length(best_joint$U)) {
      joint_regression <- rbind(
        intercept = fixed_fit$parameters$a,
        fixed_fit$parameters$beta
      )
    }
    joint_model <- list(
      S = best_joint$S,
      R = best_joint$R,
      U = best_joint$U,
      W = best_joint$W,
      criterionValue = best_joint$bic_min,
      criterion = "BIC",
      model = "joint_mnarz",
      rmodel = "joint_common_gaussian",
      imodel = "joint_common_gaussian",
      parameters = fixed_fit$parameters,
      nbcluster = best_joint$K,
      partition = fixed_fit$partition,
      proba = fixed_fit$proba,
      regparameters = joint_regression,
      imputedData = joint_imputed_original,
      clust_result = fixed_fit,
      workflow = "joint_mnarz",
      workflowRequested = workflow_requested,
      workflowEffective = workflow,
      criterionConvention = "minimize",
      criterionScope = fixed_fit$diagnostics$criterion_scope,
      selectionCriterionValue = best_joint$bic_min,
      selectionPartition = fixed_fit$partition,
      selectionResult = best_joint,
      ranking = OrderVariable,
      stoppingRule = "consecutive",
      stoppingThreshold = as.integer(hsize),
      nEvaluated = as.integer(length(best_joint$selection_trace$forward)),
      stopReason = if (isTRUE(
        best_joint$diagnostics$forward_stopped_by_rolling_c
      )) "consecutive_nonpositive" else "order_exhausted",
      wStoppingRule = "consecutive",
      wStoppingThreshold = as.integer(hsize),
      wNEvaluated = as.integer(length(best_joint$selection_trace$reverse)),
      wStopReason = if (isTRUE(
        best_joint$diagnostics$reverse_stopped_by_rolling_c
      )) "consecutive_nonpositive" else "order_exhausted",
      jointFit = best_joint,
      jointFitsByK = joint_fits,
      paperAlignment = best_joint$paper_alignment,
      preprocessing = preprocessing[setdiff(names(preprocessing), "data")],
      initialization = list(
        selected = initializer_registry[[best_index]],
        registry = initializer_registry,
        reuse = paste0(
          "the selected adapter partition initializes every joint role ",
          "candidate for its K"
        )
      )
    )
    return(ProcessModelOutput(joint_model))
  }
  
  # The supported SRUW route is unsupervised.
  supervised <- FALSE 
  knownlabels <- as.integer(1:n)
  
  # Stage B assigns SRUW roles and selects the component count.
  bestModel <- list()
  if (length(criterion) == 1) {
    if (verbose) cat("Performing variable selection with", criterion, "criterion\n")
    VariableSelectRes <- VariableSelection(
      x_imp_scaled, nbcluster, models, criterion, OrderVariable, 
      hsize, supervised, knownlabels, nbcores,
      stopping = selection_control$stopping,
      rng_seed = selection_control$seed,
      verbose = verbose
    )
    bestModel[[criterion]] <- ModelSelectionClust(
      VariableSelectRes, x_imp_scaled, rmodel, imodel, nbcores,
      rng_seed = selection_control$seed
    )
  } else {
    for (crit in criterion) {
      if (verbose) cat("Variable selection with", crit, "criterion\n")
      VariableSelectRes <- VariableSelection(
        x_imp_scaled, nbcluster, models, crit, OrderVariable, 
        hsize, supervised, knownlabels, nbcores,
        stopping = selection_control$stopping,
        rng_seed = selection_control$seed,
        verbose = verbose
      )
      
      bestModel[[crit]] <- ModelSelectionClust(
        VariableSelectRes, x_imp_scaled, rmodel, imodel, nbcores,
        rng_seed = selection_control$seed
      )
    }
  }
  
  # The decoupled result combines SRUW roles selected on the completed matrix
  # with partition, posterior probabilities, parameters, imputation, and BIC
  # from a separate role-agnostic MNARz fit.
  for (i in seq_along(bestModel)) {
    model_name <- names(bestModel)[i]
    finalModel <- bestModel[[i]]
    number_clusters <- finalModel$nbcluster
    initializer_index <- match(as.integer(number_clusters), nbcluster)
    selected_initializer <- if (!is.null(initializer_registry) &&
        !is.na(initializer_index)) {
      initializer_registry[[initializer_index]]
    } else {
      NULL
    }
    finalModel$initialization <- if (is.null(initializer_registry)) {
      NULL
    } else {
      list(
        selected = selected_initializer,
        registry = initializer_registry,
        reuse = paste0(
          "the selected adapter state initializes ranking and, when requested, ",
          "the separate MNAR fit"
        )
      )
    }
    finalModel$imputedData <- x_imp_orig
    finalModel$workflow <- "imputed_sruw"
    finalModel$workflowRequested <- workflow_requested
    finalModel$workflowEffective <- workflow
    finalModel$preprocessing <- preprocessing[setdiff(
      names(preprocessing), "data"
    )]
    # Preserve the ranking matrix together with its per-grid convergence,
    # failure, objective-trace, initialization, and warm-path attributes.
    finalModel$ranking <- OrderVariable
    finalModel$criterionConvention <- "maximize"
    finalModel$criterionScope <- paste0(
      "Additive SRUW criterion fitted to the completed working data on its ",
      "preprocessing scale; distinct from the joint observed likelihood of ",
      "(Y_obs, C)."
    )
    if (identical(workflow, "decoupled_mnarz")) {
      if (verbose) cat("Fitting final model for criterion", model_name, "\n")
      if (!exists("EMClustMNARz")) stop("EMClustMNARz function is missing")
      selection_result <- finalModel
      em_call <- c(list(x = x_scaled, K = number_clusters,
                        criterion = model_name), mnarz_control)
      shared_initializer_used <- is.null(em_call$init) &&
        !is.null(selected_initializer)
      if (is.null(em_call$init) && !is.null(selected_initializer)) {
        em_call$init <- selected_initializer
      }

      clust_result <- do.call(EMClustMNARz, em_call)
      
      # The completed matrix must preserve the observation dimension.
      if (is.null(clust_result$imputedData) || !is.matrix(clust_result$imputedData) || 
          !all(dim(clust_result$imputedData) == dim(x))) {
        stop("EMClustMNARz did not return valid imputedData")
          }

      x_imputed_final <- sweep(clust_result$imputedData, 2, sds, "*")
      x_imputed_final <- sweep(x_imputed_final,        2, centers, "+")
      dimnames(x_imputed_final) <- dimnames(x)
      finalModel$imputedData <- x_imputed_final
            
      # Hard assignments must cover the original observations.
      if (is.null(clust_result$partition) || length(clust_result$partition) != nrow(x) ||
          !all(clust_result$partition %in% 1:number_clusters)) {
        stop("EMClustMNARz returned invalid partition")
      }
      if (!isTRUE(clust_result$diagnostics$criterion_available) ||
          !is.finite(clust_result$criterionValue[[model_name]])) {
        stop(
          paste0(
            "The separate MNARz fit did not produce an available ",
            model_name, " criterion; the decoupled result is unavailable."
          ),
          call. = FALSE
        )
      }

      finalModel$workflow <- "decoupled_mnarz"
      finalModel$selectionResult <- selection_result
      finalModel$selectionCriterionValue <- selection_result$criterionValue
      finalModel$selectionPartition <- selection_result$partition
      finalModel$partition <- clust_result$partition
      finalModel$proba <- clust_result$proba
      finalModel$parameters <- clust_result$parameters
      finalModel$criterionValue <- clust_result$criterionValue[[model_name]]
      finalModel$criterionConvention <- "minimize"
      finalModel$criterionScope <- clust_result$diagnostics$criterion_scope
      finalModel$clust_result <- clust_result
      finalModel$initialization$mnar_init_source <- if (shared_initializer_used) {
        "selected shared adapter state"
      } else {
        "explicit mnarz_control$init override"
      }
    }
    bestModel[[i]] <- finalModel
  }
  # Only successfully fitted criterion candidates enter the public result.
  bestModel <- bestModel[!sapply(bestModel, is.null)]
 
  output <- PrepareOutput(bestModel)
  
  return(output)
}

# Convert criterion-specific fits to one result or a named collection.
PrepareOutput <- function(bestModel) {
  output <- list()
  for (name in names(bestModel)) {
    processed <- ProcessModelOutput(bestModel[[name]])
    if (!is.null(processed)) output[[name]] <- processed
  }
  if (length(output) == 1) {
    return(output[[1]])
  }
  if (!length(output)) return(output)
  .new_selvarmix_collection(output)
}

# Normalize a backend fit to the stable `selvarmix` schema.
ProcessModelOutput <- function(modelResult) {
  if (is.null(modelResult)) return(NULL)
  if (is.null(modelResult$imputedData)) {
    stop(
      "A successful model result must contain a completed imputedData matrix.",
      call. = FALSE
    )
  }
  if (!is.null(modelResult$regparameters)) {
    if (!is.null(modelResult$U) && length(modelResult$U) != 0)
      colnames(modelResult$regparameters) <- modelResult$U
    if (!is.null(modelResult$R) && length(modelResult$R) != 0)
      rownames(modelResult$regparameters) <- c("intercept", modelResult$R)
  }
  
  object <- list(
    S = modelResult$S,
    R = modelResult$R,
    U = modelResult$U,
    W = modelResult$W,
    criterionValue = modelResult$criterionValue,
    criterion = modelResult$criterion,
    model = modelResult$model,
    framework = modelResult$framework,
    requestedModel = modelResult$requestedModel,
    effectiveModel = modelResult$effectiveModel,
    rmodel = modelResult$rmodel,
    imodel = modelResult$imodel,
    parameters = modelResult$parameters,
    nbcluster = modelResult$nbcluster,
    partition = modelResult$partition,
    proba = modelResult$proba,
    regparameters = modelResult$regparameters,
    imputedData = modelResult$imputedData,
    parametersMNARz = modelResult$clust_result,
    stoppingRule = modelResult$stoppingRule,
    stoppingThreshold = modelResult$stoppingThreshold,
    nEvaluated = modelResult$nEvaluated,
    stopReason = modelResult$stopReason,
    wStoppingRule = modelResult$wStoppingRule,
    wStoppingThreshold = modelResult$wStoppingThreshold,
    wNEvaluated = modelResult$wNEvaluated,
    wStopReason = modelResult$wStopReason,
    workflow = modelResult$workflow,
    workflowRequested = modelResult$workflowRequested,
    workflowEffective = modelResult$workflowEffective,
    criterionConvention = modelResult$criterionConvention,
    criterionScope = modelResult$criterionScope,
    selectionCriterionValue = modelResult$selectionCriterionValue,
    selectionPartition = modelResult$selectionPartition,
    selectionResult = modelResult$selectionResult,
    ranking = modelResult$ranking,
    emptyRRefit = modelResult$emptyRRefit,
    rngSeed = modelResult$rngSeed,
    jointFit = modelResult$jointFit,
    jointFitsByK = modelResult$jointFitsByK,
    paperAlignment = modelResult$paperAlignment,
    initialization = modelResult$initialization,
    preprocessing = modelResult$preprocessing,
    schemaVersion = .selvarmix_schema_version,
    diagnostics = .selvarmix_result_diagnostics(modelResult)
  )
  .new_selvarmix_result(object)
}

is_rmixmod_model <- function(models) {
  return(grepl("mixmodGaussianModel", models))
}

is_mclust_model <- function(models){
  mclust_models <-c(
                    # Spherical models
                    "EII", "VII", 
                    
                    # Diagonal models
                    "EEI", "VEI", "EVI", "VVI", 
                    
                    # Ellipsoidal models
                    "EEE", "VEE", "EVE", "VVE", 
                    "EEV", "VEV", "EVV", "VVV"
                    )
  return(models %in% mclust_models)
}

.selvar_require_namespace <- function(package, role) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Package '", package, "' is required for ", role,
      ". Install it before using this route.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.selvar_resolve_workflow <- function(workflow, workflow_was_missing,
                                     use_missing_pattern,
                                     legacy_missingness_supplied,
                                     has_missing) {
  choices <- c(
    "auto", "decoupled_mnarz", "imputed_sruw", "joint_mnarz", "decoupled"
  )
  if (!is.logical(workflow_was_missing) || length(workflow_was_missing) != 1L ||
      is.na(workflow_was_missing) ||
      !is.logical(legacy_missingness_supplied) ||
      length(legacy_missingness_supplied) != 1L ||
      is.na(legacy_missingness_supplied)) {
    stop("Internal workflow-resolution flags are invalid.", call. = FALSE)
  }
  if (!is.logical(has_missing) || length(has_missing) != 1L ||
      is.na(has_missing)) {
    stop("has_missing must be one non-missing logical value.", call. = FALSE)
  }

  requested <- if (isTRUE(workflow_was_missing) &&
      isTRUE(legacy_missingness_supplied)) {
    "decoupled"
  } else {
    match.arg(workflow, choices)
  }
  if (!identical(requested, "decoupled") &&
      isTRUE(legacy_missingness_supplied)) {
    stop(
      paste0(
        "use_missing_pattern is a deprecated compatibility switch and cannot ",
        "be combined with an explicit modern workflow."
      ),
      call. = FALSE
    )
  }

  if (identical(requested, "decoupled")) {
    if (!is.logical(use_missing_pattern) || length(use_missing_pattern) != 1L ||
        is.na(use_missing_pattern)) {
      stop("use_missing_pattern must be TRUE or FALSE.", call. = FALSE)
    }
    warning(
      paste0(
        "workflow='decoupled' and use_missing_pattern are deprecated; use ",
        "workflow='decoupled_mnarz' or workflow='imputed_sruw'."
      ),
      call. = FALSE
    )
    effective <- if (isTRUE(use_missing_pattern)) {
      "decoupled_mnarz"
    } else {
      "imputed_sruw"
    }
  } else if (identical(requested, "auto")) {
    effective <- if (has_missing) "decoupled_mnarz" else "imputed_sruw"
  } else {
    effective <- requested
  }

  if (effective %in% c("decoupled_mnarz", "joint_mnarz") && !has_missing &&
      !identical(requested, "decoupled")) {
    stop(
      paste0(
        "workflow='", effective, "' requires a non-degenerate missingness ",
        "mask; class-only MNARz has no finite intercept on all-observed data."
      ),
      call. = FALSE
    )
  }
  list(requested = requested, effective = effective)
}

.selvar_preprocess <- function(x, scale_data = "always",
                               sd_ratio_threshold = 10,
                               cond_number_threshold = 30,
                               use = "pairwise.complete.obs") {
  x <- as.matrix(x)
  mode <- if (is.logical(scale_data) && length(scale_data) == 1L &&
      !is.na(scale_data)) {
    if (scale_data) "always" else "never"
  } else if (is.character(scale_data) && length(scale_data) == 1L &&
             !is.na(scale_data)) {
    match.arg(scale_data, c("always", "auto", "never"))
  } else {
    stop(
      "scale_data must be TRUE/FALSE or one of 'always', 'auto', and 'never'.",
      call. = FALSE
    )
  }

  observed_count <- colSums(!is.na(x))
  all_missing <- which(observed_count == 0L)
  if (length(all_missing)) {
    stop(
      "All-missing column(s) are unidentified and must be removed: ",
      paste(all_missing, collapse = ", "),
      call. = FALSE
    )
  }
  observed_sd <- apply(x, 2L, stats::sd, na.rm = TRUE)
  observed_magnitude <- vapply(seq_len(ncol(x)), function(j) {
    max(1, abs(x[!is.na(x[, j]), j]))
  }, numeric(1))
  degenerate <- observed_count < 2L | !is.finite(observed_sd) |
    observed_sd <= .Machine$double.eps * observed_magnitude
  if (any(degenerate)) {
    stop(
      paste0(
        "Non-varying or insufficiently observed column(s) cannot identify a ",
        "Gaussian covariance: ", paste(which(degenerate), collapse = ", "),
        ". Filter them before fitting."
      ),
      call. = FALSE
    )
  }

  auto_trigger <- check_scale_data(
    x, sd_ratio_threshold, cond_number_threshold, use = use
  )
  scaled <- identical(mode, "always") ||
    (identical(mode, "auto") && isTRUE(auto_trigger))
  centered <- !identical(mode, "never")
  if (identical(mode, "never") && isTRUE(auto_trigger)) {
    warning(
      paste0(
        "scale_data='never' was requested although the scale diagnostic ",
        "detected materially unequal scales or ill conditioning."
      ),
      call. = FALSE
    )
  }
  center <- if (centered) colMeans(x, na.rm = TRUE) else rep(0, ncol(x))
  scale <- if (scaled) observed_sd else rep(1, ncol(x))
  transformed <- sweep(x, 2L, center, FUN = "-")
  transformed <- sweep(transformed, 2L, scale, FUN = "/")
  dimnames(transformed) <- dimnames(x)

  list(
    data = transformed,
    mode = mode,
    centered = centered,
    scaled = scaled,
    auto_triggered = isTRUE(auto_trigger),
    center = as.numeric(center),
    scale = as.numeric(scale),
    variable_names = colnames(x),
    parameter_scale = if (scaled) "standardized" else if (centered) {
      "centered_original_units"
    } else {
      "original_units"
    },
    inverse_transform = "x_original = x_working * scale + center",
    paper_alignment = if (centered && scaled) {
      "observed-mean centering and observed-SD scaling"
    } else if (centered) {
      "observed-mean centering without scaling"
    } else {
      "no centering or scaling"
    }
  )
}

check_scale_data <- function(x,
                             sd_ratio_threshold = 10,
                             cond_number_threshold = 30,
                             use = c("pairwise.complete.obs", "median")) {
  x <- as.matrix(x)
  sds <- apply(x, 2, stats::sd, na.rm = TRUE)
  if (!length(sds) || any(!is.finite(sds)) || any(sds <= 0)) return(TRUE)
  ratio_sd <- max(sds) / min(sds)

  use <- match.arg(use)
  if (use == "pairwise.complete.obs") {
    cov_mat <- suppressWarnings(stats::cov(x, use = "pairwise.complete.obs"))
  } else {                         
    med <- matrixStats::colMedians(x, na.rm = TRUE)
    x_med <- x
    na <- is.na(x_med)
    if (any(na)) x_med[na] <- rep(med, each = nrow(x))[na]
    cov_mat <- stats::cov(x_med)
  }

  if (!is.matrix(cov_mat) || any(!is.finite(cov_mat)) ||
      min(dim(cov_mat)) == 0)
    return(TRUE)   

  singular_values <- tryCatch(
    svd(cov_mat, nu = 0L, nv = 0L)$d,
    error = function(e) numeric()
  )
  if (!length(singular_values) || any(!is.finite(singular_values)) ||
      max(singular_values) <= 0) return(TRUE)
  positive <- singular_values[
    singular_values > .Machine$double.eps * max(singular_values)
  ]
  cond_num <- if (length(positive) < ncol(cov_mat)) {
    Inf
  } else {
    max(positive) / min(positive)
  }

  do_scale <- (ratio_sd > sd_ratio_threshold) ||
              (cond_num  > cond_number_threshold)
  return(do_scale)
}

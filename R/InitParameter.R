InitParameter <- function(
    data,
    nbClust,
    init = c("kmeans", "hc"),
    n.start = 250,
    lambda_omega_0 = 50,
    epsilon_pd = sqrt(.Machine$double.eps),
    seed = NULL,
    user_init = NULL,
    previous_fit = NULL,
    multistart_methods = c("mclust_hc", "mclust_em", "kmeans", "random"),
    multistart_replicates = 1L,
    multistart_seeds = NULL,
    mclust_model = "VVV",
    return_adapter = FALSE) {
  # Normalize the default `hc` alias and all named methods before constructing
  # the numeric state.
  if (missing(init)) init <- "kmeans"
  if (!is.logical(return_adapter) || length(return_adapter) != 1L ||
      is.na(return_adapter)) {
    stop("return_adapter must be TRUE or FALSE.", call. = FALSE)
  }

  adapter <- .numeric_initialization_adapter(
    data = data,
    nbClust = nbClust,
    init = init,
    n.start = n.start,
    lambda_omega_0 = lambda_omega_0,
    epsilon_pd = epsilon_pd,
    seed = seed,
    user_init = user_init,
    previous_fit = previous_fit,
    multistart_methods = multistart_methods,
    multistart_replicates = multistart_replicates,
    multistart_seeds = multistart_seeds,
    mclust_model = mclust_model
  )
  if (isTRUE(return_adapter)) return(adapter)

  # Return the five numeric fields required by ranking; attach metadata as
  # attributes to preserve the state representation.
  state <- .initializer_state_projection(adapter)
  report_names <- setdiff(names(adapter), .numeric_initialization_state_names)
  report <- unclass(adapter[report_names])
  attr(state, "initialization_report") <- report
  attr(state, "initialization_schema_version") <-
    .numeric_initialization_schema_version
  state
}

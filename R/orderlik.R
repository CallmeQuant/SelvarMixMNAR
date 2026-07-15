orderlikC <- function(x, g, nbcores) {
  stop(
    paste0(
      "type = 'likelihood' has no fitted variable-ordering estimator. Use ",
      "type = 'lasso'."
    ),
    call. = FALSE
  )
}

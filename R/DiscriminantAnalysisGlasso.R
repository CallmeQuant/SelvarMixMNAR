.supervised_route_abort <- function() {
  condition <- structure(
    list(
      message = paste0(
        "The supervised discriminant-analysis route is unavailable in ",
        "SelvarMixMNAR 0.1.x. This package fits unsupervised clustering ",
        "models through SelvarClustLasso(); the former supervised entry ",
        "points are retained for call compatibility."
      ),
      call = NULL
    ),
    class = c(
      "selvarmix_supervised_deprecated",
      "selvarmix_supervised_lasso_unavailable",
      "selvarmix_feature_unavailable",
      "error",
      "condition"
    )
  )
  stop(condition)
}

# The former helper name returns the same condition as the exported
# supervised entry points.
.supervised_lasso_abort <- function() {
  .supervised_route_abort()
}

DiscriminantAnalysisGlasso <- function(data,
                                       nbCluster,
                                       lambda,
                                       rho,
                                       knownlabels,
                                       nbCores) {
  .supervised_route_abort()
}

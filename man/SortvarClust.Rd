\name{SortvarClust}
\alias{SortvarClust}
\title{Penalized variable ranking for model-based clustering}
\description{
Ranks variables from fitted penalized Gaussian-mixture paths. Penalty grids
and initializers may be specific to each candidate cluster count. The
warm-start mode determines whether state is reused within EM iterations,
between adjacent penalty points, or both.
}
\usage{
SortvarClust(
  x, nbcluster, type = "lasso",
  lambda = seq(20, 100, by = 10),
  rho = seq(0.1, 1, length = 5),
  group_shrinkage_method = "common", distance_method = "Euclidean",
  lambda_omega_0 = 50,
  epsilon_weighted_by_W0 = sqrt(.Machine$double.eps),
  laplacian_target_type = c("identity", "diag_Omega_hat"),
  adj_threshold = 1e-4,
  laplacian_norm_type = c("symmetric", "unsymmetric"),
  penalize_diag = FALSE, initialize = "hc",
  nbcores = min(2, parallel::detectCores(logical = FALSE)),
  n.start = 250, scale = FALSE, penalty_grids = NULL,
  warm_start = c("none", "inner", "outer", "both"),
  verbose = FALSE, min_scorable_fraction = 0.5,
  ..., init_control = list()
)
}
\arguments{
  \item{x}{Numeric matrix or data frame with observations in rows.  Missing or non-finite values are not accepted.}
  \item{nbcluster}{Positive integer vector of candidate cluster counts.}
  \item{type}{Ranking method. \code{"lasso"} fits the penalized
  Gaussian-mixture path. \code{"likelihood"} is unavailable and raises a
  classed error.}
  \item{lambda}{Non-negative mean-penalty values used when \code{penalty_grids} is \code{NULL}.}
  \item{rho}{Non-negative precision-penalty values used when \code{penalty_grids} is \code{NULL}.}
  \item{group_shrinkage_method}{Construction of precision-penalty weights.}
  \item{distance_method}{Covariance-distance method used by distance-weighted penalties.}
  \item{lambda_omega_0}{Dimensionless initial regularization multiplier. For a
  singular component covariance it is multiplied by that covariance's scale,
  so initialization is equivariant under a global change of measurement units.}
  \item{epsilon_weighted_by_W0}{Positive numerical denominator safeguard for adaptive weights.}
  \item{laplacian_target_type}{Target used by the spectral-Laplacian weighting option.}
  \item{adj_threshold}{Adjacency threshold for spectral-Laplacian weights.}
  \item{laplacian_norm_type}{Normalized-Laplacian convention.}
  \item{penalize_diag}{Logical request for diagonal precision penalization. The native ranking kernel rejects nonzero diagonal weights because the fitted objective penalizes off-diagonal entries only.}
  \item{initialize}{Initialization method. Supported values are \code{"mclust_hc"}, \code{"mclust_em"}, \code{"kmeans"}, \code{"random"}, \code{"user"}, \code{"previous_fit"}, and \code{"deterministic_multistart"}; \code{"hc"} is a compatibility alias for \code{"mclust_hc"}.}
  \item{init_control}{Strict named controls for the numeric initializer adapter, including seed, user/previous state, deterministic multistart, and per-K overrides. It cannot be combined with externally supplied \code{penalty_grids}.}
  \item{nbcores}{Number of workers.  Dependent grid points within an outer-warm path are sequential; independent cluster-count paths may be parallel.}
  \item{n.start}{Number of starts used by initialization helpers.}
  \item{scale}{Whether to center and standardize the complete input before
  direct ranking. The main \code{SelvarClustLasso} interface centers and
  scales once before supplying the transformed ranking data.}
  \item{penalty_grids}{Optional list aligned one-to-one with \code{nbcluster}; each entry supplies a K-specific grid and its pinned initialization.}
  \item{warm_start}{\code{"none"} for independent cold fits; \code{"inner"}
  for \pkg{glassoFast} continuation between EM iterations; \code{"outer"} for
  complete-state continuation along the penalty path; or \code{"both"}.}
  \item{verbose}{Logical; emit progress messages when true.}
  \item{min_scorable_fraction}{Minimum fraction of a K-specific grid that must
  converge with a valid numerical state. Failed points retain all-\code{NA}
  activity rows and are omitted from the ranking score.}
  \item{...}{Additional compatibility arguments.}
}
\details{
Outer continuation follows a canonical descending-rho backbone followed by
descending-lambda spokes. The transferred state contains proportions, means,
covariances, precisions, and responsibilities. A failed continuation leaves
its descendants unavailable and does not trigger a cold restart. The default
\code{"none"} preserves independent cold grid fits.

Initialization is performed in R. The numeric merge matrix returned by
\code{mclust::hc()} determines the hierarchical partition from which the
mixture state is constructed. Method, objective, seed/RNG, multistart, and
failure information is available in the \code{initializer_adapters} attribute.
}
\value{
A matrix with one variable ordering per candidate cluster count.  Attributes
record the K-specific grids, fitted states, initializer adapters, failure
information, and warm-path diagnostics needed to reproduce the ranking under
the same software and RNG settings.
}
\references{
Zhou, H., Pan, W., and Shen, X. (2009). Penalized model-based clustering with
unconstrained covariance matrices. \emph{Electronic Journal of Statistics},
3, 1473--1496.
}
\seealso{\code{\link{SelvarClustLasso}}}
\examples{
\dontrun{
set.seed(2025)
z <- rep(1:2, each = 20)
x <- cbind(
  x1 = rnorm(40, c(-1.5, 1.5)[z], 0.6),
  x2 = rnorm(40, c(1, -1)[z], 0.6),
  noise = rnorm(40)
)
ordering <- SortvarClust(
  x, nbcluster = 2, lambda = c(20, 100), rho = c(0.05, 0.2),
  initialize = "kmeans", n.start = 10, nbcores = 1,
  warm_start = "both"
)
}
}
\keyword{cluster}
\keyword{models}

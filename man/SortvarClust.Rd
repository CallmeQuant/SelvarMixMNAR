\name{SortvarClust}
\alias{SortvarClust}
\title{Rank Variables Using Penalized Gaussian Mixture Models}
\description{
Ranks variables separately for each candidate number of mixture components.
For every penalty-grid point, a penalized Gaussian mixture model is fitted and
a variable is counted as active when at least one component mean is nonzero.
Variables are ordered by their activity frequencies over numerically valid
fits. Penalty grids and initial states may be specified separately for each
candidate component count.
}
\usage{
SortvarClust(
  x,
  nbcluster,
  type = "lasso",
  lambda = seq(20, 100, by = 10),
  rho = seq(0.1, 1, length = 5),
  group_shrinkage_method = "common",
  distance_method = "Euclidean",
  lambda_omega_0 = 50,
  epsilon_weighted_by_W0 = sqrt(.Machine$double.eps),
  laplacian_target_type = c("identity", "diag_Omega_hat"),
  adj_threshold = 1e-4,
  laplacian_norm_type = c("symmetric", "unsymmetric"),
  penalize_diag = FALSE,
  initialize = "hc",
  nbcores = min(2, parallel::detectCores(logical = FALSE)),
  n.start = 250,
  scale = FALSE,
  penalty_grids = NULL,
  warm_start = c("none", "inner", "outer", "both"),
  verbose = FALSE,
  min_scorable_fraction = 0.5,
  ...,
  init_control = list()
)
}
\arguments{
  \item{x}{
    Numeric matrix or data frame with observations in rows and variables in
    columns. Missing and non-finite values are not accepted.
  }
  \item{nbcluster}{
    Positive integer vector giving the candidate numbers of mixture
    components.
  }
  \item{type}{
    Variable-ranking method. \code{"lasso"} fits penalized Gaussian mixture
    models along a penalty path. \code{"likelihood"} is not implemented and
    produces a classed error.
  }
  \item{lambda}{
    Vector of non-negative penalties for component means, used when
    \code{penalty_grids} is \code{NULL}.
  }
  \item{rho}{
    Vector of non-negative penalties for component precision matrices, used
    when \code{penalty_grids} is \code{NULL}.
  }
  \item{group_shrinkage_method}{
    Method used to construct the precision-penalty weight matrices.
    \code{"common"} applies a common unweighted penalty.
  }
  \item{distance_method}{
    Covariance-distance method used by distance-based weighting schemes.
    \code{"Euclidean"} uses the Frobenius distance computed internally;
    non-Euclidean methods require the optional \pkg{shapes} package.
  }
  \item{lambda_omega_0}{
    Positive, dimensionless multiplier used to regularize an initial
    component covariance matrix. If the covariance matrix is singular, the
    multiplier is combined with its scale so that initialization is
    equivariant under a common change of measurement units.
  }
  \item{epsilon_weighted_by_W0}{
    Positive numerical safeguard used in denominators of adaptive
    precision-penalty weights.
  }
  \item{laplacian_target_type}{
    Target matrix used by the spectral-Laplacian weighting scheme:
    \code{"identity"} or \code{"diag_Omega_hat"}.
  }
  \item{adj_threshold}{
    Non-negative threshold used to construct the adjacency matrix for
    spectral-Laplacian weights.
  }
  \item{laplacian_norm_type}{
    Normalization convention for the graph Laplacian:
    \code{"symmetric"} or \code{"unsymmetric"}.
  }
  \item{penalize_diag}{
    Logical indicating whether diagonal precision entries should be
    penalized. The current native objective penalizes only off-diagonal
    entries and therefore rejects \code{penalize_diag = TRUE}.
  }
  \item{initialize}{
    Initialization method. Supported values are \code{"mclust_hc"},
    \code{"mclust_em"}, \code{"kmeans"}, \code{"random"}, \code{"user"},
    \code{"previous_fit"}, and \code{"deterministic_multistart"}.
    The legacy value \code{"hc"} is an alias for \code{"mclust_hc"}.
  }
  \item{nbcores}{
    Number of worker processes. Grid points on the same outer warm-start
    path are fitted sequentially; paths for different component counts may
    be fitted in parallel.
  }
  \item{n.start}{
    Number of candidate starts used by initialization methods that perform
    multiple starts.
  }
  \item{scale}{
    Logical indicating whether the complete input is centered and
    standardized before ranking. The higher-level
    \code{\link{SelvarClustLasso}} interface performs its preprocessing
    before calling this function.
  }
  \item{penalty_grids}{
    Optional list with one entry for each element of \code{nbcluster}.
    Each entry supplies a component-count-specific penalty grid and its
    fixed initial state. When supplied, \code{lambda} and \code{rho} do not
    define the fitted grid.
  }
  \item{warm_start}{
    Numerical continuation rule. \code{"none"} fits every grid point from
    its fixed cold initialization. \code{"inner"} reuses the preceding
    \pkg{glassoFast} solution between EM iterations. \code{"outer"} transfers
    the complete fitted mixture state between successive penalty-grid
    points. \code{"both"} combines inner and outer continuation.
  }
  \item{verbose}{
    Logical; if \code{TRUE}, print progress messages.
  }
  \item{min_scorable_fraction}{
    Number in \eqn{[0,1]} giving the minimum proportion of a
    component-count-specific grid that must end in a valid numerical state.
    Invalid points remain in the grid with all-\code{NA} activity indicators
    and do not contribute to ranking scores.
  }
  \item{...}{
    Additional arguments retained for backward compatibility.
  }
  \item{init_control}{
    Named initialization controls, including the random seed, user-supplied
    or previous fitted states, deterministic multistart settings, and
    component-count-specific overrides. It cannot be combined with
    externally supplied \code{penalty_grids}.
  }
}
\details{
For each candidate component count \eqn{K}, the function fits penalized
Gaussian mixture models over the specified grid of \code{lambda} and
\code{rho} values. At a numerically valid grid point, variable \eqn{j} is
active when at least one fitted component mean for that variable is nonzero.
Its score is the number of valid grid points at which it is active. The
returned ordering is decreasing in this score. Invalid fits are excluded
rather than treated as evidence that a variable is inactive.

Outer continuation uses a deterministic traversal of the penalty grid.
It first follows a decreasing-\code{rho} sequence at a fixed
\code{lambda} value and then follows decreasing-\code{lambda} sequences
from the corresponding fitted states. The state transferred between
successive points contains mixing proportions, component means and
covariances, precision matrices, and posterior probabilities. If an outer
continuation fails, dependent downstream points are recorded as unavailable;
the function does not repeat them from a cold start. The default
\code{warm_start = "none"} fits grid points independently.

Initialization is performed in R. For hierarchical initialization, the merge
matrix returned by \code{mclust::hc()} determines the initial partition from
which the Gaussian mixture state is constructed. The
\code{initializer_adapters} attribute records the selected method, objective
value, random-number state, multistart information, and any initialization
failure.
}
\value{
A matrix containing one variable ordering for each candidate component
count. Its columns correspond to the values in \code{nbcluster}, and entries
give variable indices in decreasing ranking score. Attributes store the
component-count-specific penalty grids, fitted states, initializer records,
invalid-fit information, and warm-start diagnostics required to reproduce
the ranking under the same software and random-number settings.
}
\references{
Zhou, H., Pan, W., and Shen, X. (2009). Penalized model-based clustering
with unconstrained covariance matrices. \emph{Electronic Journal of
Statistics}, 3, 1473--1496.

Celeux, G., Maugis-Rabusseau, C., and Sedki, M. (2019). Variable selection
in model-based clustering and discriminant analysis with a regularization
approach. \emph{Advances in Data Analysis and Classification}, 13, 259--278.
}
\author{
Binh H. Ho, Long Nguyen Chi, TrungTin Nguyen, Binh T. Nguyen,
Van Ha Hoang, and Christopher Drovandi.

The implementation extends the variable-ranking code in the original
\pkg{SelvarMix} package by Mohammed Sedki, Gilles Celeux, and
Cathy Maugis-Rabusseau.
}
\seealso{
\code{\link{SelvarClustLasso}}
}
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
  x,
  nbcluster = 2,
  lambda = c(20, 100),
  rho = c(0.05, 0.2),
  initialize = "kmeans",
  n.start = 10,
  nbcores = 1,
  warm_start = "both"
)

ordering
}
}
\keyword{cluster}
\keyword{models}

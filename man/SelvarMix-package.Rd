\docType{package}
\name{SelvarMixMNAR-package}
\alias{SelvarMixMNAR-package}
\title{Variable selection in model-based clustering with missing not at random data}
\description{
  Direct R implementation accompanying the NeurIPS 2025 paper
  \emph{A Unified Framework for Variable Selection in Model-Based Clustering
  with Missing Not at Random}. The package ranks variables with penalized
  Gaussian mixtures, assigns SRUW roles, and fits a parsimonious class-only
  MNARz missingness mechanism.
}
\details{
  The SRUW model separates variables into a clustering set \eqn{S}, a
  redundant set \eqn{U} explained by predictors \eqn{R \subseteq S}, and an
  independent set \eqn{W}. A penalized Gaussian-mixture path ranks variables
  before the role-selection criterion is evaluated.

  The fitted missingness model is MNARz:
  \eqn{P(C_{ij}=1 \mid Z_i=k)=\rho_k}. Each latent component therefore has
  one missingness probability shared across coordinates. MNARy and
  variable-specific MNARzj require different likelihoods.

  In the separate decoupled fit, \code{mecha="mixed"} is a sensitivity
  specification based on a user-supplied logical \code{is_mnar} mask.
  Coordinates marked \code{TRUE} share \eqn{\rho_k}; missingness factors for
  coordinates marked \code{FALSE} are treated as ignorable and omitted from
  the modeled likelihood. The analyst fixes the mask, which introduces no
  coordinate-specific class effect.

  For incomplete input, \code{workflow="auto"} selects the default
  \code{"decoupled_mnarz"} route. SRUW roles are estimated on an imputed
  working data set, after which a separate role-agnostic MNARz mixture is
  fitted. The \code{"joint_mnarz"} route instead fits the SRUW and class-only
  MNARz terms under one observed-data criterion. Its forward and reverse role
  updates are scored by total observed-data BIC and generally require more
  computation than the decoupled route.

  \code{SelvarLearnLasso} and \code{SortvarLearn} are retained for call
  compatibility and signal that supervised-lasso fitting is unavailable in
  version 0.1.x.
}
\examples{
\dontrun{
## Seeded Gaussian example with S, U, and W variables
set.seed(2025)
n <- 90
z <- rep(1:3, each = n / 3)
x1 <- rnorm(n, c(-2, 0, 2)[z], 0.7)
x2 <- rnorm(n, c(1.5, -1.5, 0)[z], 0.7)
x <- cbind(
  x1, x2,
  u1 = 0.8 * x1 + rnorm(n, sd = 0.35),
  u2 = -0.7 * x2 + rnorm(n, sd = 0.35),
  w1 = rnorm(n), w2 = rnorm(n)
)
set.seed(2026)
obj <- SelvarClustLasso(
  x = x, nbcluster = 2:4, num_vals_penalty = 3,
  workflow = "imputed_sruw", nbcores = 1
)
summary(obj)
print(obj)
}
}
\author{
  Original SelvarMix implementation: Mohammed Sedki, Gilles Celeux, and
  Cathy Maugis-Rabusseau. See \code{citation("SelvarMixMNAR")} for the
  authors of the accompanying MNAR framework.
}
\references{
  Ho, B. H., Chi, L. N., Nguyen, T., Nguyen, B. T., Hoang, V. H., and
  Drovandi, C. (2025). A Unified Framework for Variable Selection in
  Model-Based Clustering with Missing Not at Random. \emph{Advances in Neural
  Information Processing Systems} (NeurIPS 2025).

  Maugis, C., Celeux, G., and Martin-Magniette, M. L. (2009). Variable
  selection in model-based clustering: A general variable role modeling.
  \emph{Computational Statistics and Data Analysis}, 53, 3872--3882.

  Maugis, C., Celeux, G., and Martin-Magniette, M. L. (2011). Variable
  selection in model-based discriminant analysis. \emph{Journal of
  Multivariate Analysis}, 102, 1374--1387.

  Zhou, H., Pan, W., and Shen, X. (2009). Penalized model-based clustering
  with unconstrained covariance matrices. \emph{Electronic Journal of
  Statistics}, 3, 1473--1496.

  Sedki, M., Celeux, G., and Maugis-Rabusseau, C. (2014). SelvarMix: An R
  package for variable selection in model-based clustering and discriminant
  analysis with a regularization approach. INRIA Research Report,
  \url{https://inria.hal.science/hal-01053784}.

}
\keyword{package}

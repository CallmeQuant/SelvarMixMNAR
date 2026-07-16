\docType{package}
\name{SelvarMixMNAR-package}
\alias{SelvarMixMNAR-package}
\alias{SelvarMixMNAR}
\title{Variable Selection in Model-Based Clustering with MNAR Data}
\description{
\pkg{SelvarMixMNAR} implements and extends regularized variable selection for
Gaussian model-based clustering when data may be incomplete. The package
combines penalized variable ranking, SRUW role assignment, and a parsimonious
class-only MNARz missingness model. It accompanies the NeurIPS 2025 paper
\emph{A Unified Framework for Variable Selection in Model-Based Clustering
with Missing Not at Random}.
}
\details{
The SRUW model assigns variables to three roles. The set \eqn{S} contains
variables that contribute directly to clustering. Variables in \eqn{U} are
redundant for clustering once a predictor subset \eqn{R \subseteq S} is
observed, and variables in \eqn{W} are independent noise variables. The
package first ranks variables by fitting penalized Gaussian mixture models
over a grid of mean and precision penalties and then applies an SRUW
model-selection procedure along the resulting ordering.

For incomplete data, the package supports a class-only MNARz mechanism,
\deqn{\Pr(C_{ij}=1 \mid Z_i=k)=\rho_k,}
where \eqn{C_{ij}} is the missingness indicator and \eqn{Z_i} is the latent
mixture component. Each component therefore has one missingness probability
shared across modeled coordinates. Mechanisms in which missingness depends
directly on unobserved values, or in which the class effect varies by
coordinate, require different likelihoods and are not fitted by the current
interface.

For sensitivity analysis in the separate MNARz fit,
\code{mecha = "mixed"} accepts a user-specified logical vector
\code{is_mnar}. Coordinates marked \code{TRUE} share the class-specific
probability \eqn{\rho_k}; missingness factors for coordinates marked
\code{FALSE} are treated as ignorable and omitted from the modeled
likelihood. The mask is fixed before estimation and does not introduce
coordinate-specific class effects.

Three estimation workflows are available. \code{"imputed_sruw"} performs
variable ranking and SRUW role selection on a complete working matrix and
does not model missingness. \code{"decoupled_mnarz"} performs the same role
selection and then fits a separate, role-agnostic class-only MNARz mixture to
the original incomplete data. \code{"joint_mnarz"} evaluates SRUW role
changes and missingness parameters under a common observed-data criterion.
With \code{workflow = "auto"}, complete data use \code{"imputed_sruw"} and
incomplete data use \code{"decoupled_mnarz"}.

The historical supervised functions \code{\link{SelvarLearnLasso}} and
\code{\link{SortvarLearn}} remain available as compatibility stubs. They
signal that supervised ranking and discriminant-analysis fitting are not
implemented in version 0.1.x.
}
\examples{
\dontrun{
## Gaussian example with clustering, redundant, and noise variables
set.seed(2025)
n <- 90
z <- rep(1:3, each = n / 3)

x1 <- rnorm(n, c(-2, 0, 2)[z], 0.7)
x2 <- rnorm(n, c(1.5, -1.5, 0)[z], 0.7)

x <- cbind(
  x1 = x1,
  x2 = x2,
  u1 = 0.8 * x1 + rnorm(n, sd = 0.35),
  u2 = -0.7 * x2 + rnorm(n, sd = 0.35),
  w1 = rnorm(n),
  w2 = rnorm(n)
)

set.seed(2026)
fit <- SelvarClustLasso(
  x = x,
  nbcluster = 2:4,
  models = "VVI",
  num_vals_penalty = 3,
  workflow = "imputed_sruw",
  nbcores = 1
)

print(fit)
summary(fit)
plot(fit, type = "roles")
}
}
\author{
Binh H. Ho, Long Nguyen Chi, TrungTin Nguyen, Binh T. Nguyen,
Van Ha Hoang, and Christopher Drovandi.

The package extends the original \pkg{SelvarMix} implementation by
Mohammed Sedki, Gilles Celeux, and Cathy Maugis-Rabusseau. Use
\code{citation("SelvarMixMNAR")} for citation information.
}
\references{
Ho, B. H., Nguyen Chi, L., Nguyen, T., Nguyen, B. T., Hoang, V. H., and
Drovandi, C. (2025). A unified framework for variable selection in
model-based clustering with missing not at random. In
\emph{Advances in Neural Information Processing Systems 38}.

Maugis, C., Celeux, G., and Martin-Magniette, M.-L. (2009).
Variable selection in model-based clustering: A general variable role
modeling. \emph{Computational Statistics & Data Analysis}, 53,
3872--3882.

Zhou, H., Pan, W., and Shen, X. (2009). Penalized model-based clustering
with unconstrained covariance matrices. \emph{Electronic Journal of
Statistics}, 3, 1473--1496.

Celeux, G., Maugis-Rabusseau, C., and Sedki, M. (2019). Variable selection
in model-based clustering and discriminant analysis with a regularization
approach. \emph{Advances in Data Analysis and Classification}, 13,
259--278.

Sportisse, A., Marbac, M., Laporte, F., Celeux, G., Boyer, C., Josse, J.,
and Biernacki, C. (2024). Model-based clustering with missing not at random
data. \emph{Statistics and Computing}, 34, 135.
}
\seealso{
\code{\link{SelvarClustLasso}},
\code{\link{SortvarClust}},
\code{\link{selvarmix}}
}
\keyword{package}

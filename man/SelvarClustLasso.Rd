\name{SelvarClustLasso}
\alias{SelvarClustLasso}
\title{Variable selection for model-based clustering}
\description{
Ranks variables with a penalized Gaussian-mixture fit and assigns SRUW roles.
With \code{workflow = "auto"}, incomplete data use imputed-data SRUW followed
by a separate role-agnostic class-only MNARz fit; complete data use
imputed-data SRUW alone. The \code{"imputed_sruw"}, \code{"decoupled_mnarz"}, and \code{"joint_mnarz"}
values select the corresponding estimators directly.
}
\usage{
SelvarClustLasso(
  x, nbcluster, strategy = NULL, lambda = NULL, rho = NULL,
  num_vals_penalty = 5, type = "lasso", hsize = 3,
  criterion = "BIC", models = "VVI",
  rmodel = c("LI", "LB"), imodel = c("LI", "LB"),
  nbcores = min(2, detectCores(all.tests = FALSE, logical = FALSE)),
  impute_missing = TRUE, use_copula = TRUE, scale_data = "always",
  scale_check_method = "pairwise.complete.obs",
  use_missing_pattern = FALSE, use_diag = TRUE, true_labels = NULL,
  sd_ratio_threshold = 10, cond_number_threshold = 30,
  rank = NULL, rank_control = list(), mnarz_control = list(),
  selection_control = list(),
  workflow = c("auto", "decoupled_mnarz", "imputed_sruw", "joint_mnarz",
    "decoupled"), joint_control = list(),
  init_control = list(), verbose = FALSE
)
}
\arguments{
  \item{x}{Numeric matrix or data frame; rows are observations and columns are variables.}
  \item{nbcluster}{Positive integer vector of candidate cluster counts.}
  \item{strategy}{Deprecated compatibility argument. A non-\code{NULL} value
  produces a warning and is otherwise ignored. Configure Gaussian
  initialization with \code{init_control}.}
  \item{lambda}{Optional non-negative mean-penalty path.}
  \item{rho}{Optional non-negative precision-penalty path.}
  \item{num_vals_penalty}{Number of automatically constructed penalty values per dimension.}
  \item{type}{Ranking method. \code{"lasso"} fits the penalized
  Gaussian-mixture path. \code{"likelihood"} is unavailable and raises a
  classed error.}
  \item{hsize}{Positive rolling count of consecutive nonpositive role increments used to stop a selection pass.}
  \item{criterion}{Selection criterion. The joint route accepts only \code{"BIC"}.}
  \item{models}{Mixture-model specification for SRUW selection. The default
  \code{"VVI"} requests the diagonal Mclust model. Mclust and Rmixmod models
  are supported. MixAll tokens stop before native execution because an STK
  failure may terminate the R session.}
  \item{rmodel}{Candidate covariance forms for the common regression of U on R.}
  \item{imodel}{Candidate covariance forms for independent W variables.}
  \item{nbcores}{Number of workers.  Dependent points of a warm penalty path are always sequential; warm paths may be parallelized across cluster counts.}
  \item{impute_missing}{Whether to construct the complete data set required by the ranking stage.}
  \item{use_copula}{Use \pkg{gcimputeR}; otherwise use \pkg{missRanger}.}
  \item{scale_data}{Preprocessing mode: \code{"always"} (default) centers by
  the observed mean and scales by the observed standard deviation;
  \code{"auto"} always centers and conditionally scales; \code{"never"}
  preserves the input scale. Logical \code{TRUE}/\code{FALSE} map to
  \code{"always"}/\code{"never"}.}
  \item{scale_check_method}{Missing-data covariance method used by the scale diagnostic.}
  \item{use_missing_pattern}{Deprecated compatibility switch for
  \code{workflow="decoupled"}. Do not combine it with another workflow value.}
  \item{use_diag}{Logical; use diagonal covariance matrices in the MNARz fit
  when \code{TRUE} and full covariance matrices when \code{FALSE}.}
  \item{true_labels}{Deprecated evaluation-only input.  If supplied, it is ignored during fitting and a warning is emitted.}
  \item{sd_ratio_threshold}{Scale-ratio threshold used by the scale diagnostic.}
  \item{cond_number_threshold}{Condition-number threshold used by the scale diagnostic.}
  \item{rank}{Optional permutation of \code{seq_len(ncol(x))}; if omitted or \code{NULL}, ranking is fitted.}
  \item{rank_control}{Named list passed to \code{SortvarClust}. Its
  \code{warm_start} entry is one of \code{"none"}, \code{"inner"},
  \code{"outer"}, or \code{"both"}; the default is \code{"none"}. Failed
  grid points retain their positions but are omitted from the ranking score;
  at least \code{min_scorable_fraction = 0.5} of each K-specific grid must
  remain.}
  \item{init_control}{Strict named initialization controls. Supported methods are \code{"mclust_hc"}, \code{"mclust_em"}, \code{"kmeans"}, \code{"random"}, \code{"user"}, \code{"previous_fit"}, and \code{"deterministic_multistart"}; legacy \code{"hc"} is an alias for \code{"mclust_hc"}. Optional \code{per_k} entries provide controls keyed by cluster count.}
  \item{mnarz_control}{Named controls for the separate role-agnostic
  missing-data fit in the default workflow. \code{mecha = "MNARz"} is the
  default. For a prespecified sensitivity analysis, \code{mecha = "mixed"}
  requires a logical \code{is_mnar} vector. Coordinates marked \code{TRUE}
  share the class-specific MNARz probability; factors for coordinates marked
  \code{FALSE} are treated as ignorable and omitted. The analyst fixes the
  mask, which introduces no coordinate-specific class effect.}
  \item{selection_control}{Named list controlling SRUW selection.
  \code{stopping = "consecutive"} and deterministic \code{seed = 1} are the
  defaults; \code{"legacy_block"} reproduces the former blockwise rule.}
  \item{workflow}{\code{"auto"}, \code{"imputed_sruw"},
  \code{"decoupled_mnarz"}, or \code{"joint_mnarz"}. The legacy
  \code{"decoupled"} token is deprecated. An explicit MNARz workflow requires
  a nondegenerate missingness mask.}
  \item{joint_control}{Named controls \code{diag}, \code{rmax}, \code{tol}, and \code{r_search_cap} for the joint selector and fixed-role EM.}
  \item{verbose}{Logical; print progress messages and the fitted variable
  ranking when \code{TRUE}.}
}
\details{
The four ranking start modes determine how numerical state is reused.
\code{"inner"} supplies the preceding graphical-lasso solution to the next EM
iteration. \code{"outer"} transfers the complete fitted mixture state between
adjacent successful penalty-grid points. \code{"both"} combines these forms
of continuation. A failed continuation is not retried from a cold start.

In the decoupled route, SRUW roles are selected on imputed data and the MNARz
fit is role-agnostic. The joint route uses a single observed likelihood for
the Gaussian SRUW factorization and a class-only
missingness probability shared over mask coordinates. Its criterion is
\eqn{-2\ell_{\mathrm{obs}} + q\log(n)} and is minimized. Forward and reverse
role updates are evaluated by this total observed-data BIC.

The default distance-weighted precision penalty
\code{"weighted_by_dist_to_diag_W0"} is retained. Its covariance distance is
the internal Frobenius distance for the default \code{"Euclidean"} method;
non-Euclidean distance methods are delegated to the optional \pkg{shapes}
package.

If the common regression selects an empty predictor set \eqn{R}, the fitted
\eqn{U} block is not merely relabeled as \eqn{W}. Instead, \eqn{U\cup W} is
refitted under every requested independent covariance model before cluster-
count candidates are compared.

Every initialization method is executed in R and returns one common numeric
state: proportions, means, covariance and precision arrays, responsibilities,
partition, objective, method, seed/RNG provenance, and failure status. For
\code{"mclust_em"}, the numeric merge matrix returned by \code{mclust::hc()} is
passed through \code{initialization = list(hcPairs = ...)} to
\code{mclust::Mclust()}. The same selected Gaussian state is pinned in automatic
grid construction and ranking and, by default, is adapted for the separate
MNAR fit or reused as the joint selector's partition. A supplied
\code{mnarz_control$init} replaces this default for the separate fit.
}
\value{
For one criterion, a validated object of class \code{selvarmix}. For multiple
criteria, a named \code{selvarmix_collection} whose elements are validated
\code{selvarmix} objects. Core fields are \code{S}, \code{R},
\code{U}, \code{W}, \code{criterionValue}, \code{nbcluster},
\code{partition}, \code{proba}, \code{parameters}, and \code{imputedData}.
The \code{schemaVersion} field identifies the result structure, and
\code{diagnostics} provides normalized convergence, criterion-availability,
selection-stopping, and ranking-failure information.

The result also records the requested and effective workflows, preprocessing
center and scale, inverse transformation, criterion convention and scope,
selection diagnostics, backend and covariance-model identifiers, and
initialization registry. Joint fits additionally retain \code{jointFit} and
\code{jointFitsByK}. The \code{ranking} matrix
retains the per-K variable order, penalty grids, fit status, failed points,
objective traces, initializers, and warm-path metadata from the
variable-ranking stage.
See \code{\link{selvarmix}} for the complete field description and the
standard print, summary, and plot methods.
}
\references{
Zhou, H., Pan, W., and Shen, X. (2009). Penalized model-based clustering with
unconstrained covariance matrices. \emph{Electronic Journal of Statistics},
3, 1473--1496.

Maugis, C., Celeux, G., and Martin-Magniette, M. L. (2009). Variable selection
in model-based clustering: A general variable role modeling. \emph{Computational
Statistics and Data Analysis}, 53, 3872--3882.
}
\author{Mohammed Sedki, Gilles Celeux, and Cathy Maugis-Rabusseau}
\seealso{\code{\link{SortvarClust}}}
\examples{
set.seed(2025)
z <- rep(1:2, each = 20)
x1 <- rnorm(40, c(-1.5, 1.5)[z], 0.6)
x2 <- rnorm(40, c(1, -1)[z], 0.6)
x <- cbind(x1, x2, redundant = x1 + rnorm(40, sd = 0.3),
           noise = rnorm(40))
fit <- SelvarClustLasso(
  x, nbcluster = 2, models = "VVI", rank = seq_len(ncol(x)),
  hsize = 2, nbcores = 1, workflow = "imputed_sruw"
)
summary(fit)
}
\keyword{cluster}
\keyword{models}

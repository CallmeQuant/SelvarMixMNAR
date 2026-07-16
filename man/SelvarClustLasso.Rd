\name{SelvarClustLasso}
\alias{SelvarClustLasso}
\title{Variable Selection for Gaussian Model-Based Clustering}
\description{
Fits a variable-selection procedure for Gaussian model-based clustering.
Variables are first ranked using penalized Gaussian mixture models and are
then assigned to the SRUW roles: clustering variables \eqn{S}, redundant
variables \eqn{U} explained by a subset \eqn{R \subseteq S}, and independent
variables \eqn{W}. For incomplete data, the function can supplement the SRUW
analysis with a class-only MNARz model, fitted either separately or jointly.
With \code{workflow = "auto"}, complete data use \code{"imputed_sruw"} and
incomplete data use \code{"decoupled_mnarz"}.
}
\usage{
SelvarClustLasso(
  x,
  nbcluster,
  strategy = NULL,
  lambda = NULL,
  rho = NULL,
  num_vals_penalty = 5,
  type = "lasso",
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
  sd_ratio_threshold = 10,
  cond_number_threshold = 30,
  rank = NULL,
  rank_control = list(),
  mnarz_control = list(),
  selection_control = list(),
  workflow = c(
    "auto", "decoupled_mnarz", "imputed_sruw", "joint_mnarz",
    "decoupled"
  ),
  joint_control = list(),
  init_control = list(),
  verbose = FALSE
)
}
\arguments{
  \item{x}{
    Numeric matrix or data frame with observations in rows and variables in
    columns.
  }
  \item{nbcluster}{
    Positive integer vector giving the candidate numbers of mixture
    components.
  }
  \item{strategy}{
    Deprecated compatibility argument. A non-\code{NULL} value produces a
    warning and is otherwise ignored. Use \code{init_control} to configure
    initialization.
  }
  \item{lambda}{
    Optional vector of non-negative penalties for the component means.
  }
  \item{rho}{
    Optional vector of non-negative penalties for the component precision
    matrices.
  }
  \item{num_vals_penalty}{
    Number of values generated for each penalty dimension when an explicit
    path is not supplied.
  }
  \item{type}{
    Variable-ranking method. \code{"lasso"} fits penalized Gaussian mixture
    models along a penalty path. \code{"likelihood"} is not implemented and
    produces a classed error.
  }
  \item{hsize}{
    Positive integer giving the number of consecutive non-positive role
    increments that terminates a selection scan.
  }
  \item{criterion}{
    Model-selection criterion or criteria accepted by the selected SRUW
    backend. The joint workflow supports only \code{"BIC"}.
  }
  \item{models}{
    Gaussian mixture covariance model or models used for SRUW selection.
    The default \code{"VVI"} is the diagonal \pkg{mclust} model. Model
    identifiers supported by \pkg{mclust} and \pkg{Rmixmod} may be used.
    MixAll/STK model identifiers are rejected because they cannot be
    evaluated safely by the current interface.
  }
  \item{rmodel}{
    Candidate covariance forms for the common regression of \eqn{U} on
    \eqn{R}. \code{"LI"} denotes a spherical covariance matrix and
    \code{"LB"} a diagonal covariance matrix.
  }
  \item{imodel}{
    Candidate covariance forms for the independent block \eqn{W}.
    \code{"LI"} denotes a spherical covariance matrix and \code{"LB"} a
    diagonal covariance matrix.
  }
  \item{nbcores}{
    Number of worker processes. Grid points belonging to the same outer
    warm-start path are fitted sequentially; independent paths for different
    component counts may be fitted in parallel.
  }
  \item{impute_missing}{
    Logical indicating whether an incomplete input should be completed for
    the ranking stage. Variable ranking requires a complete finite matrix.
  }
  \item{use_copula}{
    Logical indicating which imputation method is used when completion is
    required. If \code{TRUE}, imputation uses \pkg{gcimputeR}; otherwise it
    uses \pkg{missRanger}.
  }
  \item{scale_data}{
    Preprocessing rule. \code{"always"} centers each variable by its observed
    mean and scales it by its observed standard deviation. \code{"auto"}
    always centers and scales only when indicated by the scale diagnostic.
    \code{"never"} preserves the original scale. Logical values
    \code{TRUE} and \code{FALSE} are interpreted as \code{"always"} and
    \code{"never"}, respectively.
  }
  \item{scale_check_method}{
    Covariance method used by the scale diagnostic when the input contains
    missing values.
  }
  \item{use_missing_pattern}{
    Deprecated compatibility switch for \code{workflow = "decoupled"}.
    It must not be combined with another explicit workflow.
  }
  \item{use_diag}{
    Logical indicating whether the separate MNARz fit uses diagonal
    covariance matrices. If \code{FALSE}, full covariance matrices are used.
  }
  \item{true_labels}{
    Deprecated evaluation-only argument. Supplied labels are not used for
    estimation and produce a warning.
  }
  \item{sd_ratio_threshold}{
    Positive threshold for the ratio of the largest to the smallest variable
    standard deviation in the scale diagnostic.
  }
  \item{cond_number_threshold}{
    Positive condition-number threshold used by the scale diagnostic.
  }
  \item{rank}{
    Optional permutation of \code{seq_len(ncol(x))} specifying a variable
    ordering. If \code{NULL}, the ordering is estimated by
    \code{\link{SortvarClust}}.
  }
  \item{rank_control}{
    Named list of controls passed to \code{\link{SortvarClust}}.
    \code{warm_start} may be \code{"none"}, \code{"inner"},
    \code{"outer"}, or \code{"both"}; the default is \code{"none"}.
    Numerically invalid grid points retain their locations but do not
    contribute to ranking scores. By default, at least
    \code{min_scorable_fraction = 0.5} of the grid for each candidate
    component count must remain valid.
  }
  \item{mnarz_control}{
    Named controls for the separate, role-agnostic missing-data fit used by
    \code{"decoupled_mnarz"}. The default is \code{mecha = "MNARz"}.
    For a prespecified sensitivity analysis, \code{mecha = "mixed"} requires
    a logical vector \code{is_mnar}. Coordinates marked \code{TRUE} share a
    class-specific missingness probability. Missingness factors for
    coordinates marked \code{FALSE} are treated as ignorable and omitted
    from the modeled likelihood. The mask is fixed by the analyst; it does
    not introduce coordinate-specific class effects.
  }
  \item{selection_control}{
    Named list controlling SRUW role selection. The defaults are
    \code{stopping = "consecutive"} and \code{seed = 1}.
    \code{stopping = "legacy_block"} reproduces the earlier blockwise
    stopping rule.
  }
  \item{workflow}{
    Estimation workflow. \code{"imputed_sruw"} performs variable ranking and
    SRUW role selection on a complete working matrix without a missingness
    likelihood. \code{"decoupled_mnarz"} adds a separate role-agnostic MNARz
    fit to the original incomplete data. \code{"joint_mnarz"} estimates the
    SRUW and class-only MNARz terms under a common observed-data criterion.
    \code{"auto"} selects \code{"imputed_sruw"} for complete data and
    \code{"decoupled_mnarz"} for incomplete data. The legacy value
    \code{"decoupled"} is deprecated. An explicit MNARz workflow requires a
    missingness mask containing both observed and missing entries.
  }
  \item{joint_control}{
    Named controls for the joint selector and its fixed-role EM fits.
    Recognized entries are \code{diag}, \code{rmax}, \code{tol}, and
    \code{r_search_cap}.
  }
  \item{init_control}{
    Named initialization controls. Supported methods are
    \code{"mclust_hc"}, \code{"mclust_em"}, \code{"kmeans"},
    \code{"random"}, \code{"user"}, \code{"previous_fit"}, and
    \code{"deterministic_multistart"}. The legacy value \code{"hc"} is an
    alias for \code{"mclust_hc"}. Optional \code{per_k} entries specify
    controls for individual component counts.
  }
  \item{verbose}{
    Logical; if \code{TRUE}, print progress messages and the estimated
    variable ordering.
  }
}
\details{
The procedure has two principal stages. First, penalized Gaussian mixture
models are fitted over a grid of mean and precision penalties. For each
candidate component count, variables are ordered by the frequency with which
their component means remain nonzero over numerically valid grid points.
Second, the ordered variables are assigned to the SRUW roles by comparing
mixture, regression, and independent Gaussian models.

The three workflows differ in how missingness enters the analysis.
\code{"imputed_sruw"} uses a complete working matrix for both ranking and
SRUW role selection and does not include a missingness model.
\code{"decoupled_mnarz"} performs the same role selection and then fits a
separate class-only MNARz mixture to the original incomplete data; the
missing-data fit does not use the estimated SRUW roles.
\code{"joint_mnarz"} evaluates role changes using an observed-data
likelihood that combines the Gaussian SRUW factorization with one
component-specific missingness probability shared across modeled mask
coordinates.

For the joint workflow, the information criterion is
\deqn{-2\ell_{\mathrm{obs}} + q\log(n),}
where \eqn{\ell_{\mathrm{obs}}} is the observed-data log-likelihood and
\eqn{q} is the number of fitted parameters. The criterion is minimized.
Both forward and reverse role changes are evaluated using this total
observed-data BIC.

Warm starts control numerical continuation in the ranking stage.
\code{"inner"} initializes each graphical-lasso update from the preceding
EM iterate. \code{"outer"} transfers the complete fitted mixture state
between successive points of a deterministic penalty path.
\code{"both"} combines the two forms of continuation. A failed outer
continuation is recorded and is not repeated from a cold start.

The default distance-weighted precision penalty is
\code{"weighted_by_dist_to_diag_W0"}. With
\code{distance_method = "Euclidean"}, covariance distance is the Frobenius
distance computed internally. Other distance methods require the optional
\pkg{shapes} package.

If the common regression selects no predictors, so that \eqn{R = \varnothing},
the fitted \eqn{U} block is not simply reclassified as \eqn{W}. Instead,
\eqn{U \cup W} is refitted under each requested independent covariance model
before candidates with different component counts are compared.

Initializers are constructed in R and represented by a common Gaussian
mixture state containing mixing proportions, component means and
covariances, precision matrices, posterior probabilities, and hard
assignments, together with the random-number and failure information needed
for reproducibility. For \code{"mclust_em"}, the hierarchical merge matrix
returned by \code{mclust::hc()} is supplied to
\code{mclust::Mclust()} through \code{initialization = list(hcPairs = ...)}.
The selected state is used consistently in automatic grid construction and
penalized ranking. Unless replaced by \code{mnarz_control$init}, it is also
adapted to initialize the separate MNARz fit; joint selection reuses its hard
partition by default.
}
\value{
If a single criterion is requested, an object of class \code{selvarmix}.
If several criteria are requested, a named object of class
\code{selvarmix_collection} whose elements are \code{selvarmix} fits.

The principal fields of a \code{selvarmix} object are \code{S}, \code{R},
\code{U}, and \code{W}; the selected component count \code{nbcluster};
the criterion value \code{criterionValue}; hard assignments
\code{partition}; posterior component-membership probabilities
\code{proba}; fitted \code{parameters}; and, when applicable, the completed
working matrix \code{imputedData}. The \code{schemaVersion} field identifies
the result schema, and \code{diagnostics} records convergence, criterion
availability, stopping behavior, and ranking failures.

The object also records the requested and effective workflows, preprocessing
centers and scales, inverse transformations, criterion convention and scope,
backend and covariance-model identifiers, and initialization metadata. Joint
fits additionally contain \code{jointFit} and \code{jointFitsByK}. The
\code{ranking} component stores the per-component-count variable ordering,
penalty grids, fit status, objective traces, initializers, and warm-start
metadata. See \code{\link{selvarmix}} for the complete object description and
the associated print, summary, and plot methods.
}
\references{
Maugis, C., Celeux, G., and Martin-Magniette, M.-L. (2009).
Variable selection in model-based clustering: A general variable role
modeling. \emph{Computational Statistics & Data Analysis}, 53, 3872--3882.

Zhou, H., Pan, W., and Shen, X. (2009). Penalized model-based clustering
with unconstrained covariance matrices. \emph{Electronic Journal of
Statistics}, 3, 1473--1496.

Celeux, G., Maugis-Rabusseau, C., and Sedki, M. (2019). Variable selection
in model-based clustering and discriminant analysis with a regularization
approach. \emph{Advances in Data Analysis and Classification}, 13, 259--278.

Sportisse, A., Marbac, M., Laporte, F., Celeux, G., Boyer, C., Josse, J.,
and Biernacki, C. (2024). Model-based clustering with missing not at random
data. \emph{Statistics and Computing}, 34, 135.

Ho, B. H., Nguyen Chi, L., Nguyen, T., Nguyen, B. T., Hoang, V. H., and
Drovandi, C. (2025). A unified framework for variable selection in
model-based clustering with missing not at random. In
\emph{Proceedings of the 39th Conference on Neural Information Processing
Systems}.
}
\author{
Binh H. Ho, Long Nguyen Chi, TrungTin Nguyen, Binh T. Nguyen,
Van Ha Hoang, and Christopher Drovandi.

The implementation extends the original \pkg{SelvarMix} package by
Mohammed Sedki, Gilles Celeux, and Cathy Maugis-Rabusseau.
}
\seealso{
\code{\link{SortvarClust}}, \code{\link{selvarmix}}
}
\examples{
set.seed(2025)
z <- rep(1:2, each = 20)
x1 <- rnorm(40, c(-1.5, 1.5)[z], 0.6)
x2 <- rnorm(40, c(1, -1)[z], 0.6)
x <- cbind(
  x1,
  x2,
  redundant = x1 + rnorm(40, sd = 0.3),
  noise = rnorm(40)
)

fit <- SelvarClustLasso(
  x,
  nbcluster = 2,
  models = "VVI",
  rank = seq_len(ncol(x)),
  hsize = 2,
  nbcores = 1,
  workflow = "imputed_sruw"
)

summary(fit)
}
\keyword{cluster}
\keyword{models}

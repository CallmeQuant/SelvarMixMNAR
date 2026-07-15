\name{selvarmix-methods}
\alias{selvarmix-methods}
\alias{selvarmix}
\alias{selvarmix_collection}
\alias{print.selvarmix}
\alias{summary.selvarmix}
\alias{print.summary.selvarmix}
\alias{plot.selvarmix}
\alias{print.selvarmix_collection}
\alias{summary.selvarmix_collection}
\alias{print.summary.selvarmix_collection}
\alias{plot.selvarmix_collection}
\title{Results and Standard Methods for SelvarMixMNAR Fits}
\description{
Successful calls to \code{SelvarClustLasso()} return a
\code{selvarmix} object. When more than one model-selection criterion is
requested, the fitted objects are stored in a named
\code{selvarmix_collection}. Standard methods give a compact textual summary
and base-graphics displays without modifying the fitted object.
}
\usage{
\method{print}{selvarmix}(x, ..., digits = 4L, max.roles = 12L)

\method{summary}{selvarmix}(object, ...)

\method{print}{summary.selvarmix}(x, ..., digits = 4L, max.roles = 12L)

\method{plot}{selvarmix}(
  x,
  type = c("roles", "classification", "missingness"),
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  col = NULL,
  ...
)

\method{print}{selvarmix_collection}(x, ..., digits = 4L)

\method{summary}{selvarmix_collection}(object, ...)

\method{print}{summary.selvarmix_collection}(x, ..., digits = 4L)

\method{plot}{selvarmix_collection}(x, which = 1L, ...)
}
\arguments{
  \item{x, object}{A fitted \code{selvarmix} object, a named
  \code{selvarmix_collection}, or the corresponding summary object.}
  \item{digits}{Number of significant digits used to print criterion values.}
  \item{max.roles}{Maximum number of variable indices or names printed for
  each role.}
  \item{type}{Graphical summary. \code{"roles"} displays role counts;
  \code{"classification"} plots posterior probabilities after a stable sort
  by the assigned component; \code{"missingness"} displays fitted MNARz
  missingness probabilities when they are available.}
  \item{main, xlab, ylab, col}{Optional base-graphics labels and colors.}
  \item{which}{One-based position or name of the fit to plot from a
  \code{selvarmix_collection}.}
  \item{...}{Additional graphical arguments for \code{plot()}, or arguments
  reserved for extensions of the print and summary methods.}
}
\value{
\code{summary()} returns a structured \code{summary.selvarmix} object or a
named \code{summary.selvarmix_collection}. The print methods return their
input invisibly. The plot methods return the plotted counts, posterior
probability matrix, or missingness probabilities invisibly.
}
\section{Result structure}{
The result schema is identified by \code{schemaVersion}. Version 1.0 contains
the following fields:
\describe{
  \item{\code{S}, \code{U}, and \code{W}}{A partition of the variables into
  clustering, redundant, and independent roles.}
  \item{\code{R}}{The subset of \code{S} used as predictors for the
  regression model of the \code{U} variables.}
  \item{\code{criterion}, \code{criterionValue}}{The selected criterion and
  its fitted value. The direction and statistical scope are recorded in
  \code{criterionConvention} and \code{criterionScope}.}
  \item{\code{model}, \code{framework}, \code{requestedModel},
  \code{effectiveModel}}{The covariance-model request and the backend that
  fitted it. Backend-specific entries may be \code{NULL} when they do not
  apply.}
  \item{\code{nbcluster}, \code{partition}, \code{proba}}{The selected number
  of components, MAP partition, and posterior membership probabilities.}
  \item{\code{parameters}, \code{regparameters}}{Gaussian mixture and
  regression parameters in the representation returned by the selected
  backend.}
  \item{\code{imputedData}}{The completed data matrix on the original scale.}
  \item{\code{workflow}, \code{workflowRequested},
  \code{workflowEffective}}{The requested and effective fitting workflows.}
  \item{\code{diagnostics}}{A normalized diagnostic list containing
  \code{status}, \code{converged}, \code{iterations},
  \code{terminationReason}, \code{criterionAvailable}, criterion metadata,
  selection stopping information, and the number of failed ranking-grid rows
  when recorded.}
}
Advanced fields such as \code{selectionResult}, \code{parametersMNARz},
\code{jointFit}, \code{jointFitsByK}, \code{ranking}, and
\code{initialization} retain detailed backend state. They are intended for
diagnosis and reproducibility rather than as a compact interchange format.
}
\section{Validity and failed fits}{
A required ranking, role-selection, or final-model fit that is unscorable
signals a classed condition and does not produce a \code{selvarmix} object.
A returned object may record \code{diagnostics$status = "nonconverged"} or
\code{"criterion_unavailable"} when a diagnostic fit has no valid selection
criterion. Only results with an available criterion enter model comparison.
}
\details{
The primary roles \code{S}, \code{U}, and \code{W} are disjoint. The
regression-predictor role \code{R} is nested within \code{S}, so its count is
shown separately in the role plot. Posterior-probability plots use a stable
ordering by the assigned component and original row position; they do not
relabel mixture components.

The \code{"missingness"} plot is available for pure class-only MNARz fits and
for diagnostic objects that retain compatible missingness parameters. It
displays fitted probabilities; adequacy of the missingness model and
identifiability require separate assessment.
}
\examples{
fit <- structure(
  list(
    S = 1L, R = 1L, U = 2L, W = 3L,
    criterionValue = 42.1, criterion = "BIC", model = "VVI",
    framework = "Mclust", requestedModel = "VVI", effectiveModel = "VVI",
    rmodel = "LI", imodel = "LI", parameters = list(), nbcluster = 2L,
    partition = c(1L, 1L, 2L, 2L),
    proba = rbind(c(.9, .1), c(.8, .2), c(.2, .8), c(.1, .9)),
    regparameters = matrix(c(0, 1), ncol = 1L),
    imputedData = matrix(seq_len(12), nrow = 4L),
    workflow = "imputed_sruw", workflowRequested = "imputed_sruw",
    workflowEffective = "imputed_sruw", criterionConvention = "maximize",
    criterionScope = "SRUW criterion on the completed data",
    schemaVersion = "1.0",
    diagnostics = list(
      status = "completed", converged = NA, iterations = NA_integer_,
      terminationReason = NA_character_, criterionAvailable = TRUE,
      workflow = "imputed_sruw", criterionConvention = "maximize",
      criterionScope = "SRUW criterion on the completed data",
      covarianceAdjustments = NA_integer_,
      minEffectiveComponentSize = NA_real_,
      selection = list(), ranking = list()
    )
  ),
  class = "selvarmix"
)
print(fit)
summary(fit)
plot(fit, type = "roles")
}
\seealso{\code{SelvarClustLasso()}}
\keyword{models}
\keyword{cluster}

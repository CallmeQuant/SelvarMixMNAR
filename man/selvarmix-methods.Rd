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
\title{SelvarMixMNAR Fit Objects and Standard Methods}
\description{
A successful call to \code{\link{SelvarClustLasso}} returns an object of
class \code{selvarmix}. If several model-selection criteria are requested,
the corresponding fits are stored in a named
\code{selvarmix_collection}. Print, summary, and plot methods provide
compact representations of fitted roles, clustering results, diagnostics,
and, when available, class-specific missingness probabilities.
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
  \item{x, object}{
    A fitted \code{selvarmix} object, a named
    \code{selvarmix_collection}, or the corresponding summary object.
  }
  \item{digits}{
    Number of significant digits used when printing criterion values and
    other numerical summaries.
  }
  \item{max.roles}{
    Maximum number of variable indices or names printed for each SRUW role.
  }
  \item{type}{
    Type of graphical summary. \code{"roles"} displays the numbers of
    variables assigned to \eqn{S}, \eqn{R}, \eqn{U}, and \eqn{W}.
    \code{"classification"} displays posterior component-membership
    probabilities after ordering observations by their assigned component.
    \code{"missingness"} displays fitted class-specific MNARz probabilities
    when compatible missingness parameters are available.
  }
  \item{main, xlab, ylab, col}{
    Optional title, axis labels, and colors passed to the base-graphics
    method.
  }
  \item{which}{
    Name or one-based position of the fit to plot from a
    \code{selvarmix_collection}.
  }
  \item{...}{
    Additional graphical arguments passed to the selected plot method, or
    arguments reserved for extensions of the print and summary methods.
  }
}
\value{
\code{summary()} returns an object of class
\code{summary.selvarmix}, or a named
\code{summary.selvarmix_collection} for a collection of fits.

The print methods return their input invisibly. The plot methods invisibly
return the data used for the display: role counts, a posterior-probability
matrix, or fitted missingness probabilities, depending on \code{type}.
}
\section{Core result structure}{
The structure is identified by \code{schemaVersion}. Schema version 1.0
contains the following principal fields:
\describe{
  \item{\code{S}, \code{U}, and \code{W}}{
    Disjoint sets of clustering, redundant, and independent variables.
  }
  \item{\code{R}}{
    The subset of \code{S} used as predictors in the regression model for
    variables in \code{U}.
  }
  \item{\code{criterion}, \code{criterionValue}}{
    The selected model-selection criterion and its fitted value. The
    optimization direction and statistical scope are recorded in
    \code{criterionConvention} and \code{criterionScope}.
  }
  \item{\code{model}, \code{framework}, \code{requestedModel},
  \code{effectiveModel}}{
    Requested covariance model, fitting backend, and effective
    backend-specific model. Entries that do not apply to a particular
    backend may be \code{NULL}.
  }
  \item{\code{nbcluster}, \code{partition}, \code{proba}}{
    Selected number of mixture components, hard component assignments, and
    posterior component-membership probabilities.
  }
  \item{\code{parameters}, \code{regparameters}}{
    Fitted Gaussian mixture and regression parameters in the representation
    returned by the selected backend.
  }
  \item{\code{imputedData}}{
    Completed working matrix on the original scale, when a completed matrix
    is retained by the selected workflow.
  }
  \item{\code{workflow}, \code{workflowRequested},
  \code{workflowEffective}}{
    Workflow recorded for backward compatibility, together with the
    requested and effective workflows.
  }
  \item{\code{diagnostics}}{
    Normalized diagnostics including \code{status}, \code{converged},
    \code{iterations}, \code{terminationReason},
    \code{criterionAvailable}, criterion metadata, selection-stopping
    information, and the number of invalid ranking-grid rows when recorded.
  }
}

Additional fields may include \code{selectionResult},
\code{parametersMNARz}, \code{jointFit}, \code{jointFitsByK},
\code{ranking}, and \code{initialization}. These fields expose detailed
backend state for diagnosis and reproducibility. Code intended to work
across package versions should rely primarily on the documented core fields
and inspect \code{schemaVersion}.
}
\section{Validity and unsuccessful fits}{
If a required ranking, role-selection, or final-model fit cannot be assigned
a valid criterion, the fitting function signals a classed condition and
does not return a \code{selvarmix} object. A returned object may nevertheless
record unsuccessful auxiliary or diagnostic fits. Such cases are summarized
through values such as \code{diagnostics$status = "nonconverged"} or
\code{"criterion_unavailable"}. Only candidates with an available criterion
are used in model comparison.
}
\details{
The primary SRUW roles \code{S}, \code{U}, and \code{W} form a disjoint
partition. Because \code{R} is a predictor subset nested within \code{S},
the role plot reports its size separately rather than treating it as a
fourth variable role.

For \code{type = "classification"}, observations are ordered stably by hard
component assignment and then by their original row positions. The display
does not permute or relabel mixture components.

The \code{"missingness"} plot is available for class-only MNARz fits and for
diagnostic objects that retain a compatible set of missingness parameters.
It reports fitted probabilities only. The plot does not assess the adequacy
or identifiability of the missingness mechanism.
}
\examples{
fit <- structure(
  list(
    S = 1L,
    R = 1L,
    U = 2L,
    W = 3L,
    criterionValue = 42.1,
    criterion = "BIC",
    model = "VVI",
    framework = "Mclust",
    requestedModel = "VVI",
    effectiveModel = "VVI",
    rmodel = "LI",
    imodel = "LI",
    parameters = list(),
    nbcluster = 2L,
    partition = c(1L, 1L, 2L, 2L),
    proba = rbind(
      c(0.9, 0.1),
      c(0.8, 0.2),
      c(0.2, 0.8),
      c(0.1, 0.9)
    ),
    regparameters = matrix(c(0, 1), ncol = 1L),
    imputedData = matrix(seq_len(12), nrow = 4L),
    workflow = "imputed_sruw",
    workflowRequested = "imputed_sruw",
    workflowEffective = "imputed_sruw",
    criterionConvention = "maximize",
    criterionScope = "SRUW criterion on the completed data",
    schemaVersion = "1.0",
    diagnostics = list(
      status = "completed",
      converged = NA,
      iterations = NA_integer_,
      terminationReason = NA_character_,
      criterionAvailable = TRUE,
      workflow = "imputed_sruw",
      criterionConvention = "maximize",
      criterionScope = "SRUW criterion on the completed data",
      covarianceAdjustments = NA_integer_,
      minEffectiveComponentSize = NA_real_,
      selection = list(),
      ranking = list()
    )
  ),
  class = "selvarmix"
)

print(fit)
summary(fit)
plot(fit, type = "roles")
}
\seealso{
\code{\link{SelvarClustLasso}},
\code{\link{SortvarClust}}
}
\keyword{models}
\keyword{cluster}

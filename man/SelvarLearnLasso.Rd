\name{SelvarLearnLasso}
\alias{SelvarLearnLasso}
\title{Supervised SelvarMix Estimator (Currently Unavailable)}
\description{
Compatibility stub for the supervised variable-selection and discriminant-analysis
interface provided by the original \pkg{SelvarMix} package. Supervised estimation
has not yet been implemented in \pkg{SelvarMixMNAR} 0.1.x, and calling this
function signals a classed error.
}
\usage{
SelvarLearnLasso(
  x,
  z,
  lambda = seq(20, 100, by = 10),
  rho = seq(1, 2, length = 2),
  type = "lasso",
  rank,
  hsize = 3,
  models = mixmodGaussianModel(listModels = c(
    "Gaussian_pk_L_C", "Gaussian_pk_Lk_C",
    "Gaussian_pk_L_Ck", "Gaussian_pk_Lk_Ck"
  )),
  rmodel = c("LI", "LB", "LC"),
  imodel = c("LI", "LB"),
  xtest = x,
  ztest = z,
  nbcores = min(
    2,
    detectCores(all.tests = FALSE, logical = FALSE)
  )
)
}
\arguments{
  \item{x, z, lambda, rho, type, rank, hsize, models, rmodel, imodel,
  xtest, ztest, nbcores}{
    Arguments retained from the historical supervised interface. They are
    accepted only to preserve the former function signature and are not used
    for estimation in the current release.
  }
}
\details{
In the original \pkg{SelvarMix} package, \code{SelvarLearnLasso()} combined
penalized variable ranking with SRUW role assignment for model-based
discriminant analysis. The current development of \pkg{SelvarMixMNAR} is
restricted to unsupervised model-based clustering. The function is retained
so that existing code fails explicitly and informatively rather than through
an unresolved function lookup.

A future release may introduce a supervised extension, but the present
function does not fit a model, rank variables, or evaluate test data.
}
\value{
This function does not return a fitted object. It signals a condition
inheriting from \code{selvarmix_supervised_deprecated} and
\code{selvarmix_feature_unavailable}.
}
\section{Status}{
Unavailable in \pkg{SelvarMixMNAR} 0.1.x. Use
\code{\link{SelvarClustLasso}} for unsupervised variable selection in
model-based clustering.
}
\seealso{
\code{\link{SelvarClustLasso}},
\code{\link{SortvarClust}},
\code{\link{SortvarLearn}}
}
\keyword{internal}

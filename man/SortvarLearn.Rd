\name{SortvarLearn}
\alias{SortvarLearn}
\title{Supervised Variable Ranking (Currently Unavailable)}
\description{
Compatibility stub for the supervised variable-ranking interface provided by
the original \pkg{SelvarMix} package. Supervised ranking has not yet been
implemented in \pkg{SelvarMixMNAR} 0.1.x, and calling this function signals a
classed error.
}
\usage{
SortvarLearn(
  x,
  z,
  type = "lasso",
  lambda = seq(20, 100, by = 10),
  rho = seq(1, 2, length = 2),
  nbcores = min(
    2,
    detectCores(all.tests = FALSE, logical = FALSE)
  )
)
}
\arguments{
  \item{x, z, type, lambda, rho, nbcores}{
    Arguments retained from the historical supervised interface. They are
    accepted only to preserve the former function signature and are not used
    in the current release.
  }
}
\details{
In the original \pkg{SelvarMix} package, \code{SortvarLearn()} ranked
variables for model-based discriminant analysis using penalized Gaussian
models with observed class labels. The current development of
\pkg{SelvarMixMNAR} is restricted to unsupervised model-based clustering.
The function is retained so that calls written for the historical interface
fail with an explicit, classed condition.

The present function does not fit a penalized model and does not return a
variable ordering.
}
\value{
This function does not return a ranking. It signals a condition inheriting
from \code{selvarmix_supervised_deprecated} and
\code{selvarmix_feature_unavailable}.
}
\section{Status}{
Unavailable in \pkg{SelvarMixMNAR} 0.1.x. Use
\code{\link{SortvarClust}} to rank variables for unsupervised model-based
clustering.
}
\seealso{
\code{\link{SortvarClust}},
\code{\link{SelvarClustLasso}},
\code{\link{SelvarLearnLasso}}
}
\keyword{internal}

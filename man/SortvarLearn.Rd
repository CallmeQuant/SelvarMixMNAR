\name{SortvarLearn}
\alias{SortvarLearn}
\title{Unavailable supervised variable ranking}
\description{
SelvarMixMNAR 0.1.x ranks variables for unsupervised model-based clustering.
Calls to the former supervised ranking interface raise a classed condition.
}
\usage{
SortvarLearn(x, z, type = "lasso", lambda = seq(20, 100, by = 10),
             rho = seq(1, 2, length = 2),
             nbcores = min(2, detectCores(all.tests = FALSE,
                                          logical = FALSE)))
}
\arguments{
\item{x, z, type, lambda, rho, nbcores}{Former supervised arguments accepted
for call compatibility; they are not evaluated.}
}
\value{
No value. The function raises a condition inheriting from
\code{selvarmix_supervised_deprecated} and
\code{selvarmix_feature_unavailable}.
}
\section{Lifecycle}{
Deprecated and unavailable. Use \code{SortvarClust()} for unsupervised
variable ranking.
}
\seealso{\code{\link{SortvarClust}}, \code{\link{SelvarClustLasso}}}
\keyword{internal}

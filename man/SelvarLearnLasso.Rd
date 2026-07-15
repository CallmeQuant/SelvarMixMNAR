\name{SelvarLearnLasso}
\alias{SelvarLearnLasso}
\title{Unavailable supervised SelvarMix estimator}
\description{
SelvarMixMNAR 0.1.x fits unsupervised clustering models. Calls to the former
supervised discriminant-analysis/lasso interface raise a classed condition.
}
\usage{
SelvarLearnLasso(x, z, lambda = seq(20, 100, by = 10),
                 rho = seq(1, 2, length = 2), type = "lasso", rank,
                 hsize = 3,
                 models = mixmodGaussianModel(listModels = c(
                   "Gaussian_pk_L_C", "Gaussian_pk_Lk_C",
                   "Gaussian_pk_L_Ck", "Gaussian_pk_Lk_Ck")),
                 rmodel = c("LI", "LB", "LC"),
                 imodel = c("LI", "LB"), xtest = x, ztest = z,
                 nbcores = min(2, detectCores(all.tests = FALSE,
                                              logical = FALSE)))
}
\arguments{
\item{x, z, lambda, rho, type, rank, hsize, models, rmodel, imodel, xtest,
ztest, nbcores}{Former supervised arguments accepted for call compatibility;
they are not evaluated.}
}
\value{
No value. The function raises a condition inheriting from
\code{selvarmix_supervised_deprecated} and
\code{selvarmix_feature_unavailable}.
}
\section{Lifecycle}{
Deprecated and unavailable. Use \code{SelvarClustLasso()} for unsupervised
model-based clustering.
}
\seealso{\code{\link{SelvarClustLasso}}, \code{\link{SortvarClust}}}
\keyword{internal}

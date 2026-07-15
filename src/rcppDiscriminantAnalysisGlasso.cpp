#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::IntegerVector rcppDiscriminantAnalysisGlasso(
    Rcpp::NumericMatrix X_,
    Rcpp::IntegerVector labels_,
    const int nbClust,
    double l,
    double r) {
  // Preserve the registered native symbol while the supervised estimator is
  // unavailable. No supervised optimization code is entered or linked.
  (void)X_;
  (void)labels_;
  (void)nbClust;
  (void)l;
  (void)r;
  Rcpp::stop(
      "selvarmix_supervised_lasso_unavailable: the supervised discriminant "
      "kernel is disabled because its precision-penalty objective and "
      "convergence properties have not been established");
  return Rcpp::IntegerVector();
}

#include "vect.hpp"

Vect::Vect() {}

Vect::Vect(NumericMatrix Data, vector<int> experiments)
{
  this->Data = Data;
  this->experiments = experiments;
}

Vect::Vect(NumericMatrix Data)
{
  this->Data = Data;
  this->initExperiments();
}

void Vect::initExperiments()
{
  for (int i = 1; i <= Data.ncol(); ++i)
    experiments.push_back(i);
}

mat Vect::const_matrix(vector<int> variables)
{
  mat source(Data.begin(), Data.nrow(), Data.ncol(), false);
  mat result = zeros<mat>(Data.nrow(), variables.size());
  for (int i = 0; i < static_cast<int>(variables.size()); ++i) {
    if (variables[i] < 1 || variables[i] > Data.ncol())
      stop("Variable index %d is outside 1,...,%d.",
           variables[i], Data.ncol());
    result.col(i) = source.col(variables[i] - 1);
  }
  if (!result.is_finite())
    stop("SRUW regression data contain non-finite values.");
  return result;
}

// Maximization-form Schwarz criterion: 2 log L_hat - q log(n).
List Vect::bicReggen(vector<int> responseVariables,
                     vector<int> predictorVariables,
                     int covarianceModel)
{
  if (responseVariables.empty())
    stop("The SRUW regression response set must not be empty.");
  if (covarianceModel < 1 || covarianceModel > 3)
    stop("Unknown SRUW regression covariance model code %d.",
         covarianceModel);

  mat response = Vect::const_matrix(responseVariables);
  const int n = response.n_rows;
  const int v = response.n_cols;
  const int a = predictorVariables.size();
  if (n < 2)
    stop("At least two observations are required for SRUW regression.");

  mat design(n, a + 1, fill::ones);
  if (!predictorVariables.empty())
    design.cols(1, a) = Vect::const_matrix(predictorVariables);

  if (arma::rank(design) < design.n_cols)
    stop("The SRUW regression design is rank deficient.");

  mat coefficients;
  const bool solved = arma::solve(coefficients, design, response);
  if (!solved || !coefficients.is_finite())
    stop("The SRUW least-squares solve failed.");

  const mat residual = response - design * coefficients;
  if (!residual.is_finite())
    stop("The SRUW regression residuals are non-finite.");

  double criterion = NA_REAL;
  double parameterCount = NA_REAL;
  const double logN = std::log(static_cast<double>(n));

  if (covarianceModel == 1) { // LI: common spherical variance
    const double sigma2 = accu(residual % residual) /
      static_cast<double>(n * v);
    if (!std::isfinite(sigma2) || sigma2 <= 0.0)
      stop("The LI residual variance is non-positive or non-finite.");
    parameterCount = v * (a + 1) + 1;
    criterion = -n * v * (std::log(2.0 * M_PI * sigma2) + 1.0) -
      parameterCount * logN;
  } else if (covarianceModel == 2) { // LB: diagonal covariance
    const rowvec sigma2 = sum(residual % residual, 0) /
      static_cast<double>(n);
    if (!sigma2.is_finite() || any(sigma2 <= 0.0))
      stop("An LB residual variance is non-positive or non-finite.");
    parameterCount = v * (a + 1) + v;
    criterion = -n * v * (std::log(2.0 * M_PI) + 1.0) -
      n * accu(log(sigma2)) - parameterCount * logN;
  } else { // LC: full covariance
    const mat sigma = residual.t() * residual /
      static_cast<double>(n);
    double logDetSigma = NA_REAL;
    if (!log_det_sympd(logDetSigma, sigma) ||
        !std::isfinite(logDetSigma))
      stop("The LC residual covariance is not positive definite.");
    parameterCount = v * (a + 1) + 0.5 * v * (v + 1);
    criterion = -n *
      (v * (std::log(2.0 * M_PI) + 1.0) + logDetSigma) -
      parameterCount * logN;
  }

  if (!std::isfinite(criterion))
    stop("The SRUW regression criterion is non-finite.");

  return List::create(
    Named("bicvalue") = criterion,
    Named("B") = coefficients,
    Named("parameterCount") = parameterCount
  );
}

vector<int> Vect::enlever_var(vector<int>& variables,
                              vector<int>& variablesToRemove)
{
  vector<int> result = variables;
  for (int value : variablesToRemove)
    result.erase(std::remove(result.begin(), result.end(), value), result.end());
  return result;
}

vector<int> Vect::ajouter_var(vector<int>& variables,
                              vector<int>& variablesToAdd)
{
  vector<int> result = variables;
  result.insert(result.end(), variablesToAdd.begin(), variablesToAdd.end());
  sort(result.begin(), result.end());
  result.erase(unique(result.begin(), result.end()), result.end());
  return result;
}

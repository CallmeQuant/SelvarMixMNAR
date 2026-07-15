#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <map>
#include <string>
#include <vector>

using namespace Rcpp;
using namespace arma;

namespace {

struct MissingPattern {
  std::vector<int> rows;
  std::vector<int> observed;
  std::vector<int> missing;
};

arma::uvec as_uvec(const std::vector<int>& indices) {
  arma::uvec result(indices.size());
  for(size_t j = 0; j < indices.size(); ++j) {
    result[j] = static_cast<arma::uword>(indices[j]);
  }
  return result;
}

std::vector<MissingPattern> build_missing_patterns(const NumericMatrix& YNA) {
  const int n = YNA.nrow();
  const int d = YNA.ncol();
  if(n < 1 || d < 1) {
    Rcpp::stop("YNA must have at least one row and one column");
  }

  std::vector<int> observed_per_variable(d, 0);
  std::map<std::string, size_t> lookup;
  std::vector<MissingPattern> patterns;
  for(int i = 0; i < n; ++i) {
    std::string key;
    key.reserve(d);
    for(int j = 0; j < d; ++j) {
      const double value = YNA(i, j);
      if(R_IsNA(value)) {
        key.push_back('1');
      } else {
        if(R_IsNaN(value) || !std::isfinite(value)) {
          Rcpp::stop("YNA may contain NA values but not NaN or infinite values");
        }
        key.push_back('0');
        ++observed_per_variable[j];
      }
    }

    std::map<std::string, size_t>::const_iterator found = lookup.find(key);
    size_t pattern_index;
    if(found == lookup.end()) {
      pattern_index = patterns.size();
      lookup[key] = pattern_index;
      MissingPattern pattern;
      for(int j = 0; j < d; ++j) {
        if(key[j] == '0') pattern.observed.push_back(j);
        else pattern.missing.push_back(j);
      }
      patterns.push_back(pattern);
    } else {
      pattern_index = found->second;
    }
    patterns[pattern_index].rows.push_back(i);
  }

  for(int j = 0; j < d; ++j) {
    if(observed_per_variable[j] == 0) {
      Rcpp::stop(
        "Variable %d is entirely missing; its Gaussian parameters are not identifiable",
        j + 1
      );
    }
  }
  return patterns;
}

void validate_observed_state_parameters(
    int d,
    const List& mu,
    const List& sigma,
    const NumericMatrix& alpha,
    const NumericVector& prop_pi) {
  const int K = prop_pi.size();
  if(K < 1 || mu.size() != K || sigma.size() != K ||
     alpha.nrow() != K || alpha.ncol() != d) {
    Rcpp::stop("Observed-state parameter dimensions are inconsistent");
  }

  double proportion_sum = 0.0;
  for(int k = 0; k < K; ++k) {
    if(!std::isfinite(prop_pi[k]) || prop_pi[k] <= 0.0) {
      Rcpp::stop("Mixing proportions must be finite and strictly positive");
    }
    proportion_sum += prop_pi[k];
    const arma::vec mean = as<arma::vec>(mu[k]);
    const arma::mat covariance = as<arma::mat>(sigma[k]);
    if(mean.n_elem != static_cast<unsigned int>(d) ||
       covariance.n_rows != static_cast<unsigned int>(d) ||
       covariance.n_cols != static_cast<unsigned int>(d) ||
       !mean.is_finite() || !covariance.is_finite()) {
      Rcpp::stop("Mean/covariance dimensions or values are invalid");
    }
    for(int j = 0; j < d; ++j) {
      if(!std::isfinite(alpha(k, j))) {
        Rcpp::stop("Missingness intercepts must be finite");
      }
      if(j > 0) {
        const double alpha_scale = 1.0 + std::max(
          std::abs(alpha(k, 0)), std::abs(alpha(k, j))
        );
        if(std::abs(alpha(k, j) - alpha(k, 0)) >
           64.0 * std::numeric_limits<double>::epsilon() * alpha_scale) {
          Rcpp::stop(
            "Pure MNARz likelihood requires one shared missingness intercept per component"
          );
        }
      }
    }
  }
  const double sum_tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
    std::max(1, K);
  if(!std::isfinite(proportion_sum) ||
     std::abs(proportion_sum - 1.0) > sum_tolerance) {
    Rcpp::stop("Mixing proportions must sum to one");
  }
}

List loglikelihood_observed_impl(
    const NumericMatrix& YNA,
    const List& mu,
    const List& sigma,
    const NumericMatrix& alpha,
    const NumericVector& prop_pi,
    const std::vector<MissingPattern>& patterns) {
  const int n = YNA.nrow();
  const int d = YNA.ncol();
  const int K = prop_pi.size();
  validate_observed_state_parameters(d, mu, sigma, alpha, prop_pi);

  NumericMatrix log_component(n, K);
  const double log_two_pi = std::log(2.0 * std::acos(-1.0));
  for(size_t pp = 0; pp < patterns.size(); ++pp) {
    const MissingPattern& pattern = patterns[pp];
    const arma::uvec observed_index = as_uvec(pattern.observed);

    for(int k = 0; k < K; ++k) {
      const arma::vec mean = as<arma::vec>(mu[k]);
      const arma::mat covariance = as<arma::mat>(sigma[k]);
      arma::vec observed_mean;
      arma::mat chol_lower;
      double log_determinant = 0.0;
      if(!pattern.observed.empty()) {
        observed_mean = mean.elem(observed_index);
        const arma::mat observed_covariance =
          covariance.submat(observed_index, observed_index);
        if(!arma::chol(chol_lower, observed_covariance, "lower")) {
          Rcpp::stop("Observed covariance submatrix is not positive definite");
        }
        log_determinant = 2.0 * arma::sum(arma::log(chol_lower.diag()));
        if(!std::isfinite(log_determinant)) {
          Rcpp::stop("Observed covariance log-determinant is not finite");
        }
      }

      double log_mask = 0.0;
      for(int j = 0; j < d; ++j) {
        if(R_IsNA(YNA(pattern.rows[0], j))) {
          log_mask += R::pnorm(alpha(k, j), 0.0, 1.0, 1, 1);
        } else {
          log_mask += R::pnorm(alpha(k, j), 0.0, 1.0, 0, 1);
        }
      }

      for(size_t rr = 0; rr < pattern.rows.size(); ++rr) {
        const int i = pattern.rows[rr];
        double log_density = 0.0;
        if(!pattern.observed.empty()) {
          arma::vec observed_value(pattern.observed.size());
          for(size_t jj = 0; jj < pattern.observed.size(); ++jj) {
            observed_value[jj] = YNA(i, pattern.observed[jj]);
          }
          const arma::vec difference = observed_value - observed_mean;
          const arma::vec standardized = arma::solve(
            arma::trimatl(chol_lower), difference
          );
          if(!standardized.is_finite()) {
            Rcpp::stop("Observed covariance triangular solve failed");
          }
          const double quadratic = arma::dot(standardized, standardized);
          log_density = -0.5 * (
            pattern.observed.size() * log_two_pi + log_determinant + quadratic
          );
        }
        log_component(i, k) = std::log(prop_pi[k]) + log_density + log_mask;
      }
    }
  }

  NumericMatrix tik(n, K);
  double loglikelihood = 0.0;
  for(int i = 0; i < n; ++i) {
    double maximum = R_NegInf;
    for(int k = 0; k < K; ++k) maximum = std::max(maximum, log_component(i, k));
    if(!std::isfinite(maximum)) {
      Rcpp::stop("All component log-probabilities are non-finite for an observation");
    }
    double normalizer = 0.0;
    for(int k = 0; k < K; ++k) {
      tik(i, k) = std::exp(log_component(i, k) - maximum);
      normalizer += tik(i, k);
    }
    if(!std::isfinite(normalizer) || normalizer <= 0.0) {
      Rcpp::stop("Posterior normalizer is non-finite or non-positive");
    }
    loglikelihood += maximum + std::log(normalizer);
    for(int k = 0; k < K; ++k) tik(i, k) /= normalizer;
  }
  if(!std::isfinite(loglikelihood)) {
    Rcpp::stop("loglik_obs became NaN or infinite; numerical degeneracy detected");
  }

  return List::create(
    Named("loglik_obs") = loglikelihood,
    Named("tik") = tik
  );
}

} // namespace

// Shared by the pure and mixed native entry points.  A diagonal Gaussian
// coordinate needs observed variation; a full covariance additionally needs
// direct joint observations for every covariance entry.
void validate_gaussian_observation_design(
    const NumericMatrix& YNA,
    bool diagonal) {
  const int n = YNA.nrow();
  const int d = YNA.ncol();
  for(int j = 0; j < d; ++j) {
    int count = 0;
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for(int i = 0; i < n; ++i) {
      if(R_IsNA(YNA(i, j))) continue;
      ++count;
      minimum = std::min(minimum, YNA(i, j));
      maximum = std::max(maximum, YNA(i, j));
    }
    const double scale = std::max(
      std::numeric_limits<double>::min(),
      std::max(std::abs(minimum), std::abs(maximum))
    );
    if(count < 2 || !std::isfinite(minimum) || !std::isfinite(maximum) ||
       maximum - minimum <= std::sqrt(std::numeric_limits<double>::epsilon()) * scale) {
      Rcpp::stop(
        "Gaussian variable %d requires at least two distinct observed values",
        j + 1
      );
    }
  }

  if(!diagonal) {
    for(int first = 0; first < d; ++first) {
      for(int second = first + 1; second < d; ++second) {
        int jointly_observed = 0;
        for(int i = 0; i < n; ++i) {
          if(!R_IsNA(YNA(i, first)) && !R_IsNA(YNA(i, second))) {
            ++jointly_observed;
          }
        }
        if(jointly_observed < 2) {
          Rcpp::stop(
            "Full Gaussian covariance pair (%d,%d) has fewer than two jointly observed rows",
            first + 1, second + 1
          );
        }
      }
    }
  }
}

// Observed log likelihood and posterior component probabilities.
// [[Rcpp::export]]
List LoglikelihoodObsGaussian(NumericMatrix YNA, List mu, List sigma, NumericMatrix alpha, NumericVector prop_pi) {
  const std::vector<MissingPattern> patterns = build_missing_patterns(YNA);
  return loglikelihood_observed_impl(YNA, mu, sigma, alpha, prop_pi, patterns);
}

// Construct a finite Gaussian and class-only MNARz initial state.
// [[Rcpp::export]]
List InitEMGaussian(NumericMatrix YNA, int K, std::string mecha, bool diag, Nullable<List> init, Nullable<int> samplesize) {
  int n = YNA.nrow();
  int d = YNA.ncol();
  (void)build_missing_patterns(YNA);
  if(K < 1 || K > n) {
    Rcpp::stop("K must be between one and nrow(YNA)");
  }
  if(mecha != "MNARz") {
    Rcpp::stop("InitEMGaussian supports only the class-only 'MNARz' mechanism");
  }
  
  // Seed-controlled balanced allocation when no state is supplied.
  if(init.isNull()) {
    NumericMatrix Z_init(n, K);
    
    // A randomly permuted balanced allocation prevents an empty component at
    // the initial iterate while retaining seed-controlled stochastic starts.
    IntegerVector balanced_labels(n);
    for(int i = 0; i < n; ++i) balanced_labels[i] = i % K;
    IntegerVector assignments = Rcpp::sample(balanced_labels, n, false);
    
    // Convert the balanced allocation to responsibilities.
    for(int i = 0; i < n; i++) {
      if(assignments[i] < 0 || assignments[i] >= K) {
        Rcpp::stop("Internal balanced initializer produced an invalid class label");
      }
      Z_init(i, assignments[i]) = 1.0;
    }
    
    NumericVector prop_pi(K);
    for(int k = 0; k < K; k++) {
      prop_pi[k] = sum(Z_init(_, k)) / n;
    }
    
    // Global observed-data moments provide scale-adaptive substitutes when
    // componentwise moments are unavailable for a coordinate or pair.
    NumericVector global_mean(d);
    NumericMatrix global_covariance(d, d);
    for(int j = 0; j < d; ++j) {
      double total = 0.0;
      int count = 0;
      for(int i = 0; i < n; ++i) {
        if(!R_IsNA(YNA(i, j))) {
          total += YNA(i, j);
          ++count;
        }
      }
      global_mean[j] = total / count;
    }
    for(int j1 = 0; j1 < d; ++j1) {
      for(int j2 = 0; j2 <= j1; ++j2) {
        double total = 0.0;
        int count = 0;
        for(int i = 0; i < n; ++i) {
          if(!R_IsNA(YNA(i, j1)) && !R_IsNA(YNA(i, j2))) {
            total += (YNA(i, j1) - global_mean[j1]) *
              (YNA(i, j2) - global_mean[j2]);
            ++count;
          }
        }
        double value = count > 1 ? total / count :
          (j1 == j2 ? std::numeric_limits<double>::epsilon() : 0.0);
        global_covariance(j1, j2) = value;
        global_covariance(j2, j1) = value;
      }
    }

    // Compute componentwise observed-data moments.
    List mu_init(K);
    List sigma_init(K);
    
    for(int k = 0; k < K; k++) {
      IntegerVector class_obs;
      for(int i = 0; i < n; i++) {
        if(Z_init(i, k) == 1.0) class_obs.push_back(i);
      }
      
      NumericVector mu_k(d);
      for(int j = 0; j < d; j++) {
        double sum = 0.0;
        int count = 0;
        for(int idx = 0; idx < class_obs.size(); idx++) {
          int i = class_obs[idx];
          if(!R_IsNA(YNA(i, j))) {
            sum += YNA(i, j);
            count++;
          }
        }
        mu_k[j] = count > 0 ? sum / count : global_mean[j];
      }
      mu_init[k] = mu_k;
      
      NumericMatrix sigma_k(d, d);
      if(diag) {
        for(int j = 0; j < d; j++) {
          double sum_sq_diff = 0.0;
          int count = 0;
          for(int idx = 0; idx < class_obs.size(); idx++) {
            int i = class_obs[idx];
            if(!R_IsNA(YNA(i, j))) {
              double diff = YNA(i, j) - mu_k[j];
              sum_sq_diff += diff * diff;
              count++;
            }
          }
          sigma_k(j, j) = count > 1 ? sum_sq_diff / count :
            global_covariance(j, j);
        }
      } else {
        for(int j1 = 0; j1 < d; j1++) {
          for(int j2 = 0; j2 <= j1; j2++) {
            double sum_cross = 0.0;
            int count = 0;
            for(int idx = 0; idx < class_obs.size(); idx++) {
              int i = class_obs[idx];
              if(!R_IsNA(YNA(i, j1)) && !R_IsNA(YNA(i, j2))) {
                double diff1 = YNA(i, j1) - mu_k[j1];
                double diff2 = YNA(i, j2) - mu_k[j2];
                sum_cross += diff1 * diff2;
                count++;
              }
            }
            double cov_val = count > 1 ? sum_cross / count :
              global_covariance(j1, j2);
            sigma_k(j1, j2) = cov_val;
            sigma_k(j2, j1) = cov_val;
          }
        }
      }
      sigma_init[k] = sigma_k;
    }
    
    // Class-aware pure-MNARz starts pool mask entries over variables.  Jeffreys
    // smoothing is initialization-only and prevents boundary intercepts.
    NumericMatrix alpha_init(K, d);
    for(int k = 0; k < K; ++k) {
      double missing_count = 0.0;
      double component_count = 0.0;
      for(int i = 0; i < n; ++i) {
        if(Z_init(i, k) == 0.0) continue;
        ++component_count;
        for(int j = 0; j < d; ++j) {
          if(R_IsNA(YNA(i, j))) ++missing_count;
        }
      }
      const double probability = (missing_count + 0.5) /
        (d * component_count + 1.0);
      const double alpha_k = R::qnorm(probability, 0.0, 1.0, 1, 0);
      for(int j = 0; j < d; ++j) alpha_init(k, j) = alpha_k;
    }
    return List::create(
      Named("pi_init") = prop_pi,
      Named("mu_init") = mu_init,
      Named("sigma_init") = sigma_init,
      Named("alpha_init") = alpha_init
    );
  } else {
    List init_list = as<List>(init);
    return List::create(
      Named("pi_init") = as<NumericVector>(init_list["pik"]),
      Named("mu_init") = as<List>(init_list["mu"]),
      Named("sigma_init") = as<List>(init_list["sigma"]),
      Named("alpha_init") = as<NumericMatrix>(init_list["alpha"])
    );
  }
}

// [[Rcpp::export]]
NumericMatrix MechanismEMGLM(NumericMatrix YNA, NumericMatrix tik, std::string mecha) {
  if(mecha != "MNARz") {
    Rcpp::stop(
      "MechanismEMGLM supports only the identifiable class-only 'MNARz' mechanism"
    );
  }

  int n = YNA.nrow();
  int d = YNA.ncol();
  int K = tik.ncol();
  if(n < 1 || d < 1 || tik.nrow() != n || K < 1) {
    Rcpp::stop("Missingness M-step dimensions are inconsistent");
  }

  for(int i = 0; i < n; ++i) {
    double row_sum = 0.0;
    for(int k = 0; k < K; ++k) {
      const double weight = tik(i, k);
      if(!std::isfinite(weight) || weight < 0.0) {
        Rcpp::stop("Responsibilities must be finite and nonnegative");
      }
      row_sum += weight;
    }
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
      std::max(1, K);
    if(!std::isfinite(row_sum) || std::abs(row_sum - 1.0) > tolerance) {
      Rcpp::stop("Each responsibility row must sum to one");
    }
  }

  NumericMatrix alpha_new(K, d);

  // Pure MNARz has one class-specific probability shared by all variables.
  // The weighted intercept-only probit MLE is available in closed form.
  for(int k = 0; k < K; ++k) {
    double effective_size = 0.0;
    double weighted_missing = 0.0;
    for(int i = 0; i < n; ++i) {
      const double weight = tik(i, k);
      effective_size += weight;
      for(int j = 0; j < d; ++j) {
        if(R_IsNA(YNA(i, j))) {
          weighted_missing += weight;
        }
      }
    }
    const double probability = weighted_missing / (d * effective_size);
    if(!std::isfinite(probability) || probability <= 0.0 || probability >= 1.0) {
      Rcpp::stop(
        "Missingness M-step has no finite pure-MNARz intercept for component %d",
        k + 1
      );
    }
    const double alpha_k = R::qnorm(probability, 0.0, 1.0, 1, 0);
    for(int j = 0; j < d; ++j) {
      alpha_new(k, j) = alpha_k;
    }
  }

  return alpha_new;
}

// Stabilize a covariance matrix at a relative machine-scale eigenvalue floor.
// The returned flag records an estimand-changing numerical adjustment so the
// R wrapper can expose it in diagnostics.
arma::mat stabilize_covariance(const arma::mat& input, bool diagonal, bool& adjusted) {
  if(input.n_rows == 0 || input.n_rows != input.n_cols || !input.is_finite()) {
    Rcpp::stop("Covariance matrix must be finite, square, and non-empty");
  }

  arma::mat stabilized = 0.5 * (input + input.t());
  adjusted = arma::norm(stabilized - input, "inf") > 0.0;
  if(diagonal) {
    arma::mat diagonalized = arma::diagmat(stabilized.diag());
    adjusted = adjusted || arma::norm(diagonalized - stabilized, "inf") > 0.0;
    stabilized = diagonalized;
  }

  arma::vec eigenvalues;
  arma::mat eigenvectors;
  if(!arma::eig_sym(eigenvalues, eigenvectors, stabilized)) {
    Rcpp::stop("Covariance eigendecomposition failed");
  }
  double scale = eigenvalues.n_elem > 0 ? arma::abs(eigenvalues).max() : 0.0;
  double epsilon = std::numeric_limits<double>::epsilon();
  if(!std::isfinite(scale) || scale <= std::numeric_limits<double>::min()) {
    Rcpp::stop("Covariance matrix has zero numerical scale");
  }
  const double negative_tolerance = 64.0 * epsilon * scale;
  if(eigenvalues.min() < -negative_tolerance) {
    Rcpp::stop("Covariance matrix is materially indefinite");
  }
  const double eigen_floor = std::max(
    std::sqrt(epsilon) * scale,
    std::numeric_limits<double>::min()
  );

  bool requires_floor = false;
  for(unsigned int j = 0; j < eigenvalues.n_elem; ++j) {
    if(eigenvalues[j] < eigen_floor) {
      eigenvalues[j] = eigen_floor;
      requires_floor = true;
    }
  }
  if(requires_floor) {
    stabilized = eigenvectors * arma::diagmat(eigenvalues) * eigenvectors.t();
    if(diagonal) stabilized = arma::diagmat(stabilized.diag());
    adjusted = true;
  }
  return 0.5 * (stabilized + stabilized.t());
}

double validate_component_masses(const NumericMatrix& tik, double component_floor) {
  int K = tik.ncol();
  double minimum_mass = std::numeric_limits<double>::infinity();
  for(int k = 0; k < K; ++k) {
    double mass = sum(tik(_, k));
    if(!std::isfinite(mass) || mass <= component_floor) {
      Rcpp::stop(
        "Effective component size for component %d is %.17g; it must exceed %.17g",
        k + 1, mass, component_floor
      );
    }
    minimum_mass = std::min(minimum_mass, mass);
  }
  return minimum_mass;
}

namespace {

void gaussian_m_step_patterned(
    const NumericMatrix& YNA,
    const NumericMatrix& tik,
    const List& current_mu,
    const List& current_sigma,
    const std::vector<MissingPattern>& patterns,
    bool diagonal,
    List& updated_mu,
    List& updated_sigma,
    int& covariance_adjustments) {
  const int d = YNA.ncol();
  const int K = tik.ncol();

  for(int k = 0; k < K; ++k) {
    const double effective_size = sum(tik(_, k));
    const arma::vec old_mean = as<arma::vec>(current_mu[k]);
    const arma::mat old_covariance = as<arma::mat>(current_sigma[k]);
    arma::vec first_moment(d, arma::fill::zeros);
    arma::vec second_diagonal(d, arma::fill::zeros);
    arma::mat second_moment;
    if(!diagonal) second_moment.zeros(d, d);

    for(size_t pp = 0; pp < patterns.size(); ++pp) {
      const MissingPattern& pattern = patterns[pp];
      const arma::uvec observed_index = as_uvec(pattern.observed);
      const arma::uvec missing_index = as_uvec(pattern.missing);
      arma::mat gain;
      arma::mat conditional_covariance(d, d, arma::fill::zeros);

      if(!pattern.missing.empty() && pattern.observed.empty()) {
        conditional_covariance = old_covariance;
      } else if(!pattern.missing.empty()) {
        const arma::mat sigma_oo =
          old_covariance.submat(observed_index, observed_index);
        const arma::mat sigma_mo =
          old_covariance.submat(missing_index, observed_index);
        const arma::mat sigma_mm =
          old_covariance.submat(missing_index, missing_index);
        arma::mat solved;
        if(!arma::solve(
             solved, sigma_oo, sigma_mo.t(), arma::solve_opts::likely_sympd
           )) {
          Rcpp::stop("Conditional-moment covariance solve failed");
        }
        gain = solved.t();
        arma::mat conditional_mm = sigma_mm - gain * sigma_mo.t();
        conditional_mm = 0.5 * (conditional_mm + conditional_mm.t());
        if(!gain.is_finite() || !conditional_mm.is_finite()) {
          Rcpp::stop("Conditional Gaussian moments are non-finite");
        }
        conditional_covariance.submat(missing_index, missing_index) = conditional_mm;
      }

      for(size_t rr = 0; rr < pattern.rows.size(); ++rr) {
        const int i = pattern.rows[rr];
        const double weight = tik(i, k);
        arma::vec conditional_mean = old_mean;
        arma::vec observed_value(pattern.observed.size());
        for(size_t jj = 0; jj < pattern.observed.size(); ++jj) {
          observed_value[jj] = YNA(i, pattern.observed[jj]);
          conditional_mean[pattern.observed[jj]] = observed_value[jj];
        }
        if(!pattern.missing.empty() && !pattern.observed.empty()) {
          conditional_mean.elem(missing_index) = old_mean.elem(missing_index) +
            gain * (observed_value - old_mean.elem(observed_index));
        }
        if(!conditional_mean.is_finite()) {
          Rcpp::stop("Conditional Gaussian means are non-finite");
        }

        first_moment += weight * conditional_mean;
        if(diagonal) {
          second_diagonal += weight * (
            conditional_covariance.diag() + arma::square(conditional_mean)
          );
        } else {
          second_moment += weight * (
            conditional_covariance + conditional_mean * conditional_mean.t()
          );
        }
      }
    }

    const arma::vec new_mean = first_moment / effective_size;
    arma::mat new_covariance;
    if(diagonal) {
      const arma::vec new_variance = second_diagonal / effective_size -
        arma::square(new_mean);
      new_covariance = arma::diagmat(new_variance);
    } else {
      new_covariance = second_moment / effective_size - new_mean * new_mean.t();
    }
    new_covariance = 0.5 * (new_covariance + new_covariance.t());
    bool adjusted = false;
    new_covariance = stabilize_covariance(new_covariance, diagonal, adjusted);
    if(adjusted) ++covariance_adjustments;
    NumericVector new_mean_r(d);
    for(int j = 0; j < d; ++j) new_mean_r[j] = new_mean[j];
    updated_mu[k] = new_mean_r;
    updated_sigma[k] = wrap(new_covariance);
  }
}

NumericMatrix impute_from_returned_state(
    const NumericMatrix& YNA,
    const NumericMatrix& tik,
    const List& mu,
    const List& sigma,
    const std::vector<MissingPattern>& patterns) {
  const int K = tik.ncol();
  NumericMatrix imputed = clone(YNA);

  for(size_t pp = 0; pp < patterns.size(); ++pp) {
    const MissingPattern& pattern = patterns[pp];
    if(pattern.missing.empty()) continue;
    const arma::uvec observed_index = as_uvec(pattern.observed);
    const arma::uvec missing_index = as_uvec(pattern.missing);
    for(size_t rr = 0; rr < pattern.rows.size(); ++rr) {
      for(size_t jj = 0; jj < pattern.missing.size(); ++jj) {
        imputed(pattern.rows[rr], pattern.missing[jj]) = 0.0;
      }
    }

    for(int k = 0; k < K; ++k) {
      const arma::vec component_mean = as<arma::vec>(mu[k]);
      const arma::mat component_covariance = as<arma::mat>(sigma[k]);
      arma::mat gain;
      if(!pattern.observed.empty()) {
        const arma::mat sigma_oo =
          component_covariance.submat(observed_index, observed_index);
        const arma::mat sigma_mo =
          component_covariance.submat(missing_index, observed_index);
        arma::mat solved;
        if(!arma::solve(
             solved, sigma_oo, sigma_mo.t(), arma::solve_opts::likely_sympd
           )) {
          Rcpp::stop("Imputation covariance solve failed");
        }
        gain = solved.t();
      }

      for(size_t rr = 0; rr < pattern.rows.size(); ++rr) {
        const int i = pattern.rows[rr];
        arma::vec conditional_missing = component_mean.elem(missing_index);
        if(!pattern.observed.empty()) {
          arma::vec observed_value(pattern.observed.size());
          for(size_t jj = 0; jj < pattern.observed.size(); ++jj) {
            observed_value[jj] = YNA(i, pattern.observed[jj]);
          }
          conditional_missing += gain * (
            observed_value - component_mean.elem(observed_index)
          );
        }
        if(!conditional_missing.is_finite()) {
          Rcpp::stop("Conditional Gaussian imputation is non-finite");
        }
        for(size_t jj = 0; jj < pattern.missing.size(); ++jj) {
          imputed(i, pattern.missing[jj]) += tik(i, k) * conditional_missing[jj];
        }
      }
    }
  }
  return imputed;
}

} // namespace

// [[Rcpp::export]]
List EMGaussian(NumericMatrix YNA, int K, std::string mecha, bool diag, int rmax,
                Nullable<List> init = R_NilValue,
                double tol = 0.0001,
                Nullable<int> samplesize = R_NilValue) {
  const int n = YNA.nrow();
  const int d = YNA.ncol();
  const std::vector<MissingPattern> patterns = build_missing_patterns(YNA);
  validate_gaussian_observation_design(YNA, diag);

  if(K < 1 || K > n) Rcpp::stop("K must be between one and nrow(YNA)");
  if(rmax < 1) Rcpp::stop("rmax must be a positive integer");
  if(!std::isfinite(tol) || tol < 0.0) {
    Rcpp::stop("tol must be finite and non-negative");
  }
  if(mecha != "MNARz") {
    Rcpp::stop(
      "EMGaussian supports only the identifiable class-only 'MNARz' mechanism"
    );
  }

  List init_params = init.isNull() ?
    InitEMGaussian(YNA, K, mecha, diag, R_NilValue, samplesize) :
    InitEMGaussian(YNA, K, mecha, diag, init, samplesize);

  // Clone every supplied field so no native update aliases caller-owned R state.
  NumericVector pi = clone(as<NumericVector>(init_params["pi_init"]));
  List mu = clone(as<List>(init_params["mu_init"]));
  List sigma = clone(as<List>(init_params["sigma_init"]));
  NumericMatrix alpha = clone(as<NumericMatrix>(init_params["alpha_init"]));

  if(pi.size() != K || mu.size() != K || sigma.size() != K ||
     alpha.nrow() != K || alpha.ncol() != d) {
    Rcpp::stop("Initial parameter dimensions do not match K and ncol(YNA)");
  }
  double pi_sum = 0.0;
  for(int k = 0; k < K; ++k) {
    if(!std::isfinite(pi[k]) || pi[k] <= 0.0) {
      Rcpp::stop("Initial mixing proportions must be finite and strictly positive");
    }
    pi_sum += pi[k];
  }
  if(!std::isfinite(pi_sum) || pi_sum <= 0.0) {
    Rcpp::stop("Initial mixing proportions must have a positive finite sum");
  }
  for(int k = 0; k < K; ++k) pi[k] /= pi_sum;

  int covariance_adjustments = 0;
  for(int k = 0; k < K; ++k) {
    const arma::vec component_mean = as<arma::vec>(mu[k]);
    const arma::mat component_covariance = as<arma::mat>(sigma[k]);
    if(component_mean.n_elem != static_cast<unsigned int>(d) ||
       component_covariance.n_rows != static_cast<unsigned int>(d) ||
       component_covariance.n_cols != static_cast<unsigned int>(d) ||
       !component_mean.is_finite() || !component_covariance.is_finite()) {
      Rcpp::stop("Initial mean/covariance dimensions or values are invalid");
    }
    for(int j = 0; j < d; ++j) {
      if(!std::isfinite(alpha(k, j))) {
        Rcpp::stop("Initial means and missingness intercepts must be finite");
      }
      if(j > 0) {
        const double scale = 1.0 + std::max(
          std::abs(alpha(k, 0)), std::abs(alpha(k, j))
        );
        if(std::abs(alpha(k, j) - alpha(k, 0)) >
           64.0 * std::numeric_limits<double>::epsilon() * scale) {
          Rcpp::stop(
            "Pure MNARz requires one shared missingness intercept per component; row %d is not constant",
            k + 1
          );
        }
      }
    }
    bool adjusted = false;
    sigma[k] = wrap(stabilize_covariance(
      component_covariance, diag, adjusted
    ));
    if(adjusted) ++covariance_adjustments;
  }

  List observed = loglikelihood_observed_impl(
    YNA, mu, sigma, alpha, pi, patterns
  );
  double loglikelihood = as<double>(observed["loglik_obs"]);
  NumericMatrix tik = as<NumericMatrix>(observed["tik"]);
  NumericVector loglikelihood_trace = NumericVector::create(loglikelihood);

  // A complete-data full covariance requires effective mass greater than d to
  // avoid a necessarily rank-deficient empirical covariance.  The diagonal
  // model retains the scalar effective-mass threshold greater than one.
  const double component_floor = diag ? 1.0 :
    static_cast<double>(std::max(1, d));
  double minimum_component_size =
    validate_component_masses(tik, component_floor);

  int iteration = 0;
  bool converged = false;
  bool loglik_monotone = true;
  bool no_material_loglik_decrease = true;
  double final_improvement = NA_REAL;
  double final_relative_improvement = NA_REAL;
  double final_convergence_threshold = NA_REAL;
  double final_decrease_tolerance = NA_REAL;

  while(iteration < rmax && !converged) {
    ++iteration;
    NumericVector candidate_pi(K);
    for(int k = 0; k < K; ++k) {
      candidate_pi[k] = sum(tik(_, k)) / n;
    }

    List candidate_mu(K);
    List candidate_sigma(K);
    gaussian_m_step_patterned(
      YNA, tik, mu, sigma, patterns, diag,
      candidate_mu, candidate_sigma, covariance_adjustments
    );
    NumericMatrix candidate_alpha = MechanismEMGLM(YNA, tik, mecha);

    List candidate = loglikelihood_observed_impl(
      YNA, candidate_mu, candidate_sigma, candidate_alpha, candidate_pi, patterns
    );
    const double candidate_loglikelihood =
      as<double>(candidate["loglik_obs"]);
    NumericMatrix candidate_tik = as<NumericMatrix>(candidate["tik"]);
    const double candidate_minimum =
      validate_component_masses(candidate_tik, component_floor);
    minimum_component_size = std::min(
      minimum_component_size, candidate_minimum
    );

    final_improvement = candidate_loglikelihood - loglikelihood;
    const double convergence_scale = 1.0 + std::max(
      std::abs(loglikelihood), std::abs(candidate_loglikelihood)
    );
    final_relative_improvement = final_improvement / convergence_scale;
    final_convergence_threshold = tol * convergence_scale;
    final_decrease_tolerance =
      64.0 * std::numeric_limits<double>::epsilon() * convergence_scale;
    if(final_improvement < 0.0) loglik_monotone = false;
    if(final_improvement < -final_decrease_tolerance) {
      no_material_loglik_decrease = false;
      Rcpp::stop(
        "Observed log-likelihood decreased materially at iteration %d (delta=%.17g, tolerance=%.17g)",
        iteration, final_improvement, final_decrease_tolerance
      );
    }

    pi = clone(candidate_pi);
    mu = clone(candidate_mu);
    sigma = clone(candidate_sigma);
    alpha = clone(candidate_alpha);
    tik = clone(candidate_tik);
    loglikelihood = candidate_loglikelihood;
    loglikelihood_trace.push_back(loglikelihood);
    converged = std::abs(final_improvement) <= final_convergence_threshold;
  }

  NumericMatrix imputed_data = impute_from_returned_state(
    YNA, tik, mu, sigma, patterns
  );

  return List::create(
    Named("pik") = pi,
    Named("mu") = mu,
    Named("sigma") = sigma,
    Named("alpha") = alpha,
    Named("loglik_vec") = loglikelihood_trace,
    Named("tik") = tik,
    Named("imputedData") = imputed_data,
    Named("iterations") = iteration,
    Named("converged") = converged,
    Named("termination_reason") = converged ? "tolerance" : "max_iterations",
    Named("final_loglik_improvement") = final_improvement,
    Named("final_relative_loglik_improvement") = final_relative_improvement,
    Named("convergence_threshold") = final_convergence_threshold,
    Named("absolute_tolerance") = tol,
    Named("relative_tolerance") = tol,
    Named("decrease_tolerance") = final_decrease_tolerance,
    Named("loglik_monotone") = loglik_monotone,
    Named("no_material_loglik_decrease") = no_material_loglik_decrease,
    Named("min_effective_component_size") = minimum_component_size,
    Named("component_floor") = component_floor,
    Named("component_mass_rule") =
      diag ? "effective mass > 1" : "effective mass > dimension",
    Named("covariance_model_dimension") = d,
    Named("covariance_adjustments") = covariance_adjustments,
    Named("missing_pattern_count") = static_cast<int>(patterns.size()),
    Named("conditional_moment_storage") =
      "streamed sufficient statistics by missingness pattern",
    Named("sufficient_statistic_storage") =
      "streamed sufficient statistics by missingness pattern",
    Named("error") = "No error"
  );
}

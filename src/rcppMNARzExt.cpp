#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <map>
#include <string>
#include <vector>

using namespace Rcpp;
using namespace arma;

// The mixed model uses the same Gaussian initializer and covariance rule as
// pure MNARz; only the coordinates entering the missingness likelihood differ.
List InitEMGaussian(
    NumericMatrix YNA,
    int K,
    std::string mecha,
    bool diag,
    Nullable<List> init,
    Nullable<int> samplesize);
arma::mat stabilize_covariance(
    const arma::mat& input,
    bool diagonal,
    bool& adjusted);
void validate_gaussian_observation_design(
    const NumericMatrix& YNA,
    bool diagonal);

namespace {

struct MixedState {
  NumericVector pi;
  List mu;
  List sigma;
  NumericVector rho;
  NumericMatrix alpha;
  NumericMatrix beta;
};

struct MixedPattern {
  std::vector<int> rows;
  std::vector<int> observed;
  std::vector<int> missing;
};

arma::uvec mixed_uvec(const std::vector<int>& indices) {
  arma::uvec result(indices.size());
  for(size_t j = 0; j < indices.size(); ++j) {
    result[j] = static_cast<arma::uword>(indices[j]);
  }
  return result;
}

std::vector<MixedPattern> build_mixed_patterns(const NumericMatrix& YNA) {
  std::map<std::string, size_t> lookup;
  std::vector<MixedPattern> patterns;
  for(int i = 0; i < YNA.nrow(); ++i) {
    std::string key;
    key.reserve(YNA.ncol());
    for(int j = 0; j < YNA.ncol(); ++j) {
      key.push_back(R_IsNA(YNA(i, j)) ? '1' : '0');
    }
    std::map<std::string, size_t>::const_iterator found = lookup.find(key);
    size_t index;
    if(found == lookup.end()) {
      index = patterns.size();
      lookup[key] = index;
      MixedPattern pattern;
      for(int j = 0; j < YNA.ncol(); ++j) {
        if(key[j] == '0') pattern.observed.push_back(j);
        else pattern.missing.push_back(j);
      }
      patterns.push_back(pattern);
    } else {
      index = found->second;
    }
    patterns[index].rows.push_back(i);
  }
  return patterns;
}

void validate_mixed_data(const NumericMatrix& YNA) {
  const int n = YNA.nrow();
  const int d = YNA.ncol();
  if(n < 1 || d < 1) {
    Rcpp::stop("YNA must have at least one row and one column");
  }
  for(int j = 0; j < d; ++j) {
    int observed = 0;
    for(int i = 0; i < n; ++i) {
      const double value = YNA(i, j);
      if(R_IsNA(value)) continue;
      if(R_IsNaN(value) || !std::isfinite(value)) {
        Rcpp::stop("YNA may contain NA values but not NaN or infinite values");
      }
      ++observed;
    }
    if(observed == 0) {
      Rcpp::stop(
        "Variable %d is entirely missing; its Gaussian parameters are not identifiable",
        j + 1
      );
    }
  }
}

std::vector<int> mnar_indices(
    const LogicalVector& is_mnar,
    int d) {
  if(is_mnar.size() != d) {
    Rcpp::stop("Length of 'is_mnar' must equal ncol(YNA)");
  }
  std::vector<int> selected;
  for(int j = 0; j < d; ++j) {
    if(LogicalVector::is_na(is_mnar[j])) {
      Rcpp::stop("'is_mnar' must not contain NA");
    }
    if(is_mnar[j]) selected.push_back(j);
  }
  return selected;
}

void validate_zero_beta(const NumericMatrix& beta, int K, int d) {
  if(beta.nrow() != K || beta.ncol() != d) {
    Rcpp::stop("beta must have dimensions K by ncol(YNA)");
  }
  for(int k = 0; k < K; ++k) {
    for(int j = 0; j < d; ++j) {
      if(!std::isfinite(beta(k, j)) || beta(k, j) != 0.0) {
        Rcpp::stop(
          "beta is deprecated and must be zero: the mixed estimator is class-only MNARz, not self-masked MNARy"
        );
      }
    }
  }
}

NumericVector rho_from_alpha(
    const NumericMatrix& alpha,
    int K,
    int d,
    const std::vector<int>& selected) {
  if(alpha.nrow() != K || alpha.ncol() != d) {
    Rcpp::stop("alpha must have dimensions K by ncol(YNA)");
  }
  NumericVector rho(K, NA_REAL);
  std::vector<bool> is_selected(d, false);
  for(size_t jj = 0; jj < selected.size(); ++jj) {
    is_selected[selected[jj]] = true;
  }
  for(int k = 0; k < K; ++k) {
    for(int j = 0; j < d; ++j) {
      if(!is_selected[j] && !R_IsNA(alpha(k, j))) {
        Rcpp::stop(
          "alpha must be NA on coordinates not selected as MNARz"
        );
      }
    }
  }
  if(selected.empty()) return rho;

  for(int k = 0; k < K; ++k) {
    const double reference = alpha(k, selected[0]);
    if(!std::isfinite(reference)) {
      Rcpp::stop("Selected MNARz intercepts must be finite");
    }
    for(size_t jj = 1; jj < selected.size(); ++jj) {
      const double candidate = alpha(k, selected[jj]);
      if(!std::isfinite(candidate)) {
        Rcpp::stop("Selected MNARz intercepts must be finite");
      }
      const double scale = 1.0 + std::max(std::abs(reference), std::abs(candidate));
      if(std::abs(candidate - reference) >
         64.0 * std::numeric_limits<double>::epsilon() * scale) {
        Rcpp::stop(
          "Mixed class-only MNARz requires one shared intercept over selected coordinates for component %d",
          k + 1
        );
      }
    }
    rho[k] = R::pnorm(reference, 0.0, 1.0, 1, 0);
    if(!std::isfinite(rho[k]) || rho[k] <= 0.0 || rho[k] >= 1.0) {
      Rcpp::stop(
        "Initial mixed MNARz probability for component %d must be strictly between zero and one",
        k + 1
      );
    }
  }
  return rho;
}

NumericMatrix alpha_from_rho(
    const NumericVector& rho,
    int K,
    int d,
    const std::vector<int>& selected) {
  NumericMatrix alpha(K, d);
  std::fill(alpha.begin(), alpha.end(), NA_REAL);
  if(selected.empty()) return alpha;
  if(rho.size() != K) Rcpp::stop("rho must have length K");

  for(int k = 0; k < K; ++k) {
    if(!std::isfinite(rho[k]) || rho[k] <= 0.0 || rho[k] >= 1.0) {
      Rcpp::stop("rho must be finite and strictly between zero and one");
    }
    const double alpha_k = R::qnorm(rho[k], 0.0, 1.0, 1, 0);
    for(size_t jj = 0; jj < selected.size(); ++jj) {
      alpha(k, selected[jj]) = alpha_k;
    }
  }
  return alpha;
}

double validate_component_masses_mixed(
    const NumericMatrix& tik,
    double component_floor) {
  const int n = tik.nrow();
  const int K = tik.ncol();
  if(n < 1 || K < 1) {
    Rcpp::stop("tik must have at least one row and one component");
  }
  for(int i = 0; i < n; ++i) {
    double row_sum = 0.0;
    for(int k = 0; k < K; ++k) {
      const double value = tik(i, k);
      if(!std::isfinite(value) || value < 0.0 || value > 1.0) {
        Rcpp::stop("tik entries must be finite probabilities in [0, 1]");
      }
      row_sum += value;
    }
    if(!std::isfinite(row_sum) || std::abs(row_sum - 1.0) > 1e-8) {
      Rcpp::stop("Each row of tik must sum to one");
    }
  }
  double minimum_mass = std::numeric_limits<double>::infinity();
  for(int k = 0; k < K; ++k) {
    const double mass = sum(tik(_, k));
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

List mixed_observed_state(
    const NumericMatrix& YNA,
    const List& mu,
    const List& sigma,
    const NumericVector& rho,
    const NumericVector& prop_pi,
    const std::vector<int>& selected,
    const std::vector<MixedPattern>& patterns) {
  const int n = YNA.nrow();
  const int d = YNA.ncol();
  const int K = prop_pi.size();
  const int mechanism_dimension = static_cast<int>(selected.size());

  if(mu.size() != K || sigma.size() != K) {
    Rcpp::stop("mu and sigma must contain K components");
  }
  if(mechanism_dimension > 0 && rho.size() != K) {
    Rcpp::stop("rho must have length K when MNARz coordinates are selected");
  }

  double mixing_sum = 0.0;
  std::vector<arma::vec> means(K);
  std::vector<arma::mat> covariances(K);
  for(int k = 0; k < K; ++k) {
    if(!std::isfinite(prop_pi[k]) || prop_pi[k] <= 0.0) {
      Rcpp::stop("Mixing proportions must be finite and strictly positive");
    }
    mixing_sum += prop_pi[k];
    means[k] = as<arma::vec>(mu[k]);
    covariances[k] = as<arma::mat>(sigma[k]);
    if(means[k].n_elem != static_cast<unsigned int>(d) ||
       covariances[k].n_rows != static_cast<unsigned int>(d) ||
       covariances[k].n_cols != static_cast<unsigned int>(d) ||
       !means[k].is_finite() || !covariances[k].is_finite()) {
      Rcpp::stop("Mean/covariance dimensions or values are invalid");
    }
  }
  const double mixing_tolerance =
    64.0 * std::numeric_limits<double>::epsilon() * std::max(1, K);
  if(!std::isfinite(mixing_sum) ||
     std::abs(mixing_sum - 1.0) > mixing_tolerance) {
    Rcpp::stop("Mixing proportions must sum to one");
  }

  NumericMatrix log_component(n, K);
  const double log_two_pi = std::log(2.0 * std::acos(-1.0));
  for(size_t pp = 0; pp < patterns.size(); ++pp) {
    const MixedPattern& pattern = patterns[pp];
    const arma::uvec observed_index = mixed_uvec(pattern.observed);
    int selected_missing = 0;
    for(size_t jj = 0; jj < selected.size(); ++jj) {
      if(R_IsNA(YNA(pattern.rows[0], selected[jj]))) ++selected_missing;
    }

    for(int k = 0; k < K; ++k) {
      arma::vec observed_mean;
      arma::mat chol_lower;
      double log_determinant = 0.0;
      if(!pattern.observed.empty()) {
        observed_mean = means[k].elem(observed_index);
        const arma::mat observed_covariance =
          covariances[k].submat(observed_index, observed_index);
        if(!arma::chol(chol_lower, observed_covariance, "lower")) {
          Rcpp::stop("Observed covariance submatrix is not positive definite");
        }
        log_determinant =
          2.0 * arma::sum(arma::log(chol_lower.diag()));
        if(!std::isfinite(log_determinant)) {
          Rcpp::stop("Observed covariance log-determinant is not finite");
        }
      }

      double log_mask = 0.0;
      if(mechanism_dimension > 0) {
        if(!std::isfinite(rho[k]) || rho[k] <= 0.0 || rho[k] >= 1.0) {
          Rcpp::stop("MNARz probability must be strictly between zero and one");
        }
        log_mask = selected_missing * std::log(rho[k]) +
          (mechanism_dimension - selected_missing) * std::log1p(-rho[k]);
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
          log_density = -0.5 * (
            pattern.observed.size() * log_two_pi + log_determinant +
            arma::dot(standardized, standardized)
          );
        }
        log_component(i, k) =
          std::log(prop_pi[k]) + log_density + log_mask;
      }
    }
  }

  NumericMatrix tik(n, K);
  double loglikelihood = 0.0;
  for(int i = 0; i < n; ++i) {
    double maximum = R_NegInf;
    for(int k = 0; k < K; ++k) {
      maximum = std::max(maximum, log_component(i, k));
    }
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
    Rcpp::stop("Observed log-likelihood is non-finite");
  }

  return List::create(
    Named("loglik_obs") = loglikelihood,
    Named("tik") = tik
  );
}



NumericVector update_mixed_rho(
    const NumericMatrix& YNA,
    const NumericMatrix& tik,
    const std::vector<int>& selected) {
  const int n = YNA.nrow();
  const int K = tik.ncol();
  NumericVector rho(K, NA_REAL);
  if(selected.empty()) return rho;

  for(int k = 0; k < K; ++k) {
    const double effective_size = sum(tik(_, k));
    double weighted_missing = 0.0;
    for(int i = 0; i < n; ++i) {
      int missing_count = 0;
      for(size_t jj = 0; jj < selected.size(); ++jj) {
        if(R_IsNA(YNA(i, selected[jj]))) ++missing_count;
      }
      weighted_missing += tik(i, k) * missing_count;
    }
    rho[k] = weighted_missing / (selected.size() * effective_size);
    if(!std::isfinite(rho[k]) || rho[k] <= 0.0 || rho[k] >= 1.0) {
      Rcpp::stop(
        "Missingness M-step has no finite mixed class-only MNARz intercept for component %d",
        k + 1
      );
    }
  }
  return rho;
}

MixedState initialize_mixed_state(
    NumericMatrix YNA,
    int K,
    bool diag,
    const LogicalVector& is_mnar,
    const std::vector<int>& selected,
    Nullable<List> init,
    Nullable<int> samplesize,
    int& covariance_adjustments) {
  const int n = YNA.nrow();
  const int d = YNA.ncol();

  List base;
  List supplied;
  if(init.isNull()) {
    base = InitEMGaussian(YNA, K, "MNARz", diag, R_NilValue, samplesize);
  } else {
    supplied = clone(as<List>(init));
    if(!supplied.containsElementNamed("pik") ||
       !supplied.containsElementNamed("mu") ||
       !supplied.containsElementNamed("sigma")) {
      Rcpp::stop("init must contain pik, mu, and sigma");
    }
    if(!supplied.containsElementNamed("alpha")) {
      supplied["alpha"] = NumericMatrix(K, d);
    }
    Nullable<List> supplied_nullable(supplied);
    base = InitEMGaussian(YNA, K, "MNARz", diag, supplied_nullable, samplesize);
  }

  MixedState state;
  state.pi = clone(as<NumericVector>(base["pi_init"]));
  state.mu = clone(as<List>(base["mu_init"]));
  state.sigma = clone(as<List>(base["sigma_init"]));
  state.beta = NumericMatrix(K, d);
  state.beta.fill(0.0);

  if(state.pi.size() != K || state.mu.size() != K || state.sigma.size() != K) {
    Rcpp::stop("Initial parameter dimensions do not match K");
  }
  double pi_sum = 0.0;
  for(int k = 0; k < K; ++k) {
    if(!std::isfinite(state.pi[k]) || state.pi[k] <= 0.0) {
      Rcpp::stop("Initial mixing proportions must be finite and strictly positive");
    }
    pi_sum += state.pi[k];
  }
  if(!std::isfinite(pi_sum) || pi_sum <= 0.0) {
    Rcpp::stop("Initial mixing proportions must have a positive finite sum");
  }
  for(int k = 0; k < K; ++k) state.pi[k] /= pi_sum;

  for(int k = 0; k < K; ++k) {
    const arma::vec mu_k = as<arma::vec>(state.mu[k]);
    const arma::mat sigma_k = as<arma::mat>(state.sigma[k]);
    if(mu_k.n_elem != static_cast<unsigned int>(d) ||
       sigma_k.n_rows != static_cast<unsigned int>(d) ||
       sigma_k.n_cols != static_cast<unsigned int>(d) ||
       !mu_k.is_finite() || !sigma_k.is_finite()) {
      Rcpp::stop("Initial mean/covariance dimensions or values are invalid");
    }
    bool adjusted = false;
    state.sigma[k] = wrap(stabilize_covariance(sigma_k, diag, adjusted));
    if(adjusted) ++covariance_adjustments;
  }

  state.rho = NumericVector(K, NA_REAL);
  if(!init.isNull() && supplied.containsElementNamed("beta")) {
    validate_zero_beta(as<NumericMatrix>(supplied["beta"]), K, d);
  }
  if(!init.isNull() && selected.empty()) {
    if(supplied.containsElementNamed("rho")) {
      const NumericVector supplied_rho = as<NumericVector>(supplied["rho"]);
      if(supplied_rho.size() != K) Rcpp::stop("init$rho must have length K");
      for(int k = 0; k < K; ++k) {
        if(!R_IsNA(supplied_rho[k])) {
          Rcpp::stop(
            "init$rho must be all NA when no MNARz coordinate is selected"
          );
        }
      }
    }
    if(as<List>(init).containsElementNamed("alpha")) {
      const NumericMatrix supplied_alpha = as<NumericMatrix>(as<List>(init)["alpha"]);
      if(supplied_alpha.nrow() != K || supplied_alpha.ncol() != d) {
        Rcpp::stop("init$alpha must have dimensions K by ncol(YNA)");
      }
      for(int k = 0; k < K; ++k) {
        for(int j = 0; j < d; ++j) {
          if(!R_IsNA(supplied_alpha(k, j))) {
            Rcpp::stop(
              "init$alpha must be all NA when no MNARz coordinate is selected"
            );
          }
        }
      }
    }
  }
  if(!selected.empty()) {
    bool rho_supplied = false;
    bool alpha_supplied = false;
    NumericVector rho_from_supplied;
    NumericVector rho_from_supplied_alpha;
    if(!init.isNull()) {
      if(supplied.containsElementNamed("rho")) {
        rho_from_supplied = clone(as<NumericVector>(supplied["rho"]));
        if(rho_from_supplied.size() != K) Rcpp::stop("init$rho must have length K");
        for(int k = 0; k < K; ++k) {
          if(!std::isfinite(rho_from_supplied[k]) ||
             rho_from_supplied[k] <= 0.0 || rho_from_supplied[k] >= 1.0) {
            Rcpp::stop("init$rho must be finite and strictly between zero and one");
          }
        }
        rho_supplied = true;
      }
      if(as<List>(init).containsElementNamed("alpha")) {
        rho_from_supplied_alpha = rho_from_alpha(
          as<NumericMatrix>(as<List>(init)["alpha"]), K, d, selected
        );
        alpha_supplied = true;
      }
    }

    if(rho_supplied && alpha_supplied) {
      for(int k = 0; k < K; ++k) {
        const double scale = 1.0 + std::max(
          std::abs(rho_from_supplied[k]), std::abs(rho_from_supplied_alpha[k])
        );
        if(std::abs(rho_from_supplied[k] - rho_from_supplied_alpha[k]) >
           64.0 * std::numeric_limits<double>::epsilon() * scale) {
          Rcpp::stop("init$rho and init$alpha encode different MNARz probabilities");
        }
      }
    }

    if(rho_supplied) {
      state.rho = rho_from_supplied;
    } else if(alpha_supplied) {
      state.rho = rho_from_supplied_alpha;
    } else {
      double missing_count = 0.0;
      for(int i = 0; i < n; ++i) {
        for(size_t jj = 0; jj < selected.size(); ++jj) {
          if(R_IsNA(YNA(i, selected[jj]))) missing_count += 1.0;
        }
      }
      const double probability = missing_count / (n * selected.size());
      if(!std::isfinite(probability) || probability <= 0.0 || probability >= 1.0) {
        Rcpp::stop(
          "Initial selected MNARz mask is on the all-observed or all-missing boundary"
        );
      }
      state.rho.fill(probability);
    }
  }
  state.alpha = alpha_from_rho(state.rho, K, d, selected);
  return state;
}

void mixed_gaussian_m_step_patterned(
    const NumericMatrix& YNA,
    const NumericMatrix& tik,
    const List& current_mu,
    const List& current_sigma,
    const std::vector<MixedPattern>& patterns,
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
      const MixedPattern& pattern = patterns[pp];
      const arma::uvec observed_index = mixed_uvec(pattern.observed);
      const arma::uvec missing_index = mixed_uvec(pattern.missing);
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
      new_covariance = arma::diagmat(
        second_diagonal / effective_size - arma::square(new_mean)
      );
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

NumericMatrix mixed_impute_patterned(
    const NumericMatrix& YNA,
    const NumericMatrix& tik,
    const List& mu,
    const List& sigma,
    const std::vector<MixedPattern>& patterns) {
  NumericMatrix imputed = clone(YNA);
  const int K = tik.ncol();
  for(size_t pp = 0; pp < patterns.size(); ++pp) {
    const MixedPattern& pattern = patterns[pp];
    if(pattern.missing.empty()) continue;
    const arma::uvec observed_index = mixed_uvec(pattern.observed);
    const arma::uvec missing_index = mixed_uvec(pattern.missing);
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
          imputed(i, pattern.missing[jj]) +=
            tik(i, k) * conditional_missing[jj];
        }
      }
    }
  }
  return imputed;
}

} // namespace

// Observed likelihood for the semiparametric mixed specification. FALSE mask
// coordinates have an unrestricted ignorable-MAR factor common to every
// component and are omitted.  TRUE coordinates share one rho_k within each
// component.  E_mu_list/E_sig_list remain only for native ABI compatibility.
// [[Rcpp::export]]
List LoglikelihoodObsGaussianMixed(
    const NumericMatrix& YNA,
    const List& mu,
    const List& sigma,
    const NumericMatrix& alpha,
    const NumericMatrix& beta,
    const NumericVector& prop_pi,
    const LogicalVector& is_mnar,
    Nullable<List> E_mu_list = R_NilValue,
    Nullable<List> E_sig_list = R_NilValue) {
  validate_mixed_data(YNA);
  const int K = prop_pi.size();
  const int d = YNA.ncol();
  const std::vector<int> selected = mnar_indices(is_mnar, d);
  validate_zero_beta(beta, K, d);
  const NumericVector rho = rho_from_alpha(alpha, K, d, selected);
  const std::vector<MixedPattern> patterns = build_mixed_patterns(YNA);
  List result = mixed_observed_state(
    YNA, mu, sigma, rho, prop_pi, selected, patterns
  );
  result["rho"] = rho;
  result["mnar_coordinate_count"] = static_cast<int>(selected.size());
  result["conditional_arguments_ignored"] =
    E_mu_list.isNotNull() || E_sig_list.isNotNull();
  result["mechanism_specification"] =
    "semiparametric ignorable MAR plus class-only MNARz";
  return result;
}

// [[Rcpp::export]]
List InitEMGaussianMixed(
    NumericMatrix YNA,
    int K,
    std::string mecha,
    LogicalVector is_mnar,
    bool diag,
    Nullable<List> init,
    Nullable<int> samplesize) {
  if(mecha != "mixed") Rcpp::stop("InitEMGaussianMixed supports only mecha = 'mixed'");
  validate_mixed_data(YNA);
  if(K < 1) Rcpp::stop("K must be a positive integer");
  const std::vector<int> selected = mnar_indices(is_mnar, YNA.ncol());
  int covariance_adjustments = 0;
  MixedState state = initialize_mixed_state(
    YNA, K, diag, is_mnar, selected, init, samplesize, covariance_adjustments
  );
  return List::create(
    Named("pi_init") = state.pi,
    Named("mu_init") = state.mu,
    Named("sigma_init") = state.sigma,
    Named("alpha_init") = state.alpha,
    Named("beta_init") = state.beta,
    Named("rho_init") = state.rho,
    Named("covariance_adjustments") = covariance_adjustments,
    Named("mnar_coordinate_count") = static_cast<int>(selected.size()),
    Named("mechanism_specification") =
      "semiparametric ignorable MAR plus class-only MNARz"
  );
}

// Class-only MNARz M-step pooled over all selected coordinates.  Conditional
// moments and beta/current_alpha are accepted solely to preserve the native
// call signature; they do not belong to this likelihood.
// [[Rcpp::export]]
List MechanismEMGLMMixed(
    NumericMatrix YNA,
    NumericMatrix tik,
    std::string mecha,
    LogicalVector is_mnar,
    List E_y_list,
    NumericMatrix current_alpha,
    NumericMatrix current_beta) {
  if(mecha != "mixed") Rcpp::stop("MechanismEMGLMMixed supports only mecha = 'mixed'");
  validate_mixed_data(YNA);
  const int K = tik.ncol();
  const int d = YNA.ncol();
  if(tik.nrow() != YNA.nrow()) Rcpp::stop("tik must have nrow(YNA) rows");
  if(current_alpha.nrow() != K || current_alpha.ncol() != d) {
    Rcpp::stop("current_alpha must have dimensions K by ncol(YNA)");
  }
  validate_zero_beta(current_beta, K, d);
  const std::vector<int> selected = mnar_indices(is_mnar, d);
  // Although the M-step replaces alpha, an explicitly supplied current state
  // must still satisfy the active/inactive-coordinate specification used by
  // every mixed-estimator entry point.
  (void)rho_from_alpha(current_alpha, K, d, selected);
  validate_component_masses_mixed(tik, 0.0);
  const NumericVector rho = update_mixed_rho(YNA, tik, selected);
  const NumericMatrix alpha = alpha_from_rho(rho, K, d, selected);
  NumericMatrix beta(K, d);
  beta.fill(0.0);
  return List::create(
    Named("alpha_new") = alpha,
    Named("beta_new") = beta,
    Named("rho_new") = rho,
    Named("mnar_coordinate_count") = static_cast<int>(selected.size()),
    Named("conditional_arguments_ignored") = E_y_list.size() > 0,
    Named("beta_deprecated") = true,
    Named("mechanism_specification") =
      "semiparametric ignorable MAR plus class-only MNARz"
  );
}

// [[Rcpp::export]]
List EMGaussianMixed(
    NumericMatrix YNA,
    int K,
    std::string mecha,
    LogicalVector is_mnar,
    bool diag,
    int rmax,
    Nullable<List> init = R_NilValue,
    double tol = 0.0001,
    Nullable<int> samplesize = R_NilValue) {
  if(mecha != "mixed") {
    Rcpp::stop("EMGaussianMixed supports only mecha = 'mixed'");
  }
  validate_mixed_data(YNA);
  validate_gaussian_observation_design(YNA, diag);
  const int n = YNA.nrow();
  const int d = YNA.ncol();
  if(K < 1 || K > n) Rcpp::stop("K must be between one and nrow(YNA)");
  if(rmax < 1) Rcpp::stop("rmax must be a positive integer");
  if(!std::isfinite(tol) || tol < 0.0) {
    Rcpp::stop("tol must be finite and non-negative");
  }

  const std::vector<int> selected = mnar_indices(is_mnar, d);
  const std::vector<MixedPattern> patterns = build_mixed_patterns(YNA);
  const double component_floor = diag ? 1.0 :
    static_cast<double>(std::max(1, d));

  int covariance_adjustments = 0;
  MixedState state = initialize_mixed_state(
    YNA, K, diag, is_mnar, selected, init, samplesize,
    covariance_adjustments
  );
  List observed = mixed_observed_state(
    YNA, state.mu, state.sigma, state.rho, state.pi, selected, patterns
  );
  double loglikelihood = as<double>(observed["loglik_obs"]);
  NumericMatrix tik = as<NumericMatrix>(observed["tik"]);
  NumericVector loglikelihood_trace = NumericVector::create(loglikelihood);
  double minimum_component_size = validate_component_masses_mixed(
    tik, component_floor
  );

  bool converged = false;
  bool loglik_monotone = true;
  bool no_material_loglik_decrease = true;
  double final_improvement = NA_REAL;
  double final_relative_improvement = NA_REAL;
  double final_convergence_threshold = NA_REAL;
  double final_decrease_tolerance = NA_REAL;
  int iteration = 0;

  while(iteration < rmax && !converged) {
    ++iteration;
    NumericVector candidate_pi(K);
    for(int k = 0; k < K; ++k) {
      candidate_pi[k] = sum(tik(_, k)) / n;
    }
    List candidate_mu(K);
    List candidate_sigma(K);
    mixed_gaussian_m_step_patterned(
      YNA, tik, state.mu, state.sigma, patterns, diag,
      candidate_mu, candidate_sigma, covariance_adjustments
    );
    const NumericVector candidate_rho =
      update_mixed_rho(YNA, tik, selected);
    const NumericMatrix candidate_alpha =
      alpha_from_rho(candidate_rho, K, d, selected);

    List candidate = mixed_observed_state(
      YNA, candidate_mu, candidate_sigma, candidate_rho, candidate_pi,
      selected, patterns
    );
    const double candidate_loglikelihood =
      as<double>(candidate["loglik_obs"]);
    NumericMatrix candidate_tik = as<NumericMatrix>(candidate["tik"]);
    const double candidate_minimum = validate_component_masses_mixed(
      candidate_tik, component_floor
    );
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

    state.pi = clone(candidate_pi);
    state.mu = clone(candidate_mu);
    state.sigma = clone(candidate_sigma);
    state.rho = clone(candidate_rho);
    state.alpha = clone(candidate_alpha);
    tik = clone(candidate_tik);
    loglikelihood = candidate_loglikelihood;
    loglikelihood_trace.push_back(loglikelihood);
    converged = std::abs(final_improvement) <= final_convergence_threshold;
  }

  NumericMatrix imputed_data = mixed_impute_patterned(
    YNA, tik, state.mu, state.sigma, patterns
  );

  return List::create(
    Named("pik") = state.pi,
    Named("mu") = state.mu,
    Named("sigma") = state.sigma,
    Named("rho") = state.rho,
    Named("alpha") = state.alpha,
    Named("beta") = state.beta,
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
    Named("mnar_coordinate_count") = static_cast<int>(selected.size()),
    Named("ignored_mar_coordinate_count") =
      d - static_cast<int>(selected.size()),
    Named("missingness_parameter_count") = selected.empty() ? 0 : K,
    Named("mechanism_specification") =
      "semiparametric ignorable MAR plus class-only MNARz",
    Named("alpha_scope") =
      "alpha is defined only on is_mnar = TRUE coordinates; FALSE columns are NA",
    Named("beta_deprecated") = true,
    Named("error") = "No error"
  );
}

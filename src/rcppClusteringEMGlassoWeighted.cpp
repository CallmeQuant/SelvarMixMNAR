#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <vector>


using namespace std;
using namespace Rcpp;
using namespace arma;

#include "Mixture.h"
#include "Function.h"

namespace {

constexpr int kGlassoFastMaxIterations = 1000;
constexpr double kGlassoFastTolerance = 1e-7;
constexpr int kMeanMaximumCoordinateSweeps = 10000;
constexpr double kMeanKKTRelativeTolerance = 1e-8;

struct MatrixPairDiagnostics {
    double inverse_residual;
    double normalized_inverse_residual;
    double tolerance;
};

struct GlassoKKTDiagnostics {
    double residual;
    double normalized_residual;
    double tolerance;
};

bool is_symmetric_positive_definite(const arma::mat& matrix) {
    if (!matrix.is_square() || !matrix.is_finite())
        return false;
    const double scale = std::max(1.0, arma::norm(matrix, "inf"));
    if (arma::norm(matrix - matrix.t(), "inf") > 1e-8 * scale)
        return false;
    arma::mat factor;
    const arma::mat symmetric = 0.5 * (matrix + matrix.t());
    return arma::chol(factor, symmetric);
}

double inverse_validation_tolerance(arma::uword dimension) {
    return std::max(
        1e-10,
        25.0 * std::sqrt(1.0 + static_cast<double>(dimension)) *
            kGlassoFastTolerance
    );
}

MatrixPairDiagnostics covariance_precision_diagnostics(
        const arma::mat& covariance,
        const arma::mat& precision) {
    if (covariance.n_rows != precision.n_rows ||
        covariance.n_cols != precision.n_cols ||
        !covariance.is_square())
        return MatrixPairDiagnostics{
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            inverse_validation_tolerance(covariance.n_rows)
        };
    const arma::mat identity = arma::eye(covariance.n_rows,
                                         covariance.n_cols);
    const double residual = std::max(
        arma::norm(covariance * precision - identity, "inf"),
        arma::norm(precision * covariance - identity, "inf")
    );
    // This is a normwise backward error.  Unlike a raw inverse residual it
    // remains interpretable when the graphical-lasso solution is ill-scaled.
    const double denominator = 1.0 +
        arma::norm(covariance, "inf") * arma::norm(precision, "inf");
    const double normalized = residual / denominator;
    return MatrixPairDiagnostics{
        residual, normalized,
        inverse_validation_tolerance(covariance.n_rows)
    };
}

GlassoKKTDiagnostics glasso_kkt_diagnostics(
        const arma::mat& covariance,
        const arma::mat& precision,
        const arma::mat& empirical,
        const arma::mat& penalty) {
    const double precision_scale = std::max(
        1.0, arma::abs(precision).max()
    );
    const double zero_tolerance = 64.0 *
        std::numeric_limits<double>::epsilon() * precision_scale;
    double residual = 0.0;

    for (arma::uword row = 0; row < precision.n_rows; ++row) {
        for (arma::uword column = 0; column < precision.n_cols; ++column) {
            const double gradient = covariance(row, column) -
                empirical(row, column);
            const double weight = penalty(row, column);
            double entry_residual = 0.0;
            if (std::abs(precision(row, column)) > zero_tolerance) {
                entry_residual = std::abs(
                    gradient - weight *
                        (precision(row, column) > 0.0 ? 1.0 : -1.0)
                );
            } else {
                entry_residual = std::max(
                    0.0, std::abs(gradient) - weight
                );
            }
            residual = std::max(residual, entry_residual);
        }
    }

    const double scale = 1.0 + arma::abs(covariance).max() +
        arma::abs(empirical).max() + arma::abs(penalty).max();
    const double normalized = residual / scale;
    const double tolerance = std::max(
        1e-10,
        25.0 * std::sqrt(1.0 + static_cast<double>(precision.n_rows)) *
            kGlassoFastTolerance
    );
    return GlassoKKTDiagnostics{residual, normalized, tolerance};
}

MatrixPairDiagnostics require_valid_covariance_precision(
        const arma::mat& covariance,
        const arma::mat& precision,
        int component,
        const std::string& context) {
    const MatrixPairDiagnostics diagnostics =
        covariance_precision_diagnostics(covariance, precision);
    if (!is_symmetric_positive_definite(covariance) ||
        !is_symmetric_positive_definite(precision) ||
        !std::isfinite(diagnostics.normalized_inverse_residual) ||
        diagnostics.normalized_inverse_residual > diagnostics.tolerance) {
        Rcpp::stop(
            "%s returned an invalid covariance/precision pair for component %d "
            "(normalized inverse residual=%.8g, tolerance=%.8g)",
            context.c_str(), component + 1,
            diagnostics.normalized_inverse_residual, diagnostics.tolerance
        );
    }
    return diagnostics;
}

} // namespace

Mixture::Mixture(Rcpp::List     InputList,
    double         lambda_mu,
    double         lambda_omega,
    arma::cube     Pk_in,
    bool           inner_glasso_warm_start,
    bool           predecessor_fitted)
{
if (InputList.size() != 6) {
    Rcpp::stop("InputList must contain data, proportions, means, covariance, precision, and responsibilities");
}

Rcpp::NumericMatrix x_r = Rcpp::as<Rcpp::NumericMatrix>(InputList[0]);
Rcpp::NumericVector prop_r = Rcpp::as<Rcpp::NumericVector>(InputList[1]);
Rcpp::NumericMatrix mu_r = Rcpp::as<Rcpp::NumericMatrix>(InputList[2]);
Rcpp::NumericVector vecS = Rcpp::as<Rcpp::NumericVector>(InputList[3]);
Rcpp::NumericVector vecW = Rcpp::as<Rcpp::NumericVector>(InputList[4]);
Rcpp::NumericMatrix z_r = Rcpp::as<Rcpp::NumericMatrix>(InputList[5]);

n = x_r.nrow();
p = x_r.ncol();
nbClust = prop_r.size();
if (n < 1 || p < 1 || nbClust < 1) {
    Rcpp::stop("InputList dimensions must be positive");
}
if (mu_r.nrow() != p || mu_r.ncol() != nbClust ||
    z_r.nrow() != n || z_r.ncol() != nbClust ||
    vecS.size() != static_cast<R_xlen_t>(p) * p * nbClust ||
    vecW.size() != static_cast<R_xlen_t>(p) * p * nbClust) {
    Rcpp::stop("InputList contains inconsistent parameter dimensions");
}

// Each fit owns every mutable array.  In particular, covariance and
// precision updates must never write through to the caller's R objects.
Xd = arma::mat(x_r.begin(), n, p, true, true);
prop = arma::rowvec(prop_r.begin(), nbClust, true, true);
Mu = arma::mat(mu_r.begin(), p, nbClust, true, true);
CovarianceMatrix = arma::cube(vecS.begin(), p, p, nbClust, true, true);
EmpiricalCovariance = CovarianceMatrix;
PrecisionMatrix = arma::cube(vecW.begin(), p, p, nbClust, true, true);
ProbCond = arma::mat(z_r.begin(), n, nbClust, true, true);

if (!(Xd.is_finite() && prop.is_finite() && Mu.is_finite() &&
      CovarianceMatrix.is_finite() && PrecisionMatrix.is_finite() &&
      ProbCond.is_finite())) {
    Rcpp::stop("InputList contains non-finite values");
}
if (!PropOK()) {
    Rcpp::stop("Initial component proportions must be positive and sum to one");
}

lambda = lambda_mu;
rho    = lambda_omega;
innerGlassoWarmStart = inner_glasso_warm_start;
predecessorFitted = predecessor_fitted;

Pk_cube = Pk_in;   // deep copy

if (   Pk_cube.n_rows   != static_cast<uword>(p)
|| Pk_cube.n_cols   != static_cast<uword>(p)
|| Pk_cube.n_slices != static_cast<uword>(nbClust) )
Rcpp::stop("Pk_in has wrong dimension (must be p×p×K)");

for (uword k = 0; k < Pk_cube.n_slices; ++k) {
    const double penalty_scale = std::max(
        1.0, arma::norm(Pk_cube.slice(k), "inf")
    );
    if (arma::norm(
            Pk_cube.slice(k) - Pk_cube.slice(k).t(), "inf"
        ) > 1e-12 * penalty_scale) {
        Rcpp::stop("Pk_in must be symmetric for component %d",
                   static_cast<int>(k + 1));
    }
    if (arma::any(arma::abs(Pk_cube.slice(k).diag()) > 0.0))
        Rcpp::stop("Pk_in diagonal must be exactly zero for component %d",
                   static_cast<int>(k + 1));
    require_valid_covariance_precision(
        CovarianceMatrix.slice(k), PrecisionMatrix.slice(k),
        static_cast<int>(k), "Initial state"
    );
}

LastGlassoEmpiricalCovariance.zeros(p, p, nbClust);
LastGlassoPenalty.zeros(p, p, nbClust);
LastGlassoNIterations.set_size(nbClust);
LastGlassoNIterations.fill(-1);
LastGlassoInverseResiduals.zeros(nbClust);
LastGlassoNormalizedInverseResiduals.zeros(nbClust);
LastGlassoInverseTolerances.zeros(nbClust);
LastGlassoKKTResiduals.zeros(nbClust);
LastGlassoNormalizedKKTResiduals.zeros(nbClust);
LastGlassoKKTTolerances.zeros(nbClust);
LastGlassoSolvers.assign(nbClust, "none");
LastGlassoStarts.assign(nbClust, "none");
}

Mixture::Mixture(Rcpp::List InputList)
{
    if (InputList.size() != 6) {
        Rcpp::stop("InputList must contain data, proportions, means, covariance, precision, and responsibilities");
    }
    Rcpp::NumericMatrix x_r = Rcpp::as<Rcpp::NumericMatrix>(InputList[0]);
    Rcpp::NumericVector prop_r = Rcpp::as<Rcpp::NumericVector>(InputList[1]);
    Rcpp::NumericMatrix mu_r = Rcpp::as<Rcpp::NumericMatrix>(InputList[2]);
    Rcpp::NumericVector vecS = Rcpp::as<Rcpp::NumericVector>(InputList[3]);
    Rcpp::NumericVector vecW = Rcpp::as<Rcpp::NumericVector>(InputList[4]);
    Rcpp::NumericMatrix z_r = Rcpp::as<Rcpp::NumericMatrix>(InputList[5]);

    n = x_r.nrow();
    p = x_r.ncol();
    nbClust = prop_r.size();
    if (n < 1 || p < 1 || nbClust < 1 ||
        mu_r.nrow() != p || mu_r.ncol() != nbClust ||
        z_r.nrow() != n || z_r.ncol() != nbClust ||
        vecS.size() != static_cast<R_xlen_t>(p) * p * nbClust ||
        vecW.size() != static_cast<R_xlen_t>(p) * p * nbClust) {
        Rcpp::stop("InputList contains inconsistent parameter dimensions");
    }

    Xd = arma::mat(x_r.begin(), n, p, true, true);
    prop = arma::rowvec(prop_r.begin(), nbClust, true, true);
    Mu = arma::mat(mu_r.begin(), p, nbClust, true, true);
    CovarianceMatrix = arma::cube(vecS.begin(), p, p, nbClust, true, true);
    EmpiricalCovariance = CovarianceMatrix;
    PrecisionMatrix = arma::cube(vecW.begin(), p, p, nbClust, true, true);
    ProbCond = arma::mat(z_r.begin(), n, nbClust, true, true);

    if (!(Xd.is_finite() && prop.is_finite() && Mu.is_finite() &&
          CovarianceMatrix.is_finite() && PrecisionMatrix.is_finite() &&
          ProbCond.is_finite()) || !PropOK()) {
        Rcpp::stop("InputList contains invalid initial parameters");
    }

    lambda = 0.0;        
    rho    = 0.0;
    innerGlassoWarmStart = false;
    predecessorFitted = false;

    Pk_cube.ones(p, p, nbClust);
    for (int k = 0; k < nbClust; ++k) {
        Pk_cube.slice(k).diag().zeros();  // no penalty on diagonal
        require_valid_covariance_precision(
            CovarianceMatrix.slice(k), PrecisionMatrix.slice(k), k,
            "Initial state"
        );
    }
    LastGlassoEmpiricalCovariance.zeros(p, p, nbClust);
    LastGlassoPenalty.zeros(p, p, nbClust);
    LastGlassoNIterations.set_size(nbClust);
    LastGlassoNIterations.fill(-1);
    LastGlassoInverseResiduals.zeros(nbClust);
    LastGlassoNormalizedInverseResiduals.zeros(nbClust);
    LastGlassoInverseTolerances.zeros(nbClust);
    LastGlassoKKTResiduals.zeros(nbClust);
    LastGlassoNormalizedKKTResiduals.zeros(nbClust);
    LastGlassoKKTTolerances.zeros(nbClust);
    LastGlassoSolvers.assign(nbClust, "none");
    LastGlassoStarts.assign(nbClust, "none");
}

// Cholesky log determinant for a finite positive-definite matrix.
double Mixture::safe_log_det(const mat& M, bool& success) const
{
    double log_det_val = -std::numeric_limits<double>::infinity();
    success = false;
    if (M.is_finite() && M.is_square()) {
        mat R;
        if (chol(R, M)) {
            // M = R'R and R is triangular, hence log det(M) is twice the
            // sum of the logarithms of the diagonal entries of R.
            log_det_val = 2.0 * accu(log(R.diag()));
            if (std::isfinite(log_det_val)) {
                 success = true;
            } else {
                log_det_val = -std::numeric_limits<double>::infinity();
            }
        }
    }
    return log_det_val;
}

double Mixture::effective_component_size(int component) const
{
    const colvec weights = ProbCond.col(component);
    const double mass = arma::accu(weights);
    const double squared_mass = arma::dot(weights, weights);
    if (!std::isfinite(mass) || !std::isfinite(squared_mass) ||
        mass <= 0.0 || squared_mass <= 0.0) {
        return 0.0;
    }
    return mass * mass / squared_mass;
}

bool Mixture::component_fully_regularized(int component) const
{
    if (!(rho > 0.0)) return false;
    for (int row = 0; row < p; ++row) {
        for (int column = row + 1; column < p; ++column) {
            if (!(Pk_cube(row, column, component) > 0.0) ||
                !(Pk_cube(column, row, component) > 0.0)) {
                return false;
            }
        }
    }
    return true;
}

double Mixture::minimum_effective_component_size(int component) const
{
    // An unpenalized full covariance needs more than p effective
    // observations.  With every off-diagonal direction regularized, the
    // graphical-lasso subproblem is allowed once there is more than one
    // effective observation; positive coordinate variances and the solver
    // KKT conditions are checked separately below.
    const double base = component_fully_regularized(component) ?
        1.0 : static_cast<double>(p);
    return base + 64.0 * std::numeric_limits<double>::epsilon() *
        (1.0 + static_cast<double>(p));
}

bool Mixture::ComponentSupportOK() const
{
    for (int component = 0; component < nbClust; ++component) {
        if (!(effective_component_size(component) >
              minimum_effective_component_size(component))) {
            return false;
        }
    }
    return true;
}

// Observed Gaussian-mixture log likelihood minus the mean and precision
// penalties used by the variable-ranking stage.
double Mixture::PenLogLik(void)
{
    mat lD = zeros<mat>(n, nbClust);
    double SLogDet_Sigma = 0.0;
    bool log_det_ok;

    for(int k = 0; k < nbClust; k++)
    {
        // log det(Sigma_k) = -log det(Omega_k).
        double log_det_Omega = this->safe_log_det(PrecisionMatrix.slice(k), log_det_ok);
        if (!log_det_ok) {
             lD.col(k).fill(-std::numeric_limits<double>::infinity());
             continue;
        }
        SLogDet_Sigma = -log_det_Omega;

        for(int i = 0; i < n; i++) {
            lD(i,k) = log(prop(k)) + ldcppmvt(trans(Xd.row(i)), Mu.col(k), PrecisionMatrix.slice(k), SLogDet_Sigma);
        }
    }

    // Normalize the component contributions by log-sum-exp.
    double loglik_obs = 0.0;
    for (int i = 0; i < n; ++i) {
        rowvec log_dens_i = lD.row(i);
        uvec finite_components = find_finite(log_dens_i);
        if (finite_components.n_elem == 0) {
             Rcpp::stop("All component log-densities are invalid for observation %d", i + 1);
        }
        vec finite_log_dens = log_dens_i.elem(finite_components);
        double max_log_d = finite_log_dens.max();
        loglik_obs += max_log_d + log(accu(exp(finite_log_dens - max_log_d)));
    }

    double mu_penalty = lambda * accu(abs(Mu));

    double precision_penalty = 0.0;
    for(int k = 0; k < nbClust; k++) {
        precision_penalty += accu(abs(Pk_cube.slice(k) % PrecisionMatrix.slice(k)));
    }

    double penloglik = loglik_obs - mu_penalty - (rho * precision_penalty);

    return std::isfinite(penloglik) ? penloglik : -std::numeric_limits<double>::infinity();
};

// E-step responsibilities.
void Mixture::GetProbCond(void){
    mat T = zeros<mat>(n, nbClust);
    mat lD = zeros<mat>(n, nbClust);
    double SLogDet_Sigma = 0.0;
    bool log_det_ok;

    for(int k = 0; k < nbClust; k++)
    {
        double log_det_Omega = this->safe_log_det(PrecisionMatrix.slice(k), log_det_ok);
        if (!log_det_ok) {
             lD.col(k).fill(-std::numeric_limits<double>::infinity());
             continue;
        }
        SLogDet_Sigma = -log_det_Omega;

        for(int i = 0; i < n; i++)
        {
            lD(i,k) = log(prop(k)) + ldcppmvt(trans(Xd.row(i)), Mu.col(k), PrecisionMatrix.slice(k), SLogDet_Sigma);
        }
    }

    // Normalize probabilities using log-sum-exp for numerical stability
    for(int i = 0; i < n; i++) {
        rowvec log_p = lD.row(i);
        uvec finite_components = find_finite(log_p);
        if (finite_components.n_elem == 0) {
            Rcpp::stop("All component log-densities are invalid for observation %d", i + 1);
        }
        vec finite_log_p = log_p.elem(finite_components);
        double max_log_p = finite_log_p.max();
        vec p_norm = exp(finite_log_p - max_log_p);
        double sum_p_norm = accu(p_norm);
        if (!std::isfinite(sum_p_norm) ||
            sum_p_norm <= std::numeric_limits<double>::epsilon()) {
            Rcpp::stop("Responsibility normalization failed for observation %d", i + 1);
        }
        T.row(i).zeros();
        for (uword index = 0; index < finite_components.n_elem; ++index) {
            T(i, finite_components(index)) = p_norm(index) / sum_p_norm;
        }
    }
    ProbCond = T;
};

// --- M-Step: Update Proportions ---
void Mixture::GetClassesSizes(void){
    rowvec updated_prop = mean(ProbCond, 0);
    const double mass_floor = std::numeric_limits<double>::epsilon();
    if (!updated_prop.is_finite() || arma::any(updated_prop <= mass_floor)) {
        Rcpp::stop("A component has zero or near-zero responsibility mass");
    }
    if (std::abs(accu(updated_prop) - 1.0) > 1e-8) {
        Rcpp::stop("Updated component proportions do not sum to one");
    }
    for (int component = 0; component < nbClust; ++component) {
        const double effective_size = effective_component_size(component);
        const double required_size =
            minimum_effective_component_size(component);
        if (!(effective_size > required_size)) {
            Rcpp::stop(
                "Component %d has insufficient effective support for a "
                "%d-dimensional covariance (n_eff=%.8g must exceed %.8g; %s)",
                component + 1, p, effective_size, required_size,
                component_fully_regularized(component) ?
                    "all off-diagonal directions penalized" :
                    "an unpenalized precision direction remains"
            );
        }
    }
    prop = updated_prop;
}

// --- M-Step: exact coordinate descent for the penalized mean block ---
// For fixed responsibilities and Omega_k, this maximizes
//   -n_k/2 (mu_k-xbar_k)' Omega_k (mu_k-xbar_k) - lambda ||mu_k||_1.
// Gauss--Seidel updates are iterated to an explicit subgradient KKT bound;
// a single simultaneous/Jacobi sweep is not a valid ECM update in general.
void Mixture::UpdateMeans(int em_iteration)
{
    const rowvec component_masses = sum(ProbCond, 0);
    const double mass_floor = std::numeric_limits<double>::epsilon();
    if (!component_masses.is_finite() ||
        arma::any(component_masses <= mass_floor)) {
        Rcpp::stop(
            "A component has zero or near-zero responsibility mass in the mean update"
        );
    }

    for (int k = 0; k < nbClust; ++k) {
        const double nk = component_masses(k);
        const mat& precision = PrecisionMatrix.slice(k);
        if (arma::any(precision.diag() <= mass_floor)) {
            Rcpp::stop(
                "Precision matrix has a non-positive diagonal for component %d",
                k + 1
            );
        }

        const colvec weighted_mean =
            Xd.t() * ProbCond.col(k) / nk;
        colvec updated_mean = Mu.col(k);
        const auto block_objective = [&](const colvec& candidate) {
            const colvec difference = candidate - weighted_mean;
            return -0.5 * nk *
                arma::as_scalar(difference.t() * precision * difference) -
                lambda * arma::accu(arma::abs(candidate));
        };
        const double objective_before = block_objective(updated_mean);
        const double score_scale = arma::abs(
            nk * precision * weighted_mean
        ).max();
        const double kkt_tolerance = std::max(
            1e-10,
            kMeanKKTRelativeTolerance *
                (1.0 + lambda + score_scale)
        );

        int sweeps = 0;
        double max_kkt_residual = std::numeric_limits<double>::infinity();
        bool converged = false;

        // With no mean penalty the unique block maximizer is the weighted
        // mean, independently of the current precision matrix.
        if (lambda == 0.0) {
            updated_mean = weighted_mean;
            max_kkt_residual = arma::abs(
                nk * precision * (weighted_mean - updated_mean)
            ).max();
            sweeps = 1;
            converged = max_kkt_residual <= kkt_tolerance;
        } else {
            colvec score = nk * precision *
                (weighted_mean - updated_mean);
            for (sweeps = 1;
                 sweeps <= kMeanMaximumCoordinateSweeps;
                 ++sweeps) {
                for (int j = 0; j < p; ++j) {
                    const double curvature = nk * precision(j, j);
                    const double numerator =
                        score(j) + curvature * updated_mean(j);
                    double new_value = 0.0;
                    if (numerator > lambda) {
                        new_value = (numerator - lambda) / curvature;
                    } else if (numerator < -lambda) {
                        new_value = (numerator + lambda) / curvature;
                    }
                    const double change = new_value - updated_mean(j);
                    if (change != 0.0) {
                        updated_mean(j) = new_value;
                        score -= nk * precision.col(j) * change;
                    }
                }

                // Recompute rather than accumulate the score across sweeps so
                // the reported KKT residual is not contaminated by drift.
                score = nk * precision * (weighted_mean - updated_mean);
                const double zero_tolerance = 64.0 *
                    std::numeric_limits<double>::epsilon() *
                    std::max(1.0, arma::abs(updated_mean).max());
                max_kkt_residual = 0.0;
                for (int j = 0; j < p; ++j) {
                    double residual = 0.0;
                    if (std::abs(updated_mean(j)) > zero_tolerance) {
                        residual = std::abs(
                            score(j) - lambda *
                                (updated_mean(j) > 0.0 ? 1.0 : -1.0)
                        );
                    } else {
                        residual = std::max(
                            0.0, std::abs(score(j)) - lambda
                        );
                    }
                    max_kkt_residual = std::max(
                        max_kkt_residual, residual
                    );
                }
                if (max_kkt_residual <= kkt_tolerance) {
                    converged = true;
                    break;
                }
            }
        }

        if (!converged) {
            Rcpp::stop(
                "Penalized mean coordinate descent failed its KKT check for "
                "component %d after %d sweeps (residual=%.8g, tolerance=%.8g)",
                k + 1, kMeanMaximumCoordinateSweeps,
                max_kkt_residual, kkt_tolerance
            );
        }

        const double objective_after = block_objective(updated_mean);
        const double objective_change = objective_after - objective_before;
        const double objective_tolerance = 128.0 *
            std::numeric_limits<double>::epsilon() *
            (1.0 + std::abs(objective_before) + std::abs(objective_after));
        if (!std::isfinite(objective_after) ||
            objective_change < -objective_tolerance) {
            Rcpp::stop(
                "Penalized mean block decreased its conditional objective for "
                "component %d (change=%.8g, tolerance=%.8g)",
                k + 1, objective_change, objective_tolerance
            );
        }

        Mu.col(k) = updated_mean;
        meanEMIterations.push_back(em_iteration + 1);
        meanComponents.push_back(k + 1);
        meanCoordinateSweeps.push_back(sweeps);
        meanKKTResiduals.push_back(max_kkt_residual);
        meanKKTTolerances.push_back(kkt_tolerance);
        meanBlockObjectiveChanges.push_back(objective_change);
        meanBlockObjectiveTolerances.push_back(objective_tolerance);
    }
};

// M-step covariance moments
void Mixture::GetEmpiricalCovariance(void){
    rowvec ProbCond_cols_sums = sum(ProbCond, 0);
    const double mass_floor = std::numeric_limits<double>::epsilon();
    if (!ProbCond_cols_sums.is_finite() ||
        arma::any(ProbCond_cols_sums <= mass_floor)) {
        Rcpp::stop("A component has zero or near-zero responsibility mass in the covariance update");
    }

    mat centered_Xd = Xd;

    for(int k = 0; k < nbClust; k++)
    {
        // BLAS-backed weighted crossproduct avoids an R/C++ level loop over
        // n rank-one matrices.
        centered_Xd = Xd;
        centered_Xd.each_row() -= trans(Mu.col(k));
        centered_Xd.each_col() %= arma::sqrt(ProbCond.col(k));
        const mat empirical =
            centered_Xd.t() * centered_Xd / ProbCond_cols_sums(k);
        const vec variances = empirical.diag();
        const double maximum_variance = variances.max();
        const double variance_floor = 64.0 *
            std::numeric_limits<double>::epsilon() *
            (1.0 + static_cast<double>(p)) * maximum_variance;
        if (!empirical.is_finite() || !(maximum_variance > 0.0) ||
            arma::any(variances <= variance_floor)) {
            Rcpp::stop(
                "Component %d has a zero or numerically unresolved weighted "
                "variance in its %d-dimensional covariance update",
                k + 1, p
            );
        }
        EmpiricalCovariance.slice(k) = 0.5 *
            (empirical + empirical.t());
    }
}

// --- M-Step: Update Covariance/Precision Matrices (Weighted Penalty) ---
void Mixture::UpdateCovarianceMatrices(int em_iteration){
    Environment glassoFast_pkg = Environment::namespace_env("glassoFast");
    Function RglassoFast = glassoFast_pkg["glassoFast"];

    rowvec ProbCond_cols_sums = sum(ProbCond, 0);
    const double mass_floor = std::numeric_limits<double>::epsilon();
    if (!ProbCond_cols_sums.is_finite() ||
        arma::any(ProbCond_cols_sums <= mass_floor)) {
        Rcpp::stop("A component has zero or near-zero responsibility mass in the precision update");
    }

    for(int k = 0; k < nbClust; k++)
    {
        long double nk = ProbCond_cols_sums(k);

        // The graphical-lasso subproblem uses L_k=(2*rho/n_k)P_k.
        mat rho_k_matrix = (2.0 * rho / nk) * Pk_cube.slice(k);
        mat empirical = EmpiricalCovariance.slice(k);
        mat W_new;
        mat Wi_new;
        int errflag = 0;
        int niter = 0;
        std::string solver = "glassoFast";
        std::string start = "cold";
        std::string seed_source = "none";

        if (!(empirical.is_finite() && rho_k_matrix.is_finite())) {
            Rcpp::stop("Non-finite graphical-lasso input for cluster %d", k + 1);
        }

        const mat off_diagonal = empirical - arma::diagmat(empirical.diag());
        const bool exactly_diagonal = arma::accu(arma::abs(off_diagonal)) == 0.0;

        if (exactly_diagonal) {
            // glassoFast 1.0.1's exact-diagonal early return omits diag(S).
            // Solve that convex subproblem analytically and label it explicitly;
            // this is neither a warm/cold glassoFast call nor a retry.
            const vec covariance_diagonal = empirical.diag() + rho_k_matrix.diag();
            if (!covariance_diagonal.is_finite() ||
                arma::any(covariance_diagonal <= 0.0)) {
                Rcpp::stop(
                    "Analytic diagonal graphical-lasso update has a non-positive diagonal for cluster %d",
                    k + 1
                );
            }
            W_new.zeros(p, p);
            Wi_new.zeros(p, p);
            W_new.diag() = covariance_diagonal;
            Wi_new.diag() = 1.0 / covariance_diagonal;
            solver = "analytic_diagonal";
            start = "analytic";
            seed_source = "none";
        } else {
            const bool use_predecessor_seed = predecessorFitted && em_iteration == 0;
            const bool use_previous_em_seed = innerGlassoWarmStart && em_iteration > 0;
            const bool use_warm = use_predecessor_seed || use_previous_em_seed;

            List GlassoResult;
            if (use_warm) {
                require_valid_covariance_precision(
                    CovarianceMatrix.slice(k), PrecisionMatrix.slice(k), k,
                    "Warm-start seed"
                );
                start = "warm";
                seed_source = use_predecessor_seed ?
                    "path_predecessor" : "previous_em";
                GlassoResult = RglassoFast(
                    Named("S") = empirical,
                    Named("rho") = rho_k_matrix,
                    Named("start") = "warm",
                    Named("w.init") = CovarianceMatrix.slice(k),
                    Named("wi.init") = PrecisionMatrix.slice(k),
                    Named("thr") = kGlassoFastTolerance,
                    Named("maxIt") = kGlassoFastMaxIterations
                );
            } else {
                GlassoResult = RglassoFast(
                    Named("S") = empirical,
                    Named("rho") = rho_k_matrix,
                    Named("start") = "cold",
                    Named("thr") = kGlassoFastTolerance,
                    Named("maxIt") = kGlassoFastMaxIterations
                );
            }

            if (!GlassoResult.containsElementNamed("w") ||
                !GlassoResult.containsElementNamed("wi") ||
                !GlassoResult.containsElementNamed("errflag") ||
                !GlassoResult.containsElementNamed("niter")) {
                Rcpp::stop("glassoFast returned an incomplete result for cluster %d",
                           k + 1);
            }
            W_new = as<mat>(GlassoResult["w"]);
            Wi_new = as<mat>(GlassoResult["wi"]);
            errflag = as<int>(GlassoResult["errflag"]);
            niter = as<int>(GlassoResult["niter"]);

            if (errflag != 0) {
                Rcpp::stop("glassoFast returned errflag=%d for cluster %d", errflag, k + 1);
            }
            if (!(W_new.is_finite() && Wi_new.is_finite())) {
                Rcpp::stop("glassoFast returned non-finite matrices for cluster %d", k + 1);
            }
            if (niter < 1) {
                Rcpp::stop("glassoFast returned invalid niter=%d for cluster %d",
                           niter, k + 1);
            }
            if (niter >= kGlassoFastMaxIterations) {
                Rcpp::stop(
                    "glassoFast reached the maximum iteration cap (%d) "
                    "for cluster %d; the precision update is unscorable",
                    kGlassoFastMaxIterations, k + 1
                );
            }
        }

        const MatrixPairDiagnostics pair_diagnostics =
            require_valid_covariance_precision(W_new, Wi_new, k, solver);
        const GlassoKKTDiagnostics kkt_diagnostics =
            glasso_kkt_diagnostics(
                W_new, Wi_new, empirical, rho_k_matrix
            );
        if (!std::isfinite(kkt_diagnostics.normalized_residual) ||
            kkt_diagnostics.normalized_residual >
                kkt_diagnostics.tolerance) {
            Rcpp::stop(
                "%s failed the graphical-lasso KKT check for component %d "
                "(normalized residual=%.8g, tolerance=%.8g)",
                solver.c_str(), k + 1,
                kkt_diagnostics.normalized_residual,
                kkt_diagnostics.tolerance
            );
        }

        CovarianceMatrix.slice(k) = W_new;
        PrecisionMatrix.slice(k) = Wi_new;
        LastGlassoEmpiricalCovariance.slice(k) = empirical;
        LastGlassoPenalty.slice(k) = rho_k_matrix;
        LastGlassoNIterations(k) = niter;
        LastGlassoInverseResiduals(k) =
            pair_diagnostics.inverse_residual;
        LastGlassoNormalizedInverseResiduals(k) =
            pair_diagnostics.normalized_inverse_residual;
        LastGlassoInverseTolerances(k) = pair_diagnostics.tolerance;
        LastGlassoKKTResiduals(k) = kkt_diagnostics.residual;
        LastGlassoNormalizedKKTResiduals(k) =
            kkt_diagnostics.normalized_residual;
        LastGlassoKKTTolerances(k) = kkt_diagnostics.tolerance;
        LastGlassoSolvers[k] = solver;
        LastGlassoStarts[k] = start;
        glassoEMIterations.push_back(em_iteration + 1);
        glassoComponents.push_back(k + 1);
        glassoNIterations.push_back(niter);
        glassoErrorFlags.push_back(errflag);
        glassoSolvers.push_back(solver);
        glassoStarts.push_back(start);
        glassoSeedSources.push_back(seed_source);
        glassoInverseResiduals.push_back(
            pair_diagnostics.inverse_residual
        );
        glassoNormalizedInverseResiduals.push_back(
            pair_diagnostics.normalized_inverse_residual
        );
        glassoInverseTolerances.push_back(pair_diagnostics.tolerance);
        glassoKKTResiduals.push_back(kkt_diagnostics.residual);
        glassoNormalizedKKTResiduals.push_back(
            kkt_diagnostics.normalized_residual
        );
        glassoKKTTolerances.push_back(kkt_diagnostics.tolerance);
    } // End loop k
}

// A variable is active when at least one component mean is nonzero.
rowvec Mixture::VarRole(void){
    rowvec MuSum = zeros<rowvec>(p);
    rowvec alive = ones<rowvec>(p);

    MuSum = trans(sum(abs(Mu), 1));

    for(int j = 0; j < p; ++j) {
        if(MuSum(j) < std::numeric_limits<double>::epsilon()) {
            alive(j) = 0;
        }
    }
    return alive;
}

Rcpp::List Mixture::FitState(void) const {
    Rcpp::NumericVector prop_out(prop.begin(), prop.end());
    return Rcpp::List::create(
        Rcpp::Named("prop") = prop_out,
        Rcpp::Named("Mu") = Rcpp::wrap(Mu),
        Rcpp::Named("SigmaCube") = Rcpp::wrap(CovarianceMatrix),
        Rcpp::Named("OmegaCube") = Rcpp::wrap(PrecisionMatrix),
        Rcpp::Named("Z") = Rcpp::wrap(ProbCond)
    );
}

Rcpp::List Mixture::GlassoDiagnostics(void) const {
    const int n_calls = static_cast<int>(glassoSolvers.size());
    Rcpp::List calls(n_calls);
    int cold_calls = 0;
    int warm_calls = 0;
    int analytic_calls = 0;
    int total_iterations = 0;
    bool inverse_validated = true;
    bool kkt_validated = true;
    double maximum_normalized_inverse_residual = 0.0;
    double maximum_normalized_kkt_residual = 0.0;

    for (int index = 0; index < n_calls; ++index) {
        if (glassoSolvers[index] == "analytic_diagonal") {
            analytic_calls++;
        } else if (glassoStarts[index] == "warm") {
            warm_calls++;
        } else if (glassoStarts[index] == "cold") {
            cold_calls++;
        }
        total_iterations += glassoNIterations[index];
        inverse_validated = inverse_validated &&
            glassoNormalizedInverseResiduals[index] <=
                glassoInverseTolerances[index];
        kkt_validated = kkt_validated &&
            glassoNormalizedKKTResiduals[index] <=
                glassoKKTTolerances[index];
        maximum_normalized_inverse_residual = std::max(
            maximum_normalized_inverse_residual,
            glassoNormalizedInverseResiduals[index]
        );
        maximum_normalized_kkt_residual = std::max(
            maximum_normalized_kkt_residual,
            glassoNormalizedKKTResiduals[index]
        );
        calls[index] = Rcpp::List::create(
            Rcpp::Named("em_iteration") = glassoEMIterations[index],
            Rcpp::Named("component") = glassoComponents[index],
            Rcpp::Named("solver") = glassoSolvers[index],
            Rcpp::Named("start") = glassoStarts[index],
            Rcpp::Named("seed_source") = glassoSeedSources[index],
            Rcpp::Named("niter") = glassoNIterations[index],
            Rcpp::Named("errflag") = glassoErrorFlags[index],
            Rcpp::Named("solver_tolerance") = kGlassoFastTolerance,
            Rcpp::Named("inverse_residual") =
                glassoInverseResiduals[index],
            Rcpp::Named("normalized_inverse_residual") =
                glassoNormalizedInverseResiduals[index],
            Rcpp::Named("inverse_tolerance") =
                glassoInverseTolerances[index],
            Rcpp::Named("kkt_residual") = glassoKKTResiduals[index],
            Rcpp::Named("normalized_kkt_residual") =
                glassoNormalizedKKTResiduals[index],
            Rcpp::Named("kkt_tolerance") =
                glassoKKTTolerances[index]
        );
    }

    Rcpp::List final_inputs(nbClust);
    for (int k = 0; k < nbClust; ++k) {
        final_inputs[k] = Rcpp::List::create(
            Rcpp::Named("component") = k + 1,
            Rcpp::Named("empirical_covariance") =
                Rcpp::wrap(LastGlassoEmpiricalCovariance.slice(k)),
            Rcpp::Named("penalty") = Rcpp::wrap(LastGlassoPenalty.slice(k)),
            Rcpp::Named("covariance") =
                Rcpp::wrap(CovarianceMatrix.slice(k)),
            Rcpp::Named("precision") =
                Rcpp::wrap(PrecisionMatrix.slice(k)),
            Rcpp::Named("solver") = LastGlassoSolvers[k],
            Rcpp::Named("start") = LastGlassoStarts[k],
            Rcpp::Named("niter") = LastGlassoNIterations(k),
            Rcpp::Named("solver_tolerance") = kGlassoFastTolerance,
            Rcpp::Named("inverse_residual") =
                LastGlassoInverseResiduals(k),
            Rcpp::Named("normalized_inverse_residual") =
                LastGlassoNormalizedInverseResiduals(k),
            Rcpp::Named("inverse_tolerance") =
                LastGlassoInverseTolerances(k),
            Rcpp::Named("kkt_residual") = LastGlassoKKTResiduals(k),
            Rcpp::Named("normalized_kkt_residual") =
                LastGlassoNormalizedKKTResiduals(k),
            Rcpp::Named("kkt_tolerance") = LastGlassoKKTTolerances(k)
        );
    }

    return Rcpp::List::create(
        Rcpp::Named("calls") = calls,
        Rcpp::Named("final_inputs") = final_inputs,
        Rcpp::Named("n_calls") = n_calls,
        Rcpp::Named("glassoFast_calls") = cold_calls + warm_calls,
        Rcpp::Named("cold_calls") = cold_calls,
        Rcpp::Named("warm_calls") = warm_calls,
        Rcpp::Named("analytic_diagonal_calls") = analytic_calls,
        Rcpp::Named("total_iterations") = total_iterations,
        Rcpp::Named("solver_tolerance") = kGlassoFastTolerance,
        Rcpp::Named("inverse_validated") = inverse_validated,
        Rcpp::Named("kkt_validated") = kkt_validated,
        Rcpp::Named("maximum_normalized_inverse_residual") =
            maximum_normalized_inverse_residual,
        Rcpp::Named("maximum_normalized_kkt_residual") =
            maximum_normalized_kkt_residual
    );
}

Rcpp::List Mixture::MeanDiagnostics(void) const {
    const int n_calls = static_cast<int>(meanComponents.size());
    Rcpp::List calls(n_calls);
    bool kkt_validated = true;
    bool block_nondecreasing = true;
    int total_sweeps = 0;
    double maximum_residual = 0.0;
    double minimum_objective_change =
        std::numeric_limits<double>::infinity();

    for (int index = 0; index < n_calls; ++index) {
        kkt_validated = kkt_validated &&
            meanKKTResiduals[index] <= meanKKTTolerances[index];
        block_nondecreasing = block_nondecreasing &&
            meanBlockObjectiveChanges[index] >=
                -meanBlockObjectiveTolerances[index];
        total_sweeps += meanCoordinateSweeps[index];
        maximum_residual = std::max(
            maximum_residual, meanKKTResiduals[index]
        );
        minimum_objective_change = std::min(
            minimum_objective_change, meanBlockObjectiveChanges[index]
        );
        calls[index] = Rcpp::List::create(
            Rcpp::Named("em_iteration") = meanEMIterations[index],
            Rcpp::Named("component") = meanComponents[index],
            Rcpp::Named("coordinate_sweeps") =
                meanCoordinateSweeps[index],
            Rcpp::Named("kkt_residual") = meanKKTResiduals[index],
            Rcpp::Named("kkt_tolerance") = meanKKTTolerances[index],
            Rcpp::Named("block_objective_change") =
                meanBlockObjectiveChanges[index],
            Rcpp::Named("block_objective_tolerance") =
                meanBlockObjectiveTolerances[index]
        );
    }
    if (n_calls == 0) minimum_objective_change = 0.0;

    return Rcpp::List::create(
        Rcpp::Named("calls") = calls,
        Rcpp::Named("n_calls") = n_calls,
        Rcpp::Named("kkt_validated") = kkt_validated,
        Rcpp::Named("block_nondecreasing") = block_nondecreasing,
        Rcpp::Named("maximum_kkt_residual") = maximum_residual,
        Rcpp::Named("minimum_block_objective_change") =
            minimum_objective_change,
        Rcpp::Named("total_coordinate_sweeps") = total_sweeps,
        Rcpp::Named("maximum_coordinate_sweeps") =
            kMeanMaximumCoordinateSweeps,
        Rcpp::Named("relative_kkt_tolerance") =
            kMeanKKTRelativeTolerance
    );
}

Rcpp::List Mixture::ComponentDiagnostics(void) const {
    Rcpp::NumericVector effective_size(nbClust);
    Rcpp::NumericVector minimum_size(nbClust);
    Rcpp::LogicalVector fully_regularized(nbClust);
    for (int component = 0; component < nbClust; ++component) {
        effective_size[component] = effective_component_size(component);
        minimum_size[component] =
            minimum_effective_component_size(component);
        fully_regularized[component] =
            component_fully_regularized(component);
    }
    return Rcpp::List::create(
        Rcpp::Named("effective_sample_size") = effective_size,
        Rcpp::Named("minimum_exclusive") = minimum_size,
        Rcpp::Named("all_off_diagonal_directions_penalized") =
            fully_regularized,
        Rcpp::Named("dimension") = p,
        Rcpp::Named("valid") = ComponentSupportOK()
    );
}

//[[Rcpp::export]]
IntegerVector rcppClusteringEMGlassoWeighted(
    List InputList, double l, double r, arma::cube Pk_in,
    double tol = 1e-3, int max_iter = 250,
    bool inner_warm_start = false, bool predecessor_fitted = false
){
    if (!std::isfinite(l) || l < 0.0 || !std::isfinite(r) || r < 0.0) {
        Rcpp::stop("lambda and rho must be finite and non-negative");
    }
    if (!std::isfinite(tol) || tol <= 0.0) {
        Rcpp::stop("tol must be finite and strictly positive");
    }
    if (max_iter < 0) {
        Rcpp::stop("max_iter must be non-negative");
    }
    if (!Pk_in.is_finite() || arma::any(arma::vectorise(Pk_in) < 0.0)) {
        Rcpp::stop("Pk_in must contain finite, non-negative weights");
    }

    Mixture MyMixture(
        InputList, l, r, Pk_in, inner_warm_start, predecessor_fitted
    );

    double PenLogLik_1 = 0.0, PenLogLik_0 = -std::numeric_limits<double>::infinity();
    int itr = 0;
    double relative_diff = std::numeric_limits<double>::infinity();
    bool converged = false;
    bool monotone_non_decreasing = true;
    bool no_material_decrease = true;
    double largest_decrease = 0.0;
    std::vector<double> objective_trace;

    PenLogLik_1 = MyMixture.PenLogLik();
    if (!std::isfinite(PenLogLik_1)) {
        Rcpp::stop("Initial penalized log-likelihood is non-finite (lambda=%.4f, rho=%.4f)", l, r);
    }
    const double initial_penloglik = PenLogLik_1;
    objective_trace.push_back(PenLogLik_1);

    while(itr < max_iter)
    {
        PenLogLik_0 = PenLogLik_1;

        MyMixture.GetProbCond();
        if (!MyMixture.ProbOK()) {
            Rcpp::stop("Non-finite responsibilities at iteration %d (lambda=%.4f, rho=%.4f)", itr + 1, l, r);
        }
        MyMixture.GetClassesSizes();
        if (!MyMixture.PropOK()) {
            Rcpp::stop("Invalid component proportions at iteration %d (lambda=%.4f, rho=%.4f)", itr + 1, l, r);
        }

        MyMixture.UpdateMeans(itr);
        if (!MyMixture.MuOK()) {
            Rcpp::stop("Non-finite means at iteration %d (lambda=%.4f, rho=%.4f)", itr + 1, l, r);
        }

        MyMixture.GetEmpiricalCovariance();
        if (!MyMixture.EmpCovOK()) {
            Rcpp::stop("Non-finite empirical covariance at iteration %d (lambda=%.4f, rho=%.4f)", itr + 1, l, r);
        }

        MyMixture.UpdateCovarianceMatrices(itr);
        if (!(MyMixture.OmegaOK() && MyMixture.SigmaOK())) {
            Rcpp::stop("Invalid covariance or precision matrix at iteration %d (lambda=%.4f, rho=%.4f)", itr + 1, l, r);
        }

        PenLogLik_1 = MyMixture.PenLogLik();
        if (!std::isfinite(PenLogLik_1)) {
            Rcpp::stop("Penalized log-likelihood is non-finite at iteration %d (lambda=%.4f, rho=%.4f)", itr + 1, l, r);
        }
        objective_trace.push_back(PenLogLik_1);

        const double objective_change = PenLogLik_1 - PenLogLik_0;
        const double decrease_tolerance = 1e-8 * (1.0 + std::abs(PenLogLik_0));
        if (objective_change < 0.0) {
            monotone_non_decreasing = false;
            largest_decrease = std::min(largest_decrease, objective_change);
        }
        if (objective_change < -decrease_tolerance) {
            no_material_decrease = false;
            Rcpp::stop(
                "Penalized objective decreased materially at iteration %d "
                "(change=%.8g, lambda=%.4f, rho=%.4f)",
                itr + 1, objective_change, l, r
            );
        }

        if (std::abs(PenLogLik_0) < tol || !std::isfinite(PenLogLik_0)) {
             relative_diff = (std::abs(PenLogLik_1 - PenLogLik_0) < tol) ? 0.0 : 1.0;
        } else {
             relative_diff = std::abs(PenLogLik_1 - PenLogLik_0) / (1.0 + std::abs(PenLogLik_0));
        }

        itr++;
        if (relative_diff <= tol) {
            converged = true;
            break;
        }
    }

    // Synchronize the returned responsibilities with the final parameters.
    // This does not modify the parameters or the activity decision.
    MyMixture.GetProbCond();
    if (!MyMixture.ProbOK()) {
        Rcpp::stop("Final responsibilities are non-finite (lambda=%.4f, rho=%.4f)", l, r);
    }
    if (!MyMixture.ComponentMassOK()) {
        Rcpp::stop("A final component has zero responsibility mass (lambda=%.4f, rho=%.4f)", l, r);
    }
    if (!MyMixture.ComponentSupportOK()) {
        Rcpp::stop(
            "A final component has insufficient dimension-aware covariance "
            "support (lambda=%.4f, rho=%.4f)", l, r
        );
    }

    const bool state_valid = MyMixture.PropOK() && MyMixture.ProbOK() &&
        MyMixture.ComponentMassOK() && MyMixture.ComponentSupportOK() &&
        MyMixture.MuOK() &&
        MyMixture.SigmaOK() && MyMixture.OmegaOK();
    List glasso_diagnostics = MyMixture.GlassoDiagnostics();
    List mean_diagnostics = MyMixture.MeanDiagnostics();
    List component_diagnostics = MyMixture.ComponentDiagnostics();
    const bool glasso_validated =
        as<bool>(glasso_diagnostics["inverse_validated"]) &&
        as<bool>(glasso_diagnostics["kkt_validated"]);
    const bool mean_kkt_validated =
        as<bool>(mean_diagnostics["kkt_validated"]) &&
        as<bool>(mean_diagnostics["block_nondecreasing"]);
    const bool scorable = converged && no_material_decrease && state_valid &&
        glasso_validated && mean_kkt_validated;

    IntegerVector activity(MyMixture.getDimP(), NA_INTEGER);
    if (scorable) {
        rowvec role = MyMixture.VarRole();
        for (int j = 0; j < MyMixture.getDimP(); ++j) {
            activity[j] = static_cast<int>(role(j));
        }
    }

    NumericVector trace(objective_trace.begin(), objective_trace.end());
    List status = List::create(
        Named("converged") = converged,
        Named("scorable") = scorable,
        Named("reason") = converged ? "tolerance" : "max_iterations",
        Named("iterations") = itr,
        Named("max_iter") = max_iter,
        Named("tolerance") = tol,
        Named("relative_difference") = relative_diff,
        Named("initial_penalized_loglik") = initial_penloglik,
        Named("final_penalized_loglik") = PenLogLik_1,
        Named("monotone_non_decreasing") = monotone_non_decreasing,
        Named("no_material_decrease") = no_material_decrease,
        Named("largest_decrease") = largest_decrease,
        Named("state_valid") = state_valid,
        Named("glasso_validated") = glasso_validated,
        Named("glasso_inverse_validated") =
            glasso_diagnostics["inverse_validated"],
        Named("glasso_kkt_validated") =
            glasso_diagnostics["kkt_validated"],
        Named("glasso_maximum_normalized_inverse_residual") =
            glasso_diagnostics["maximum_normalized_inverse_residual"],
        Named("glasso_maximum_normalized_kkt_residual") =
            glasso_diagnostics["maximum_normalized_kkt_residual"],
        Named("mean_kkt_validated") = mean_kkt_validated,
        Named("mean_maximum_kkt_residual") =
            mean_diagnostics["maximum_kkt_residual"],
        Named("mean_minimum_block_objective_change") =
            mean_diagnostics["minimum_block_objective_change"],
        Named("mean_coordinate_sweeps_total") =
            mean_diagnostics["total_coordinate_sweeps"],
        Named("component_support_validated") =
            component_diagnostics["valid"],
        Named("glasso_calls") = glasso_diagnostics["n_calls"],
        Named("glassoFast_calls") = glasso_diagnostics["glassoFast_calls"],
        Named("glasso_cold_calls") = glasso_diagnostics["cold_calls"],
        Named("glasso_warm_calls") = glasso_diagnostics["warm_calls"],
        Named("analytic_diagonal_calls") =
            glasso_diagnostics["analytic_diagonal_calls"],
        Named("glasso_iterations_total") =
            glasso_diagnostics["total_iterations"]
    );

    activity.attr("fit_status") = status;
    activity.attr("fit_state") = MyMixture.FitState();
    activity.attr("objective_trace") = trace;
    activity.attr("glasso_diagnostics") = glasso_diagnostics;
    activity.attr("mean_diagnostics") = mean_diagnostics;
    activity.attr("component_diagnostics") = component_diagnostics;
    std::string glasso_start_label = "cold";
    if (predecessor_fitted && inner_warm_start) {
        glasso_start_label = "predecessor_then_previous_em_warm";
    } else if (predecessor_fitted) {
        glasso_start_label = "predecessor_warm_then_cold";
    } else if (inner_warm_start) {
        glasso_start_label = "cold_then_previous_em_warm";
    }
    activity.attr("fit_metadata") = List::create(
        Named("state_ownership") = "deep_copy",
        Named("glasso_start") = glasso_start_label,
        Named("warm_start") = inner_warm_start || predecessor_fitted,
        Named("inner_warm_start") = inner_warm_start,
        Named("predecessor_fitted") = predecessor_fitted,
        Named("mean_solver") =
            "gauss_seidel_coordinate_descent_to_subgradient_kkt",
        Named("mean_relative_kkt_tolerance") =
            kMeanKKTRelativeTolerance,
        Named("glasso_solver_tolerance") = kGlassoFastTolerance,
        Named("glasso_validation") =
            "spd_plus_normwise_inverse_backward_error_plus_kkt",
        Named("component_support_rule") =
            "p_effective_unpenalized_or_positive_variance_fully_penalized"
    );
    return activity;
};

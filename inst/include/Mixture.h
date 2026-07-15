#ifndef MIXTURE_H
#define MIXTURE_H
#include <RcppArmadillo.h>
#include <limits> 
#include <string>
#include <vector>

class Mixture
{
private:
    // Data and dimensions.
    arma::mat        Xd;                // n x p data matrix
    int              n;                 // number of observations
    int              p;                 // number of variables
    int              nbClust;           // number of mixture components

    // Mixture parameters.
    arma::rowvec     prop;              // length-K mixing proportions
    arma::mat        Mu;                // p x K component means
    arma::cube       CovarianceMatrix;  // p x p x K covariance matrices
    arma::cube       PrecisionMatrix;   // p x p x K precision matrices

    // EM working state.
    arma::cube       EmpiricalCovariance; // weighted component covariances
    arma::cube       Pk_cube;             // adaptive precision penalties
    arma::mat        ProbCond;            // posterior responsibilities

    // Graphical-lasso continuation and diagnostics.
    bool             innerGlassoWarmStart;
    bool             predecessorFitted;
    std::vector<int> glassoEMIterations;
    std::vector<int> glassoComponents;
    std::vector<int> glassoNIterations;
    std::vector<int> glassoErrorFlags;
    std::vector<std::string> glassoSolvers;
    std::vector<std::string> glassoStarts;
    std::vector<std::string> glassoSeedSources;
    std::vector<double> glassoInverseResiduals;
    std::vector<double> glassoNormalizedInverseResiduals;
    std::vector<double> glassoInverseTolerances;
    std::vector<double> glassoKKTResiduals;
    std::vector<double> glassoNormalizedKKTResiduals;
    std::vector<double> glassoKKTTolerances;
    arma::cube       LastGlassoEmpiricalCovariance;
    arma::cube       LastGlassoPenalty;
    arma::ivec       LastGlassoNIterations;
    arma::vec        LastGlassoInverseResiduals;
    arma::vec        LastGlassoNormalizedInverseResiduals;
    arma::vec        LastGlassoInverseTolerances;
    arma::vec        LastGlassoKKTResiduals;
    arma::vec        LastGlassoNormalizedKKTResiduals;
    arma::vec        LastGlassoKKTTolerances;
    std::vector<std::string> LastGlassoSolvers;
    std::vector<std::string> LastGlassoStarts;

    // Coordinate-descent diagnostics for the penalized mean block.
    std::vector<int> meanEMIterations;
    std::vector<int> meanComponents;
    std::vector<int> meanCoordinateSweeps;
    std::vector<double> meanKKTResiduals;
    std::vector<double> meanKKTTolerances;
    std::vector<double> meanBlockObjectiveChanges;
    std::vector<double> meanBlockObjectiveTolerances;

    // Penalty parameters.
    double           lambda;            // Lambda_mu  (mean penalty)
    double           rho;               // Lambda_Omega  (precision penalty)

    // Cholesky log determinant; returns -Inf when M is invalid.
    double safe_log_det(const arma::mat& M, bool& ok) const;
    double effective_component_size(int component) const;
    bool component_fully_regularized(int component) const;
    double minimum_effective_component_size(int component) const;
public:
    // Common-shrinkage formulation of Zhou, Pan, and Shen (2009).
    explicit Mixture(Rcpp::List InputList);

    // Adaptive formulation with component-specific penalty matrices.
    Mixture(Rcpp::List    InputList,
            double        lambda_mu,
            double        lambda_omega,
            arma::cube    Pk_in,
            bool          inner_glasso_warm_start = false,
            bool          predecessor_fitted = false);

    ~Mixture(){}

    // ECM updates.
    double PenLogLik( void );
    void   GetProbCond( void );
    void   GetClassesSizes( void );
    void   UpdateMeans( int em_iteration );
    void   GetEmpiricalCovariance( void );
    void   UpdateCovarianceMatrices( int em_iteration );

    // Variable activity and fitted-state diagnostics.
    arma::rowvec VarRole( void );
    Rcpp::List FitState( void ) const;
    Rcpp::List GlassoDiagnostics( void ) const;
    Rcpp::List MeanDiagnostics( void ) const;
    Rcpp::List ComponentDiagnostics( void ) const;

    // State checks used by the ECM driver.
    inline int  getDimP()      const { return p; }
    inline bool ProbOK()       const {
        return ProbCond.is_finite() &&
            arma::all(arma::vectorise(ProbCond) >= 0.0) &&
            arma::all(arma::abs(arma::sum(ProbCond, 1) - 1.0) <= 1e-8);
    }
    inline bool ComponentMassOK() const {
        const arma::rowvec masses = arma::sum(ProbCond, 0);
        return masses.is_finite() && arma::all(masses > 0.0);
    }
    bool ComponentSupportOK() const;
    inline bool MuOK()         const { return Mu.is_finite(); }
    inline bool SigmaOK()      const { return CovarianceMatrix.is_finite(); }
    inline bool OmegaOK()      const { return PrecisionMatrix.is_finite(); }
    inline bool EmpCovOK()     const { return EmpiricalCovariance.is_finite(); }
    inline bool PropOK()       const {
        return prop.is_finite() && arma::all(prop > 0.0) &&
            std::abs(arma::accu(prop) - 1.0) <= 1e-8;
    }

    inline int  GetNbClust()   const { return nbClust; }
};

#endif /* MIXTURE_H */

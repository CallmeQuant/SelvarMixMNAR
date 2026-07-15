#include "select.hpp"

namespace {

DataFrame empty_missing_values() {
    return DataFrame::create(
        Named("row") = IntegerVector(0),
        Named("col") = IntegerVector(0),
        Named("value") = NumericVector(0)
    );
}

void validate_stopping_rule(const std::string& stoppingRule, int threshold) {
    if (threshold < 1)
        stop("The stopping threshold must be a positive integer.");
    if (stoppingRule != "consecutive" && stoppingRule != "legacy_block")
        stop("Unknown stopping rule '%s'. Use 'consecutive' or 'legacy_block'.",
             stoppingRule.c_str());
}

void validate_model_result(const List& fit, const std::string& context) {
    if (!fit.containsElementNamed("error"))
        stop("%s did not return an error-status field.", context.c_str());

    const std::string fitError = as<std::string>(fit["error"]);
    if (fitError != "No error") {
        const std::string errorClass = fit.containsElementNamed("errorClass") &&
                !Rf_isNull(fit["errorClass"])
            ? as<std::string>(fit["errorClass"])
            : "selvarmix_sruw_backend_fit_error";
        throw SRUWStatusException(
            context + " failed: " + fitError, errorClass
        );
    }

    const char* required[] = {
        "criterionValue", "model", "criterion", "nbcluster",
        "parameters", "proba", "partition"
    };
    for (const char* name : required) {
        if (!fit.containsElementNamed(name))
            stop("%s did not return required field '%s'.", context.c_str(), name);
    }

    const double criterionValue = as<double>(fit["criterionValue"]);
    if (!std::isfinite(criterionValue))
        stop("%s returned a non-finite criterion value.", context.c_str());
}

} // namespace

Select::Select()
    : packSize(1), stoppingRule("consecutive"), lastWEvaluated(0),
      lastWStopReason("not_run") {}

Select::Select(Vect v, CritClust b, SelectReg sReg, int packSize,
               string stoppingRule)
{
    validate_stopping_rule(stoppingRule, packSize);
    this->v = v;
    this->b = b;
    this->sReg = sReg;
    this->packSize = packSize;
    this->stoppingRule = stoppingRule;
    this->lastWEvaluated = 0;
    this->lastWStopReason = "not_run";
}

Select::Select(Vect v, SelectReg sReg, int packSize,
               string stoppingRule)
{
    validate_stopping_rule(stoppingRule, packSize);
    this->v = v;
    this->sReg = sReg;
    this->packSize = packSize;
    this->stoppingRule = stoppingRule;
    this->lastWEvaluated = 0;
    this->lastWStopReason = "not_run";
}

List Select::selectS(std::vector<int> Order) {
    if (Order.empty())
        stop("The SRUW ordering must contain at least one variable.");

    int InitialProjectsNb = this->v.experiments.size();
    const int regressionModel = 1;
    std::vector<int> varSelectClust;
    varSelectClust.push_back(Order.at(0));

    List currentFit = b.ClustBestModel(varSelectClust);
    validate_model_result(currentFit, "Initial clustering fit");
    double currentCriterion = as<double>(currentFit["criterionValue"]);
    DataFrame missingVals = currentFit.containsElementNamed("missingValues")
        ? as<DataFrame>(currentFit["missingValues"])
        : empty_missing_values();

    int nEvaluated = 0;
    std::string stopReason = "order_exhausted";

    auto evaluateCandidate = [&](int orderIndex) -> bool {
        ++nEvaluated;
        std::vector<int> candidate;
        candidate.push_back(Order.at(orderIndex));
        std::vector<int> regressionVariables =
            this->sReg.selectReg(varSelectClust, candidate,
                                 InitialProjectsNb);
        std::vector<int> augmented =
            this->v.ajouter_var(varSelectClust, candidate);

        List candidateFit = b.ClustBestModel(augmented);
        validate_model_result(candidateFit, "Candidate clustering fit");

        List regressionFit =
            this->v.bicReggen(candidate, regressionVariables, regressionModel);
        if (!regressionFit.containsElementNamed("bicvalue"))
            stop("Regression comparison did not return 'bicvalue'.");

        const double candidateCriterion =
            as<double>(candidateFit["criterionValue"]);
        const double regressionCriterion =
            as<double>(regressionFit["bicvalue"]);
        const double criterionDifference =
            candidateCriterion - currentCriterion - regressionCriterion;
        if (!std::isfinite(regressionCriterion) ||
            !std::isfinite(criterionDifference))
            stop("SRUW candidate comparison returned a non-finite criterion.");

        const bool improvement = criterionDifference > 0.0;
        if (improvement) {
            varSelectClust = augmented;
            currentFit = candidateFit;
            currentCriterion = candidateCriterion;
            missingVals = candidateFit.containsElementNamed("missingValues")
                ? as<DataFrame>(candidateFit["missingValues"])
                : empty_missing_values();
        }
        return improvement;
    };

    if (stoppingRule == "consecutive") {
        ConsecutiveStopper stopper(packSize);
        for (int idx = 1; idx < static_cast<int>(Order.size()); ++idx) {
            const bool improvement = evaluateCandidate(idx);
            if (stopper.observe(improvement)) {
                stopReason = "consecutive_nonpositive";
                break;
            }
        }
    } else {
        int firstIndex = 1;
        while (firstIndex < static_cast<int>(Order.size())) {
            const int lastIndex = std::min(
                firstIndex + packSize, static_cast<int>(Order.size()));
            int blockImprovements = 0;
            for (int idx = firstIndex; idx < lastIndex; ++idx)
                blockImprovements += evaluateCandidate(idx) ? 1 : 0;
            if (blockImprovements == 0) {
                stopReason = "legacy_empty_block";
                break;
            }
            firstIndex = lastIndex;
        }
    }

    return List::create(
        Named("S") = wrap(varSelectClust),
        Named("model") = currentFit["model"],
        Named("criterionValue") = currentFit["criterionValue"],
        Named("criterion") = currentFit["criterion"],
        Named("nbcluster") = currentFit["nbcluster"],
        Named("parameters") = currentFit["parameters"],
        Named("proba") = currentFit["proba"],
        Named("partition") = currentFit["partition"],
        Named("missingValues") = missingVals,
        Named("framework") = currentFit.containsElementNamed("framework") ? currentFit["framework"] : R_NilValue,
        Named("requestedModel") = currentFit.containsElementNamed("requestedModel") ? currentFit["requestedModel"] : R_NilValue,
        Named("effectiveModel") = currentFit.containsElementNamed("effectiveModel") ? currentFit["effectiveModel"] : R_NilValue,
        Named("stoppingRule") = stoppingRule,
        Named("stoppingThreshold") = packSize,
        Named("nEvaluated") = nEvaluated,
        Named("stopReason") = stopReason
    );
}

vector<int> Select::selectW(vector<int> Order, vector<int> OtherVar)
{
    vector<int> varIndep;
    lastWEvaluated = 0;
    lastWStopReason = "order_exhausted";
    if (Order.empty())
        return varIndep;

    int InitialProjectsNb = this->v.experiments.size();
    auto evaluateCandidate = [&](int orderIndex) -> bool {
        ++lastWEvaluated;
        vector<int> candidate;
        candidate.push_back(Order.at(orderIndex));
        vector<int> regressionVariables =
            sReg.selectReg(OtherVar, candidate, InitialProjectsNb);
        const bool improvement = regressionVariables.empty();
        if (improvement)
            varIndep.push_back(Order[orderIndex]);
        return improvement;
    };

    if (stoppingRule == "consecutive") {
        ConsecutiveStopper stopper(packSize);
        for (int idx = static_cast<int>(Order.size()) - 1; idx >= 0; --idx) {
            const bool improvement = evaluateCandidate(idx);
            if (stopper.observe(improvement)) {
                lastWStopReason = "consecutive_nonpositive";
                break;
            }
        }
    } else {
        int lastIndex = static_cast<int>(Order.size());
        while (lastIndex > 0) {
            const int firstIndex = std::max(0, lastIndex - packSize);
            int blockImprovements = 0;
            for (int idx = lastIndex - 1; idx >= firstIndex; --idx)
                blockImprovements += evaluateCandidate(idx) ? 1 : 0;
            if (blockImprovements == 0) {
                lastWStopReason = "legacy_empty_block";
                break;
            }
            lastIndex = firstIndex;
        }
    }
    return varIndep;
}

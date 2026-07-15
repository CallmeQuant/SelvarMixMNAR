#include "critClust.hpp"

#include <cctype>

namespace {

bool contains_token(const std::vector<std::string>& tokens,
                    const std::string& value) {
    return std::find(tokens.begin(), tokens.end(), value) != tokens.end();
}

const std::vector<std::string>& mclust_models() {
    static const std::vector<std::string> models = {
        "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE",
        "EVE", "VVE", "EEV", "VEV", "EVV", "VVV"
    };
    return models;
}

const std::vector<std::string>& mixall_models() {
    static const std::vector<std::string> models = {
        "gaussian_pk_sjk", "gaussian_pk_sj", "gaussian_pk_sk",
        "gaussian_pk_s", "gaussian_p_sjk", "gaussian_p_sj",
        "gaussian_p_sk", "gaussian_p_s"
    };
    return models;
}

std::string compact_token(const std::string& token) {
    std::string compact;
    compact.reserve(token.size());
    for (unsigned char value : token) {
        if (!std::isspace(value))
            compact.push_back(static_cast<char>(value));
    }
    return compact;
}

bool rmixmod_family_from_token(const std::string& token,
                               std::string& family) {
    const std::string compact = compact_token(token);
    if (compact == "mixmodGaussianModel()") {
        family = "general";
        return true;
    }

    const std::vector<std::string> families = {
        "general", "diagonal", "spherical", "all"
    };
    for (const std::string& candidate : families) {
        const std::string doubleQuoted =
            "mixmodGaussianModel(family=\"" + candidate + "\")";
        const std::string singleQuoted =
            "mixmodGaussianModel(family='" + candidate + "')";
        if (compact == doubleQuoted || compact == singleQuoted) {
            family = candidate;
            return true;
        }
    }
    return false;
}

std::string canonical_rmixmod_model(const std::string& family) {
    return "mixmodGaussianModel(family=\"" + family + "\")";
}

DataFrame empty_missing_values() {
    return DataFrame::create(
        Named("row") = IntegerVector(0),
        Named("col") = IntegerVector(0),
        Named("value") = NumericVector(0)
    );
}

bool labels_are_valid(const IntegerVector& labels, int n, int k) {
    if (labels.size() != n)
        return false;
    for (int label : labels) {
        if (IntegerVector::is_na(label) || label < 1 || label > k)
            return false;
    }
    return true;
}

} // namespace

// Store an exact backend/model request for one clustering criterion fit.
CritClust::CritClust() {}

CritClust::CritClust(int k, std::string framework, std::string model_name,
                     NumericMatrix data, std::string crit,
                     IntegerVector knownlabels, bool DA)
{
    this->crit = crit;
    this->framework = framework;
    this->model_name = model_name;
    this->k = k;
    this->data = data;
    this->knownlabels = knownlabels;
    this->DA = DA;
}

List CritClust::ClustBestModel(std::vector<int> numExp)
{
    List defaultResult = List::create(
        Named("criterionValue") = NA_REAL,
        Named("criterion") = crit,
        Named("nbcluster") = k,
        Named("model") = model_name,
        Named("parameters") = R_NilValue,
        Named("proba") = R_NilValue,
        Named("partition") = R_NilValue,
        Named("error") = "Model fitting failed",
        Named("errorClass") = "selvarmix_sruw_backend_fit_error",
        Named("framework") = framework,
        Named("requestedModel") = model_name,
        Named("effectiveModel") = R_NilValue,
        Named("missingValues") = empty_missing_values(),
        Named("S") = IntegerVector(0)
    );

    auto fail = [&](const std::string& message,
                    const std::string& errorClass) -> List {
        defaultResult["error"] = message;
        defaultResult["errorClass"] = errorClass;
        return defaultResult;
    };

    try {
        if (k < 1)
            return fail("The SRUW cluster count must be positive.",
                        "selvarmix_sruw_backend_error");
        if (numExp.empty())
            return fail("The SRUW clustering variable set must be non-empty.",
                        "selvarmix_sruw_backend_error");
        std::vector<bool> seen(data.ncol() + 1, false);
        for (int variable : numExp) {
            if (variable < 1 || variable > data.ncol() || seen[variable])
                return fail("The SRUW clustering variable indices are invalid.",
                            "selvarmix_sruw_backend_error");
            seen[variable] = true;
        }
        if (framework != "Mclust" && framework != "Rmixmod" &&
            framework != "MixAll") {
            return fail(
                "Unsupported SRUW framework token '" + framework +
                    "'. Use Mclust, Rmixmod, or MixAll.",
                "selvarmix_sruw_framework_unavailable"
            );
        }

        Environment base = Environment::namespace_env("base");
        Function dataframe = base["data.frame"];
        NumericMatrix dataAux(data.nrow(), numExp.size());
        for (size_t j = 0; j < numExp.size(); ++j)
            dataAux(_, j) = data(_, numExp[j] - 1);
        DataFrame df = dataframe(dataAux);

        if (framework == "Mclust") {
            if (!contains_token(mclust_models(), model_name)) {
                return fail(
                    "Unsupported Mclust model token '" + model_name + "'.",
                    "selvarmix_sruw_model_unavailable"
                );
            }
            if (crit != "BIC") {
                return fail(
                    "Mclust SRUW backend does not implement requested criterion '" +
                        crit + "'; supported exact token: BIC.",
                    "selvarmix_sruw_criterion_unavailable"
                );
            }
            if (DA) {
                return fail(
                    "Mclust SRUW does not implement supervised fitting; the request was not downgraded to unsupervised fitting.",
                    "selvarmix_sruw_supervision_unavailable"
                );
            }

            // In one dimension, mclust's exact covariance constraint is E or V.
            // Map the requested volume constraint deterministically; do not
            // broaden the request to a search over both univariate models.
            const std::string effectiveModel =
                dataAux.ncol() == 1 ? model_name.substr(0, 1) : model_name;

            Environment mclustEnv = Environment::namespace_env("mclust");
            Function MclustBICFunc = mclustEnv["mclustBIC"];
            Function SummaryMclustBICFunc = mclustEnv["summaryMclustBIC"];
            Function doCall = base["do.call"];

            List mclustArgs;
            mclustArgs["data"] = df;
            mclustArgs["G"] = k;
            mclustArgs["modelNames"] = CharacterVector::create(effectiveModel);
            mclustArgs["verbose"] = false;

            SEXP bicResult = doCall(MclustBICFunc, mclustArgs);
            List summaryArgs = List::create(
                Named("object") = bicResult,
                Named("data") = df,
                Named("G") = k,
                Named("modelNames") = CharacterVector::create(effectiveModel)
            );
            List mclustResult = doCall(SummaryMclustBICFunc, summaryArgs);

            const std::string returnedModel =
                as<std::string>(mclustResult["modelName"]);
            if (returnedModel != effectiveModel) {
                return fail(
                    "Mclust returned model '" + returnedModel +
                        "' instead of the exact effective model '" +
                        effectiveModel + "'.",
                    "selvarmix_sruw_backend_result_error"
                );
            }

            return List::create(
                Named("criterionValue") = as<double>(mclustResult["bic"]),
                Named("criterion") = "BIC",
                Named("nbcluster") = as<int>(mclustResult["G"]),
                Named("model") = returnedModel,
                Named("parameters") = mclustResult["parameters"],
                Named("proba") = mclustResult["z"],
                Named("partition") = mclustResult["classification"],
                Named("error") = "No error",
                Named("errorClass") = R_NilValue,
                Named("framework") = "Mclust",
                Named("requestedModel") = model_name,
                Named("effectiveModel") = effectiveModel,
                Named("missingValues") = empty_missing_values()
            );
        }

        if (framework == "Rmixmod") {
            std::string family;
            if (!rmixmod_family_from_token(model_name, family)) {
                return fail(
                    "Unsupported Rmixmod model token '" + model_name +
                        "'. Only a validated mixmodGaussianModel family token is accepted.",
                    "selvarmix_sruw_model_unavailable"
                );
            }
            const bool validCriterion = DA
                ? (crit == "BIC" || crit == "CV")
                : (crit == "BIC" || crit == "ICL");
            if (!validCriterion) {
                return fail(
                    "Rmixmod SRUW backend does not implement requested criterion '" +
                        crit + "' for the requested supervision mode.",
                    "selvarmix_sruw_criterion_unavailable"
                );
            }
            if (DA && !labels_are_valid(knownlabels, data.nrow(), k)) {
                return fail(
                    "Rmixmod supervised labels must have one value in 1,...,K for every row.",
                    "selvarmix_sruw_supervision_unavailable"
                );
            }

            Environment Rmixmod = Environment::namespace_env("Rmixmod");
            Function RmixmodLearn = Rmixmod["mixmodLearn"];
            Function RmixmodCluster = Rmixmod["mixmodCluster"];
            Function RmixmodStrategy = Rmixmod["mixmodStrategy"];
            Function RmixmodGaussianModel = Rmixmod["mixmodGaussianModel"];
            SEXP modelObject = RmixmodGaussianModel(Named("family") = family);
            const std::string effectiveModel = canonical_rmixmod_model(family);

            SEXP bestResultObject = R_NilValue;
            if (!DA) {
                S4 mixmodstrategy = RmixmodStrategy(
                    Named("nbTry") = 4,
                    Named("nbIterationInAlgo") = 200,
                    Named("nbTryInInit") = 100,
                    Named("nbIterationInInit") = 10,
                    Named("epsilonInAlgo") = 1e-8,
                    Named("initMethod") = "SEMMax"
                );
                S4 xem = RmixmodCluster(
                    Named("data") = df,
                    Named("nbCluster") = k,
                    Named("models") = modelObject,
                    Named("strategy") = mixmodstrategy,
                    Named("criterion") = CharacterVector::create(crit)
                );
                bestResultObject = xem.slot("bestResult");
            } else {
                S4 xem = RmixmodLearn(
                    Named("data") = df,
                    Named("knownLabels") = knownlabels,
                    Named("models") = modelObject,
                    Named("criterion") = CharacterVector::create(crit)
                );
                bestResultObject = xem.slot("bestResult");
            }
            S4 bestResult(bestResultObject);

            if (bestResult.hasSlot("error")) {
                SEXP errorSlot = bestResult.slot("error");
                if (!Rf_isNull(errorSlot)) {
                    const std::string backendError =
                        as<std::string>(errorSlot);
                    if (!backendError.empty() && backendError != "No error") {
                        return fail(
                            "Rmixmod fit failed: " + backendError,
                            "selvarmix_sruw_backend_fit_error"
                        );
                    }
                }
            }
            const std::string returnedCriterion =
                as<std::string>(bestResult.slot("criterion"));
            if (returnedCriterion != crit) {
                return fail(
                    "Rmixmod returned criterion '" + returnedCriterion +
                        "' instead of requested criterion '" + crit + "'.",
                    "selvarmix_sruw_backend_result_error"
                );
            }

            NumericVector criterionValues = bestResult.slot("criterionValue");
            if (criterionValues.size() < 1 ||
                !std::isfinite(criterionValues[0])) {
                return fail(
                    "Rmixmod returned no finite criterion value.",
                    "selvarmix_sruw_backend_result_error"
                );
            }

            return List::create(
                Named("criterionValue") = -criterionValues[0],
                Named("criterion") = returnedCriterion,
                Named("nbcluster") = bestResult.slot("nbCluster"),
                Named("model") = bestResult.slot("model"),
                Named("parameters") = bestResult.slot("parameters"),
                Named("proba") = DA ? R_NilValue : bestResult.slot("proba"),
                Named("partition") = bestResult.slot("partition"),
                Named("error") = "No error",
                Named("errorClass") = R_NilValue,
                Named("framework") = "Rmixmod",
                Named("requestedModel") = model_name,
                Named("effectiveModel") = effectiveModel,
                Named("missingValues") = empty_missing_values()
            );
        }

        if (!contains_token(mixall_models(), model_name)) {
            return fail(
                "Unsupported MixAll model token '" + model_name + "'.",
                "selvarmix_sruw_model_unavailable"
            );
        }
        if (crit != "BIC" && crit != "ICL") {
            return fail(
                "MixAll SRUW backend does not implement requested criterion '" +
                    crit + "'; supported exact tokens: BIC, ICL.",
                "selvarmix_sruw_criterion_unavailable"
            );
        }
        if (DA && !labels_are_valid(knownlabels, data.nrow(), k)) {
            return fail(
                "MixAll supervised labels must have one value in 1,...,K for every row.",
                "selvarmix_sruw_supervision_unavailable"
            );
        }

        Environment MixAll = Environment::namespace_env("MixAll");
        Function clusterDiagGaussian = MixAll["clusterDiagGaussian"];
        Function clusterStrategy = MixAll["clusterStrategy"];
        Function learnDiagGaussian = MixAll["learnDiagGaussian"];
        Function missingValuesFunc = MixAll["missingValues"];

        S4 strategy = clusterStrategy(
            Named("nbTry") = 1,
            Named("nbInit") = 100,
            Named("initMethod") = "random",
            Named("initAlgo") = "SEM",
            Named("nbInitIteration") = 10,
            Named("initEpsilon") = 1e-4,
            Named("nbShortRun") = 10,
            Named("shortRunAlgo") = "EM",
            Named("nbShortIteration") = 100,
            Named("shortEpsilon") = 1e-7,
            Named("longRunAlgo") = "EM",
            Named("nbLongIteration") = 200,
            Named("longEpsilon") = 1e-8
        );

        S4 xem;
        if (!DA) {
            xem = clusterDiagGaussian(
                Named("data") = df,
                Named("nbCluster") = k,
                Named("models") = model_name,
                Named("strategy") = strategy,
                Named("criterion") = crit,
                Named("nbCore") = 1
            );
        } else {
            xem = learnDiagGaussian(
                Named("data") = df,
                Named("labels") = knownlabels,
                Named("models") = model_name,
                Named("algo") = "simul",
                Named("nbIter") = 100,
                Named("epsilon") = 1e-8,
                Named("criterion") = crit,
                Named("nbCore") = 1
            );
        }

        if (xem.hasSlot("error")) {
            SEXP errorSlot = xem.slot("error");
            if (!Rf_isNull(errorSlot)) {
                const std::string backendError = as<std::string>(errorSlot);
                if (!backendError.empty() && backendError != "No error") {
                    return fail(
                        "MixAll fit failed: " + backendError,
                        "selvarmix_sruw_backend_fit_error"
                    );
                }
            }
        }
        const std::string returnedCriterion =
            as<std::string>(xem.slot("criterionName"));
        if (returnedCriterion != crit) {
            return fail(
                "MixAll returned criterion '" + returnedCriterion +
                    "' instead of requested criterion '" + crit + "'.",
                "selvarmix_sruw_backend_result_error"
            );
        }

        NumericMatrix imputedData = missingValuesFunc(xem);
        if (imputedData.ncol() < 3) {
            return fail(
                "MixAll returned an invalid missing-value table.",
                "selvarmix_sruw_backend_result_error"
            );
        }
        DataFrame missingVals = DataFrame::create(
            Named("row") = imputedData(_, 0),
            Named("col") = imputedData(_, 1),
            Named("value") = imputedData(_, 2)
        );

        const double criterionValue = -as<double>(xem.slot("criterion"));
        if (!std::isfinite(criterionValue)) {
            return fail(
                "MixAll returned a non-finite criterion value.",
                "selvarmix_sruw_backend_result_error"
            );
        }
        return List::create(
            Named("criterionValue") = criterionValue,
            Named("criterion") = returnedCriterion,
            Named("nbcluster") = xem.slot("nbCluster"),
            Named("model") = model_name,
            Named("parameters") = xem.slot("component"),
            Named("proba") = xem.slot("tik"),
            Named("partition") = xem.slot("zi"),
            Named("error") = "No error",
            Named("errorClass") = R_NilValue,
            Named("framework") = "MixAll",
            Named("requestedModel") = model_name,
            Named("effectiveModel") = model_name,
            Named("missingValues") = missingVals
        );
    }
    catch (std::exception &ex) {
        return fail(
            framework + " SRUW backend raised an exception: " + ex.what(),
            "selvarmix_sruw_backend_fit_error"
        );
    }
    catch (...) {
        return fail(
            framework + " SRUW backend raised an unknown native exception.",
            "selvarmix_sruw_backend_fit_error"
        );
    }
}

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <Rdefines.h>
#include "select.hpp"
#include "selectRegGen.hpp"
#include "critClust.hpp"
#include "sruwSearchControl.hpp"

namespace {

void validate_variable_indices(const std::vector<int>& variables,
                               int nVariables,
                               const std::string& label) {
    std::vector<bool> seen(nVariables + 1, false);
    for (int variable : variables) {
        if (variable < 1 || variable > nVariables)
            stop("%s contains variable index %d outside 1,...,%d.",
                 label.c_str(), variable, nVariables);
        if (seen[variable])
            stop("%s contains duplicate variable index %d.",
                 label.c_str(), variable);
        seen[variable] = true;
    }
}

} // namespace

//[[Rcpp::export]]
List rcppSelectS(NumericMatrix X, std::vector<int> Order, const int nbCluster, 
                 std::string framework, std::string model_name, const int packSize, 
                 std::string Crit, IntegerVector knownlabels, IntegerVector DA,
                 std::string stoppingRule = "consecutive"){
    try {
        if (X.nrow() < 2 || X.ncol() < 1)
            stop("SRUW selection requires at least two rows and one variable.");
        if (Order.empty())
            stop("The SRUW ordering must contain at least one variable.");
        validate_variable_indices(Order, X.ncol(), "Order");
        if (packSize < 1)
            stop("The stopping threshold must be a positive integer.");
        Vect v(X);
        SelectReg sReg(v);
        CritClust b(nbCluster, framework, model_name, X, Crit, knownlabels, as<bool>(DA));
        Select s(v, b, sReg, packSize, stoppingRule);
        return wrap(s.selectS(Order));
    } catch (const SRUWStatusException &ex) {
        return List::create(
            Named("error") = ex.what(),
            Named("errorClass") = ex.errorClass()
        );
    } catch (std::exception &ex) {
        return List::create(Named("error") = ex.what());
    } catch (...) {
        return List::create(Named("error") = "Unknown error occurred in rcppSelectS");
    }
}

//[[Rcpp::export]]
IntegerVector rcppSelectW(NumericMatrix X, std::vector<int> Order,
                          std::vector<int> OtherVar, const int packSize,
                          std::string stoppingRule = "consecutive"){
    validate_variable_indices(Order, X.ncol(), "Order");
    validate_variable_indices(OtherVar, X.ncol(), "OtherVar");
    Vect v(X);
    SelectReg sReg(v);
    Select s(v, sReg, packSize, stoppingRule);
    IntegerVector result = wrap(s.selectW(Order, OtherVar));
    result.attr("stoppingRule") = stoppingRule;
    result.attr("stoppingThreshold") = packSize;
    result.attr("nEvaluated") = s.getLastWEvaluated();
    result.attr("stopReason") = s.getLastWStopReason();
    return result;
}


//[[Rcpp::export]]
IntegerVector rcppSelectR(NumericMatrix X, std::vector<int> S, std::vector<int> U, std::string regmodel){
    validate_variable_indices(S, X.ncol(), "S");
    validate_variable_indices(U, X.ncol(), "U");
    int nummodel = 0;
    if(regmodel == "LI")
        nummodel = 1;
    else
        if(regmodel == "LB")
            nummodel = 2;
        else if(regmodel == "LC")
            nummodel = 3;
        else
            stop("Unknown SRUW regression model '%s'. Use LI, LB, or LC.",
                 regmodel.c_str());
    
    Vect v(X);
    int InitialProjectsNb = v.experiments.size();
    SelectRegGen sRegGen(v);
    try {
        vector<int> varReg =
            sRegGen.selectReggen(S, U, nummodel, InitialProjectsNb);
        IntegerVector result = wrap(varReg);
        result.attr("selectionStatus") = "bounded_converged";
        result.attr("iterationCap") =
            sruw_regression_iteration_cap(S.size());
        result.attr("cycleHandling") = "fail_on_repeated_full_state";
        return result;
    } catch (const SRUWStatusException &ex) {
        stop("%s: %s", ex.errorClass().c_str(), ex.what());
    }

    return IntegerVector(0);
}

// Deterministic reference for consecutive and blockwise stopping rules.
// A non-finite token represents a failed fit and is never counted as a
// non-positive criterion comparison.
//[[Rcpp::export]]
List rcppSelectionTrace(NumericVector differences, const int c,
                        std::string stoppingRule = "consecutive",
                        std::string direction = "forward") {
    if (c < 1)
        stop("The stopping threshold must be a positive integer.");
    if (stoppingRule != "consecutive" && stoppingRule != "legacy_block")
        stop("Unknown stopping rule '%s'.", stoppingRule.c_str());
    if (direction != "forward" && direction != "reverse")
        stop("Unknown traversal direction '%s'.", direction.c_str());

    std::vector<int> traversal;
    traversal.reserve(differences.size());
    if (direction == "forward") {
        for (int i = 0; i < differences.size(); ++i)
            traversal.push_back(i);
    } else {
        for (int i = static_cast<int>(differences.size()) - 1; i >= 0; --i)
            traversal.push_back(i);
    }

    std::vector<int> processed, accepted, countTrace;
    int stopIndex = NA_INTEGER;
    std::string reason = "end_of_order";
    bool valid = true;

    if (stoppingRule == "consecutive") {
        ConsecutiveStopper stopper(c);
        for (int originalIndex : traversal) {
            processed.push_back(originalIndex + 1);
            const double difference = differences[originalIndex];
            if (!std::isfinite(difference)) {
                valid = false;
                reason = "invalid_fit";
                stopIndex = originalIndex + 1;
                countTrace.push_back(stopper.count());
                break;
            }
            const bool improvement = difference > 0.0;
            if (improvement)
                accepted.push_back(originalIndex + 1);
            const bool stopNow = stopper.observe(improvement);
            countTrace.push_back(stopper.count());
            if (stopNow) {
                reason = "c_nonpositive";
                stopIndex = originalIndex + 1;
                break;
            }
        }
    } else {
        int traversalIndex = 0;
        int withinBlockNonPositive = 0;
        while (traversalIndex < static_cast<int>(traversal.size())) {
            const int blockEnd = std::min(
                traversalIndex + c, static_cast<int>(traversal.size()));
            bool blockImprovement = false;
            for (; traversalIndex < blockEnd; ++traversalIndex) {
                const int originalIndex = traversal[traversalIndex];
                processed.push_back(originalIndex + 1);
                const double difference = differences[originalIndex];
                if (!std::isfinite(difference)) {
                    valid = false;
                    reason = "invalid_fit";
                    stopIndex = originalIndex + 1;
                    countTrace.push_back(withinBlockNonPositive);
                    break;
                }
                const bool improvement = difference > 0.0;
                if (improvement) {
                    accepted.push_back(originalIndex + 1);
                    blockImprovement = true;
                    withinBlockNonPositive = 0;
                } else {
                    ++withinBlockNonPositive;
                }
                countTrace.push_back(withinBlockNonPositive);
            }
            if (!valid)
                break;
            if (!blockImprovement) {
                reason = "legacy_empty_block";
                stopIndex = processed.back();
                break;
            }
        }
    }

    return List::create(
        Named("processed") = wrap(processed),
        Named("accepted") = wrap(accepted),
        Named("stop_index") = stopIndex,
        Named("count_trace") = wrap(countTrace),
        Named("termination_reason") = reason,
        Named("valid") = valid
    );
}

// Regression-criterion entry point used by independent numerical checks.
//[[Rcpp::export]]
List rcppRegressionBIC(NumericMatrix X, std::vector<int> response,
                       std::vector<int> predictors, std::string model) {
    int modelCode = 0;
    if (model == "LI") modelCode = 1;
    else if (model == "LB") modelCode = 2;
    else if (model == "LC") modelCode = 3;
    else stop("Unknown SRUW regression model '%s'. Use LI, LB, or LC.",
              model.c_str());
    Vect v(X);
    return v.bicReggen(response, predictors, modelCode);
}

//[[Rcpp::export]]
List rcppCrit(NumericMatrix X, List MyList, std::vector<std::string> rgm, std::vector<std::string> idm){
    typedef vector<int> stdivec;
    typedef vector<double> stddvec;
    Vect v(X);
    int InitialProjectsNb = v.experiments.size();
    SelectReg sReg(v);
    SelectRegGen sRegGen(v);
    int rhat = 0, lhat = 0, initsave = 0;
    mat reg;
    stdivec varSelectClust, varIndep, varNonIndep, varReg, SFinal, RFinal, UFinal, WFinal, Empty, regmodel, indepmodel;
    long double critClustFinal, BicRegFinal, crit, Lmax;
    stddvec BicIndepFinal;
    BicIndepFinal.clear(); BicRegFinal=0.0; crit=0.0; Lmax = 0.0;

    // Validate and unpack the selected role partition.
    try {
        if (MyList.containsElementNamed("error") &&
            !Rf_isNull(MyList["error"]))
            stop("Cannot evaluate a failed SRUW selection result: %s",
                 as<std::string>(MyList["error"]).c_str());
        varSelectClust = as<stdivec>(MyList["S"]);
        varNonIndep = as<stdivec>(MyList["U"]);
        varIndep = as<stdivec>(MyList["W"]);     
        critClustFinal = as<double>(MyList["criterionValue"]);
        if (!std::isfinite(static_cast<double>(critClustFinal)))
            stop("The clustering criterion is non-finite.");
    } catch (std::exception& e) {
        Rcerr << "Error accessing List elements: " << e.what() << std::endl;
        return List::create(Named("error") = e.what());
    }
  
    for(int p = 0; p < (int)rgm.size(); ++p)
    {
        if(rgm[p] == "LI")
            regmodel.push_back(1);
        else if(rgm[p] == "LB")
            regmodel.push_back(2);
        else if(rgm[p] == "LC")
            regmodel.push_back(3);
        else
            stop("Unknown SRUW regression model '%s'.", rgm[p].c_str());
    }
    
    for(int p = 0; p < (int)idm.size(); ++p)
    {
        if(idm[p] == "LI")
            indepmodel.push_back(1);
        else if(idm[p] == "LB")
            indepmodel.push_back(2);
        else
            stop("Unknown SRUW independent model '%s'.", idm[p].c_str());
    }
    if (!varNonIndep.empty() && regmodel.empty())
        stop("At least one regression model is required when U is non-empty.");
    if (!varIndep.empty() && indepmodel.empty())
        stop("At least one independent model is required when W is non-empty.");
    if (varIndep.size()==0)                 // W is empty
        if (varNonIndep.size()==0)          // U is empty
        {
            crit = critClustFinal;
            if ((initsave==0) || ((initsave==1) & (crit>Lmax)))
            {
                initsave=1;
                SFinal=varSelectClust; RFinal=Empty; UFinal=Empty; WFinal=Empty;
                rhat=0; lhat=0;
                Lmax = crit;
            }
        }
        else                                // U is nonempty and W is empty
            if (varNonIndep.size()==1)      // univariate redundant block
            {
                varReg=sReg.selectReg(varSelectClust,varNonIndep,InitialProjectsNb);
                List mylist = v.bicReggen(varNonIndep,varReg,regmodel[0]);
                BicRegFinal= mylist["bicvalue"];
                crit = critClustFinal + BicRegFinal;
                if ((initsave==0) || ((initsave==1) & (crit>Lmax)))
                {
                    initsave=1;
                    SFinal=varSelectClust; RFinal=varReg; UFinal=varNonIndep; WFinal=Empty;
                    rhat=1; 
                    lhat=0;
                    reg = as<mat>(mylist("B"));
                    Lmax=crit;
                }
            }
            else                              // multivariate redundant block
            {   
                for (int p=0; p < (int)regmodel.size();++p)
                {
                    varReg=sRegGen.selectReggen(varSelectClust,varNonIndep,regmodel[p],InitialProjectsNb);
                   List mylist =v.bicReggen(varNonIndep,varReg,regmodel[p]); 
                    BicRegFinal= mylist["bicvalue"];
                    crit = critClustFinal + BicRegFinal;
                    if ((initsave==0) || ((initsave==1) & (crit>Lmax)))
                    {
                        initsave=1;
                        SFinal=varSelectClust; RFinal=varReg; UFinal=varNonIndep; WFinal=Empty;
                        reg = as<mat>(mylist("B"));
                        rhat=regmodel[p]; lhat=0;
                        Lmax=crit;
                    }
                }
            }
            else                                // W is nonempty
            {
                // Evaluate each independent-block covariance model once.
                for (int l=0; l < (int)indepmodel.size();++l)
                {   List mylist = v.bicReggen(varIndep,Empty,indepmodel[l]);
                    BicIndepFinal.push_back(mylist["bicvalue"]);
                }  
                if (varNonIndep.size()==0)         // U is empty
                {
                    for (int l=0; l < (int)indepmodel.size();++l)
                    {
                        crit = critClustFinal + BicIndepFinal[l];
                        if ((initsave==0) || ((initsave==1) & (crit>Lmax)))
                        {
                            initsave=1;
                            SFinal=varSelectClust; RFinal=Empty; UFinal=Empty; WFinal=varIndep;
                            rhat=0; lhat=indepmodel[l];
                            Lmax=crit;
                            
                        }
                    }
                }
                else                                 // U is nonempty
                {
                    if (varNonIndep.size()==1)
                    {
                        varReg=sReg.selectReg(varSelectClust,varNonIndep,InitialProjectsNb);
                        List mylist = v.bicReggen(varNonIndep,varReg,regmodel[0]); 
                        BicRegFinal= mylist["bicvalue"];
                        for (int l=0; l < (int)indepmodel.size();++l)
                        {
                            crit=critClustFinal + BicRegFinal + BicIndepFinal[l];
                            if ((initsave==0) || ((initsave==1) & (crit>Lmax)))
                            {
                                initsave=1;
                                SFinal=varSelectClust; RFinal=varReg; UFinal=varNonIndep; WFinal=varIndep;
                                reg = as<mat>(mylist("B"));
                                rhat=1; lhat=indepmodel[l];
                                Lmax=crit;
                                
                            }
                        }
                    }
                    else
                    {
                        for (int p=0; p < (int)regmodel.size();++p)
                        {
                            varReg=sRegGen.selectReggen(varSelectClust,varNonIndep,regmodel[p],InitialProjectsNb);
                            List mylist = v.bicReggen(varNonIndep,varReg,regmodel[p]); 
                            BicRegFinal = mylist["bicvalue"];
                            for (int l=0; l < (int)indepmodel.size();++l)
                            {
                                crit = critClustFinal + BicRegFinal+ BicIndepFinal[l]; 
                                if ((initsave==0) || ((initsave==1) & (crit>Lmax)))
                                {
                                    initsave=1;                            
                                    SFinal=varSelectClust; RFinal=varReg; UFinal=varNonIndep; WFinal=varIndep; 
                                    rhat=regmodel[p]; lhat=indepmodel[l];
                                    reg = as<mat>(mylist("B"));
                                    Lmax=crit;
                                } 
                            }
                        }
                    }
                }
            }
    
    string rhats = "", lhats = ""; 
    if(rhat == 1)
        rhats = "LI";
    if(rhat == 2)
        rhats = "LB";
    if(rhat == 3)
        rhats = "LC";
    if(lhat == 1)
        lhats = "LI";
    if(lhat == 2)
        lhats = "LB";
    
    return List::create(Named("S") = wrap(SFinal), 
                        Named("R") = wrap(RFinal), 
                        Named("U") = wrap(UFinal), 
                        Named("W") = wrap(WFinal),
                        Named("criterionValue") = Lmax, 
                        Named("criterion") = MyList["criterion"],
                        Named("nbcluster") = MyList["nbcluster"],
                        Named("model") = MyList["model"],
                        Named("framework") = MyList.containsElementNamed("framework") ? MyList["framework"] : R_NilValue,
                        Named("requestedModel") = MyList.containsElementNamed("requestedModel") ? MyList["requestedModel"] : R_NilValue,
                        Named("effectiveModel") = MyList.containsElementNamed("effectiveModel") ? MyList["effectiveModel"] : R_NilValue,
                        Named("rmodel") = rhats, 
                        Named("imodel") = lhats,
                        Named("parameters") = MyList["parameters"],
                        Named("proba") = MyList["proba"],
                        Named("partition") = MyList["partition"],
                        Named("missingValues") = MyList.containsElementNamed("missingValues") ? MyList["missingValues"] : R_NilValue,
                        Named("regparameters")= wrap(reg),
                        Named("stoppingRule") = MyList.containsElementNamed("stoppingRule") ? MyList["stoppingRule"] : R_NilValue,
                        Named("stoppingThreshold") = MyList.containsElementNamed("stoppingThreshold") ? MyList["stoppingThreshold"] : R_NilValue,
                        Named("nEvaluated") = MyList.containsElementNamed("nEvaluated") ? MyList["nEvaluated"] : R_NilValue,
                        Named("stopReason") = MyList.containsElementNamed("stopReason") ? MyList["stopReason"] : R_NilValue,
                        Named("wStoppingRule") = MyList.containsElementNamed("wStoppingRule") ? MyList["wStoppingRule"] : R_NilValue,
                        Named("wStoppingThreshold") = MyList.containsElementNamed("wStoppingThreshold") ? MyList["wStoppingThreshold"] : R_NilValue,
                        Named("wNEvaluated") = MyList.containsElementNamed("wNEvaluated") ? MyList["wNEvaluated"] : R_NilValue,
                        Named("wStopReason") = MyList.containsElementNamed("wStopReason") ? MyList["wStopReason"] : R_NilValue);  
}

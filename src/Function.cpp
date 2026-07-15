#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <Rdefines.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

#include "Function.h"

long double Quad_Form(const colvec& x, const colvec& mu,
                      const mat& S_Inv){
    const colvec centered = x - mu;
    const colvec transformed = S_Inv * centered;
    return as_scalar(trans(centered) * transformed);
};

long double ldcppmvt(const colvec& x, const colvec& mu,
                     const mat& SInv, double SLogDet){
    const long double log2pi = log(2.0 * M_PI);
    int xdim = x.size();
    long double constants = -0.5 * xdim * log2pi ;
    long double Qf =  Quad_Form(x, mu, SInv);
    long double lret = constants - (0.5 * SLogDet) - (0.5 * Qf);
    return(lret);
    
};

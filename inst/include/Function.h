#ifndef ____Function__
#define ____Function__

#include <iostream>

long double Quad_Form(const colvec& x, const colvec& mu, const mat& S_Inv);

long double ldcppmvt(const colvec& x, const colvec& mu,
                     const mat& SInv, double SLogDet);
#endif /* defined(____Function__) */

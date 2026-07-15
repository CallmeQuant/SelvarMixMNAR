#ifndef SELECTREG_H 
#define SELECTREG_H 

#include "vect.hpp"

class SelectReg{

  Vect v;    


public:
  SelectReg();  
  SelectReg(Vect v); 

  // One backward and one forward update for the LI regression model.
  void exclusion_reg(vector<int>& varSelectReg, vector<int>& varNonSig,vector<int>& jE, vector<int>& jI, int& stop, int& InitialProjectsNb);

  void inclusion_reg(vector<int> varSelect, vector<int>& varSelectReg, vector<int>& varNonSig,vector<int>& jE, vector<int>& jI, int& arret, int& InitialProjectsNb);

  // Alternating subset search for the common spherical residual model.
  vector<int> selectReg(vector<int> varSelect,vector<int>& varNonSig, int& InitialProjectsNb);

};

#endif

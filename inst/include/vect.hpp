#ifndef VECT_H
#define VECT_H

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <Rdefines.h>

#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace std;
using namespace Rcpp;
using namespace arma;


class Vect {
 
private:
  NumericMatrix Data;

public:
  vector<int> experiments;  // zero-based variable indices
 
  Vect();
  Vect(NumericMatrix Data,vector<int> experiments);
  Vect(NumericMatrix Data);

  // Initialize the full variable index set.
  void initExperiments();                
  mat const_matrix(vector<int> vecteur);
  
  // Maximization-form regression BIC for response and predictor blocks.
  List bicReggen(vector<int> vectH, vector<int> vectY, int numr);
  // Set difference and ordered union on variable-index vectors.
  vector<int> enlever_var(vector<int>& vecteur, vector<int>& varenv);
  vector<int> ajouter_var(vector<int>& vecteur, vector<int>& varajout);
  
};
#endif

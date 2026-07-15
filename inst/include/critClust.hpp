#ifndef BICCLUST_H 
#define BICCLUST_H 

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
#include "sruwException.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

class CritClust {
  
  // Attributes
  string crit;
  int k; 
  string framework;          // "Mclust", "Rmixmod", "MixAll"
  string model_name;         // "VVV", "gaussian_pk_sjk", etc.
  NumericMatrix data;        // data matrix
  IntegerVector knownlabels; // labels for supervised classification
  bool DA;  
  
public:
  // Constructors
  CritClust();
  CritClust(int k, string framework, string model_name, NumericMatrix data, string crit, IntegerVector knownlabels, bool DA);
  
  // Method
  List ClustBestModel(vector<int> numExp);
};

#endif

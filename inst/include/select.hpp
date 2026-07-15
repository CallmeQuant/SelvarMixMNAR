#ifndef SELECT_H 
#define SELECT_H 
#include "vect.hpp" 
#include "critClust.hpp"
#include "selectReg.hpp"

class ConsecutiveStopper {
  int threshold;
  int consecutiveNonPositive;

public:
  explicit ConsecutiveStopper(int threshold)
    : threshold(threshold), consecutiveNonPositive(0) {
    if (threshold < 1)
      Rcpp::stop("The stopping threshold must be a positive integer.");
  }

  bool observe(bool improvement) {
    if (improvement) {
      consecutiveNonPositive = 0;
      return false;
    }
    ++consecutiveNonPositive;
    return consecutiveNonPositive >= threshold;
  }

  int count() const { return consecutiveNonPositive; }
};

class Select{
  
  Vect v;
  CritClust b;
  SelectReg sReg;
  int packSize;
  string stoppingRule;
  int lastWEvaluated;
  string lastWStopReason;
public:

  /**** Constructors ****/
  Select();
  Select(Vect v, CritClust b, SelectReg sReg, int packSize,
         string stoppingRule = "consecutive");
  Select(Vect v, SelectReg sReg, int packSize,
         string stoppingRule = "consecutive");

  /**** Methods ****/
  List selectS(vector<int> Order);
  vector<int> selectW(vector<int> Order, vector<int> OtherVar);
  int getLastWEvaluated() const { return lastWEvaluated; }
  string getLastWStopReason() const { return lastWStopReason; }
  
};

#endif

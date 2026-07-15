#include "selectRegGen.hpp"
#include "sruwException.hpp"
#include "sruwSearchControl.hpp"

#include <set>

SelectRegGen::SelectRegGen(){}
 
SelectRegGen::SelectRegGen(Vect v)
{ 
  this->v = v;
}

// Remove the predictor whose deletion gives the smallest loss under the
// selected multivariate regression covariance model.
void SelectRegGen::exclusion_reggen(vector<int>& varSelectReg, vector<int>& varNonSig,vector<int>& jE, vector<int>& jI, int& stopreg, int& nummodel, int& InitialProjectsNb)
{
  List mylist = (this->v).bicReggen(varNonSig,varSelectReg,nummodel);  
  double bicRegTotal = as<double>(mylist["bicvalue"]);
  double bicDiffReg = 0.0;
  vector<int> aux;
  aux.push_back(varSelectReg[0]); 
  vector<int> numExpAux = (this->v).enlever_var(varSelectReg,aux); 
  vector<int> jEmin;
  jEmin.push_back(varSelectReg[0]);
  mylist = (this->v).bicReggen(varNonSig,numExpAux,nummodel);  
  bicDiffReg = bicRegTotal - as<double>(mylist["bicvalue"]);   
  
  aux.clear(); numExpAux.clear();
  double bicDiffReg_aux = 0.0;
  for (int j=1; j < (int)varSelectReg.size();++j)
     {
        aux.push_back(varSelectReg[j]);       
        numExpAux = (this->v).enlever_var(varSelectReg,aux);    
        List mylist = (this->v).bicReggen(varNonSig,numExpAux,nummodel);
        bicDiffReg_aux = bicRegTotal -  as<double>(mylist["bicvalue"]);
    
        if (bicDiffReg_aux<=bicDiffReg)
          {
             bicDiffReg = bicDiffReg_aux;
             jEmin.clear();
             jEmin.push_back(varSelectReg[j]);
          }
          
        aux.clear(); numExpAux.clear();
     }
 
  if (bicDiffReg<=0)
    {       
       varSelectReg = (this->v).enlever_var(varSelectReg,jEmin);  
       jE.clear();
       jE.push_back(jEmin[0]);    
       if (jE==jI)
         stopreg = 1; 
       else     
         stopreg = 0; 
    }
  else
    {
       jE.clear();
       if (jI.empty())
         stopreg = 1; 
       else
         stopreg = 0;   
    }
}


// Add the excluded predictor that gives the largest positive BIC increment.
void SelectRegGen::inclusion_reggen(vector<int> varSelect, vector<int>& varSelectReg, vector<int>& varNonSig,vector<int>& jE, vector<int>& jI, int& stopreg, int& nummodel, int& InitialProjectsNb)
{
  List mylist = (this->v).bicReggen(varNonSig,varSelectReg,nummodel);
  double bicRegTotal = as<double>(mylist["bicvalue"]);
  vector<int> varSelectRegBis = (this->v).enlever_var(varSelect,varSelectReg);   
  if (varSelectRegBis.empty())
    {
      // No admissible inclusion remains. This is a converged boundary state,
      // not an indexable candidate set.
      jI.clear();
      stopreg = 1;
      return;
    }
  double bicDiffReg = 0.0;
  vector<int> aux;
  aux.push_back(varSelectRegBis[0]);   
  vector<int> numExpAux = (this->v).ajouter_var(varSelectReg,aux);  
  vector<int> jImax;
  jImax.push_back(varSelectRegBis[0]);
  mylist = (this->v).bicReggen(varNonSig,numExpAux,nummodel);
  bicDiffReg = -bicRegTotal + as<double>(mylist["bicvalue"]);   
 
  aux.clear(); numExpAux.clear();
  double bicDiffReg_aux = 0.0;
  for (int j=1; j < (int)varSelectRegBis.size();++j)
     {       
        aux.push_back(varSelectRegBis[j]);
        numExpAux = (this->v).ajouter_var(varSelectReg,aux);
        mylist = (this->v).bicReggen(varNonSig,numExpAux,nummodel);
        bicDiffReg_aux = -bicRegTotal + as<double>(mylist["bicvalue"]); 
        if (bicDiffReg_aux>bicDiffReg)
          {
             bicDiffReg = bicDiffReg_aux;
             jImax.clear();
             jImax.push_back(varSelectRegBis[j]);
          }
     aux.clear();
   }

   if (bicDiffReg>0)
     {
       if (jImax==jE)
         stopreg = 1;
       else
         {        
           varSelectReg = (this->v).ajouter_var(varSelectReg,jImax);   
           jI.clear();             
           jI.push_back(jImax[0]);      
           stopreg = 0;    
         }
     }  
   else
     {
        stopreg = 0;  
        jI.clear();
     }
}


// Alternate exclusion and inclusion until the regression subset is stable.
vector<int> SelectRegGen::selectReggen(vector<int> varSelect, vector<int>& varNonSig, int nummodel, int& InitialProjectsNb)
{
  vector<int> varSelectReg;
  varSelectReg = varSelect;
  vector<int> jI,jE;         
  int stopreg = 0;            
  const int iterationCap = sruw_regression_iteration_cap(varSelect.size());
  int iterations = 0;
  std::set<std::vector<int> > visitedStates;
  while (stopreg==0 && !varSelectReg.empty())
       {
         if (iterations >= iterationCap)
           throw SRUWStatusException(
             "SRUW multivariate regression search exceeded its deterministic iteration cap of " +
               std::to_string(iterationCap) + ".",
             "selvarmix_sruw_regression_iteration_cap"
           );
         const std::vector<int> state =
           sruw_regression_state_key(varSelectReg, jI, jE);
         if (!visitedStates.insert(state).second)
           throw SRUWStatusException(
             "SRUW multivariate regression search entered a repeated selection state.",
             "selvarmix_sruw_regression_cycle"
           );
         ++iterations;
         SelectRegGen::exclusion_reggen(varSelectReg,varNonSig,jE,jI,stopreg,nummodel,InitialProjectsNb); 
         if (stopreg==0)
            SelectRegGen::inclusion_reggen(varSelect,varSelectReg,varNonSig,jE,jI,stopreg,nummodel,InitialProjectsNb); 
       }
return varSelectReg;
}

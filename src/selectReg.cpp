#include "selectReg.hpp"
#include "sruwException.hpp"
#include "sruwSearchControl.hpp"

#include <set>

SelectReg::SelectReg(){}
 
SelectReg::SelectReg(Vect v)
{ 
  this->v = v;
}

// Remove the predictor whose deletion gives the smallest loss in the
// maximization-form regression BIC.
void SelectReg::exclusion_reg(vector<int>& varSelectReg, vector<int>& varNonSig, vector<int>& jE, vector<int>& jI, int& stop, int& InitialProjectsNb)
{
 
  const int numeromodeleaux=1;
  List mylist = v.bicReggen(varNonSig, varSelectReg, numeromodeleaux);
  double bicRegTotal = as<double>(mylist["bicvalue"]);
   
  double bicDiffReg = 0.0;
  vector<int> aux;
  aux.push_back(varSelectReg[0]); 
  vector<int> numProjets_aux = (this->v).enlever_var(varSelectReg,aux);
  vector<int> jEmin;
  jEmin.push_back(varSelectReg[0]);
  mylist = v.bicReggen(varNonSig, numProjets_aux, numeromodeleaux);
  bicDiffReg = bicRegTotal - as<double>(mylist["bicvalue"]);
  aux.clear(); numProjets_aux.clear();
  
  double bicDiffReg_aux = 0.0;
  for (int j=1; j < (int)varSelectReg.size();++j)
     {
       aux.push_back(varSelectReg[j]);       
       numProjets_aux = (this->v).enlever_var(varSelectReg,aux);
       List mylist = v.bicReggen(varNonSig, numProjets_aux, numeromodeleaux);
       bicDiffReg_aux = bicRegTotal - as<double>(mylist["bicvalue"]);

       if (bicDiffReg_aux<=bicDiffReg)
         {
            bicDiffReg = bicDiffReg_aux;
            jEmin.clear();
            jEmin.push_back(varSelectReg[j]);
         }
    
       aux.clear(); numProjets_aux.clear();
     }
  if (bicDiffReg<=0)
    {       
      varSelectReg = (this->v).enlever_var(varSelectReg,jEmin);  
      jE.clear();
      jE.push_back(jEmin[0]);    
      if (jE==jI)
         stop = 1; 
      else     
         stop = 0;    
    }
  else
    {
      jE.clear();
      if (jI.empty())
         stop = 1; 
      else
         stop = 0;     
    }

}


// Add the excluded predictor that gives the largest positive BIC increment.
void SelectReg::inclusion_reg(vector<int> varSelect, vector<int>& varSelectReg, vector<int>& varNonSig,vector<int>& jE,vector<int>& jI,int& stop, int& InitialProjectsNb)
{
  const int numeromodeleaux=1; 
  List mylist = v.bicReggen(varNonSig, varSelectReg, numeromodeleaux);
  double bicRegTotal = as<double>(mylist["bicvalue"]);   
  
  vector<int> varSelectRegBis = (this->v).enlever_var(varSelect,varSelectReg);   
  if (varSelectRegBis.empty())
    {
      // No admissible inclusion remains. This is a converged boundary state,
      // not an indexable candidate set.
      jI.clear();
      stop = 1;
      return;
    }
  double bicDiffReg = 0.0;
  vector<int> aux;
  aux.push_back(varSelectRegBis[0]);   
  vector<int> numProjets_aux = (this->v).ajouter_var(varSelectReg,aux);  
 
  vector<int> jImax;
  jImax.push_back(varSelectRegBis[0]);
  mylist = v.bicReggen(varNonSig, numProjets_aux, numeromodeleaux);
  bicDiffReg = -bicRegTotal + as<double>(mylist["bicvalue"]); 

  aux.clear(); numProjets_aux.clear();
  
  double bicDiffReg_aux = 0.0;
  for (int j=1; j < (int)varSelectRegBis.size();++j)
     {       
         aux.push_back(varSelectRegBis[j]);
         numProjets_aux = (this->v).ajouter_var(varSelectReg,aux);
         mylist = v.bicReggen(varNonSig, numProjets_aux, numeromodeleaux);
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
         stop = 1;
       else
       {        
         varSelectReg =  (this->v).ajouter_var(varSelectReg,jImax);   
         jI.clear();             
         jI.push_back(jImax[0]);      
         stop = 0;    
       } 
     }  
   else
     {
        stop = 0;  
        jI.clear();
     } 
}


// Alternate exclusion and inclusion until the LI regression subset is stable.
vector<int> SelectReg::selectReg(vector<int> varSelect,vector<int>& varNonSig, int& InitialProjectsNb)
{
  vector<int> varSelectReg;
  varSelectReg = varSelect;
  vector<int> jI,jE;         
  int stop = 0;           
  const int iterationCap = sruw_regression_iteration_cap(varSelect.size());
  int iterations = 0;
  std::set<std::vector<int> > visitedStates;
  
  while (stop==0 && !varSelectReg.empty())
       {
         if (iterations >= iterationCap)
           throw SRUWStatusException(
             "SRUW LI regression search exceeded its deterministic iteration cap of " +
               std::to_string(iterationCap) + ".",
             "selvarmix_sruw_regression_iteration_cap"
           );
         const std::vector<int> state =
           sruw_regression_state_key(varSelectReg, jI, jE);
         if (!visitedStates.insert(state).second)
           throw SRUWStatusException(
             "SRUW LI regression search entered a repeated selection state.",
             "selvarmix_sruw_regression_cycle"
           );
         ++iterations;
         SelectReg::exclusion_reg(varSelectReg,varNonSig,jE,jI,stop,InitialProjectsNb); 
         if (stop==0)
            SelectReg::inclusion_reg(varSelect,varSelectReg,varNonSig,jE,jI,stop,InitialProjectsNb); 
       }
       
return varSelectReg;
}

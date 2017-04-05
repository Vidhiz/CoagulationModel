
#include <iostream>
#include <map>
#include <string>
#include <vector>

/// ODESolver
class ODESolver
{
public:

  enum evars { PBz2ba, PBe2ba, PBz5ba, PBe5ba, PBz8ba, PBe8ba, PBz9ba, PBe9ba,
               PBz10ba, PBe10ba, PBten, PBpro, PBz2ba_pro, PBz5ba_e2ba, 
               PBz5ba_e10ba, PBz8ba_e2ba, PBz8ba_e10ba, PBz10ba_ten, 
               PBapc_e5ba, PBapc_e8ba, PBe9starba, PBtenstar, PBz10ba_tenstar, 
               PBPltsADP, PBPltse2, e2mtot, z2mtot, z5mtot, e5mtot, z8mtot, 
               e8mtot, z9mtot, e9mtot, e10mtot, z10mtot  
  };

  enum econsts {
    /*A09*/ k2on = 0, N2b, N2se, k2off, kz2promplus, kz2promminus,
    /*A10*/ kz2promcat, kz5e2mcat, kz5e2mminus, kz5e2mplus, kz8e2mcat,
    kz8e2mminus, kz8e2mplus,
    /*A11*/ k5on, N5b, N5se, k5off, kz5e10mplus, kz5e10mminus,
    /*A12*/ kz5e10mcat, kprominus, kproplus, kapce5mplus, kapce5minus,
    /*A13*/ k8on, N8b, N8se, k8off, kz8e10mplus, kz8e10mminus,
    /*A14*/ kz8e10cat, ktenminus, ktenplus, kapce8mplus, kapce8mminus,
    /*A15*/ k9on, N9b, N9se, k9off,
    /*A16*/
    /*A17*/ k10on, N10b, N10se, k10off, k10tenmminus, kz10tenmplus,
    /*A18*/ kz2e10mminus, kz10tenmcat,
    /*A19*/ kz10tenmminus,
    /*A20*/
    /*A21*/
    /*A22*/
    /*A23*/
    /*A24*/
    /*A25*/
    /*A26*/
    /*A27*/ kapce5mcat,
    /*A28*/ kapce8mcat,
    /*A29*/ N9starb, N9starse
    /*A30*/
    /*A31*/
  };

  typedef std::map<econsts, const double> ConstMap;
  typedef std::map<evars, std::vector<double> > VarMap;

  const static int numV = 35;
  const static int numC = 52;

  /// Constructor
  ODESolver()
  {
    std::cout<<"Created ODESolver.\n";
    //this->populateConstMap();
  }

  /// Destructor
  ~ODESolver()
  {
    std::cout<<"Destroyed ODESolver.\n";
  }

  /// getConstMap
  ConstMap& getConsts()
  {
    return this->consts;
  }

  /// getVarMap
  VarMap& getVars()
  {
    return this->vars;
  }


private:

  //consts
  ConstMap consts = { {k2on, 1.1}, {k2off, 4.4} };
  
  //vars
  VarMap vars;

};

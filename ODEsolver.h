
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
			   e8mtot, z9mtot, e9mtot, e10mtot, z10mtot, eLASTV  
  };

  enum econsts {
	/*A09*/ k2on, N2b, N2se, k2off, kz2promplus, kz2promminus,
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


private:

  //consts
  ConstMap consts = { {k2on, 1.0e7}, {N2b, 1000.0}, {N2se, 1000.0}, {k2off, 5.9},
					  {kz2promplus, 1.03e8}, {kz2promminus, 1.0},
		/*A10*/ {kz2promcat,30.0}, {kz5e2mcat,0.23}, {kz5e2mminus,1.0}, 
				{kz5e2mplus,1.73e7}, /**/{kz8e2mcat,11.11}, {kz8e2mminus,12.12},
				{kz8e2mplus, 13.13},
		/*A11*/ {k5on,14.14}, {N5b,15.15},{N5se,16.16}, {k5off,17.17}, 
				{kz5e10mplus,18.18}, {kz5e10mminus, 19.19},
		/*A12*/ {kz5e10mcat,20.20}, {kprominus,21.21}, {kproplus,1.0e8}, 
				{kapce5mplus,23.23}, {kapce5minus, 24.24},
		/*A13*/ {k8on, 25.25}, {N8b,26.26}, {N8se,27.27}, {k8off,28.28}, 
				{kz8e10mplus,29.29}, {kz8e10mminus,30.30},
		/*A14*/ {kz8e10cat,31.31}, {ktenminus,32.32}, {ktenplus,33.33}, 
				{kapce8mplus,34.34}, {kapce8mminus,35.35},
		/*A15*/ {k9on,36.36}, {N9b,37.37}, {N9se,38.38}, {k9off,39.39},
		/*A16*/
		/*A17*/ {k10on,40.40}, {N10b,41.41}, {N10se,42.42}, {k10off,43.43}, 
				{k10tenmminus,44.44}, {kz10tenmplus,45.45},
		/*A18*/ {kz2e10mminus,46.46}, {kz10tenmcat,47.47},
		/*A19*/ {kz10tenmminus,48.48}
   };
  

};

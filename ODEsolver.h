#ifndef _REACTIONS_ODESOLVER_H_
#define _REACTIONS_ODESOLVER_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "myEnums.h"
#include "Domain.h"
#include "Point.h"

namespace Reactions{


	typedef const std::vector<double> ConstMap;

	/// ODESolver
	class ODESolver
	{
	public:

	  const static int numV = 32;
	  const static int numC = 52;

	  // Final time and time-step
	  double T;
	  double dt, tchar; 

	  // Constants provided by user to solver
	  double Pba, Psea;

	  /// Default Constructor
	  ODESolver()
	  {
		std::cout<<"Created deafult ODESolver.\n";
		T = 5.0;
		dt = 0.001;
		tchar = 1e-5; // what should this be?
	  }

	  /// Parametrized constructor
	  ODESolver(double T, double dt)
	  {
		std::cout<<"Created parameterized ODESolver.\n";
		T = T;
		ODESolver::dt = dt;
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

	  void setPba(double pba)
	  {
	  	Pba = pba; 
	  }

	  void setPsea(double psea)
	  {
	  	Psea = psea;
	  }

	  double RK2solve(Domain &d);

	  // for RK2 co-efficient computations
	  void computeK(Mesh &m, Point &kp);



	private:

	  //consts
	  ConstMap consts = { 
	  	1.0e7, 1000.0, 1000.0, 5.9, 1.03e8, 1.0, // /*A09*/ k2on, N2b, N2se, k2off, kz2promplus, kz2promminus,
		30.0, 0.23, 1.0, 1.73e7, 11.11, 12.12, 13.13, // /*A10*/ kz2promcat, kz5e2mcat, kz5e2mminus, kz5e2mplus, kz8e2mcat, kz8e2mminus, kz8e2mplus
		14.14, 15.15, 16.16, 17.17, 18.18, 19.19, // /*A11*/ k5on, N5b, N5se, k5off, kz5e10mplus, kz5e10mminus,
		20.20, 21.21, 1.0e8, 23.23, 24.24, // /*A12*/ kz5e10mcat, kprominus, kproplus, kapce5mplus, kapce5minus,
	    25.25, 26.26, 27.27, 28.28, 29.29, 30.30, // /*A13*/ k8on, N8b, N8se, k8off, kz8e10mplus, kz8e10mminus,
		31.31, 32.32, 33.33, 34.34, 35.35, // /*A14*/ kz8e10cat, ktenminus, ktenplus, kapce8mplus, kapce8mminus,
		36.36, 37.37, 38.38, 39.39, // /*A15*/ k9on, N9b, N9se, k9off,
									// /*A16*/
		40.40, 41.41, 42.42, 43.43, 44.44, 45.45, // /*A17*/ k10on, N10b, N10se, k10off, k10tenmminus, kz10tenmplus,
		46.46, 47.47, // /*A18*/ kz2e10mminus, kz10tenmcat,
		48.48, // /*A19*/ kz10tenmminus,
		49.49,  // /*A27*/ kapce5mcat,
		50.50, // /*A28*/ kapce8mcat,
	  	51.51 // /*A29*/ N9starb, N9starse
	   };


	};
}
#endif //_REACTIONS_ODESOLVER_H_
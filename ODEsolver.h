#ifndef _REACTIONS_ODESOLVER_H_
#define _REACTIONS_ODESOLVER_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "Consts.h"
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
		//double Pba, Psea;

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
			ODESolver::T = T;
			ODESolver::dt = dt;
			tchar = dt/100;
		}


		/// Destructor
		~ODESolver()
		{
			std::cout<<"Destroyed ODESolver.\n";
		}

		/*void setPba(double pba)
		{
			Pba = pba; 
		}

		void setPsea(double psea)
		{
			Psea = psea;
		}*/

		double RK2solve(Domain &d);

		// for RK2 co-efficient computations
		void computeK(Point &p, Point &kp);



	private:

	  //consts
	  ConstMap consts = { 
	  	1.0e7, 1000.0, 1000.0, 5.9, 1.03e8, 1.0, // /*A09*/ k2on, N2b, N2se, k2off, kz2promplus, kz2promminus,
		30.0, 0.23, 1.0, 1.73e7, 0.9, 1.0, 2.64e7, // /*A10*/ kz2promcat, kz5e2mcat, kz5e2mminus, kz5e2mplus, kz8e2mcat, kz8e2mminus, kz8e2mplus
		5.7e7, 3000, 3000, 0.17, 1.0e8, 1.0, // /*A11*/ k5on, N5b, N5se, k5off, kz5e10mplus, kz5e10mminus,
		4.6e-2, 0.01, 1.0e8, 1.2e8, 1.0, // /*A12*/ kz5e10mcat, kprominus, kproplus, kapce5mplus, kapce5minus,
	    5.0e7, 450.0, 450.0, 0.17, 5.1e7, 1.0, // /*A13*/ k8on, N8b, N8se, k8off, kz8e10mplus, kz8e10mminus,
		2.3e-2, 0.01, 1.0e8, 1.2e8, 1.0, // /*A14*/ kz8e10cat, ktenminus, ktenplus, kapce8mplus, kapce8mminus,
		1.0e7, 250.0, 250.0, 2.5e-2, // /*A15*/ k9on, N9b, N9se, k9off,
									// /*A16*/
		1.0e7, 2700.0, 2700.0, 2.5e-2, 1.0, 1.31e8, // /*A17*/ k10on, N10b, N10se, k10off, kz10tenmminus, kz10tenmplus,
		20.0, // /*A18*/ kz10tenmcat,
	    // /*A19*/
		0.5,  // /*A27*/ kapce5mcat,
		0.5, // /*A28*/ kapce8mcat,
	  	250.0, 250.0 // /*A29*/ N9starb, N9starse
	   };

	};
}
#endif //_REACTIONS_ODESOLVER_H_

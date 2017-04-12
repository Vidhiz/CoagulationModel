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
			ODESolver::T = T;
			ODESolver::dt = dt;
			tchar = dt/100;
		}


		/// Destructor
		~ODESolver()
		{
		std::cout<<"Destroyed ODESolver.\n";
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
		void computeK(Point &p, Point &kp);



	private:

	};
}
#endif //_REACTIONS_ODESOLVER_H_

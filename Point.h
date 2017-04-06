
#include <iostream>
#include <map>
#include <string>
#include <vector>

class Point
{
	public:
		
		/// Default Constructor
		Point()
		{
			std::cout<<"Point object !";
  			Cvals.resize(ODESolver::eLASTV);
		    for (int i=0;i<ODESolver::eLASTV;i++)
		    {
		      Cvals[i]=0.0;
		    }
		}

		/// Destructor
		~Point()
		{

		}

		/// vector which contains the chemical concentrations
		std::vector<double> Cvals;
};
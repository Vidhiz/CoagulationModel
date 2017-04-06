
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace ODESolver{
class Point
{
	public:
		
		/// Default Constructor
		Point()
		{
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
}
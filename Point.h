
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
		}

		/// Destructor
		~Point()
		{

		}

		/// vector which contains the chemical concentrations
		std::vector<double> Cval;
};
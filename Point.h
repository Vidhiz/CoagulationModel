
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

		struct Pt
		{
	  		double x;
	  		double y;
	  		std::vector<double> C; // information of all the chemicals at this point
		};	

		/// vector which contains x coordinate, y coordinate and chemical concentrations
		std::vector<Pt> Ptdata;
};
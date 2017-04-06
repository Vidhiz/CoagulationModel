
#include <iostream>
#include <map>
#include <string>
#include <vector>


/// Domain
class Domain{
	public:
		/// Constructor
		Domain()
		{
			std::cout<<"Created ODESolver.\n";
			xlen = 2.0;
			ylen = 1.0;
			nx = 40;
			ny = 20;
		}

		/// Destructor
		~Domain()
		{
			std::cout<<"Destroyed ODESolver.\n";
		}	

		/// x - length
		double xlen;

		/// y - length
		double ylen;

		/// points in x
		double nx;

		/// points in y
		double ny;

		/// get and set variables:
		double get_nx()
		{
			return nx;
		}
		double get_ny()
		{
			return ny;
		}
		double get_xlen()
		{
			return xlen;
		}
		double get_ylen()
		{
			return ylen;
		}
		
};
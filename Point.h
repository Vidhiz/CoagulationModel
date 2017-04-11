#ifndef _REACTIONS_POINT_H_
#define _REACTIONS_POINT_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "ODEsolver.h"

namespace Reactions{
class Point
{
	public:
		
		/// Default Constructor
		Point();
		

		/// Destructor
		~Point();
		
		/// vector which contains the chemical concentrations
		std::vector<double> Cvals;
};
}

#endif //_REACTIONS_POINT_H_
#ifndef _REACTIONS_POINT_H_
#define _REACTIONS_POINT_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace Reactions
{

enum evars { PBz2ba, PBe2ba, PBz5ba, PBe5ba, PBz8ba, PBe8ba, PBz9ba, PBe9ba,
					   PBz10ba, PBe10ba, PBten, PBpro, PBz2ba_pro, PBz5ba_e2ba, 
					   PBz5ba_e10ba, PBz8ba_e2ba, PBz8ba_e10ba, PBz10ba_ten, PBapc,
					   PBapc_e5ba, PBapc_e8ba, PBe9starba, PBtenstar, PBz10ba_tenstar, 
					   PBPltsADP, PBPltse2, e2mtot, z2mtot, ze5mtot, ze8mtot, 
					   ze9mtot, ze10mtot, eLASTV  
};

class Point
{
public:
		
  /// Default Constructor
  Point();
		
  /// Destructor
  ~Point();
		
  ///chemical concentrations
  std::vector<double> Cvals;
};

}

#endif //_REACTIONS_POINT_H_

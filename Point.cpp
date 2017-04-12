#include "Consts.h"
#include "Point.h"

namespace Reactions
{

Point::Point()
{
  Cvals.resize(eLASTV);
  for (int i=0;i<eLASTV;i++)
  {
    Cvals[i]=0.0;
  }
}

/// Destructor
Point::~Point()
{

}

}

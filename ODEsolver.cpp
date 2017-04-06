// Platelet-bound chemical specied : Purely reaction

#include "ODEsolver.h"
#include "Domain.h"

using namespace std;

/////////////////////////////
int main(int argc, char *argv[])
{
  	ODESolver *solver = new ODESolver();
  	Domain *domain = new Domain(50,20,100,80);

  	ODESolver::ConstMap &consts = solver->getConsts();

	cout<<"Constants assigned as follows:"<<endl;
	for (const auto &kv : consts)
	{
		cout << kv.first<<" = "<< kv.second << endl;
	}



	return 0;
}

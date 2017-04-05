// Platelet-bound chemical specied : Purely reaction

#include "ODEsolver.h"

using namespace std;

/////////////////////////////
int main(int argc, char *argv[])
{
  ODESolver *solver = new ODESolver();
  ODESolver::ConstMap &consts = solver->getConsts();

	cout<<"Constants assigned as follows:"<<endl;
	for (const auto &kv : consts)
	{
		cout << kv.first<<" = "<< kv.second << endl;
	}

	return 0;
}
